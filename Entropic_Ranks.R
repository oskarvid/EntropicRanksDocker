arguments <- commandArgs(TRUE)
data_under_analysis <- read.table(arguments[1],sep="\t",dec=".",header=TRUE)
population_vector <- read.table(arguments[2],sep="\t",dec=".",header=TRUE)[,1]
if (arguments[3]!="null"){data_origin <- read.table(arguments[3],sep="\t",dec=".",header=TRUE)[,1]} else {data_origin <- NULL}
granularity <- as.integer(arguments[4])
supervised <- as.logical(arguments[5])
process_log <- as.logical(arguments[6])
export_plots <- as.logical(arguments[7])
create_output_files <- as.logical(arguments[8])
is_logged <- as.logical(arguments[9])
logbase <- as.integer(arguments[10])
huge_feature_list <- as.logical(arguments[11])

library("RankProd")
library("entropy")
library("factoextra")
setwd("/data")

entropic_ranks <- function(data_under_analysis,population_vector,data_origin,granularity,supervised,process_log,export_plots,create_output_files,is_logged,logbase,huge_feature_list)
{
  if (is.null(data_origin))
    data_origin <- rep(1,length(population_vector))
  
  message("Calculating Rank Products. May take a long time, depending on data set size.")
  comparison <- RPadvance(data_under_analysis,cl=population_vector,origin=data_origin,logged=is_logged,na.rm=FALSE,gene.names=rownames(data_under_analysis),plot=process_log,huge=TRUE)
  if (huge_feature_list){
    message("Investigating only the first 20000 features.")
    rank_product_lists <- topGene(comparison,num.gene=20000,logged=is_logged,logbase=logbase,gene.names=rownames(data_under_analysis))
  }else{
    rank_product_lists <- topGene(comparison,cutoff=0.99,method="pfp",logged=is_logged,logbase=logbase,gene.names=rownames(data_under_analysis))
  }
  
  if (export_plots)
  {
    if(!file.exists(paste(getwd(),"Entropic Ranks plots",sep="/")))
      dir.create(paste(getwd(),"Entropic Ranks plots",sep="/"))
    path_down <- paste(getwd(),"Entropic Ranks plots","Downregulated",sep="/")
    path_up <- paste(getwd(),"Entropic Ranks plots","Upregulated",sep="/")
  }
  if (!supervised)
  {
    if (create_output_files)
    {
      if (!is.null(rownames(rank_product_lists$Table1)))
        write.table(file="Downregulated list [original].txt",rank_product_lists$Table1[,3:5],sep="\t",quote=FALSE,row.names=TRUE,col.names=TRUE)
      if (!is.null(rownames(rank_product_lists$Table2)))
        write.table(file="Upregulated list [original].txt",rank_product_lists$Table2[,3:5],sep="\t",quote=FALSE,row.names=TRUE,col.names=TRUE)
    }
    
    message("Trimming down the list of downregulated features in unsupervised mode...")
    rank_product_lists$Table1 <- rank_product_lists$Table1[1:isolate_significant_elements(rank_product_lists$Table1[,2],granularity,supervised,process_log,export_plots,path=path_down),]
    message("Trimming down the list of upregulated features in unsupervised mode...")
    rank_product_lists$Table2 <- rank_product_lists$Table2[1:isolate_significant_elements(rank_product_lists$Table2[,2],granularity,supervised,process_log,export_plots,path=path_up),]
    
    if (create_output_files)
    {
      if (!is.null(rownames(rank_product_lists$Table1)))
        write.table(file="Downregulated list [information-dense].txt",rank_product_lists$Table1[,3:5],sep="\t",quote=FALSE,row.names=TRUE,col.names=TRUE)
      if (!is.null(rownames(rank_product_lists$Table2)))
        write.table(file="Upregulated list [information-dense].txt",rank_product_lists$Table2[,3:5],sep="\t",quote=FALSE,row.names=TRUE,col.names=TRUE)
      message("Output files created successfully.")
    }
  } else {
    message("Calculating suggested cutoff points for the downregulated feature list...")
    isolate_significant_elements(rank_product_lists$Table1[,2],granularity,supervised,process_log,export_plots,path=path_down)
    message("Calculating suggested cutoff points for the upregulated feature list...")
    isolate_significant_elements(rank_product_lists$Table2[,2],granularity,supervised,process_log,export_plots,path=path_up)
  }
  return(rank_product_lists)
}

isolate_significant_elements <- function(ordered_vector,granularity=1,supervised=FALSE,process_log=FALSE,export_plots=FALSE,path=NULL)
{
  bin_min <- 15
  bin_max <- 35
  bin_increment <- 5
  window_min <- 60
  window_max <- 250
  window_increment <- 10
  
  suggested_surfaces <- c()
  if (process_log)
  {
    par(mfrow=c(2,1))
  } else {
    progress <- txtProgressBar(max=((bin_max-bin_min)/bin_increment +1) * ((window_max-window_min)/window_increment +1),char="=",style=3)
  }
  
  for (i in seq(from=bin_min, to=bin_max, by=bin_increment))
    for (j in seq(from=window_min, to=window_max, by=window_increment))
    {
      if (process_log)
      {
        cat("Calculating using ",i," bins and a sliding window of ",j," features","\n",sep="")
      } else {
        setTxtProgressBar(progress,(i-bin_min)/bin_increment*((window_max-window_min)/window_increment + 1) + (j-window_min)/window_increment + 1)
      }
      suggested_surfaces <- c(suggested_surfaces,entropic_analysis(ordered_vector,step_up=granularity,bins=i,window_size=j,verbose=process_log,export_plots,path))
    }
  if (process_log)
  {
    par(mfrow=c(1,1))
  } else {
    close(progress)
  }
  
  if (supervised||process_log)
  {
    print(table(suggested_surfaces))
    message("Most consistent cutoff point: feature no ",as.integer(rownames(table(suggested_surfaces))[table(suggested_surfaces) == max(table(suggested_surfaces))])[1],".")
  }
  if (!supervised)
    return(as.integer(rownames(table(suggested_surfaces)))[table(suggested_surfaces) == max(table(suggested_surfaces))][1])
}

entropic_analysis <- function(ordered_vector,step_up=1,window_size,bins,verbose=FALSE,export_plots=FALSE,path=NULL)
{
  if (export_plots)
  {
    if (is.null(path))
      path <- paste(getwd(),"Entropic Ranks plots",sep="/")
    if (!file.exists(path))
      dir.create(path)
    write_file <- file.exists(paste(path,"Sliding window distributions.txt",sep="/"))
    if (write_file)
      write_file <- length(read.table(paste(path,"Sliding window distributions.txt",sep="/"),nrows=1,header=FALSE))<42
  }
  
  differences <- ordered_vector[seq(2,length(ordered_vector))]-ordered_vector[seq(1,length(ordered_vector)-1)]
  entropy_plotter <- vector(length=floor((length(differences)-window_size)/step_up))
  if (export_plots)
    if (write_file || !file.exists(paste(path,"Sliding window distributions.txt",sep="/")))
      mean_differences <- c()
  for (i in 0:(length(entropy_plotter)-1))
  {
    entropy_plotter[i+1] <- entropy(discretize(differences[(i*step_up+1):(i*step_up+window_size)],numBins=bins),method="Laplace")
    if (export_plots)
      if (write_file || !file.exists(paste(path,"Sliding window distributions.txt",sep="/")))
        mean_differences <- c(mean_differences,mean(differences[(i*step_up+1):(i*step_up+window_size)]))
  }
  entropy_clusters <- eclust(entropy_plotter, "kmeans", k=2, nstart=200, graph=FALSE)
  
  if (verbose)
  {
    cat("Calculating entropies of ",length(entropy_plotter)," ovelapping windows","\n","Suggested cutoff at feature no ",seq(length(entropy_clusters$cluster))[entropy_clusters$cluster!=as.integer(entropy_clusters$centers==max(entropy_clusters$centers))[1]+1][1]-1,", at a mean entropy of ",mean(entropy_plotter[1:seq(length(entropy_clusters$cluster))[entropy_clusters$cluster!=as.integer(entropy_clusters$centers==max(entropy_clusters$centers))[1]+1][1]-1]),"\n","Last 1/3 minimum entropy: ",min(entropy_plotter[floor(length(entropy_plotter)*2/3):length(entropy_plotter)]),"\n",sep="")
    barplot(entropy_plotter,border=c("gold1","dodgerblue3")[as.vector(as.integer(entropy_clusters$cluster==(as.integer(entropy_clusters$centers==max(entropy_clusters$centers))[1]+1))+1)],col=c("gold1","dodgerblue3")[as.vector(as.integer(entropy_clusters$cluster==(as.integer(entropy_clusters$centers==max(entropy_clusters$centers))[1]+1))+1)],xlab=c("Granularity = ",step_up),ylab=c("Window size = ",window_size),main=c(bins," bins"),names.arg=seq(length(entropy_plotter))*step_up)
    barplot(entropy_clusters$silinfo$widths$sil_width,border=c("gold1","dodgerblue3")[as.integer(entropy_clusters$silinfo$widths$cluster==(as.integer(entropy_clusters$centers==max(entropy_clusters$centers))[1]+1))+1],col=c("gold1","dodgerblue3")[as.integer(entropy_clusters$silinfo$widths$cluster==(as.integer(entropy_clusters$centers==max(entropy_clusters$centers))[1]+1))+1],ylim=c(-0.3,1),ylab="Silhouette width",xlab="Clustered elements",main="K-means clustering quality")
  }
      
  if (export_plots)
  {
    if (!file.exists(paste(path,"Sliding window distributions.txt",sep="/")))
    {
      write.table(file=paste(path,"Sliding window distributions.txt",sep="/"),cbind(differences,differences/max(differences)),row.names=seq(length(differences)),col.names=c("Differences","Ratio_to_maximum"),sep="\t",quote=FALSE)
      write_file <- TRUE
    }
    if (write_file)
    {
      distributions <- read.table(paste(path,"Sliding window distributions.txt",sep="/"),header=TRUE)
      write.table(file=paste(path,"Sliding window distributions.txt",sep="/"),cbind(distributions,c(mean_differences,rep(NA,dim(distributions)[1]-length(mean_differences))),c(mean_differences/max(mean_differences),rep(NA,dim(distributions)[1]-length(mean_differences)))),row.names=seq(dim(distributions)[1]),col.names=c(colnames(distributions),paste("Mean_difference_window_size",window_size,sep="_"),paste("Ratio_to_maximum_mean_window_size",window_size,sep="_")),sep="\t",quote=FALSE)
      rm(distributions,mean_differences)
    }
    png(file=file.path(path,paste("Entropy - ",bins," bins - ",window_size," window size",".png", sep = "")),width=640,height=640)
    barplot(entropy_plotter,border=c("gold1","dodgerblue3")[as.vector(as.integer(entropy_clusters$cluster==(as.integer(entropy_clusters$centers==max(entropy_clusters$centers))[1]+1))+1)],col=c("gold1","dodgerblue3")[as.vector(as.integer(entropy_clusters$cluster==(as.integer(entropy_clusters$centers==max(entropy_clusters$centers))[1]+1))+1)],xlab=c("Granularity = ",step_up),ylab=c("Window size = ",window_size),main=c(bins," bins"),names.arg=seq(length(entropy_plotter))*step_up)
    dev.off()
    png(file=file.path(path,paste("Clustering quality - ",bins," bins - ",window_size," window size",".png", sep = "")),width=640,height=640)
    barplot(entropy_clusters$silinfo$widths$sil_width,border=c("gold1","dodgerblue3")[as.integer(entropy_clusters$silinfo$widths$cluster==(as.integer(entropy_clusters$centers==max(entropy_clusters$centers))[1]+1))+1],col=c("gold1","dodgerblue3")[as.integer(entropy_clusters$silinfo$widths$cluster==(as.integer(entropy_clusters$centers==max(entropy_clusters$centers))[1]+1))+1],ylim=c(-0.3,1),ylab="Silhouette width",xlab="Clustered elements",main="K-means clustering quality")
    dev.off()
  }
  return(seq(length(entropy_clusters$cluster))[entropy_clusters$cluster!=as.integer(entropy_clusters$centers==max(entropy_clusters$centers))[1]+1][1]-1)
}

entropic_ranks(data_under_analysis,population_vector,data_origin,granularity,supervised,process_log,export_plots,create_output_files,is_logged,logbase,huge_feature_list)

#Code written by Hector-Xavier de Lastic
#Development & testing by Hector-Xavier de Lastic & Irene Liampa
#Contact:
#hector.xavier.de.lastic@gmail.com
#irini.liampa@gmail.com
