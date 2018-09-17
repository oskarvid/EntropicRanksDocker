require("RankProd")
#reliable compatibility with RankProd versions up to 2.44.0
require("entropy")
require("factoextra")

entropic_ranks <- function(data_under_analysis,population_vector,data_origin=NULL,granularity=1,supervised=FALSE,process_log=FALSE,export_plots=FALSE,create_output_files=FALSE,is_logged=TRUE,logbase=2,huge_feature_list=FALSE)
{
  if (is.null(data_origin))
    data_origin <- rep(1,length(population_vector))
  
  message("Calculating Rank Products. May take a lot of time, depending on data size.")
  comparison <- RPadvance(data_under_analysis,cl=population_vector,origin=data_origin,logged=is_logged,na.rm=FALSE,gene.names=rownames(data_under_analysis),plot=process_log,huge=TRUE)
  if (huge_feature_list){
    message("Investigating only the first 5000 features.")
    rank_product_lists <- topGene(comparison,num.gene=5000,logged=is_logged,logbase=logbase,gene.names=rownames(data_under_analysis))
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
  differences <- ordered_vector[seq(2,length(ordered_vector))]-ordered_vector[seq(1,length(ordered_vector)-1)]
  entropy_plotter <- vector(length=floor((length(differences)-window_size)/step_up))
  for (i in 0:(length(entropy_plotter)-1))
    entropy_plotter[i+1] <- entropy(discretize(differences[(i*step_up+1):(i*step_up+window_size)],numBins=bins),method="Laplace")
  entropy_clusters <- eclust(entropy_plotter, "kmeans", k=2, nstart=200, graph=FALSE)
  
  if (verbose)
  {
    cat("Calculating entropies of ",length(entropy_plotter)," ovelapping windows","\n","Suggested cutoff at feature no ",seq(length(entropy_clusters$cluster))[entropy_clusters$cluster!=as.integer(entropy_clusters$centers==max(entropy_clusters$centers))[1]+1][1]-1,", at a mean entropy of ",mean(entropy_plotter[1:seq(length(entropy_clusters$cluster))[entropy_clusters$cluster!=as.integer(entropy_clusters$centers==max(entropy_clusters$centers))[1]+1][1]-1]),"\n","Last 1/3 minimum entropy: ",min(entropy_plotter[floor(length(entropy_plotter)*2/3):length(entropy_plotter)]),"\n",sep="")
    barplot(entropy_plotter,border=c("gold1","dodgerblue3")[as.vector(as.integer(entropy_clusters$cluster==(as.integer(entropy_clusters$centers==max(entropy_clusters$centers))[1]+1))+1)],col=c("gold1","dodgerblue3")[as.vector(as.integer(entropy_clusters$cluster==(as.integer(entropy_clusters$centers==max(entropy_clusters$centers))[1]+1))+1)],xlab=c("Granularity = ",step_up),ylab=c("Window size = ",window_size),main=c(bins," bins"),names.arg=seq(length(entropy_plotter))*step_up)
    barplot(entropy_clusters$silinfo$widths$sil_width,border=c("gold1","dodgerblue3")[as.integer(entropy_clusters$silinfo$widths$cluster==(as.integer(entropy_clusters$centers==max(entropy_clusters$centers))[1]+1))+1],col=c("gold1","dodgerblue3")[as.integer(entropy_clusters$silinfo$widths$cluster==(as.integer(entropy_clusters$centers==max(entropy_clusters$centers))[1]+1))+1],ylim=c(-0.3,1),ylab="Silhouette width",xlab="Clustered elements",main="K-means clustering quality")
  }
  if (export_plots)
  {
    if (is.null(path))
      path <- paste(getwd(),"Entropic Ranks plots",sep="/")
    if(!file.exists(path))
      dir.create(path)
    png(file=file.path(path,paste("Entropy - ",bins," bins - ",window_size," window size",".png", sep = "")),width=640,height=640)
    barplot(entropy_plotter,border=c("gold1","dodgerblue3")[as.vector(as.integer(entropy_clusters$cluster==(as.integer(entropy_clusters$centers==max(entropy_clusters$centers))[1]+1))+1)],col=c("gold1","dodgerblue3")[as.vector(as.integer(entropy_clusters$cluster==(as.integer(entropy_clusters$centers==max(entropy_clusters$centers))[1]+1))+1)],xlab=c("Granularity = ",step_up),ylab=c("Window size = ",window_size),main=c(bins," bins"),names.arg=seq(length(entropy_plotter))*step_up)
    dev.off()
    png(file=file.path(path,paste("Clustering quality - ",bins," bins - ",window_size," window size",".png", sep = "")),width=640,height=640)
    barplot(entropy_clusters$silinfo$widths$sil_width,border=c("gold1","dodgerblue3")[as.integer(entropy_clusters$silinfo$widths$cluster==(as.integer(entropy_clusters$centers==max(entropy_clusters$centers))[1]+1))+1],col=c("gold1","dodgerblue3")[as.integer(entropy_clusters$silinfo$widths$cluster==(as.integer(entropy_clusters$centers==max(entropy_clusters$centers))[1]+1))+1],ylim=c(-0.3,1),ylab="Silhouette width",xlab="Clustered elements",main="K-means clustering quality")
    dev.off()
  }
  
  return(seq(length(entropy_clusters$cluster))[entropy_clusters$cluster!=as.integer(entropy_clusters$centers==max(entropy_clusters$centers))[1]+1][1]-1)
}
