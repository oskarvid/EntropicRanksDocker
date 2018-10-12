# Entropic Ranks (docker) readme

entropic_ranks

Description: Performs an Entropic Ranks analysis on a data set, returning a list containing downregulated and upregulated features. May be used supervised, returning the full feature list and printing the suggested cutoff points for later manual trimming, or unsupervised, returning only the information-rich feature list. In the unsupervised mode, the lists of information-rich features may be exported as tab-delimited .txt files automatically.

Usage: entropic_ranks (data_under_analysis,population_vector,data_origin=NULL,granularity=1,supervised=FALSE,process_log=FALSE,export_plots=FALSE,create_output_files=FALSE,is_logged=TRUE,logbase=2,huge_feature_list=FALSE)

Arguments:
data_under_analysis - Table with rows representing features, columns representing samples and cells containing the values to be compared. Rownames and column names must be unique.

population_vector - Binary integer vector (0 or 1), of length equal to the number of columns of data_under_analysis. Denotes the two sample subpopulations to be compared.

data_origin - A vector containing the origin labels of the samples. Must be of length equal to the number of columns of data_under_analysis. If NULL, it defaults to assuming that data are of the same origin. To be used only if the data are from different experiments or publications.

granularity - The sliding window step, corresponding to the granularity of the partitioning process (feature-by-feature, or partitioning by 5-feature steps).

supervised - If TRUE, the full list of differentially behaving features is returned and the tables of suggested cutoff points are printed. If FALSE, only the list of information-rich features is returned.

process_log - If TRUE, statistics of the entropic_ranks execution will be printed and plots of the entropy distributions and clustering qualities will be generated.

export_plots - If TRUE, png plots of the entropy distributions and clustering qualities will be exported as files in a folder system created in the current working directory.

create_output_files - If TRUE, the feature lists of information-rich features will be automatically exported in the working directory as tab-delimited .txt files. Ignored if supervised is set to TRUE.

is_logged - Set to TRUE if the values were log-transformed and you want to export the Fold Change instead of the Log Fold Change in .txt files. Ignored if supervised is set to TRUE.

logbase - The base of the log transformation. Ignored if supervised is set to TRUE or if create_output_files is set to FALSE.

huge_feature_list - Only set to TRUE if the entropic_ranks fails to run due to huge feature lists returned by RankProd (e.g. more than 20000-30000 features) and you are sure that less than 500 are differentially expressed and information-rich. If TRUE, entropic_analysis will only investigate the first 5000 features and isolare information-rich features from among them.
