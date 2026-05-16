#' Convert summary_statistic file output from Inspector into a dataframe
#'
#'@author Gen Morinaga
#' The python program Inspector (found at: https://github.com/Maggi-Chen/Inspector)
#' outputs summary statistics for genome assemblies. Unfortunately, the output is not
#' a format that is convenient for reading into R. This function attempts to take this
#' output and convert it into a dataframe. It should be noted that the structural errors found
#' in this file is a tally (i.e., the number of times an error type was found), and thus
#' denoted in the dataframe with the prefix num_* (as in: number of error_type).
#'
#' As of this writing (Feb. 3, 2023), Inspector does not name its output files with identifying information (by design).
#' For use with this function, the file must be copied and renamed with identifying information separated by underscores (_).
#' For example: specimen_assembly_summary_statistics
#' Note that if you need to use inspector-correct, then the files cannot be renamed (hence copying and renaming).
#' @param path_to_insp_summ The absolute path to the *_summary_statistics output from
#' Inspector. The \* should be some sort of ID string
#' @details This function is best used in lapply to import and load multiple summary_statistic
#' files. This can be done by creating a vector of absolute paths to the summary_statistic files
#' or by aggregating the files into a single directory, setting that directory as the working directory,
#' and passing the list.files() argument into lapply. See examples for details.
#' @examples
#' #write example output of a summary_statistics output from Inspector to the current working directory
#' cat(`Nmex03_ipa-np_insp_summary_statistics`, "Nmex03_ipa-np_insp_summary_statistics")
#' cat(`Nmex03_ipa-np_pd_insp_summary_statistics`, "Nmex03_ipa-np_pd_insp_summary_statistics")
#'
#' #for a single file
#' insp_sum2df("Nmex03_ipa-np_insp_summary_statistics")
#'
#' #for multiple files
#' do.call('rbind', lapply(list.files(pattern = 'summary_statistics' ), insp_sum2df))
#'
#' @export
insp_sum2df<-function(path_to_insp_summ) {
  txt<-readLines(path_to_insp_summ)
  vec<-gsub('>','',gsub('_/%','', gsub(':', '', gsub(' ', '_', txt))))
  vec<-vec[vec != '']
  vec<-vec[c(-1, -12)]
  df<-data.frame(t(matrix(ncol = 2, unlist((strsplit(vec, '\t'))), byrow = T)))
  names(df)<-df[1,]
  df<-df[-1,]
  df<-data.frame(t(apply(df, 2, as.numeric)))
  df$assembly<-gsub('_summary_statistics', '', fs::path_file(path_to_insp_summ))
  colnames(df)[17:21]<-paste0('num_', colnames(df)[17:21])
  return(df)
}
