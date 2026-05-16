#' Read in short summary outputs from BUSCO and compleasm
#'
#' @param file path to the output you want to upload into R
#' @param format Can be either "short" or "full". "short" format are the short summaries outputs from either BUSCO or compleasm. "full" is the full output from BUSCO.
#' @param program Can be either "BUSCO" or "compleasm"
#' @param assembly.name OPTIONAL name of the assembly
#' @details
#' A series of simple wrapper functions to read BUSCO and compleasm summary outputs. It also reads in the BUSCO full output, and adds sensible column names.
#' Reading full table output from compleasm is trivial, but this function will do it for sake of completeness.
#'

#' @export
read.busco<-function(file, format = c("short", "full"), program = c("BUSCO", "compleasm"), assembly.name = NULL) {
  if(format == "short" & program == "BUSCO") {
    tmp<-readLines(file)
    tmp.dat<-data.frame(
      category = c("Single", "Duplicated", "Fragmented", "Missing", "Total"),
      count = as.numeric(gsub("(\t)([0-9]+)(\t.+)", "\\2", tmp[11:15])),
      percent = c(gsub("%", "",unlist(regmatches(tmp[9], gregexpr("\\d+(\\.\\d+){0,1}%", tmp[9]))))[2:5], 100)
    )
  }
  if(format == "long" & program == "BUSCO") {
    tmp.dat<-read.table(file, sep = "\t", quote = "", header = F, fill = T)
    names(tmp.dat)<-c("busco_id", "status", "sequence", "gene_start",
                  "gene_end", "strand", "score", "length", "orthodb_url", "description")
  }
  if(format == "long" & program == "compleasm") {
    tmp.dat<-read.table(file, sep = "\t", quote = "", header = TRUE, fill = T)
  }
  if(format == "short" & program == "compleasm") {
    tmp<-readLines(file)[2:7]
    tmp.dat<-data.frame(
      category = c("Single", "Duplicated", "Fragmented", "I_Fragmented", "Missing", "Total"),
      count = as.numeric(gsub("N\\:", "",gsub("(.+, )([0-9]+)","\\2", tmp))),
      percent = as.numeric(gsub("%", "", c(unlist(regmatches(tmp, gregexpr("\\d+(\\.\\d+){0,1}%", tmp))), 100)))
    )
  }
  if(!is.null(assembly.name)) {
    tmp.dat$assembly.name<-assembly.name
  }
  return(tmp.dat)
}




