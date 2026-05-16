#' BUSCO contig filter
#'
#' Use information of BUSCO gene locations in the genome assembly to create a filter
#' list of contigs to keep
#' @param contig.info A two column dataframe containing contig names (first column) and lengths (second column) (REQUIRED)
#' @param compleasm.out A dataframe output from compleasm, normally called full_table.tsv (https://github.com/huangnengCSU/compleasm)
#' @param busco.out A dataframe output from BUSCO, normally called *_full_table.tsv (https://busco.ezlab.org/)
#' @param out.prefix A prefix for the output file. If no prefix is provided, outputs the list of contigs to keep to console (default = NULL; optional)
#'
#' @details
#' This function creates a list that can be used by an external program (like seqtk subseq) to subset a
#' fasta file, keeping the listed sequences. Rules for which contigs are kept are as follows:
#' 1) All contigs that have a BUSCO gene considered to be complete and single (S) copy
#' 2) All contigs that have a BUSCO gene considered to be fragmented (F) or interspaced (I) (only output from compleasm has this)
#' 3) All contigs that do not have ANY BUSCO genes on them
#' 4) The longest contig with duplicated (D) BUSCO genes
#'
#' The user should supply EITHER the BUSCO output *_full_table.tsv OR a compleasm full_table.tsv. Note that compleasm outputs both formats
#' and calls the BUSCO format full_table_busco_format.tsv

#' @export
busco.contig.filter <- function(
    contig.info = NULL,
    compleasm.out = NULL,
    busco.out = NULL,
    out.prefix = NULL) {
  if (is.null(compleasm.out) & is.null(busco.out)) {
    stop('Please supply EITHER a TSV file output from compleasm or BUSCO')
    }
  if (!(is.null(compleasm.out))) {
    busco.dat <- compleasm.out
  }
  if (!(is.null(busco.out))) {
    names(busco.out) <- c(
      'Gene',
      'Status',
      'Sequence',
      'Gene.Start',
      'Gene.End',
      'Strand',
      'Score',
      'Length',
      'OrthoDB_url',
      'Description'
    )
    busco.dat <- busco.out
    busco.dat$Status<-gsub('Complete', 'Single', busco.dat$Status)
  }
  busco.dat <- busco.dat[which(busco.dat$Status != 'Missing'), ]
  names(contig.info) <- c('sequence.id', 'sequence.length')
  comp.seqs <- unique(busco.dat$Sequence)
  keep <- contig.info$sequence.id[!(contig.info$sequence.id %in% comp.seqs)] #keep all sequences that don't have BUSCOs on them

  contig.status.counts <- as.data.frame.matrix(table(busco.dat[c('Sequence', 'Status')]))
  contig.status.counts$sequence <- row.names(contig.status.counts)
  row.names(contig.status.counts) <- NULL

  if(is.null(compleasm.out)) {
    keep <- c(keep, contig.status.counts[which(contig.status.counts$Single > 0 |
                                                 contig.status.counts$Fragmented > 0), 'sequence']) #keep all sequences that have single copy or fragmented BUSCOs
  }
  if(is.null(busco.out)) {
    if(!('Interspaced' %in% names(contig.status.counts))) {
      contig.status.counts$Interspaced<-0
    }
    keep <- c(keep, contig.status.counts[which(contig.status.counts$Single > 0 |
                                                 contig.status.counts$Fragmented > 0) |
                                           contig.status.counts$Interspaced > 0, 'sequence']) #keep all sequences that have single copy, fragmented or interspaced BUSCOs
  }
  dups <- busco.dat[which(busco.dat$Status == 'Duplicated'), c('Gene', 'Status', 'Sequence')]
  names(contig.info)<-c("Sequence", "Length")

  tmp <- by(dups, list(dups$Gene), function(x)
    c(seq = x$Sequence), simplify = F)
  tmp <- lapply(tmp, function(x)
    data.frame(t(x)))
  tmp.names <- unique(unlist(lapply(tmp, names)))
  dup.seqs <- do.call(rbind, c(lapply(tmp, function(x)
    data.frame(c(
      x, sapply(setdiff(tmp.names, names(x)), function(y)
        NA)
    )))))
  rm(tmp, tmp.names)
  dup.seqs$n.seqs.keep <- sapply(1:nrow(dup.seqs), function(x)
    sum(dup.seqs[x, ] %in% keep), simplify = T)
  tmp <- dup.seqs[which(dup.seqs$n.seqs.keep == 0), -ncol(dup.seqs)]
  if(nrow(tmp) > 0) {
    contig.lengths <- setNames(contig.info$sequence.length, contig.info$sequence.id)
    length.mat <- as.matrix((sapply(1:ncol(tmp), function(x)
      contig.lengths[tmp[, x]], simplify = T)))
    row.names(length.mat) <- NULL
  ind <- apply(length.mat, 1, function(x)
    which(max(x, na.rm = T) == x, arr.ind = T))
  dup.keep <- unique(unlist(mapply(
    function(x, y)
      tmp[x, y],
    x = 1:nrow(tmp),
    y = ind,
    SIMPLIFY = T
  )))} else {
    dup.keep<-unique(dups$Sequence)
  }

  keep <- c(keep, dup.keep)
  keep <- unique(keep)
  if (is.null(out.prefix)) {
    return(keep)
  } else {
    write.table(
      keep,
      file = paste(out.prefix, 'busco.contig.filter.list.txt', sep = '_'),
      row.names = F,
      col.names = F,
      quote = F,
      sep = '\t'
    )
  }
}
