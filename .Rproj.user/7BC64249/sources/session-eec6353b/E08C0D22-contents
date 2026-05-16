#' Look for missing BUSCO in the haplotig assembly
#'
#'
#' @param asm.busco1 dataframe containing the full_table.tsv or full_table.txt output from BUSCO or compelasm for the main/primary assembly (required)
#' @param asm.busco2 dataframe containing the full_table.tsv or full_table.txt output from BUSCO or compelasm for the secondary/haplotig assembly (required)
#' @param summary summarise the number of BUSCO genes that are found in the subset of contigs from the secondary/haplotig assembly (optional; logical)
#' @param out supplying a output prefix will write a headerless, single-column file of contig names from the secondar/haplotig assembly. File name will be in the format [out].txt (optional, string)
#'
#' @details Determine if there are BUSCO genes missing in the main assembly that are
#' are found in the alternate assembly. BUSCO should be run on two sister
#' assemblies—that is they should be output from a single set of reads from the
#' same genome assembler as main and alternate assemblies OR output from a
#' haplotig/duplicate contig purging program such as purge_haplotigs or
#' purge_dups. Function requires two BUSCO (or compleasm) full_table.tsv files.
#' @export
asmBUSCOcheck<-function(asm1.busco = NULL, asm2.busco = NULL, summary = FALSE, out = NULL) {

  try(
    if (is.null(asm1.busco) == TRUE | is.null(asm2.busco) == TRUE)
      stop("please specify asm1.busco and/or asm2.busco in dataframe format")
  )

  check1<-sum(grepl("Duplicated|Single|Fragmented|Missing|Interspaced", asm1.busco[1, 2]))
  check2<-sum(grepl("Duplicated|Single|Fragmented|Missing|Interspaced", asm2.busco[1, 2]))

  try(
    if(check1 != 1)
      stop("asm1.busco is not the full table output from BUSCO or compleasm")
    )

  try(
    if(check2 != 1)
      stop("asm2.busco is not the full table output from BUSCO or compleasm")
    )

  if(ncol(asm1.busco) == 10) {
    names(asm1.busco)<-c("Gene", "Status", "Sequence", "Gene.Start",
                         "Gene.End", "Strand", "Score", "Length")
  }

  sub.col<-c("Gene", "Status", "Sequence")

  asm1.sub<-asm1.busco[, sub.col]
  asm2.sub<-asm2.busco[, sub.col]

  asm1.missing<-asm1.sub[which(asm1.sub[, "Status"] == "Missing"), 1]
  asm2.not.missing<-asm2.sub[which(asm2.sub[, "Status"] != "Missing"), ]

  asm2.not.in.asm1<-asm2.not.missing[asm2.not.missing[, "Gene"] %in% asm1.missing, ]

  # hap<-unique(asm2.not.in.asm1[, "Sequence"])
  # busco<-unique(asm2.not.in.asm1[, "Gene"])

  ####something to consider to recover fragmented BUSCOs from asm2
  # asm1.frag.in.asm2<-asm2.sub[which(asm2.sub[, 1] %in% asm1.sub[which(asm1.sub[,2] == "Fragmented"),1]),]
  #
  # asm1.frag.asm2.single<-asm1.frag.in.asm2[which(asm1.frag.in.asm2[, 2] == "Single"), 3]

  # single<-asm2.not.in.asm1[which(asm2.not.in.asm1[, "Status"] == "Single"), "Sequence"]
  # keep<-unique(c(single, asm1.frag.asm2.single))
  ####
  keep<-unique(asm2.not.in.asm1[, "Sequence"])
  if(length(keep) == 0)
    warning("all missing BUSCOs in asm1 are also missing in asm2", call. = FALSE)

  if(is.null(out) == FALSE) {
    write.table(keep, file = paste(out, "txt", sep = "."),
                quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
  }
  if(summary == TRUE & length(keep) != 0) {
    asm2.keep.summary<-aggregate(Gene ~ Sequence + Status,
                                 asm2.not.in.asm1[which(asm2.not.in.asm1[, 3] %in% keep),],
                                 FUN = length)
    print(asm2.keep.summary)
    } else {
      print(keep)
  }
  }
