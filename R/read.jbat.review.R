#' Reads the *.review.assembly output from Juicebox HiC visualization software
#'
#' @details
#' Reads the *.review.assembly output from Juicebox (https://github.com/aidenlab/Juicebox) and
#' returns a list of two dataframes, each comprising the information contained in *.review.assembly
#'
#' *.review.assembly contains two sets of information: the first half is a list of contigs, its line number,
#' and length. The name may also contain a label that allows Juciebox to designate how a contig may have been
#' curated by the user (e.g., misjoined contigs might be broken, or duplicate stretches sequence might be placed in 'debris').
#' Juicebox takes this information and organizes the contigs into scaffolds (or superscaffolds) for visualization, which is shown in the second half.
#' The second half of these files contain scaffold information. Each scaffold is represented by a row, and each row contains a space-delimited list of
#' contigs represented by the line number shown in the first half. The orientation of the contig becomes flipped (inverse) if it is shown with a minus preceding
#' it (e.g., -42 will refer to the contig in line 42 in the first half, with inverted orientation).
#'
#' We create a a list of dataframes: one representing the contig data, and the other the scaffolding data, each with
#' three columns. The first element contains contig data, (contig name, line number, and length). Name here also contains the status
#' of the contig (e.g., debris, fragment, or simply unlabeled). The second element contains scaffold
#' data (scaffold name, line from element 1, and orientation of the contig).
#' This list can be fed into scaffold.filter() along with BUSCO information to create a filter list of scaffolds.

#' @export
read.jbat.review<-function(path) {
   review<-readLines(path) #read in review.assembly from JBAT
   review.top<-noquote(
    matrix(
      unlist(
        strsplit(
          review[grepl(">", review)],
          split = ' ')
        ),
      byrow = T,
      ncol = 3
      )
    )

  review.top<-data.frame(name = review.top[, 1], line = review.top[, 2], length = review.top[, 3])

  review.bottom<-data.frame(line = noquote(review[!grepl(">", review)]))
  review.bottom$scaffold<-paste('scaffold', 1:nrow(review.bottom), sep = '_')
  line<-strsplit(review.bottom$line, split = ' ')
  review.bottom<-data.frame(scaffold = rep(review.bottom$scaffold, sapply(line, length)),
                            line  = as.numeric(unlist(line)))
  review.bottom$orientation<-ifelse(review.bottom$line > 0, "+", "-")
  review.bottom$line<-abs(review.bottom$line)
  scaffold.dat<-merge(review.bottom, review.top)
  scaffold.dat$jbat_status<-ifelse(grepl("debris", scaffold.dat$name), "debris", "scaffold")

 x<-list(review.top, review.bottom)
 return(x)
}
