#" Generate a filter list of scaffolds
#"
#" @param jbat.list a list object generated from the read.jbat.review() function (REQUIRED)
#" @param busco.dat a dataframe generated from the read.busco() function (REQUIRED)
#" @param out a string that will be used as output file prefixes (default = NULL; optional)
#" @details
#" This creates a filter list and a summary table of scaffolds which will output
#" tab-delimited txt files. The main goal of this function is to create two lists which will help to subset
#" a scaffold assembly to contain only main scaffolds and unplaceable contigs/scaffolds. The generated lists
#" and their contents are:
#"  (out_)scaffold_summary.txt
#"    sequence: name of the sequence, scaffold_1:n where n is the total number of scaffolds
#"    jbat_status: either "scaffold" or "debris". when editing HiC contact maps in Juicebox
#"                 some sequences may be designated as "debris", which are typically
#"                 sequenced duplicates. fragments are considered scaffolds.
#"    length: length of scaffold derived from taking the sum of all contigs that comprise
#"             that scaffold
#"    complete: the number of complete single (S) copy BUSCOs found on a sequence
#"    dup: the number of complete duplicated (D) BUSCOs found on a sequence
#"    frag: the number of fragmented BUSCOs (F) found on a sequence
#"    keep: TRUE/FALSE was the scaffold placed in the keep list? this determination
#"          is made by a set of criteria:
#"          for TRUE:
#"             1) the sequence contain any complete or fragmented BUSCOs (regardless of JBAT status)
#"             2) the sequence is considered a scaffold and does not contain any BUSCOs
#"             3) the sequence contains duplicate BUSCOs that cannot be found on any
#"                other sequence
#"          for FALSE:
#"             1) the sequence contains only duplicate BUSCOs and these BUSCOs are
#"                also located on sequences in TRUE
#"             2) the jbat_status is considered "debris" and does not contain any BUSCOs
#"
#" (out_)list_of_scaff_to_keep.txt
#"    a single column list without headers containing sequence names of scaffolds
#"    that should be kept. these are sequences where keep is TRUE in (out_)scaffold.summary.txt
#"
#" (out_)debris_and_dups_scaffolds.txt
#"    a single column list without headers containing sequence names of scaffolds
#"    that contain only BUSCO duplicates and debris. these are sequences where
#"    keep is FALSE in (out_)scaffold.summary.txt
#"
#" The user can feed either list into their fasta handling program such as
#" seqtk subseq or samtools faidx to create a subsetted fasta containing only
#" the sequences included in either list.

#" @export
busco.scaffold.filter<-function(jbat.list = jbat.list, busco.dat = busco.dat, out = NULL) {
  jbat.list[[1]]$line<-as.numeric(jbat.list[[1]]$line)
  scaffold.dat<-merge(jbat.list[[2]], jbat.list[[1]])
  # scaffold.dat$jbat_status<-ifelse(grepl("debris", scaffold.dat$name), "debris", "scaffold")
  scaffold.dat$scaf_num<-as.numeric(gsub("scaffold_", "", scaffold.dat$scaffold))
  last.scaff.num<-scaffold.dat[scaffold.dat$line == max(jbat.list[[1]]$line), "scaf_num"]
  scaffold.dat$jbat_status<-ifelse(scaffold.dat$scaf_num > last.scaff.num, "debris", "scaffold")


  scaffold.dat$length<-as.numeric(scaffold.dat$length)
  busco.dat$sequence<-gsub("(scaffold_[0-9]+)(\\:.+)", "\\1", busco.dat$sequence)
  busco.summarised<-by(busco.dat, list(busco.dat$sequence), function(x)
    c(
      sequence = unique(x$sequence),
      complete = sum(x$status == "Complete"),
      dup = sum(x$status == "Duplicated"),
      frag = sum(x$status == "Fragmented")
    )
  )
  busco.summarised<-data.frame(do.call(rbind, busco.summarised))
  rownames(busco.summarised)<-NULL
  busco.summarised[, c(2:4)]<-apply(busco.summarised[, c(2:4)], 2, as.numeric)
  busco.summarised<-busco.summarised[which(busco.summarised$sequence != ""),]

  scaffold.summarised<-by(scaffold.dat, list(scaffold.dat$scaffold), function(x)
    c(
      sequence = unique(x$scaffold),
      jbat_status = unique(x$jbat_status),
      length = sum(x$length)
    )
  )
  scaffold.summarised<-data.frame(do.call(rbind, scaffold.summarised))
  rownames(scaffold.summarised)<-NULL
  scaffold.summarised$length<-as.numeric(scaffold.summarised$length)
  combined.summary<-merge(scaffold.summarised, busco.summarised, all.x = T)

  scaff.keep<-rbind(
    combined.summary[which(combined.summary$complete > 0 | combined.summary$frag > 0),],
    combined.summary[which(combined.summary$jbat_status == "debris" & combined.summary$complete > 0),],
    combined.summary[which(combined.summary$jbat_status == "debris" & combined.summary$frag > 0),],
    combined.summary[which(combined.summary$jbat_status == "scaffold" &
                             is.na(combined.summary$complete) == T &
                             is.na(combined.summary$dup) == T &
                             is.na(combined.summary$frag) == T),]
  )
  debris.dups<-rbind(
    combined.summary[which(combined.summary$dup > 0 & combined.summary$complete == 0 & combined.summary$frag == 0), ],
    combined.summary[which(combined.summary$jbat_status == "debris" &
                             is.na(combined.summary$complete == T) &
                             is.na(combined.summary$dup == T) &
                             is.na(combined.summary$frag == T)
    ), ]
  )
  dup.busco<-busco.dat[which(busco.dat$status == "Duplicated"), c("busco_id", "status", "sequence")]
  tmp<-by(dup.busco, list(dup.busco$busco_id), function(x)
    c(
      seq = x$sequence
    ), simplify = F
  )
  tmp<-lapply(tmp, function(x) data.frame(t(x)))
  tmp.names<-unique(unlist(lapply(tmp, names)))
  dup.scaff<-do.call(rbind,
                     lapply(tmp, function(x)
                       data.frame(c(x,
                                    sapply(setdiff(tmp.names, names(x)), function(y)
                                      NA
                                    )
                       )
                       )
                     )
  )
  dup.scaff$remove<-sapply(1:nrow(dup.scaff), function(x) sum(dup.scaff[x, ] %in% scaff.keep$sequence), simplify = T)
  dup.keep<-dup.scaff[dup.scaff$remove<1, ]
  dup.keep<-dup.keep[, -ncol(dup.keep)]
  dup.keep<-dup.keep[, !apply(is.na(dup.keep), 2, all)]
  dup.keep<-merge(merge(dup.keep, data.frame(seq1 = scaffold.summarised[which(scaffold.summarised$sequence %in% dup.keep$seq1),"sequence"],
                                             seq1_length = scaffold.summarised[which(scaffold.summarised$sequence %in% dup.keep$seq1),"length"])),
                  data.frame(seq2 = scaffold.summarised[which(scaffold.summarised$sequence %in% dup.keep$seq2),"sequence"],
                             seq2_length = scaffold.summarised[which(scaffold.summarised$sequence %in% dup.keep$seq2),"length"]))
  if(nrow(dup.keep) > 0) {
    dup.keep$keep<-NA} else {
      dup.keep
    }
  for(i in 1:nrow(dup.keep)){
    dup.keep$keep[i]<-dup.keep[i, gsub("_length", "", names(dup.keep)[grepl("_length", names(dup.keep))][apply(dup.keep[grepl("_length", names(dup.keep))], 1, function(x) which.max(x), simplify = T)])[i]]
  }
  scaff.keep<-unique(rbind(scaff.keep, combined.summary[which(combined.summary$sequence %in% dup.keep$keep),]))
  debris.dups<-unique(debris.dups[!(debris.dups$sequence %in% dup.keep$keep),])
  combined.summary$keep<-ifelse(combined.summary$sequence %in% scaff.keep$sequence, T, F)
  if(is.null(out)) {
    write.table(combined.summary, file = "scaffold_summary.txt", row.names = F, quote = F, sep = "\t")
    write.table(scaff.keep$sequence, file = "list_of_scaff_to_keep.txt", row.names = F, col.names = F, quote = F, sep = "\t")
    write.table(debris.dups$sequence, file = "debris_and_dup_scaffolds.txt", row.names = F, col.names = F, quote = F, sep = "\t")
  } else {
    write.table(combined.summary, file = paste(out, "scaffold_summary.txt", sep = "_"), row.names = F, quote = F, sep = "\t")
    write.table(scaff.keep$sequence, file = paste(out, "list_of_scaff_to_keep.txt", sep = "_"), row.names = F, col.names = F, quote = F, sep = "\t")
    write.table(debris.dups$sequence, file = paste(out, "debris_and_dup_scaffolds.txt", sep ="_"), row.names = F, col.names = F, quote = F, sep = "\t")
  }
}

