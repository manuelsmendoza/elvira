#' Reverse translate
#'
#' @description
#' `rev_translate()` accepts a protein sequence as input and uses the standard genetic code
#'  to generate a DNA sequence representing the most likely non-degenerate coding sequence.
#'
#' @param aaseq Protein sequence
#' @return A DNA sequence.
#'
#' @examples
#' sq <- "MACDEFGHIKL*"
#' rev_translate(sq)
rev_translate <- function(aaseq) {
  aasq <- as.character(aaseq)
  ntsq <- ""

  # Convert the aminoacids into the correspondent codon
  # If there are more than one option (genetic code degeneration)
  # the codon is choose randomly between the different options
  genetic_code <- Biostrings::GENETIC_CODE
  for (I in 1:nchar(aasq)) {
    cods <- names(genetic_code)[which(genetic_code == str_sub(string = aasq, start = I, end = I))]

    if (length(cods) > 1) {
      ntsq <- paste0(ntsq, sample(x = cods, size = 1, replace = FALSE))
    } else {
      ntsq <- paste0(ntsq, cods)
    }
  }

  # Add the sequence name (if available)
  if (is.null(names(aasq))) {
    ntsq <- DNAString(ntsq)
  } else {
    ntsq <- DNAStringSet(ntsq)
    names(ntsq) <- names(aasq)
  }

  return(ntsq)
}


#' Transcriptome completeness checking
#'
#' @description
#' `check_completeness()`verify if the transcripts were assembled completely and the
#' coverage of those sequences
#'
#' @param prot_seqs Transcriptome sequences translated to proteins
#' @param prot_blast Proteins BLAST output
#' @return A dataset containing the completeness calssification and coverage
check_completeness <- function(prot_seqs, prot_blast) {
  # Load the protein sequences and the blast output
  ptseqs <- suppressMessages(
    readAAStringSet(filepath = prot_seqs, use.names = TRUE)
  )

  # Load the blast result
  colnames <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "slen", "qlen", "qcovs")
  pthomo <- suppressMessages(
    read_delim(file = prot_blast, delim = "\t", col_names = colnames) %>%
      mutate("txcov" = 100 * abs(qend - qstart) / qlen,
             "ptcov" = 100 * abs(send - sstart) / slen) %>%
      select(qseqid, qlen, txcov, ptcov) %>%
      rename("txname" = qseqid, "txlen" = qlen)
  )

  # Check the completeness
  compl_tbl <- tibble("txname" = NA, "txcompl" = NA)

  for (I in 1:length(ptseqs)) {
    # Extract the sequence name
    ptname <- names(ptseqs[I])

    # Extract the protein sequence
    ptseq  <- as.character(ptseqs[I])

    # Extract the first and last amino acid
    start_pos <- str_sub(string = ptseq, start =  1, end =  1)
    end_pos   <- str_sub(string = ptseq, start = -1, end = -1)

    # Check protein completeness
    if (start_pos == "M" & end_pos == "*") {
      sqcomp <- "Complete"
    } else if (start_pos == "M" & end_pos != "*") {
      sqcomp <- "N-Terminal"
    } else if (start_pos != "M" & end_pos == "*") {
      sqcomp <- "C-Terminal"
    } else {
      sqcomp <- "Internal"
    }

    # Add this information to the data.frame
    compl_tbl[I, "txname"]  <- ptname
    compl_tbl[I, "txcompl"] <- sqcomp
  }

  # Merge both dataset
  completeness_tbl <- inner_join(pthomo, compl_tbl)

  return(completeness_tbl)
}
