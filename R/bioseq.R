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
