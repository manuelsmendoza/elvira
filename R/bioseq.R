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

  genetic_code <- Biostrings::GENETIC_CODE
  for (I in 1:width(aasq)) {
    cods <- genetic_code[which(names(genetic_code) == str_sub(string = aasq, start = I, end = I))]

    if (length(cods) > 1) {
      ntsq <- paste0(ntsq, sample(x = cods, size = 1, replace = FALSE))
    } else {
      ntsq <- paste0(ntsq, cods)
    }
  }

  return(ntsq)
}
