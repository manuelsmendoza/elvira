#' Transcripts functional annotation
#'
#' @description
#' `go_annotation()` Uses the information about homologous proteins and protein domains
#' detected in the transcripts to predict the Gene Ontologies associated with the sequences
#' assembled
#'
#'
#' @param trans_blast Result of BLAST alignment if should be output format 6.
#' @param trans_hmm Result of protein domain identification
#' @param pfam_goid Conversion table between protein domains and GO IDs
#' @param unip_goid Conversion table between protein accession in UniProt DB and GO IDS
#' @return Conversion table between the transcripts name and GO IDs
go_annotation <- function(trans_blast, trans_hmm, pfam_goid = NULL, unip_goid = NULL) {
  # Check the conversion data
  if (is.null(unip_goid)) {
    stop("UniProt to Gene Ontology conversion is required")
  } else {
    unip_go <- suppressMessages(
      read_delim(file = unip_goid, delim = "\t", col_names = c("ptacc", "ptgoid"))
    )
  }

  # Download the conversion between protein domains and Gene Ontologies
  if (is.null(pfam_goid)) {
    pfam_url <- "http://current.geneontology.org/ontology/external2go/pfam2go"
    pfam_go  <- suppressMessages(
      read_delim(pfam_url, comment = "!", delim = ";", col_names = FALSE) %>%
        rename("ptdom" = X1, "pdgoid" = X2) %>%
        mutate(ptdom = str_split(string = ptdom, pattern = " ", simplify = TRUE)[, 1],
               ptdom = str_replace(string = ptdom, pattern = "Pfam:", replacement = ""),
               pdgoid  = str_replace(string = pdgoid, pattern = " ", replacement = "")) %>%
        group_by(ptdom) %>%
        summarise(pdgoid = paste(pdgoid, collapse = ";"))
    )
  } else {
    pfam_go <- suppressMessages(
      read_delim(file = pfam_goid, delim = "\t", col_names = c("ptdom", "pdgoid"))
    )
  }

  # Load the transcriptome annotation
  txunip_tbl <- suppressMessages(
    read_delim(file = trans_blast, delim = "\t", col_names = FALSE) %>%
      select(X1, X2) %>%
      rename("txname" = X1, "ptacc" = X2)
  )

  txpfam_txt <- system(paste("grep -v \"#\" ", trans_hmm, " |  awk \'{print $1\"\t\"$5}\'"), intern = TRUE)
  txpfam_tbl <- read.table(text = txpfam_txt, col.names = c("txname", "ptdom")) %>%
    as_tibble() %>%
    mutate(ptdom = str_replace(string = ptdom, pattern = "\\..*", replacement = ""))

  # Convert the protein homologous match into GO annotation
  txunip_tbl <- suppressMessages(
    txunip_tbl %>%
    inner_join(unip_go) %>%
    group_by(txname) %>%
    summarise(ptgoid = paste(ptgoid, collapse = ";"))
  )

  # Convert the protein domains into GO annotation
  txpfam_tbl <- suppressMessages(
    txpfam_tbl %>%
    inner_join(pfam_go) %>%
    group_by(txname) %>%
    summarise(pdgoid = paste(pdgoid, collapse = ";"))
  )

  # Convert the protein domains and homologous into GO IDs
  txgoid_tbl <- suppressMessages(
    full_join(txunip_tbl, txpfam_tbl) %>%
    group_by(txname) %>%
    mutate(txgoid = paste(ptgoid, pdgoid, sep = ";", collapse = ";")) %>%
    select(txname, txgoid)
  )

  # Simplify the GO annotations
  for (I in 1:nrow(txgoid_tbl)) {
    txgoid_tbl[I, "txgoid"] <- paste(unique(str_split(string = txgoid_tbl[I, "txgoid"], pattern = ";")[[1]]), collapse = ";")
  }

  # Remove NA
  txgoid_tbl <- txgoid_tbl %>%
    mutate(txgoid = str_replace(string = txgoid, pattern = "NA;|;NA$", replacement = ""))

  return(txgoid_tbl)
}
