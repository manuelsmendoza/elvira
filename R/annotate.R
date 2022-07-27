#' Check if the homologous proteins identified by local alignment come from a certain phylum
#'
#' @param nr_matches The BLASTP result against the NCBI non-redundant protein database (BLAST output format 6)
#' @param tr_matches The BLASTP result against the UniPrto/TrEMBL database (BLAST output format 6)
#' @param nr_taxonomy Dataset with information to convert NCBI accession number to taxonomy information
#' @param tr_taxonomy Dataset with information to convert UniProt accession number to taxonomy information
#' @param phylum_name Phylum name to check
#' @param evalmax Maximum e-value allowed to qualify a match
#' @return If the matches correspond to the phylum specified
check_taxonomy <- function(nr_matches, nr_taxonomy = NULL, tr_matches, tr_taxonomy = NULL, phylum_name = "Mollusca", evalmax = 1e-06) {
    # Attach the dataset to convert proteins ID to taxonomy information
    if (is.null(nr_taxonomy)) {
        load("data/protnr2taxa.Rda")
    } else {
        nr2taxa <- read_delim(nr_taxonomy)
    }

    if (is.null(tr_taxonomy)) {
        load("data/prottr2taxa.Rda")
    } else {
        tr2taxa <- read_delim(tr_taxonomy)
    }


    # Read the blast output using both databases
    header <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue",
        "bitscore")
    nr_tbl <- read_delim(file = nr_matches, delim = "\t", col_names = header) %>%
        select(qseqid, sseqid, evalue) %>%
        filter(evalue <= evalmax) %>%
        rename(accession = sseqid)
    tr_tbl <- read_delim(file = tr_matches, delim = "\t", col_names = header) %>%
        select(qseqid, sseqid, evalue) %>%
        filter(evalue <= evalmax) %>%
        rename(accession = sseqid)

    # Check the phyla
    tx_tbl <- nr_tbl %>%
        select(qseqid, accession) %>%
        rename(nr_acc = accession) %>%
        full_join(tr_tbl %>%
            select(qseqid, accession) %>%
            rename(tr_acc = accession)) %>%
        full_join(nr2taxa %>%
            select(accession, phylum) %>%
            rename(nr_acc = accession, nr_phylum = phylum)) %>%
        full_join(tr2taxa %>%
            select(accession, phylum) %>%
            rename(tr_acc = accession, tr_phylum = phylum)) %>%
        select(qseqid, nr_phylum, tr_phylum) %>%
        mutate(is_phylum = case_when(is.na(nr_phylum) & is.na(tr_phylum) ~ FALSE,
                                     is.na(nr_phylum) & tr_phylum == phylum_name ~ TRUE,
                                     nr_phylum == phylum_name & is.na(tr_phylum) ~ TRUE,
                                     nr_phylum == phylum_name & tr_phylum == phylum_name ~ TRUE,
                                     TRUE ~ FALSE)) %>%
        filter(!is.na(qseqid))

    # Return if the match was from the correct phylum
    in_phylum <- tx_tbl %>%
        pull(is_phylum)
    names(in_phylum) <- tx_tbl %>%
        pull(qseqid)

    return(in_phylum)
}









