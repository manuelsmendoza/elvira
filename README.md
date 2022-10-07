
<h1 align="center">
<p>
Transcriptome landscape of kleptoplastic sea slug <i>Elysia viridis</i>
</p>
</h1>
<p align="center">
<a href="https://scholar.google.com/citations?user=sLWYp9QAAAAJ&hl=en&oi=ao"><b>Manuel
Mendoza</b></a>,
<a href="https://scholar.google.com/citations?user=RWCFV-IAAAAJ&hl=en&oi=sra">Sara
Rocha</a>,
<a href="https://scholar.google.com/citations?user=BmYFiWEAAAAJ&hl=en&oi=ao">Jesus
Troncoso</a>,
<a href="https://scholar.google.com/citations?user=sXJWNsYAAAAJ&hl=en">David
Posada</a> and
<a href="https://scholar.google.com/citations?hl=en&user=X-KLBboAAAAJ">Carlos
A. Canchaya</a>
</p>

Sea slugs from the superorder
[Sacoglossa](https://www.marinespecies.org/aphia.php?p=taxdetails&id=167)
can sequester functional chloroplast through feeding and keep them
photosynthetically active inside their digestive tubules (de Vries,
Christa, and Gould 2014). A small polyphyletic group of sacoglossan
species can maintain the stolen plastids
([kleptoplasts](https://en.wikipedia.org/wiki/Kleptoplasty)) functional
for more than a month (Händeler et al. 2009). Despite an extensive
research record (de Vries et al. 2014), some questions remain unsolved:
How are the plastids recognised from the other components of prey’s
cells? Why do the sea slugs remain alive after weeks of starvation if
the photosynthates are not essential? Even if the ability to sequester
the plastids has multiple independent origins along the evolution, can
we find orthologs related to the time plastids remain active?

We tried to help to answer the three questions sequencing a novel
transcriptome from a cosmopolitan species present along the European
Atlantic shore, *Elysia viridis* (Montagu, 1804) (Jensen 2007) and
comparing it with other species of sea slugs. The sample collection and
the methodology used to find the answers are described below. Although
we report new findings, we are sure that further studies will be
required to find a clear answer to the questions. We here describe the
full pipeline of our analysis.

    @article{
      title={Transcriptome landscape of kleptoplastic sea slug \textit{Elysia viridis}},
      author={Mendoza, Manuel and Rocha, Sara and Troncoso, Jes\'us and Posada, David and Canchaya, Carlos A.},
      journal={bioRxiv},
      year={2022},
      publisher={Cold Spring Harbor Laboratory}
    }

# Working environment set up

In our analysis we used miniconda3 v4.11.0. This analysis was done using
a high-performance computer ([CESGA Finisterrae
II](https://www.cesga.es/infraestructuras/computacion/finisterrae-ii/))
using the queue system [SLURM](https://slurm.schedmd.com), so we added
some SLURM-specific options in the different chunks.

``` bash
# Crete the environment
conda create --yes --name elvira_env

# Add the channels to the list to download the required tools
conda config --add channels conda-forge 
conda config --add channels bioconda
conda config --add channels r

# Install all the tools required for the analysis
conda install --yes --name elvira_env python=3 r-base=4 fastp trinity transrate transdecoder busco blast hmmer fastqc
```

## Reference generation

**STEP 1**: Reads quality assessment. We checked the quality of our
reads using the modules in FastQC.

``` bash
# Reads quality assessment 
fastqc \
  --extract \
  --nogroup \
  --threads $SLURM_NTASKS \
  --outdir qc_dir \
  sample_1.fq.gz sample_2.fq.gz 
```

**STEP 2**: Reads trimming using fastp (Chen et al. 2018). We removed
remaining sequencing adapters (detected by reads-pairs overlapping),
poly-X tails (we did not allow more than 15bp consecutive position with
the nucleotide), and low complexity reads (reads with many repetitions).
After we trimmed the low-quality positions, requiring a mean quality of
30 phred-score. The final length required after trimming was 70bp.

``` bash
# Reads quality control
fastp \
  --thread $SLURM_NTASKS \
  --in1 sample_1.fq.gz \
  --in2 sample_2.fq.gz \
  --out1 sample_1.trim.fq.gz \
  --out2 sample_2.trim.fq.gz \
  --json sample.trim.json\
  --html sample.trim.html \
  --detect_adapter_for_pe \
  --low_complexity_filter \
  --complexity_threshold 50 \
  --cut_right \
  --cut_window_size 10 \
  --cut_mean_quality 30 \
  --qualified_quality_phred 30 \
  --unqualified_percent_limit 25 \
  --average_qual 30 \
  --length_required 70 \
  --correction \
  --overrepresentation_analysis \
  --overrepresentation_sampling 5
```

**STEP 3**: Transcriptome *de novo* assembly. We used the clean reads
obtained in the previous step to reconstruct the transcripts using
Trinity (Haas et al. 2013).

``` bash
# Transcriptome assembly
Trinity \
  --CPU $SLURM_NTASKS \
  --max_memory 64G \
  --seqType fq \
  --min_contig_length 300 \
  --KMER_SIZE 20 \
  --min_glue 10 \
  --min_per_id_same_path 84 \
  --max_diffs_same_path 16 \
  --left sample_1.trim.fq.gz \
  --right sample_2.trim.fq.gz \
  --output sample-Trinity \
  --full_cleanup
```

**STEP 4**: Transcriptome assembly evaluation. We checked the number of
molluscan orthologs assembled using BUSCO (Manni et al. 2021)
(transcriptome completeness); and we checked the number of transcripts
that were assembled correctly using TransRate (Smith-Unna et al. 2016)
(transcriptome correctness).

``` bash
# Check the transcriptome completeness 
busco \
  --force \
  --cpu $SLURM_NTASKS \
  --evalue 1e-6 \
  --in sample_transcripts.fna \
  --out mollusca_ort \
  --mode transcriptome \
  --lineage_dataset mollusca_odb10  

# Check the transcriptome correctness 
transrate \
  --threads $SLURM_NTASKS \
  --assembly sample_transcripts.fna \
  --left sample_1.trim.fq.gz \
  --right sample_1.trim.fq.gz \
  --output sample_cor
```

**STEP 5**: Protein-coding sequences identification. We extracted the
protein coding sequences from the transcript assembled correctly using
TransDecoder. After that, we removed the redundant sequences using
seqkit (Shen et al. 2016).

``` bash
# Create the transcript to gene map
${TRINITY_HOME}/util/support_scripts/get_Trinity_gene_to_trans_map.pl \
  sample.good.fasta > sample.good.gtm

# Extract the protein-coding sequences
TransDecoder.LongOrfs \
  -t sample.good.fasta \
  --gene_trans_map sample.good.gtm \
  --output_dir sample.good-cds

TransDecoder.Predict \
  --single_best_only \
  -t sample.good.fasta \
  --output_dir sample.good-cds

# Remove redundant sequences
seqkit rmdup \
  --threads $SLURM_NTASKS \
  --ignore-case \
  --by-seq \
  --out-file sample.good.rmdup.faa \
  sample.good.faa
```

**STEP 6**: Possible biological contamination removal. After removing
the possible miss-assemblies, we pulled the potential contaminants from
other organisms (bacteria, algae, etc.). We aligned the reads to the
NCBI non-redundant protein database and from that we extracted the
taxonomic information from these matches using the [taxonomizr R
package](https://cran.r-project.org/web/packages/taxonomizr/index.html).
We only kept the sequences that matched with a molluscan protein from
the genera *Elysia*, *Aplysia* and *Plakobranchus*. Previously to define
the filtering condition, we explored the taxonomy of all the matches to
find the most appropriate limits.

``` bash
# Local alignment to identify the product
blastp \
  -num_threads $SLURM_NTASKS \ 
  -task blastp-fast \
  -evalue 1e-3 \
  -max_hsps 1 \
  -max_target_seqs 1 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs" \
  -db nr \
  -query sample.good.rmdup.faa \
  -out sample.good.rmdup.nr.tsv
```

``` r
# Attach the packages
library(taxonomizr)
library(Biostrings)
library(stringr)
library(readr)
library(dplyr)

# Prepare the database (do only the 1st time)
getNamesAndNodes()
getAccession2taxid(types = c("nucl_wgs", "nucl_gb", "prot"))
read.names.sql("names.dmp", "accessionTaxa.sql")
read.nodes.sql("nodes.dmp", "accessionTaxa.sql")
read.accession2taxid(list.files(".", "accession2taxid.gz$"), "accessionTaxa.sql")

# Load the transcripts sequences and the output from STEP 5A
prot_good <- readAAStringSet("sample.good.rmdup.faa")
blast_out <- read_delim(file = "sample.good.rmdup.nr.tsv", col_names = FALSE) %>%
  select("X1", "X2") %>%
  rename("tx_name" = "X1", "nr_accs" = "X2")

# Extract the taxonomic information from the proteins ID
prot_taxa <- tibble(
  "nr_accs" = unique(sort(blast_out %>% pull("nr_accs"))),
  "nr_egid" = accessionToTaxa(nr_accs, "accessionTaxa.sql")
  ) %>%
  bind_cols(getTaxonomy(nr_accs %>% pull(nr_egid), "accessionTaxa.sql"))

# Filter by Genus to keep only sea slugs
prot_filt <- prot_taxa %>% 
  filter(genus %in% c("Aplysia", "Elysia", "Plakobranchus")) %>%
  inner_join(blast_out)
prot_filt <- prot_good[prot_filt %>% pull(tx_name)]

# Rename the final transcripts to handle simple names
new_sqaccs <- paste0("SPNAME", str_pad(1:length(prot_good), width = nchar(length(prot_good)), pad = "0"))
old_sqaccs <- names(prot_filt)
accs_convr <- tibble("qseqid" = old_sqaccs, "qseqid_new" = new_sqaccs)
names(prot_filt) <- new_sqaccs

# Change the sequence accession in the alignment
col_names <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovs")
res_blast <- read_delim(file = "sample.good.rmdup.nr.tsv", col_names = col_names) %>%
  inner_join(accs_convr) %>%
  select(!qseqid) %>%
  rename("qseqid" = "qseqid_new")

# Export the final transcriptome
writeXStringSet(x = prot_filt, filepath = "sample.good.rmdup.slugs.faa", append = FALSE)
write_delim(x = res_blast, file = "sample.good.rmdup.slugs.nr.tsv", delim = "\t", col_names = FALSE)
```

**STEP 7**: Protein homology prediction. We tried to predict the product
of the transcripts by local alignment using BLASTP (Altschul et al.
1990; **shiryev2007improved?**; Camacho et al. 2009). We annotated the
transcriptome using four databases: NCBI Nr, NCBI RefSeq Protein (Sayers
et al. 2021), UniProt TrEMBL and UniProt TrEMBL limited to molluscan
entries (UniProt Consortium 2021).

``` bash
# Identify the homologous proteins 
blastp \
  -num_threads $SLURM_NTASKS \
  -task blastp \
  -evalue 1e-6 \
  -max_hsps 1 \
  -max_target_seqs 1 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs" \
  -db dbname \
  -query sample.good.rmdup.slugs.faa \
  -out sample.good.rmdup.slugs.dbname.tsv
```

**STEP 8**: Protein domain identification. We identified the different
protein domains from the filtered transcripts using HMMER (Eddy 2011)
and the database Pfam (Mistry et al. 2021).

``` bash
# Identify the different protein domains
hmmsearch \
  --cpu $SLURM_NTASKS \
  --tblout sample.good.rmdup.slugs.seqdom.tsv \
  --domtblout sample.good.rmdup.slugs.ptrdom.tsv \
  -E 1e-6 \
  --domE 1e-6 \
  --incE 1e-6 \
  --incdomE 1e-6 \
  pfam \
  sample.good.rmdup.slugs.faa
```

**STEP 9**: Gene Ontology annotation. We predicted the GO associated
with each transcript using the homologous proteins (STEP 7 using TrEMBL
DB) and the protein domains (STEP 8). Additionally, we required an extra
file containing the conversion between protein accession and the GO IDs;
this information can be obtained using the ID mapping tool.

``` bash
# Extract the protein accessions to convert into GO IDs
awk '{print $1}' sample.good.rmdup.slugs.trembl.tsv | sort -u > sample.good.rmdup.slugs.prot2go.tsv
```

``` r
# Attach the packages
library(elvira)

# Add GO annotations to the transcripts 
txann <- go_annotation( 
  trans_blast = "sample.good.rmdup.slugs.trembl.tsv", 
  trans_hmm   = "sample.good.rmdup.slugs.seqdom.tsv", 
  unip_goid   = "sample.good.rmdup.slugs.prot2go.tsv"
  )

# Export the annotation
write_delim(txann, file = "sample.good.rmdup.slugs.trans2go.tsv", delim = "\t", col_names = FALSE)
```

**STEP 10**: Additional functional annotation using an orthologs
database. After filtering the sequences and obtaining the final
transcriptome (`sample.good.rmdup.slugs.faa`) we sent it for functional
annotation using the eggNOG webportal (Huerta-Cepas et al. 2019;
Cantalapiedra et al. 2021) with the following options:

- protein sequences
- minimum hit e-value = `1e-6`
- percentage identity = `60`
- minimum % of query coverage = `50`
- Gene Ontology evidence = `Transfer all annotations`
- PFAM refinement = `realign queries to the whole PFAM DB`
- SMART annotation = `Perform SMART annotation`.

The other options were set by default.

## Additional steps

**Over-represented sequences analysis**: We extracted the
overrepresented sequences identified using FastQC and aligned them to
the Nucleotide database. The result is showed in the Table S2.

``` r
# Attach the packages
library(fastqcr)
library(dplyr)
library(Biostrings)

# Extract the overrepresented sequences and their abundance
overseq_tbl <- full_join(
  qc_plot(qc_read("sample_1_fastqc.zip"), "Overrepresented sequences") %>% select(Sequence, Count) %>% rename("frcts" = Count),
  qc_plot(qc_read("sample_2_fastqc.zip"), "Overrepresented sequences") %>% select(Sequence, Count) %>% rename("rvcts" = Count)
  )

# Export overrepresented sequences
overseq <- DNAStringSet(overseq_tbl %>% pull(Sequence))
names(overseq) <- paste0("ELVIRA_OVERSEQ-", 1:nrow(overseq_tbl))
writeXStringSet(x = overseq, filepath = "sample.overseq.fna", append = FALSE)
```

``` bash
# Align the overrepresented sequences to the Nucleotide database 
blastn \
  -num_threads $SLURM_NTASKS \
  -task blastn-short \
  -db nt \
  -query sample.overseq.fna
  -out sample.overseq.nt.tsv \
  -evalue 1e-3 \
  -outfmt 6 \
  -perc_identity 75 \
  -max_target_seqs 1 \
  -max_hsps 1
```

**Transcriptome fragmentation**: We checked if the transcript were
assembled completely and the proteins coverage (we required the
alignment obtained in the STEP 7).

``` r
# Calculate fragmentation and protein coverage
trans_frag <- check_completeness(
  prot_seqs  = "sample.good.rmdup.slugs.faa",
  prot_blast = "sample.good.rmdup.slugs.nr.tsv"
)

# Export the result
write_delim(trans_frag, file = "sample.fragmentation.tsv", delim = "\t", append = FALSE, col_names = TRUE)
```

**Orthogroups detection**: We detected the different orthogroups using
OrthoFinder (Emms and Kelly 2019). After, these sequences were aligned
using MAFFT (Katoh et al. 2002). OrthoFinder is not reporting the final
species tree correctly (Issue
\#732)\[<https://github.com/davidemms/OrthoFinder/issues/732>\] so, we
had to build it using an alternaltive approach.

``` bash
# Create a dump directory to store the different transcriptomes
mkdir ref_sequences

# Run OrthoFinder
orthofinder \
  -f ref_sequences \
  -a $SLURM_NTASKS \
  -M msa \
  -A mafft \
  -T raxml-ng \
  -o slugs_ortphy
```

**Building the species tree**:

``` bash
# Estimate the substitution model
modeltest-ng \
  --force \
  --processes $SLURM_NTASKS \ 
  --datatype aa \
  --input slugs_ortphy/*/MultipleSequenceAlignments/SpeciesTreeAlignment.fa \
  --output SpeciesTreeAlignment.model \
  --topology ml \
  --frequencies e \
  --model-het f \
  --template raxml

# Build the ML tree
raxml-ng \
  --redo \
  --threads $SLURM_NTASKS \
  --all \
  --check \
  --msa slugs_ortphy/*/MultipleSequenceAlignments/SpeciesTreeAlignment.fa \
  --model MODEL \
  --blopt nr_safe \
  --bs-trees 10000 
```

**Orthogroups gene-ontology enrichment**

``` r
# Attach the libraries
library(topGO)
library(readr)
library(dplyr)
library(stringr)

# Load the Orthogroups information
cnames      <- c("Orthogroup", "ACA", "ECH", "ECO", "ECR", "EOR", "ETI", "EVI", "OVI", "POC")
orthogroups <- read_delim(file = "orthofinder_dir/Orthogroups.tsv", delim = "\t", col_names = cnames, skip = 1) 

# Load ELVIRA annotation
txannotation <- readMappings(file = "sample.good.rmdup.slugs.trans2go.tsv", sep = "\t", IDsep = ";")  
```

``` r
# Pick up the genes of interest
# KLEPTOPLASTY
txsubset <- orthogroups %>% 
  filter(is.na(ACA), is.na(OVI), !is.na(ECO), !is.na(EOR), !is.na(ECR), !is.na(ETI), !is.na(ECH), !is.na(EVI), !is.na(POC)) %>%
  pull(EVI) 

# LONG-TERM RETAINERS
txsubset <- orthogroups %>% 
  filter(is.na(ACA), is.na(OVI), is.na(ECO), is.na(EOR), !is.na(ECR), !is.na(ETI), !is.na(ECH), !is.na(EVI), !is.na(POC)) %>%
  pull(EVI) 

# LTR ULVOPHYCEAE-FEEDER
txsubset <- orthogroups %>% 
  filter(is.na(ACA), is.na(OVI), is.na(ECO), is.na(EOR), !is.na(ECR), !is.na(ETI), is.na(ECH), !is.na(EVI), !is.na(POC)) %>%
  pull(EVI) 

# Select the genes of interest 
txsubset <- unlist(str_split(string = txsubset, pattern = ", "))
txsubset <- factor(as.integer(names(txannotation) %in% txsubset))
names(txsubset) <- names(txannotation)
```

``` r
# Perform GO-term enrichment
# Biological Process
bp_godata <- new("topGOdata", ontology = "BP", allGenes = txsubset, annot = annFUN.gene2GO, gene2GO = txannotation)
bp_enrichment <- runTest(object = bp_godata, algorithm = "weight01", statistic = "t")
bp_enrichment <- GenTable(object = bp_godata, pvalue = bp_enrichment, topNodes = 2583)

# Cellular Component
cc_godata <- new("topGOdata", ontology = "CC", allGenes = txsubset, annot = annFUN.gene2GO, gene2GO = txannotation)
cc_enrichment <- runTest(object = cc_godata, algorithm = "weight01", statistic = "t")
cc_enrichment <- GenTable(object = cc_godata, pvalue = cc_enrichment, topNodes = 683)

# Molecular Function
mf_godata <- new("topGOdata", ontology = "MF", allGenes = txsubset, annot = annFUN.gene2GO, gene2GO = txannotation)
mf_enrichment <- runTest(object = mf_godata, algorithm = "weight01", statistic = "t")
mf_enrichment <- GenTable(object = mf_godata, pvalue = mf_enrichment, topNodes = 1606)

# Binding enrichment result
go_enrichment <- tibble(bind_rows(
  bp_enrichment %>% mutate("go_class" = "BP"),
  cc_enrichment %>% mutate("go_class" = "CC"),
  mf_enrichment %>% mutate("go_class" = "MF")
)) %>%
  mutate(pvalue   = str_replace(string = pvalue, pattern = "< ", replacement = ""),
         pvalue   = as.numeric(pvalue),
         go_class = factor(go_class, levels = c("BP", "CC", "MF")))

# Export the result
write_delim(x = go_enrichment, file = "go_enrichment.tsv", delim = "\t", append = FALSE, col_names = TRUE)
```

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-altschul1990blast" class="csl-entry">

Altschul, Stephen F, Warren Gish, Webb Miller, Eugene W Myers, and David
J Lipman. 1990. “Basic Local Alignment Search Tool.” *Journal of
Molecular Biology* 215 (3): 403–10.

</div>

<div id="ref-camacho2009blastplus" class="csl-entry">

Camacho, Christiam, George Coulouris, Vahram Avagyan, Ning Ma, Jason
Papadopoulos, Kevin Bealer, and Thomas L Madden. 2009. “BLAST+:
Architecture and Applications.” *BMC Bioinformatics* 10 (1): 1–9.

</div>

<div id="ref-cantalapiedra2021eggnogmapper" class="csl-entry">

Cantalapiedra, Carlos P, Ana Hernández-Plaza, Ivica Letunic, Peer Bork,
and Jaime Huerta-Cepas. 2021. “eggNOG-Mapper V2: Functional Annotation,
Orthology Assignments, and Domain Prediction at the Metagenomic Scale.”
*Molecular Biology and Evolution* 38 (12): 5825–29.

</div>

<div id="ref-chen2018fastp" class="csl-entry">

Chen, Shifu, Yanqing Zhou, Yaru Chen, and Jia Gu. 2018. “Fastp: An
Ultra-Fast All-in-One FASTQ Preprocessor.” *Bioinformatics* 34 (17):
i884–90.

</div>

<div id="ref-deVries2014plastid" class="csl-entry">

de Vries, Jan, Gregor Christa, and Sven B Gould. 2014. “Plastid Survival
in the Cytosol of Animal Cells.” *Trends in Plant Science* 19 (6):
347–50.

</div>

<div id="ref-deVries2014guide" class="csl-entry">

de Vries, Jan, Cessa Rauch, Gregor Christa, and Sven B Gould. 2014. “A
Sea Slug’s Guide to Plastid Symbiosis.” *Acta Societatis Botanicorum
Poloniae* 83 (4).

</div>

<div id="ref-eddy2011accelerated" class="csl-entry">

Eddy, Sean R. 2011. “Accelerated Profile HMM Searches.” *PLoS
Computational Biology* 7 (10): e1002195.

</div>

<div id="ref-emms2019orthofinder" class="csl-entry">

Emms, David M, and Steven Kelly. 2019. “OrthoFinder: Phylogenetic
Orthology Inference for Comparative Genomics.” *Genome Biology* 20 (1):
1–14.

</div>

<div id="ref-haas2013Trinity" class="csl-entry">

Haas, Brian J, Alexie Papanicolaou, Moran Yassour, Manfred Grabherr,
Philip D Blood, Joshua Bowden, Matthew Brian Couger, et al. 2013. “De
Novo Transcript Sequence Reconstruction from RNA-Seq Using the Trinity
Platform for Reference Generation and Analysis.” *Nature Protocols* 8
(8): 1494–1512.

</div>

<div id="ref-handeler2009functional" class="csl-entry">

Händeler, Katharina, Yvonne P Grzymbowski, Patrick J Krug, and Heike
Wägele. 2009. “Functional Chloroplasts in Metazoan Cells-a Unique
Evolutionary Strategy in Animal Life.” *Frontiers in Zoology* 6 (1):
1–18.

</div>

<div id="ref-huerta2019eggnogort" class="csl-entry">

Huerta-Cepas, Jaime, Damian Szklarczyk, Davide Heller, Ana
Hernández-Plaza, Sofia K Forslund, Helen Cook, Daniel R Mende, et al.
2019. “eggNOG 5.0: A Hierarchical, Functionally and Phylogenetically
Annotated Orthology Resource Based on 5090 Organisms and 2502 Viruses.”
*Nucleic Acids Research* 47 (D1): D309–14.

</div>

<div id="ref-jensen2007biogeography" class="csl-entry">

Jensen, Kathe R. 2007. “Biogeography of the Sacoglossa (Mollusca,
Opisthobranchia).” *Bonner Zoologische Beiträge* 55 (3/4): 255–81.

</div>

<div id="ref-katoh2002mafft" class="csl-entry">

Katoh, Kazutaka, Kazuharu Misawa, Kei-ichi Kuma, and Takashi Miyata.
2002. “MAFFT: A Novel Method for Rapid Multiple Sequence Alignment Based
on Fast Fourier Transform.” *Nucleic Acids Research* 30 (14): 3059–66.

</div>

<div id="ref-manni2021busco" class="csl-entry">

Manni, Mosè, Matthew R Berkeley, Mathieu Seppey, Felipe A Simão, and
Evgeny M Zdobnov. 2021. “BUSCO Update: Novel and Streamlined Workflows
Along with Broader and Deeper Phylogenetic Coverage for Scoring of
Eukaryotic, Prokaryotic, and Viral Genomes.” *Molecular Biology and
Evolution* 38 (10): 4647–54.

</div>

<div id="ref-mistry2021pfam" class="csl-entry">

Mistry, Jaina, Sara Chuguransky, Lowri Williams, Matloob Qureshi,
Gustavo A Salazar, Erik LL Sonnhammer, Silvio CE Tosatto, et al. 2021.
“Pfam: The Protein Families Database in 2021.” *Nucleic Acids Research*
49 (D1): D412–19.

</div>

<div id="ref-sayers2021ncbidb" class="csl-entry">

Sayers, Eric W, Jeffrey Beck, Evan E Bolton, Devon Bourexis, James R
Brister, Kathi Canese, Donald C Comeau, et al. 2021. “Database Resources
of the National Center for Biotechnology Information.” *Nucleic Acids
Research* 49 (D1): D10.

</div>

<div id="ref-smith2016transrate" class="csl-entry">

Smith-Unna, Richard, Chris Boursnell, Rob Patro, Julian M Hibberd, and
Steven Kelly. 2016. “TransRate: Reference-Free Quality Assessment of de
Novo Transcriptome Assemblies.” *Genome Research* 26 (8): 1134–44.

</div>

<div id="ref-uniprot2021uniprot" class="csl-entry">

UniProt Consortium. 2021. “UniProt: The Universal Protein Knowledgebase
in 2021.” *Nucleic Acids Research* 49 (D1): D480–89.

</div>

</div>
