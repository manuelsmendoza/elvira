<h1 align=center>
  <p>Transcriptome landscape of kleptoplastic sea slug <i>Elysia viridis</i></p>
</h1>

<p align=center>
  <a href="https://scholar.google.com/citations?user=sLWYp9QAAAAJ&hl=en&oi=ao"><b>Manuel Mendoza</b></a>, 
  <a href="https://scholar.google.com/citations?user=RWCFV-IAAAAJ&hl=en&oi=sra">Sara Rocha</a>, 
  <a href="https://scholar.google.com/citations?user=BmYFiWEAAAAJ&hl=en&oi=ao">Jesus Troncoso</a>,
  <a href="https://scholar.google.com/citations?user=sXJWNsYAAAAJ&hl=en">David Posada</a> and 
  <a href="https://scholar.google.com/citations?hl=en&user=X-KLBboAAAAJ">Carlos A. Canchaya</a>
</p>

Certain sacoglossan sea slugs can sequester and maintain photosynthetically active chloroplasts through algae feeding, a phenomenon called kleptoplasty ([Cruz *et al*., 2013](https://academic.oup.com/jxb/article/64/13/3999/436339); [Rumpho *et al*., 2001](https://www.sciencedirect.com/user/identity/landing?code=Hth8LrKrfJuZBXgyrF2UD68jSyVVVzzZWMlXiELR&state=retryCounter%3D0%26csrfToken%3De3dbbab8-4b7d-4c90-bbb8-e3847797f01a%26idpPolicy%3Durn%253Acom%253Aelsevier%253Aidp%253Apolicy%253Aproduct%253Ainst_assoc%26returnUrl%3D%252Fscience%252Farticle%252Fpii%252FS0944200604700355%26prompt%3Dnone%26cid%3Darp-7e0b9f60-7591-41fa-919c-93101b032496)). Adding to the available data to shed some light on these processes, here we describe for the first time the transcriptomic landscape of the sea slug [*Elysia viridis* (Montagu, 1804)](https://www.marinespecies.org/aphia.php?p=taxdetails&id=139686) assembled *de novo* from a pool of ten individual. The reads and assembled transcriptome were deposited in the NCBI database under the BioProject accession number [**PRJNA549923**](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA549923). 

```
@article{
  title={Transcriptome landscape of kleptoplastic sea slug \textit{Elysia viridis}},
  author={Mendoza, Manuel and Rocha, Sara and Troncoso, Jesus and Posada, David and Canchaya, Carlos A},
  journal={bioRxiv},
  year={2022},
  publisher={Cold Spring Harbor Laboratory}
}
```

## Analysis pipeline
1. Set up the conda environment.
```bash
conda create --yes --name elvira_env
conda activate elvira_env
conda config --add channels bioconda
conda install fastqc==0.11.8 fastp==0.23.2  
```

2. Reads quality assessment.
```bash
fastq \
  --threads $(nproc --all) \
  --outdir reada_qa \
  --noextract \
  --nogroup \
  sample_*.fastq.gz
```

3. Reads trimming and biological contamination filtering.
```bash
fastp \
  --thread $(nproc --all) \
  --in1 sample_1.fastq.gz \
  --in2 sample_2.fastq.gz \
  --out1 sample_1.trim.fastq.gz \
  --out2 sample_2.trim.fastq.gz \
  --json sample_qc.fastp.json \
  --html sample_qc.fastp.html \
  --compression 9 \
  --dedup \
  --detect_adapter_for_pe \
  --trim_poly_g \
  --trim_poly_x \
  --qualified_quality_phred 30 \
  --unqualified_percent_limit 25 \
  --average_qual 30 \
  --low_complexity_filter \
  --length_required 70 \
  --correction \
  --overrepresentation_analysis
```



























