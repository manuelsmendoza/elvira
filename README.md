
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

# References

<div id="refs" class="references csl-bib-body hanging-indent">

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

<div id="ref-handeler2009functional" class="csl-entry">

Händeler, Katharina, Yvonne P Grzymbowski, Patrick J Krug, and Heike
Wägele. 2009. “Functional Chloroplasts in Metazoan Cells-a Unique
Evolutionary Strategy in Animal Life.” *Frontiers in Zoology* 6 (1):
1–18.

</div>

<div id="ref-jensen2007biogeography" class="csl-entry">

Jensen, Kathe R. 2007. “Biogeography of the Sacoglossa (Mollusca,
Opisthobranchia).” *Bonner Zoologische Beiträge* 55 (3/4): 255–81.

</div>

</div>
