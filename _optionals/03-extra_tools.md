---
title: "Another bioinformatic tool: PyANI" 
objectives:
- To calculate average nucleotide identity of a genome with related genomes. 

---

In the regular lessons, we implemented three bioinformatic tools: 
`blast` for homology search,
`maaft` for sequence alignment, and
`raxml` for phylogenetic analysis.

In this section, we will discuss two more tools.

## PyANI
PyANI is an open-source python-based tool for calculating 
Average Nucleotide Identity (ANI) between two or more sequences.
When comparing two genomes, first syntenic regions are identified
using tools such as `mummer` or `blast`.
Then the nucleotide identity is calculated in the syntenic regions.

The source code for **PyANI** is available at 
[widdowquinn/pyani](https://github.com/widdowquinn/pyani){: target="_blank"}.
The documentation for basic usage is available 
[here](https://github.com/widdowquinn/pyani/blob/master/README_v_0_2_x.md){: target="_blank"}.

**PyANI** v2 is available in Hipergator, but has to be loaded.
The dependencies, `mummer` and `blast+` will be loaded together with `pyani`.

~~~
$ ml pyani
~~~
{: .language-bash}

~~~
Lmod is automatically replacing "python/3.8" with "pyani/0.2.10".
~~~
{: .output}

We will be using the genomes present in `files/ani` for computing ANI.
The file `UXhortspp.fasta` contains genome of a unknown *X. hortorum* species.
The other sequences are genome of some *X. hortorum* pathovars 
downloaded from NCBI.

> ## Getting genome sequences from NCBI
> PyANI has a script called `genbank_get_genomes_by_taxon.py` to download 
> all genomes for a taxon from NCBI.
> For usage, check the documentation linked above.
{: .tips}

The objective now is to perform pairwise comparisons of all reference genomes
and calculate ANI. This can be performed with following command.

~~~
average_nucleotide_identity.py -i files/ani -o ani -m ANIm -g --gformat png,pdf
~~~
{: .language-bash}

> - `average_nucleotide_identity.py` is the name of the script
> - `-i` is used to specify directory containing input genomes/sequences.
> - `-o` is used to specify output directory.
> Note that the program will exit if this directory preexists.
> - `-m` is used to specify mode for alignment of syntenic region. 
> `ANIm` specifies `mummer` and `ANIb` specifies `blast+`.
> - `-g` is used to generate graphic output, i.e., heatmap.
> `--gformat` specifies the graphic output formats.
{: .notes}

~~~
$ ls ani
~~~
{: .language-bash}

~~~
ANIm_alignment_coverage.pdf  ANIm_hadamard.pdf             ANIm_similarity_errors.pdf
ANIm_alignment_coverage.png  ANIm_hadamard.png             ANIm_similarity_errors.png
ANIm_alignment_coverage.tab  ANIm_hadamard.tab             ANIm_similarity_errors.tab
ANIm_alignment_lengths.pdf   ANIm_percentage_identity.pdf  nucmer_output.tar.gz
ANIm_alignment_lengths.png   ANIm_percentage_identity.png
ANIm_alignment_lengths.tab   ANIm_percentage_identity.tab
~~~
{: .output}

You can now transfer `ANIm_percentage_identity.png` 
to your computer to view the heatmap.
For numeric values, you can use `ANIm_percentage_identity.tab` table.

<img src="/fig/ANIm_percentage_identity.png" height="500px">

Based on ANI, the unknown strain seems to be *X. hortorum pv. gardneri*.
