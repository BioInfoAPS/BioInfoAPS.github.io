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

The first stage of PyANI is 

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

<img src="/assets/img/ANIm_percentage_identity.png" height="500px">

You can also get numeric ANI values from the `ANIm_percentage_identity.tab` file.

~~~
$ awk 'NR==1{print "-"$0; next}{for (i=2; i<=NF; i++) {$i=substr($i,1,5)}; print $0}' ani8/ANIm_percentage_identity.tab | column -t
~~~
{: .language-bash}

~~~
-             Xpopuli  Xhgardneri  Xhcynarae  Xhhederae  Xhvitians  Xhunknown  Xhcarotae  Xhtaraxaci  Xhpelargonii
Xpopuli       1.0      0.913       0.913      0.912      0.912      0.913      0.913      0.912       0.913
Xhgardneri    0.913    1.0         0.993      0.960      0.982      0.999      0.961      0.974       0.958
Xhcynarae     0.913    0.993       1.0        0.960      0.984      0.993      0.961      0.975       0.958
Xhhederae     0.912    0.960       0.960      1.0        0.960      0.960      0.965      0.955       0.961
Xhvitians     0.912    0.982       0.984      0.960      1.0        0.982      0.961      0.976       0.958
Xhunknown     0.913    0.999       0.993      0.960      0.982      1.0        0.961      0.974       0.958
Xhcarotae     0.913    0.961       0.961      0.965      0.961      0.961      1.0        0.956       0.961
Xhtaraxaci    0.912    0.974       0.975      0.955      0.976      0.974      0.956      1.0         0.954
Xhpelargonii  0.913    0.958       0.958      0.961      0.958      0.958      0.961      0.954       1.0
~~~
{: .output}

Based on ANI, the unknown strain seems closest to *X. hortorum pv. gardneri*.
