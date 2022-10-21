---
title: "ML trees with bootstraps"
objectives:
- Creating maximum likelihood trees with RAxML
- Creating bootstrapped trees with RAxML
- Creating bipartition tree (ML + bootstrap)
---

We showed how to generate a simple tree in the lesson "Phylogenetic tree".
However, for publication, you are expected to show more evidence, 
such as the stability of tree topology.
One way to show this is using bootstrap values, 
which is the probability that a particular node
(i.e. a dichotomous branching with a particular set of samples in each branch)
appears among a large number of trees generated 
while resampling within the sequence alignment.
Bootstrap values are often used in the same fashion as confidence intervals.

A common method for creating bootstrapped trees using RAxML consists of 
a 3 step approach:
- Generate a few maximum likelihood (ML) trees.
- Generate many bootstrapped trees.
- Apply the bootstrap information to the best ML tree.

### Generating an ML tree

Maximum likelihood tree generation is computationally expensive, 
but the resulting tree is considered superior to other rapid methods.
Thus, we' usually generate a small number of ML trees (16 in the following example).

~~~
$ raxmlHPC -T 8 -m GTRGAMMA -p 144 -# 16 -s input.fasta -n treeML -w outdir

$ ls -t outdir/*treeML*
~~~
{: .language-bash}

> Output directory specified by `-w` must be an absolute path in RAxML.
> `$(pwd)/outdir` may be used if **outdir** is a relative path.
> This directory needs to be created prior to running the command.
> Alternatively, you can use the current directory as output directory
> by not including this `-w` parameter.
{: .caution}

~~~
RAxML_parsimonyTree.treeML.RUN.0  RAxML_parsimonyTree.treeML.RUN.1
RAxML_log.treeML.RUN.0            RAxML_parsimonyTree.treeML.RUN.2
RAxML_log.treeML.RUN.1            RAxML_parsimonyTree.treeML.RUN.3
...
...
RAxML_result.treeML.RUN.0         RAxML_result.treeML.RUN.1
RAxML_result.treeML.RUN.2         RAxML_result.treeML.RUN.3
...
...
RAxML_info.treeML                  RAxML_bestTree.treeML
~~~
{: .output}

RAxML will automatically select the best tree among the outputs 
and store it in the file **RAxML_bestTree.xxx**.

In the command above,
- `-T` specifies number of CPU threads to be used.
- `-m` specified the substitution model.
- `-p` specifies random seed for starting parsimony tree.
- `-#` specifies the number of trees to generate using unique starting tree.
- `-s` specifies input file containing sequence alignment.
- `-n` specifies suffix for output files.
- `-w` specifies output directory.

### Generating bootstraps

Next, we can generate a large number of computationally permissive trees
for calculating bootstrap values.

~~~
$ raxmlHPC -T 8 -m GTRGAMMA -p 144 -b 144 -# 1000 -s input.fasta -n treeML -w outdir

$ ls outdir/*treeBS*
~~~
{: .language-bash}

~~~
RAxML_info.treeBS    RAxML_bootstrap.treeBS
~~~
{: .output}

In the command above, `-b` specifies bootstrapping with supplied random seed, 
and `-#` specifies the number of bootstraps.

> A newer rapid bootstrap method
> [↗](https://doi.org/10.1080/10635150802429642){:target="_blank"}
> can be employed in place of standard bootstraping
> [↗](https://doi.org/10.1111/j.1558-5646.1985.tb00420.x){:target="_blank"}
> by using the argument `-x` instead of `-b`.
{: .tips}

RAxML can also perform posterior bootstrap convergence analysis to determine
if the number of bootstraps is adequate.

~~~
$ raxmlHPC -m GTRGAMMA -p144 -z outdir/RAxML_bootstrap.treeBS -I autoMRE -n BStest -w outdir

$ tail -n1 outdir/RAxML_info.BStest
~~~
{: .language-bash}

~~~
Converged after 900 replicates
~~~
{: .output}

In the command above, `-I` initiates convergence testing and 
specifies which criterion to use for the test.
`-z` specifies input bootstrap tree file to test.

### Applying bootstrap values to the best ML tree

The final step is to apply the bootstrap values to the best ML tree.

~~~
$ raxmlHPC -T 8 -m GTRGAMMA -p 144 -f b -t outdir/RAxML_bestTree.treeML -z outdir/RAxML_bootstrap.treeBS -n treeBP -w outdir

$ ls outdir/*treeBP*
~~~
{: .language-bash}

~~~
RAxML_bipartitionsBranchLabels.treeBP    RAxML_bipartitions.treeBP
~~~
{: .output}

In the command above, `-f b` instructs creation of bipartition tree from 
best ML tree (supplied with `-t`) and 
the bootstrap trees (specified with `-z`).

The output files can be used to visualize the trees.
The two output files have similar information except 
the branch support information is supplied in a 
slightly different format (node label vs branch label).
Select the file that is correctly interpreted by your visualization 
program.

### A single-step approach

RAxML can perform all three steps above with a single line of code.
However, only the newer rapid approach can be used for bootstraping.
By default, 20 ML trees are generated.

~~~
$ raxmlHPC -T 8 -m GTRGAMMA -p 144 -f a -x 144 -# 1000 -s input.fasta -n treeALL -w outdir

$ ls outdir/*treeALL*
~~~
{: .language-bash}

~~~
RAxML_bestTree.treeALL                    RAxML_bootstrap.treeALL
RAxML_bipartitionsBranchLabels.treeALL    RAxML_info.treeALL
RAxML_bipartitions.treeALL
~~~
{: .language-bash}

> ## RAxML resources
> - [RAxML v8 manual](https://cme.h-its.org/exelixis/resource/download/NewManual.pdf){:target="_blank"}  
> - [A more extensive RAxML guide by the authors](https://cme.h-its.org/exelixis/web/software/raxml/hands_on.html){:target="_blank"}  
> - [ExaML - parallelized approach for whole genomes datasets](https://github.com/stamatak/ExaML){:target="_blank"}
{: .notes}