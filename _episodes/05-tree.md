---
title: "Phylogenetic tree"
teaching: 0
exercises: 0
questions:
- "Key question (FIXME)"
objectives:
- "First learning objective. (FIXME)"
keypoints:
- "First key point. Brief Answer to questions. (FIXME)"
---

## Sequence alignment

[Sequence alignment](https://en.wikipedia.org/wiki/Sequence_alignment) is the basic and the most important step in phylogenetic analysis. 
A wrong sequence alignment will lead to wrong result. 
Various alignment tools are used to identify similarity regions that indicate fuctional, 
structural, and/or evolutionary relationships sequences. 
Multiple sequence alignment tools such as ClustalOmega, Muscle, MAFFT are commonly used. 
They vary in terms of algorithms used; [ClustalOmega](https://www.ebi.ac.uk/Tools/msa/clustalo/) uses HMM profile-profile techniques, 
and [MAFFT](https://mafft.cbrc.jp/alignment/software/) uses Fast Fourier Transforms and 
both are considered good for medium to large alignments. 
[Muscle](https://www.ebi.ac.uk/Tools/msa/muscle/) is regarded good with proteins for medium alignments.

We will use MAFFT sequence alignment tool.

~~~
$ module load mafft

$ mafft avrBs2_all_genomes.fas > avrBs2_all_genomes_aligned.fas
~~~
{: .language-bash}

~~~
nthread = 0
generating a scoring matrix for nucleotide (dist=200) ... done
Gap Penalty = -1.53, +0.00, +0.00

Making a distance matrix ..
    1 / 20
done.
...
...
~~~
{: .output}

~~~
$ more avrBs2_all_genomes_aligned.fas
~~~
{: .language-bash}

> Make sure you are still in phylogeny folder. You can use `pwd` to check.
{: .caution}

The two files avrBs2_all_genomes.fas and avrBs2_all_genomes_aligned.fas 
may not look very different at the moment considering the gene was present 
in all the sample genomes used here as seen from blast results. 
However, aligning the sequence prior to phylogenetic reconstruction is always valuable.
Now, we have an input file ready to use for a software that can run a 
[maximum likelihood phylogentic analysis](https://www.ncbi.nlm.nih.gov/Class/NAWBIS/Modules/Phylogenetics/phylo15.html). 

## Phylogenetic tree

There are several phylogenetic analyses software and tools that can be used. 
We will run a Maximum likelhood phylogenetic analysis. 
There are other models/methods such as Neighbor-Joining, Maximum parsimony, and Bayesian. 
However, our focus is on application of command line, 
we will only focus on running a specific analysis here.
Maximum Likelihood uses probabilistic values to model the evolution and 
the tree with the highest probability is shared.

We will use [RAxML - Randomized Axelerated Maximum Likelihood](https://cme.h-its.org/exelixis/web/software/raxml/) tool to 
construct the tree that we will submit through a batch file.

The SLURM submission script is present in `/blue/general_workshop/share/scripts/slurm_tree.sh.sh

~~~
$ cp ../../share/scripts/slurm_tree.sh ./

$ tail -n8 slurm_tree.sh
~~~
{: .language-bash}

~~~
# Load RAxML
module load raxml

# Make tree with RAxML
raxmlHPC -d -p 12345 -m GTRGAMMAI -s avrBs2_all_genomes_aligned.fas -n avrBs2_tree

# Success message
echo 'ML tree output in RAxML_bestTree.avrBs2_tree'
~~~
{: .output}

The argument `-n` allows us to name the output suffix. In our case, the output will
be name `RAxML_bestTree.avrBs2_tree`.

> ## Help with arguments
> To understand what other arguments mean, you can check the help file `ml raxml; raxmlHPC -h`.
{: .tips}

Edit the email address in the SLURM submission script and submit the job.

> ## Running bootstraps
> To bootstrap the tree, add arguments `-b -#1000` 
> where `1000` is the number of bootstraps.
> If you run bootstraps, the output will be in a file named `RAxML_bootstrap.<suffix>`.
> To understand about the arguments, you can run `ml raxml; raxmlHPC -h`
{: .tips}

We can download this tree only and visualize in our own computer using ‘Figtree’. 
Please use file transfer tool such as FileZilla or Cyberduck to download the file. 
The tree will open as follows:

![Phylogenentic tree](/fig/tree.png)

> ## (Optional) Trees in ASCII!
> To quickly visualize the tree, try this:
> ~~~
> $ ml newick_utils
>
> $ nw_display RAxML_bestTree.avrBs2_tree
> ~~~
> {: .language-bash}
{: .challenge}
---

## Phylogenetic tree pipeline

We have now successfully extracted a gene of interest from our genomes and 
generated a phylogenetic tree. 
However, we completed the entire process using several distinct steps where
we evaluated output of one step before starting to run the next step.
In real world scenario, most of the intermediate results are not of 
interest, so those disjoint steps are often chained
into one or few longer workflow, referred to as **pipeline**. 

Question: Would it be possible to develop a script where we only need to 
provide external (generated outside the pipeline usch as gene and
genome sequences) inputs and get the final result (the phylogennetic tree) 
in return directly?

YES*

For this, we will merge all the steps together into one script.

Warning: *Chaining multiple steps is not always straightforward. 
We have to be vigilant for what might change when the steps are automated and 
intermediate  results are not evlauated.
For example, we have to make sure that the genes of interest are actually in the genome. 
Otherwise, the BLAST may output match to different closely related gene, 
gene sequence alignment may add gaps and give us a wrong output and so on. 
Understanding what might go wrong comes from experience in bioinformatics.

Let's put everything together into one script. The input files are
available in `/blue/general_workshop/share/all_in_one/` and the 
script is available as `/blue/general_workshop/share/scripts/slurm_pipeline.sh`.

~~~
$ cd /blue/general_workshop/<username>

$ mkdir pipeline

$ cd pipeline

$ cp -r ../../share/phylogeny/* ./

$ cp ../../share/scripts/slurm_pipeline.sh ./

$ cp ../../share/scripts/blast2fasta.sh ./

$ ls
~~~
{: .language-bash}

~~~
avrBs2.fas         GEV1054.fasta     GEV909.fasta     GEV968.fasta          Xeu.fasta
blast2fasta.sh     GEV1063.fasta     GEV915.fasta     GEV993.fasta          Xg.fasta
GEV1001.fasta      GEV839.fasta      GEV917.fasta     slurm_pipeline.sh     Xp.fasta
GEV1026.fasta      GEV893.fasta      GEV936.fasta     TB15.fasta
GEV1044.fasta      GEV904.fasta      GEV940.fasta     Xc.fasta
~~~
{: .output}

~~~
$ tail -n23 slurm_pipeline.sh
~~~
{: .language-bash}

~~~
# Load programs
module load ncbi_blast
module load mafft
module load raxml

# Loop over each genome
for genome in `ls *.fasta | sed 's/.fasta//g'`
do
  # Make database
   makeblastdb -in "$genome.fasta" -dbtype nucl -out "$genome"

  # Run blastn on that database
   blastn -query avrBs2.fas -db "$genome" -out $genome"_avrBs2.out" -outfmt 5 -evalue 0.001
done

# Combine all outputs and parse into multiFASTA file
cat *_avrBs2.out | ./blast2fasta.sh > avrBs2_all_genomes.fas

# Align sequences with MAAFT
mafft avrBs2_all_genomes.fas > avrBs2_all_genomes_aligned.fas

# Make tree with RAxML
raxmlHPC -d -p 12345 -m GTRGAMMAI -s avrBs2_all_genomes_aligned.fas -n avrBs2_tree
~~~
{: .output}

After editing the email address in `nano`, you are ready to submit the submission script.

~~~
$ sbatch slurm_pipeline.sh
~~~
{: .language-bash}


> ## Visualize tree
> Once the job finishes, you can transfer all the outputs to your personal computer. 
> Then, you can visualize the ‘RAxML_bestTree.avrBs2_tree’ using Figtree software. 
> [Setup page](/setup.html) has instructions on how to install Figtree.
{: .tips}
