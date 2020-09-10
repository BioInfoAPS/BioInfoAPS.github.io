---
title: "BLAST search"
teaching: 20
exercises: 30
questions:
- "Key question (FIXME)"
objectives:
- "First learning objective. (FIXME)"
keypoints:
- "First key point. Brief Answer to questions. (FIXME)"
---

**B**asic **L**ocal **A**lignment **S**earch **T**ool ([Altschul S et al.](https://doi.org/10.1016/s0022-2836(05)80360-2)) is a set of programs 
that search sequence database for statistically significant similarities. 
There are five traditional BLAST programs available in BLAST+ tools in commandline:

| Program | Description |
| ------- | ----------- |
| BLASTN | compares a nucleotide query sequence against a nucleotide sequence database. |
| BLASTP | compares an amino acid query sequence against a protein sequence database. |
| BLASTX | compares a nucleotide query sequence translated in a ll reading frames against a protein sequence database. |
| TBLASTN | compares a protein query sequence against a nucleotide sequence database dynamically translated in all reading frames. |
| TBLASTX | compares the six-frame translations of a nucleotide query sequence against the six-frame translations of a nucleotide sequence database. |

These are the same blast applications what you can see in the [NCBI blast website](https://blast.ncbi.nlm.nih.gov/Blast.cgi). 
Read more about command line [BLAST 
here](https://open.oregonstate.education/computationalbiology/chapter/command-line-blast/).

## BLAST terminology

- Query
  - The sequence for which the search is performed.

- Database
  - A set of sequences on which search is performed to find mathces for the query sequence.

- Subject
  - A single sequence within the database.

-	Score (S)
    -	A value representing quality of each alignment.
    - It is calculated as the sum of substitution and gap scores.

-	Expectation value (E-value or E).
    -	A value representing probabilistic significance of each alignment.
    - A number of different alignments with scores equivalent to or better than ‘S’ that are expected to occur in database search by chance.
    -	The lower the E value, the more significant the score.


## Application of BLAST

BLAST can be run online (e.g. NCBI blast) as well as locally. Few applications of local BLAST are:

1. to compare against multiple genomes of interest,
2. to output in different formats that can be used in downstream applications.

## Outline of sequences to BLAST

We will use the local BLAST to identify if the bacterial genomes carry a gene of interest. 
We will use bacterial spot causing *Xanthomonas* and one of the conserved effector ‘AvrBs2’ sequence as an example today. More information on *Xanthomonas* is [available here](https://www.nature.com/articles/s41579-020-0361-8). 

Note: Effectors are proteins secreted by pathogenic bacteria into the host cells. 
Note: *avrBs2* is one of the most studied and conserved effector gene found in multiple *Xanthomonas* species.

In the folder `/blue/share/xanthomonas/`, there are few *Xanthomonas* species genomes.

- *Xanthomonas euvesicatoria* &rarr; `Xeu.fasta`
- *Xanthomonas perforans* &rarr; `Xp.fasta`
- *Xanthomonas gardneri* &rarr; `Xg.fasta`
- *Xanthomonas citri* &rarr; `Xc.fasta`

The folder also contains nucleotide sequence for *avrBs2* gene as `avrBs2.fas`.

Our objective would be to see if the bacterial genomes carry the effector gene.

Before starting, copy all the required files to your working directory.

~~~
$ cd /blue/general_workshop/<username>

$ cp -r ../share/xanthomonas ./

$ ls
~~~
{: .language-bash}

~~~
demo     strep     slurm     xanthomonas
~~~
{: .output}

~~~
$ cd xanthomonas

$ ls
~~~
{: .language-bash}

~~~
avrBs2.fas     Xc.fasta     Xeu.fasta     Xg.fasta     Xp.fasta
~~~
{: .output}

### Loading BLAST in Hipergator cluster

First, we need the BLAST program. 
BLAST tools is already installed in HiperGator, but not available by default.
We can use `module load` command to load the program.
~~~
$ module load ncbi_blast
~~~
{: .language-bash}

### Creating BLAST database

Next, we need to specify the database to blast against, 
which in our case are the genome files `Xeu.fasta`, `Xp.fasta`, `Xg.fasta`, and `Xc.fasta`.
`makeblastdb` command, which is a part of BLAST+ package, 
is used for creating BLAST databases from FASTA files.

~~~
$ makeblastdb -in Xeu.fasta -out Xeu -dbtype nucl
~~~
{: .language-bash}

~~~
Building a new DB, current time: 09/15/2020 02:25:22
New DB name:   /blue/general_workshop/<username>/xanthomonas/Xeu
New DB title:  Xeu.fasta
Sequence type: Nucleotide
Keep MBits: T
Maximum file size: 1000000000B
Adding sequences from FASTA; added 3 sequences in 0.0926349 seconds.
~~~
{: .output}

~~~
$ ls
~~~
{: .language-bash}

~~~
avrBs2.fas     Xeu.ndb     Xeu.not     Xeu.ntf     Xg.fasta
Xc.fasta       Xeu.nhr     Xeu.nsq     Xeu.nto     Xp.fasta
Xeu.fasta      Xeu.nin
~~~
{: .output}

> New database files with extensions `.ndb`, `.nhr`, `.nsq` etc. have been created.
{: .notes}

> ## Understanding command arguments
> For most unix commandline commands, you can use `-h` or `-help` argument to see help or options.
> This will include list of arguments and their meanings.
> In this case, `makeblastdb -h` command displays options avaiable for `makeblastdb`.
{: .tips}

### Creating BLAST database in batch

Making database individually for all genomes will be manually laborious. 
Instead, we can identify common patterns in the name of the genome files and 
execute `makeblastdb` in a loop.

> What is one pattern that is common in all the genomes file??
> <details>
> <summary></summary>
> The extension - ‘.fasta’.
> </details>
{: .challenge}

We can specify all files in current path with fasta extension by `./*.fasta`.

~~~
$ genomes=`ls *.fasta | sed 's/.fasta//g'`

$ for genome in $genomes; do makeblastdb -in "$genome.fasta" -out $genome -dbtype nucl; done
~~~
{: .language-bash}

~~~
...
...
~~~
{: .output}

> `sed` command is used to remove `.fasta` extension from list of names.
{: .notes}

> For maually selecting the databases, you can use `for` loop like this:
> `for genome in Xp Xg Xc; do makeblastdb -in "$genome.fasta" -out $genome -dbtype nucl; done`
{: .tips}

> ## One liners
> Short loops can be written in a same line by separating commands with `;`. 
>`;` is equivalent to pressing <kbd>Enter</kbd>.
{: .tips}

> ## Merging database with `blastdb_aliastool`
> BLAST+ also includes a command `blastdb_aliastool` for combining databases; 
> however, it is outside the scope of this workshop.  
> Usage:  
> `blastdb_aliastool -dblist "Xeu Xp Xg Xc" -dbtype nucl -out Xspp -title "Xanthomonas genomic"`
{: .notes}

### Performing BLAST search

We are now ready to do a BLAST search. Since both the query (`avrBs2.fasta`) 
and the database are nucleotide sequences, we will perform `blastn`.

~~~
$ blastn -query avrBs2.fas -db Xeu -out Xeu_avrBs2.out -evalue 0.001

$ ls *.out
~~~
{: .language-bash}

~~~
Xeu_avrBs2.out
~~~
{: .output}

### Performing BLAST search in multiple databases in batch

We can blast multiple databases in a loop as well. 
Let's do this in a different way from previous example.
First lets create a list of all databases to BLAST against and save it into a file.
We can use `nano` for writing to a file.

~~~
$ nano
~~~
{: .language-bash}

A basic text editor will open. Type in all the databases to search against. 
~~~
-----------------------------------------------------------------------------------------------
 GNU nano 3.3 beta 02                      New Buffer
-----------------------------------------------------------------------------------------------
Xp
Xg
Xc

-----------------------------------------------------------------------------------------------
^G Get Help     ^O WriteOut     ^R Read File     ^Y Prev Page     ^K Cut Text       ^C Cur Pos
^X Exit         ^J Justify      ^W Where Is      ^V Next Page     ^U UnCut Text     ^T To Spell
-----------------------------------------------------------------------------------------------
~~~
{: .terminal}

You do not have to press enter and go to the next line after typing 'Xc'. The next line will be read as an empty line and give error message. Although the additional line and error message due to the empty line does not affect our outputs here, it might be relevant for other situations. 

Press <kbd>Ctrl</kbd>+<kbd>o</kbd> (<kbd>Cmd</kbd>+<kbd>o</kbd> in MacOS) to save the file.
Give it a name `dblist.txt` and press <kbd>Enter</kbd>.

Press <kbd>Ctrl</kbd>+<kbd>x</kbd> (<kbd>Cmd</kbd>+<kbd>x</kbd> in MacOS) to return to bash prompt.

Now we can run the `blastn` in loop.

~~~
$ while read -r dbname
> do
>   blastn -query avrBs2.fas -db "$dbname" -out $dbname"_avrBs2.out" -evalue 0.001
> done < dblist.txt
~~~
{: .language-bash}

> After copying the code block, do not forget to remove the `$` or `>` sign.
> Or click [here](/code/blastn_loop.txt){: target="_blank"} to open a separate page to copy the code.
{: .caution}

> ## Exercise: Performing blast search in SLURM
> 
> We can also run `blast` as a SLURM job, which is useful for long and resource intensive search.
> 
> The SLURM submission has been prepared for your and is available as 
> `/blue/general_workshop/share/scripts/slurm_blast.sh`. 
> Genome and query files are available in `/blue/general_workshop/share/xanthomonas`.
> You can use the checklist to mark progress.
> 
> 1. Change your location to your working directory `/blue/general_workshop/'username';` <input type="checkbox">
> 2. Make a folder in your working directory called `slurm_blast` and enter that directory. <input type="checkbox">
> 3. Copy all genome and query fasta files to current directory. <input type="checkbox">
> 4. Copy the submission script from to the current directory. <input type="checkbox">
> 5. Open the script in nano and edit the email address. <input type="checkbox">
> 6. Submit the job to SLURM. <input type="checkbox">
> 8. After job is completed, check if output files exist in current directory. <input type="checkbox">
> 
> <details markdown="1">
>   <summary></summary>
> 
> ~~~
> #1
> $ cd /blue/general_workshop/<username>
> 
> #2
> $ mkdir slurm_blast
> 
> $ cd slurm_blast
> 
> #3
> $ cp ../../share/xanthomonas/* ./
> 
> #4
> $ cp ../../share/scripts/slurm_blast.sh ./
> 
> #5
> $ nano slurm_blast.sh
> → edit email address → ctrl+x → y
> 
> #6
> $ sbatch slurm_blast.sh
> 
> #7
> $ ls *.out
> ~~~
> {: .language-bash}
> 
> </details>
{: .challenge}
