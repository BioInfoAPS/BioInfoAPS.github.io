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

#  BLAST

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

These are the same blast applications what you can see in the NCBI blast website. 
Read more about command line BLAST 
[here](https://open.oregonstate.education/computationalbiology/chapter/command-line-blast/).

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

## BLAST exercise

We will use the local BLAST to identify if the bacterial genomes carry a gene of interest. 
We will use bacterial spot causing *Xanthomonas* and one of the conserved effector ‘AvrBs2’ sequence as an example today. 

Note: Effectors are proteins secreted by pathogenic bacteria into the host cells. 
Note: *avrBs2* is one of the most studied and conserved effector gene found in multiple *Xanthomonas* species.

In the folder `/blue/share/bacteria_genomes/xanthomonas/`, there are few *Xanthomonas* species genomes.

- *Xanthomonas euvesicatoria* &rarr; `Xeu.fasta`
- *Xanthomonas perforans* &rarr; `Xp.fasta`
- *Xanthomonas gardneri* &rarr; `Xg.fasta`
- *Xanthomonas citri* &rarr; `Xc.fasta`

The folder also contains nucleotide sequence for *avrBs2* gene as `avrBs2.fas`.

Our objective would be to see if the bacterial genomes carry the effector gene.

Before starting, copy all the required files to your working directory.

```sh
$ cd /blue/general_workshop/<username>

$ cp ../share/bacteria_genomes/xanthomonas ./

$ ls
strep     slurm     xanthomonas

$ cd xanthomonas
```

### Loading BLAST in Hipergator cluster

First, we need the BLAST program. 
BLAST tools is already installed in HiperGator, but not available by default.
We can use `module load` command to load the program.
```sh
$ module load ncbi_blast
```
Tip: `ml` is shortcut for `module load`.

### Creating BLAST database

Next, we need to specify the database to blast against, 
which in our case are the genome files `Xeu.fasta`, `Xp.fasta`, `Xg.fasta`, and `Xc.fasta`.
`makeblastdb` command, which is a part of BLAST+ package, 
is used for creating BLAST databases from FASTA files.

```sh
$ makeblastdb -in Xeu.fasta -dbtype nucl

Building a new DB, current time: 09/06/2020 23:17:22
New DB name:   /blue/general_workshop/<username>/xanthomonas/Xeu.fasta
New DB title:  Xeu.fasta
Sequence type: Nucleotide
Keep MBits: T
Maximum file size: 1000000000B
Adding sequences from FASTA; added 3 sequences in 0.0926349 seconds.

$ ls
avrBs2.fasta      Xc.fasta          Xeu.fasta         Xeu.fasta.ndb     Xeu.fasta.nhr 
Xeu.fasta.nin     Xeu.fasta.not     Xeu.fasta.nsq     Xeu.fasta.ntf     Xeu.fasta.nto
Xg.fasta          Xp.fasta
```
Tip: `makeblastdb -h` command displays options avaialbe for `makeblstdb`.

Note that new database files with extensions `.ndb`, `.nhr`, `.nsq` etc. have been created.

### Creating BLAST database in batch

Makeing database individually for all genomes will be manually laborious. 
Instead, we can identify common patterns in the name of the genome files nad 
execute `makeblastdb` in a loop.

Exercise:
What is one pattern that is common in all the genomes file??
Extension - ‘.fasta’.

We can specify all files in current path with fasta extension by `./*.fasta`.

```sh
$ for genome in ./*.fasta; do makeblastdb -in "${genome}" --dbtype nucl; done
```
Tip: Short loops can be written in a same line by separating commands with `;`. 
`;` is equivalent to pressing <kbd>Enter</kbd>.

### Performing BLAST search

We are now ready to do a BLAST search. Since both the query (`avrBs2.fasta`) 
and the database are nucleotide sequences, we will perform `blastn`.

```sh
$ blastn -query avrBs2.fas -db Xeu.fasta -out Xeu_avrBs2.out -outfmt 0 -evalue 0.001

$ ls
avrBs2.fasta      Xeu_avrBs2.out     Xc.fasta          Xeu.fasta         Xeu.fasta.ndb
...
...
```

Tip: Run `blastn -h` for explanation of the arguments used.

### Performing BLAST search in multiple databases in batch

We can blast multiple databases in a lop as well. 
First lets create a list of all databases to BLAST against and save it into a file.
We can use a small text editor program called `nano` for writing to a file.

```sh
$ nano
```

A basic text editor will open. Type in all the databases to search against. 
```
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
```

Note, we are not putting `.fasta` part of the database. This is come handy later in naming outputs.

Press <kbd>Ctrl</kbd>+<kbd>o</kbd> (<kbd>Cmd</kbd>+<kbd>o</kbd> in MacOS) to save the file.
Give it a name `dblist.txt` and press <kbd>Enter</kbd>.

Press <kbd>Ctrl</kbd>+<kbd>x</kbd> (<kbd>Cmd</kbd>+<kbd>x</kbd> in MacOS) to return to bash prompt.

Now we can run the `blastn` in loop.

```sh
$ while read -r dbname
$ do
$   blastn -query avrBs2.fas -db "$dbname.fasta" -out $dbname"_avrBs2.out" -outfmt 0 -evalue 0.001
$ done < dblist.txt
```

Note: BLAST+ also includes a command `blastdb_aliastool` for combining databases; 
however, it is outside the scope of this workshop.

```sh
$ blastdb_aliastool -dblist "Xeu.fasta Xp.fasta Xg.fasta Xc.fasta" -dbtype nucl -out xanthomonas_all -title "Xanthomonas genomic"
```

## Exercise: Performing blast search in SLURM

We can also run `blast` as a SLURM job, which is useful for long and resource intensive search.

The SLURM submission script is available in `/blue/general_workshop/share/scripts/slurm_ncbi.sh`. 
Genome and query files are available in `/blue/general_workshop/share/bacteria_genomes/xanthomonas`

1. Change your location to your working directory `/blue/general_workshop/&lt;username&gt;`
2. Make a folder in your working directory called `slurm_blast` and enter that directory.
3. Copy all genome and query fasta files to current directory.
4. Copy the submission script from to the current directory.
5. Open the script in nano and edit the email address.
6. Submit the job to SLURM.
7. Check status of the job as it is running.
8. After job is completed, check the list of files in current directory.

SLURM script: (/blue/general_workshop/share/scripts/slurm_ncbi.sh)
```
#!/bin/bash
#SBATCH --job-name=serial_job_test    # Job name
#SBATCH --account=general_workshop    # Account to run the computational task
#SBATCH --qos=general_workshop        # Account allocation
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<email_address>   # Where to send mail  
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=1gb                     # Job memory request
#SBATCH --time=00:05:00               # Time limit hrs:min:sec
#SBATCH --output=serial_test_%j.log   # Standard output and error log
pwd; hostname; date

# Load blast
module load ncbi_blast

# Loop over each genome
for genome in `ls *.fasta | sed 's/.fasta//g'`
do
  # Make database
  makeblastdb -in "$genome.fasta" --dbtype nucl

  # Run blastn on that database
  blastn -query avrBs2.fas -db "$genome.fasta" -out $genome"_avrBs2.out" -outfmt 5 -evalue 0.001
done
```

Answer:
```sh
cd /blue/general_workshop/<username>
mkdir slurm_blast
cd slurm_blast
cp /blue/general_workshop/share/bacteria_genomes/xanthomonas/* ./
cp /blue/general_workshop/share/scripts/slurm_ncbi.sh ./
nano slurm_ncbi.sh > edit email > Ctrl+x > y
sbatch slurm_ncbi.sh
squeue -u <username>
ls
```