---
title: "Parsing BLAST output"
teaching: 30
exercises: 15
questions:
- How can we modify files in UNIX?
objectives:
- Modify file formats 
- Running bash script
keypoints:
- You can run bash command as a script in .sh file format. 
---

## File parsing

In a bioinformatics project, it is common to use several programs in series
i.e. output of a program is used as input to another program. 
There are many standard data formats such as 
[FASTA](https://en.wikipedia.org/wiki/FASTA_format){: target="_blank"}, 
[SAM and BAM](https://samtools.github.io/hts-specs/SAMv1.pdf){: target="_blank"} which are 
commonly accepted by most programs. 
This makes transition from one program to another program easier.
Unfortunately, sometimes programs output non-standard or 
proprietary formats, causing disconnect in using downstream programs. 
This is where file parsing is important. 
File parsing simply means reading the file and facilitating the data 
transformation for subsequent analysis.

Note: While file parsing does not sound important, a bioinformatician has 
to spend significant time correcting file formats.

## Running BLASTn for phylogenetic tree

In next lesson, we will phylogenetically compare ‘avrBs2’ gene 
from various strains of *Xanthomonas* species.
The ouline of activity will be as follows:

BLAST &rarr; parse to fasta format &rarr; align sequence &rarr; phylogenetic analysis
 
To get started, we will run BLASTn in SLURM just like last section, 
however, we will work with a larger set (20) of *Xanthomonas* genomes.
You can find the genomes in the directory 
`/blue/general_workshop/share/phylogeny`. 

~~~
$ cd /blue/general_workshop/<username>

$ cp -r ../share/phylogeny ./

$ cd phylogeny

$ ls
~~~
{: .language-bash}

~~~
avrBs2.fas     GEV1063.fasta     GEV915.fasta     GEV993.fasta     Xp.fasta
...
...
~~~
{: .output}

The SLURM submission script is located in `/blue/general_workshop/share/scripts/slurm_blast_xml.sh`

> This is not the same script as last section. Use this script, not the last one.
{: .caution}

~~~
$ tail -n9 ../../share/scripts/slurm_blast_xml.sh
~~~
{: .language-bash}

~~~
# Loop over each genome
for genome in `ls *.fasta | sed 's/.fasta//g'`
do
  # Make database
  makeblastdb -in "$genome.fasta" -dbtype nucl -out "$genome"

  # Run blastn on that database
  blastn -query avrBs2.fas -db "$genome" -out $genome"_avrBs2.out" -𝗼𝘂𝘁𝗳𝗺𝘁 𝟱 -evalue 0.001
done
~~~
{: .output}

> This script has a new argument `-outfmt 5` for `blastn` command. This will output the result in XML format.
> The `.out` files you generated in previous section is not in XML format.
{: .notes}



> ## Exercise: Submit a job with BLAST script
> 
> The obejctive of this exercise is to be able to run SLURM script  
> You can use checklist to track progress.
> 1. Copy script `slurm_blast_xml.sh` to your working directory (`phylogeny`). <input type="checkbox">
> 2. Edit email address in SLURM submission script. <input type="checkbox">
> 3. Submit the SLURM script. <input type="checkbox">
> 4. Make sure the all output files are present at the end of the job. <input type="checkbox">
> 
> <details markdown="1">
>   <summary></summary>
> ~~~
> #1
> $ cp ../../share/scripts/slurm_blast_xml.sh ./
> 
> #2
> $ nano slurm_blast_xml.sh
> → edit email address → ctrl+x → y → press enter
> 
> #3
> $ sbatch slurm_blast_xml.sh
>
> #4
> $ ls *_avrBs2.out | wc -l
> ~~~
> {: .language-bash}
>
> The output of answer #4 should be 20 if everything went fine.  
>
> </details>
{: .challenge}

## BLAST XML output Vs FASTA format

Now, we can convert the BLAST `XML` output in a more common format
that can be read by multiple software to run analyses. 
Specifically, we will convert it to `FASTA` format, which we
will subsequently use to generate a phylogenetic tree.

BLAST XML format includes a lot of information on query and database name, 
sequence identity, statistics, alignments and so on.
Below, you can see what BLAST XML output looks like:

~~~
$ cat Xeu_avrBs2.out
~~~
{: .language-bash}

~~~
<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_program>blastn</BlastOutput_program>
  <BlastOutput_version>BLASTN 2.10.1+</BlastOutput_version>
  ...
  ...
  <𝗕𝗹𝗮𝘀𝘁𝗢𝘂𝘁𝗽𝘂𝘁_𝗱𝗯>𝗫𝗲𝘂</𝗕𝗹𝗮𝘀𝘁𝗢𝘂𝘁𝗽𝘂𝘁_𝗱𝗯>
  ...
  ...
</BlastOutput>
...
...
<Hit>
  <Hit_num>1</Hit_num>
  <Hit_hsps>
    <Hsp>
      <Hsp_bit-score>4170.86</Hsp_bit-score>
      <Hsp_score>2258</Hsp_score>
      <Hsp_evalue>0</Hsp_evalue>
      ...
      ...
      <Hsp_qseq>AAGTGCTGGCAACGCGTCCAAACACGAGCAGGCCAGGCAGACCGAGACGGATTGA</Hsp_qseq>
      <𝗛𝘀𝗽_𝗵𝘀𝗲𝗾>𝗔𝗔𝗚𝗧𝗚𝗖𝗧𝗚𝗚𝗖𝗔𝗔𝗖𝗚𝗖𝗚𝗧𝗖𝗖𝗔𝗔𝗔𝗖𝗔𝗖𝗖𝗔𝗚𝗖𝗔𝗚𝗚𝗖𝗖𝗔𝗚𝗚𝗖𝗔𝗚𝗔𝗖𝗖𝗚𝗔𝗚𝗔𝗖𝗚𝗚𝗔𝗧𝗧𝗚𝗔</𝗛𝘀𝗽_𝗵𝘀𝗲𝗾>
      <Hsp_midline>||||||||||||||||||||||||| |||||||||||||||||||||||||||||</Hsp_midline>
    </Hsp>
  </Hit_hsps>
</Hit>
~~~
{: .output}

> The output above is for demonstration purpose and the real output may look different.
{: .notes}

We will have to convert it into FASTA format which only needs two line from the XML:
`BlastOutput_db` and `Hsp_hseq`. The FASTA format looks like:

~~~
> Xeu
AAGTGCTGGCAACGCGTCCAAACACCAGCAGGCCAGGCAGACCGAGACGGATTGA
~~~
{: .terminal}

### Bash script to parse XML to FASTA

We will now use a bash script to parse the BLAST output to 
extract the file as a multiFASTA file.
The script is located in path `/blue/general_workshop/share/scripts/blast2fasta.sh`

~~~
$ cp ../../share/scripts/blast2fasta.sh ./   

$ cat blast2fasta.sh
~~~
{: .language-bash}

~~~
#!/bin/bash
while read line
do
        [[ `grep "BlastOutput_db" <<< $line` ]] && echo -n ">" <<< $line
        [[ `grep -E "BlastOutput_db|Hsp_hseq" <<< $line` ]] && sed -n "s:.*>\(.*\)</.*:\1:p" <<< $line
done
~~~
{: .output}

What this script does is loop line by line (`while`) over the file, 
and look for (`grep`) keywords `BlastOutput_db` and `Hsp_hseq`.
If it find the keywords, then it will extract (`sed`) the contents
between those tags and output in appropriate format.

> Considering the course objective, we won't cover the exact syntax of the script at this time. 
{: .tips}

### Running a bash script

We can call the bash script in two ways

- `./<script_name>`
- `sh <scriptname>`

~~~
$ cat *_avrBs2.out | ./blast2fasta.sh > avrBs2_all_genomes.fas
~~~
{: .language-bash}

The FASTA output is stored in file `avrBs2_all_genomes.fas`.

~~~
$ cat avrBs2_all_genomes.fas
~~~
{: .language-bash}

~~~
>GEV1026
AAGTGCTGGCAACGCGTCCAAACAAGGCCTGCGCCGCACGCCTGCCAGCGCGCGCAACGCAGGCATCGTTTCGCATCCGGG
CGGTACTTTTCGCCTAATTTGCCAATTGTCATATGCCACGCGCTTTACTGGCCGCCCGCCGCGTTTTCGAGGTCATCATGC
GCATCGGTCCTCTGCAACCTTCTATCGCGCACACTGC
...
...
>GEV1001
GTGCTGGCAACGCGTCCAAACAAGGCCTGCGCCGCACGCCTGCCAGCGCGCGTAACGCAGGCATCGTTTCGCATCGGCGGG
CGGTACTTTTCGCCTAATTTGCCAATTGTCAGATGCCACGCGCTTTACTGGCCGCCCGACGCGTTTTCGAGGTCATCATGC
GCATCGGTCCTCTGCAACCTTCTATCGCGCA
...
...
~~~
{: .output}

> The output above is for demonstration purpose and the real output may look different.
{: .notes}

> ## When to use SLURM?
> We are not using SLURM for this job because it is fast and not resource intensive.
> Use SLURM for resource intensive or time-consuming tasks.
{: .tips}

Now we have a multiFASTA file with sequence id for each genome database 
we used along with their respective sequence. 
We can now use this FASTA formatted file to run sequence alignment and 
subsequently construct a tree.
