---
title: "FASTA manipulation with awk"
teaching: 
exercises: 
questions:
- 
objectives:
- 
keypoints:
- 

---

## The FASTA format

One of the most common file format when working in bioinformatics is
the FASTA file. FASTA format holds a nucleotide or amino acid sequences,
following a (unique) identifier, called a description line.
A FASTA file is a text file, often with extension `.fasta`, `.fas`, or `.fa`
for both nucleotide and proteins, `.fna` for nucleotide only,
`.faa` for amino acid sequence only.
A multiFASTA file is the concatenation of multiple
FASTA sequences of one type (nucleotide or amino acid).

The description line precedes the sequence.
It starts with a `>` followed by the sequence identifier and
is delimited within a single line for each sequence.
`>` sign, in general, does not appear anywhere else
except in the beginning of the description line.
The sequence starts from immediately next line, and long sequences are
often into many lines often with ~80 characters/bases per line.
Both uppercase and lowercase characters are acceptable,
although some programs assign meaning to the case.
Degenerate bases are accepted, so are `-` for gaps and `*` for stop codon.
Any other characters, including spaces, are ignored.

Here is an example of a valid FASTA file.

~~~
$ cat files/fasta/example.fasta
~~~
{: .language-bash}

~~~
>Xanthomonas perforans strain GV872 contig_50
TGCGGGCGGTGAGGTCGCCCGCCGCGATCGCCTGCAATACCGTCGAGAGCGCGTGCAGATTCGCATCGGAGGTGGCCATC
AATTGATTGAGGCTGTCCACCATCGCGCGGAAGTCGTACTGGAAGCGCTGTGCGTCGCCGCGCACGCTGAAATCGCCATT
GGCTGCGGCCTGCGCAAGCAGCTTGATCTCGCTGTTCATCGCCGACAGGTTG
>Xanthomonas perforans strain GV872 contig_51
CCCTCCACAGGAAGCGGACCGGACAATGACGCTTGATTGCCGCGCAGGAAATGCCCCCGTTCCCGCAGGCATACACCCCA
TTTGCCGCGGACATACGCCGCACATTTCTCGGAGCCTGGACAAGTCGCCACAGCTTTGATGGAACGACGGAAAAAACGCC
CCGGAGTTGACCAGAACAACCGCCAGTACTGAAATCTGGACCTTAAGTTCCAGTTCCTGCAGCTCTGGACCTGGTTCATC
AGCGAGCTTCTATCGCTTCGCGGCGGCTTCGATCGCAGCGTGTGCGCGCCATCACTTGGTCGTGCGAAATCGTCGGCCGC
GTGTCGGCCACGGCGCGCTCGATCTTCGCCCGGAACTACCAGTCGTAGTCCTCGGCGGCTTCGGTCGACTCAAATTCGGA
AACAGGTGGTCGAGTGATACGCCCAGGGTTC
>Xanthomonas perforans strain GV872 contig_52
CGCCAGCGCAGCCGACCACGACGCCCGAATCCAGCGTGCCGGCCTCGTTGCCGCCCGCGTCCACGTCGGGCAATTAAGCC
GCGCCCCGCGTTTCCGGTCATGGAACTTCGCCCAGAGCAGCTGGCCGGCCAGTCCAGCCAGCCATTGCAGCCGGTGTATC
TGATCGCCGGCCCCGAGACCTTGCGTGTGCTTGAGGCCGCCGATGCGGTGCGCGCGCGCGCGCGCGC
>Xanthomonas perforans strain GV872 plasmid_contig_53
GGCGTTGTTCCAGCCCAGCACCATTTCCAACAGATTCAGGGCCTGCGCCGGCGACACCGGATTGGGCGTGGCCAAGGTCA
CGGTGCCCTGCACGCCCGGCGCGATGACGTAGTTCTGGCCGAGCATGTCGCCCAGGATCGCCTTGACCACCGCCTGCAGC
GATTCGCCTTCGAAATTGAAGGTGGCGCTGCCGCTGCTGGCCATGCCCAGCGTG
>Xanthomonas perforans strain GV872 plasmid_contig_54
AAATCCTGGGCGTATCAGACCGCTGCCGAGCGTTTGCTTTGGTACCGCCTGCGTGATCGGCGTCCCCCC/àùbä';ç3cC
¦Å¼OÈ’LäÖ’sëmôzÒÒ${–-ÃêFY;+ã}:ãž0}%ÂZDÅì9D‹Rž÷D×ç°:âñÆó²9h‰ÐžæK¼ë™XÕJK¥·—(ü7¢]ÃÆ
Ño5®o&ÇDçJ^ÃŒZòsÀš‹Ü´†ü¿ÄŒ—xÎ—¸¬µ»?à=Of±)2?îP3ñLapÈëY_{DŸm)ã}ž[Êº„ódMÅbóSl¤ÌŽÑ‚¼
Óh+¶\Éz½/­Z¹â²kgfSû%äÖÿzvŽ#?69éhÎÔøK×‰Êí¨/4ŠSVfÍõÔûÁÒ&½ûE÷kÑ²:k ÷yxÙgH2/ÂžQ,[@ûIß
€|Ë²Ï2,ôCÆ
~~~
{: .output}

### FASTA file processing

The simplicity of FASTA file makes it a good target for learning text processing
for beginner bioinformaticians. To be clear, there are plethora of tools and 
programs that can be used to achieve all kinds of processing
(such as filtration, editing, metadata computation, etc.).
Rather, the objective of this section is to learn to do those processing 
by writing your own code.

There are many tools which can be achieve same processing steps.
In this section, we will use a powerful text processing utility called `awk`, 
which is available by default in all modern UNIX systems.

## Processing a file in `awk`

The syntax for awk to process a file looks like follows.

~~~
$ awk '
    commands to process ...
    ...
  ' input_file

$ cat input_file | awk '
    commands to process ...
    ...
  '
~~~
{: .language-bash}

### Processing one FASTA sequence at a time

Despite its simple structure, processing a FASTA using `awk`
(or most other UNIX tools) poses one important challenge.
`Awk` processes one line at a time, but ideally,
we'd like to process one FASTA sequence at a time.
Fortunately, there are ways to handle this in `awk`.

Consider a simple table:

~~~
abc def ghi
jkl mno pqr
stu vwx yz
~~~
{: .language-bash}

In awk, columns are called fields and rows are called records.
In the table above, space character ` ` separates different columns (or fields)
and is called a field separator.
Similarly, newline character, represented by `\n`, separates different rows
(or records) and is called a record separator.
Input and output files can have different field and record separators in `awk`.
They are stored in following variables within `awk`:

| Variables | Meaning |
|-----------|---------|
| FS | Field separator (input, and output too if OFS not specified) |
| OFS | Output field separator |
| RS | Record separator (input, and output too if ORS not specified) |
| ORS | Output record separator |

A individual FASTA record starts with a `>` character, so if we use `>` as 
the record separator, then `awk` would process one entire FASTA at a time.

> ## Accessing fields in `awk`
> 
> Field values can be accessed in `awk` using numbered variables.
> For the record currently being processed,
> - `$0` stores the entire record
> - `$1` stores first field
> - `$2` stores second field and so on.
> - `NF` stores the total number of fields.
{: .notes}

Similarly, if we assign newline character as the field separator, 
then each line for that FASTA record will act as different columns.
The first column would be the description line, and the remaining columns
would be the sequences for that FASTA record.

This can be accomplished with option `-v` which allows setting a variable.

~~~
$ awk -vRS=">" -vORS="\n" -vFS="\n" -vOFS="\t" '
    {$1=$1; print $0}
  ' files/fasta/example.fasta
~~~
{: .language-bash}

~~~

Xanthomonas perforans strain GV872 contig_50    TGCGGGCGGTGAGGTCGCCCGCCGCGATCGCCTGCAATACCGTCGAGAGCGCGTGCAGATTCGCATCGGAGGTGGCCATC        AATTGATTGAGGCTGTCCACCATCGCGCGGAAGTCGTACTGGAAGCGCTGTGCGTCGCCGCGCACGCTGAAATCGCCATT    GGCTGCGGCCTGCGCAAGCAGCTTGATCTCGCTGTTCATCGCCGACAGGTTG
Xanthomonas perforans strain GV872 contig_51    CCCTCCACAGGAAGCGGACCGGACAATGACGCTTGATTGCCGCGCAGGAAATGCCCCCGTTCCCGCAGGCATACACCCCA        TTTGCCGCGGACATACGCCGCACATTTCTCGGAGCCTGGACAAGTCGCCACAGCTTTGATGGAACGACGGAAAAAACGCC    CCGGAGTTGACCAGAACAACCGCCAGTACTGAAATCTGGACCTTAAGTTCCAGTTCCTGCAGCTCTGGACCTGGTTCATC        AGCGAGCTTCTATCGCTTCGCGGCGGCTTCGATCGCAGCGTGTGCGCGCCATCACTTGGTCGTGCGAAATCGTCGGCCGC   GTGTCGGCCACGGCGCGCTCGATCTTCGCCCGGAACTACCAGTCGTAGTCCTCGGCGGCTTCGGTCGACTCAAATTCGGA AACAGGTGGTCGAGTGATACGCCCAGGGTTC
Xanthomonas perforans strain GV872 contig_52    CGCCAGCGCAGCCGACCACGACGCCCGAATCCAGCGTGCCGGCCTCGTTGCCGCCCGCGTCCACGTCGGGCAATTAAGCC        GCGCCCCGCGTTTCCGGTCATGGAACTTCGCCCAGAGCAGCTGGCCGGCCAGTCCAGCCAGCCATTGCAGCCGGTGTATC    TGATCGCCGGCCCCGAGACCTTGCGTGTGCTTGAGGCCGCCGATGCGGTGCGCGCGCGCGCGCGCGC
Xanthomonas perforans strain GV872 plasmid_contig_53    GGCGTTGTTCCAGCCCAGCACCATTTCCAACAGATTCAGGGCCTGCGCCGGCGACACCGGATTGGGCGTGGCCAAGGTCA        CGGTGCCCTGCACGCCCGGCGCGATGACGTAGTTCTGGCCGAGCATGTCGCCCAGGATCGCCTTGACCACCGCCTGCAGC    GATTCGCCTTCGAAATTGAAGGTGGCGCTGCCGCTGCTGGCCATGCCCAGCGTG
Xanthomonas perforans strain GV872 plasmid_contig_54    AAATCCTGGGCGTATCAGACCGCTGCCGAGCGTTTGCTTTGGTACCGCCTGCGTGATCGGCGTCCCCCC/àùbä';ç3cC
¦Å¼OÈ’LäÖ’sëmôzÒÒ${–-ÃêFY;+ã}:ãž0}%ÂZDÅì9D‹Rž÷D×ç°:âñÆó²9h‰ÐžæK¼ë™XÕJK¥·—(ü7¢]ÃÆ    Ño5®o&ÇDçJ^ÃŒZòsÀš‹Ü´†ü¿ÄŒ—xÎ—¸¬µ»?à=Of±)2?îP3ñLapÈëY_{DŸm)ã}ž[Êº„ódMÅbóSl¤ÌŽÑ‚¼    Óh+¶\Éz½/­Z¹â²kgfSû%äÖÿzvŽ#?69éhÎÔøK×‰Êí¨/4ŠSVfÍõÔûÁÒ&½ûE÷kÑ²:k ÷yxÙgH2/ÂžQ,[@ûIß    €|Ë²Ï2,ôCÆ
~~~
{: .output}

> `$1=$1` is used to recompute the records, so `OFS` can be applied.
{: .notes}

Notice the empty line in the first line of output. 
It exists because `awk` thinks there is a record before the first `>`,
as it is a record separator, not record starter.
To avoid this empty line, we can ask `awk` to not process the first record.
This can be done with a another `awk` variable, `NR`. 
`NR` stands for number of records, and is roughly equivalent to 
row number or line number being processed.

Conditional execution of code can be done my writing the condition 
before the braces.

~~~
$ awk -vRS=">" -vORS="\n" -vFS="\n" -vOFS="\t" '
    NR>1 {$1=$1; print $0} # Only execute if record number is greater than 1.
  ' files/fasta/example.fasta

$ awk -vRS=">" -vORS="\n" -vFS="\n" -vOFS="\t" '
    NR==1 {next} # Skip to next record if record number is equal to 1.
    {$1=$1; print $0}
  ' files/fasta/example.fasta
~~~
{: .language-bash}

## FASTA filtration

If is often necessary to remove one or more FASTA record from a multiFASTA file.
Small files can be edited manually,
but large files or files with a large number of records
can be better handled programmatically.

FASTA filtration is often done using two criteria.
1. Sequence ID, i.e., the description line
2. Length of the sequence

### FASTA filtration by record ID

FASTA filtration by ID is necessary in various cases such as removal of 
plasmid, plastid, or mitochondrial sequence from chromosomal sequence, 
removal of aberrant sequences before alignment etc.

In the code above, the sequence IDs can be accessed within `awk` with 
variable `$1`. Here are a few conditions for string matching.
- `$1=="keyword"`: sequence ID exactly matches `keyword` (except `>` symbol).
- `$1~/keyword/`: sequence ID contains `keyword`.
- `$1~/^keyword/`: sequence ID starts with `keyword`.
    - `^` is resolved to start anchor.
- `$1~/keyword$/`: sequence ID ends with `keyword`.
    - `$` is resolved to end anchor.
- `$1!~/keyword/`: sequence ID does not contain `keyword`.
    - `!` is resolved to **NOT**.
- `^.a`: sequence ID with `a` as the second character.
    - `.` can be resolved into any single character.
- `^\.a`: sequence ID starts with `.a`.
    - `\` is an escape character, thus `.` is used in literal fashion.
- `$1~/[a-c]/`: sequence ID contains `a`, `b`, or `c`.
    - `[...]` is resolved to character group.
- `$1~/^[2-4]/`: sequence ID starts with `2`, `3`, or `4`.
    - `-` is resolved to character range.
- `$1~/2[0-9]$/`: sequence ID ends between `20` and `29`.
- `$1~/19$/ || $1~/2[0-1]$/`: sequence ID ends with `19`, `20`, or `21`.
    - `||` is resolved to **OR**.
- `$1~/key1/ && $1~/key2/`: sequence ID contains both `key1` and `key2` in no order.
    - `&&` is resolved to **AND**.
- `$1~/key1.*key2/`: sequence ID contains both `key1` and `key2` with `key1` before `key2`.
    - `.*` is resolved to any characters, including nothing.
- `$1~/^key1.*key2$/`: sequence ID starts with `key1` and ends with `key2`.

The example FASTA above has two plasmid contigs, 
which can be separated with the term "plasmid" in description line.
To remove plasmid sequences, we can do:

~~~
$ awk -vRS=">" -vORS="" -vFS="\n" -vOFS="\n" '
    NR>1 && $1!~/plasmid/ {print ">"$0}
  ' files/fasta/example.fasta
~~~
{: .language-bash}

~~~
>Xanthomonas perforans strain GV872 contig_50
TGCGGGCGGTGAGGTCGCCCGCCGCGATCGCCTGCAATACCGTCGAGAGCGCGTGCAGATTCGCATCGGAGGTGGCCATC
AATTGATTGAGGCTGTCCACCATCGCGCGGAAGTCGTACTGGAAGCGCTGTGCGTCGCCGCGCACGCTGAAATCGCCATT
GGCTGCGGCCTGCGCAAGCAGCTTGATCTCGCTGTTCATCGCCGACAGGTTG
>Xanthomonas perforans strain GV872 contig_51
CCCTCCACAGGAAGCGGACCGGACAATGACGCTTGATTGCCGCGCAGGAAATGCCCCCGTTCCCGCAGGCATACACCCCA
TTTGCCGCGGACATACGCCGCACATTTCTCGGAGCCTGGACAAGTCGCCACAGCTTTGATGGAACGACGGAAAAAACGCC
CCGGAGTTGACCAGAACAACCGCCAGTACTGAAATCTGGACCTTAAGTTCCAGTTCCTGCAGCTCTGGACCTGGTTCATC
AGCGAGCTTCTATCGCTTCGCGGCGGCTTCGATCGCAGCGTGTGCGCGCCATCACTTGGTCGTGCGAAATCGTCGGCCGC
GTGTCGGCCACGGCGCGCTCGATCTTCGCCCGGAACTACCAGTCGTAGTCCTCGGCGGCTTCGGTCGACTCAAATTCGGA
AACAGGTGGTCGAGTGATACGCCCAGGGTTC
>Xanthomonas perforans strain GV872 contig_52
CGCCAGCGCAGCCGACCACGACGCCCGAATCCAGCGTGCCGGCCTCGTTGCCGCCCGCGTCCACGTCGGGCAATTAAGCC
GCGCCCCGCGTTTCCGGTCATGGAACTTCGCCCAGAGCAGCTGGCCGGCCAGTCCAGCCAGCCATTGCAGCCGGTGTATC
TGATCGCCGGCCCCGAGACCTTGCGTGTGCTTGAGGCCGCCGATGCGGTGCGCGCGCGCGCGCGCGC
~~~
{: .output}

> A `>` was added to print statement to recover the `>` in the front of
> description line.
{: .notes}

### FASTA filtration by length

FASTA filtration by length is necessary cases such as removal of contigs
smaller than a threshold length.

The numbers of characters in a FASTA record can be derived by `awk` function, 
`length`. 
Of course, this way, all characters, 
including those not standard in FASTA are counted.

~~~
$ awk -vRS=">" -vORS="\n" -vFS="\n" -vOFS="\t" '
    NR>1 {a=$1; $1 = ""; print a, length($0)-NF+1}
  ' files/fasta/example.fasta
~~~
{: .language-bash}

~~~
Xanthomonas perforans strain GV872 contig_50    212
Xanthomonas perforans strain GV872 contig_51    431
Xanthomonas perforans strain GV872 contig_52    227
Xanthomonas perforans strain GV872 plasmid_contig_53    214
Xanthomonas perforans strain GV872 plasmid_contig_54    331
~~~
{: .output}

> ## How it works
> `a=$1` stores description line in variable `a`.
> `$1=""` empties description line, so `$0` only holds sequence only.
> `length($0)` calculates length of the sequence.
> `-NF+1` adjusts character length to string as `awk` also counts field separators. 
{: .notes}

Next, we can do filtration.
To remove contigs smaller than 300 base pairs, we can do:

~~~
$ awk -vRS=">" -vORS="" -vFS="\n" -vOFS="\n" '
    NR>1 {a=$0; $1 = ""; if(length-NF+1>=300) {print ">"a}}
  ' files/fasta/example.fasta
~~~
{: .language-bash}

~~~
>Xanthomonas perforans strain GV872 contig_51
CCCTCCACAGGAAGCGGACCGGACAATGACGCTTGATTGCCGCGCAGGAAATGCCCCCGTTCCCGCAGGCATACACCCCA
TTTGCCGCGGACATACGCCGCACATTTCTCGGAGCCTGGACAAGTCGCCACAGCTTTGATGGAACGACGGAAAAAACGCC
CCGGAGTTGACCAGAACAACCGCCAGTACTGAAATCTGGACCTTAAGTTCCAGTTCCTGCAGCTCTGGACCTGGTTCATC
AGCGAGCTTCTATCGCTTCGCGGCGGCTTCGATCGCAGCGTGTGCGCGCCATCACTTGGTCGTGCGAAATCGTCGGCCGC
GTGTCGGCCACGGCGCGCTCGATCTTCGCCCGGAACTACCAGTCGTAGTCCTCGGCGGCTTCGGTCGACTCAAATTCGGA
AACAGGTGGTCGAGTGATACGCCCAGGGTTC
>Xanthomonas perforans strain GV872 plasmid_contig_54
AAATCCTGGGCGTATCAGACCGCTGCCGAGCGTTTGCTTTGGTACCGCCTGCGTGATCGGCGTCCCCCC/àùbä';ç3cC
¦Å¼OÈ’LäÖ’sëmôzÒÒ${–-ÃêFY;+ã}:ãž0}%ÂZDÅì9D‹Rž÷D×ç°:âñÆó²9h‰ÐžæK¼ë™XÕJK¥·—(ü7¢]ÃÆ
Ño5®o&ÇDçJ^ÃŒZòsÀš‹Ü´†ü¿ÄŒ—xÎ—¸¬µ»?à=Of±)2?îP3ñLapÈëY_{DŸm)ã}ž[Êº„ódMÅbóSl¤ÌŽÑ‚¼
Óh+¶\Éz½/­Z¹â²kgfSû%äÖÿzvŽ#?69éhÎÔøK×‰Êí¨/4ŠSVfÍõÔûÁÒ&½ûE÷kÑ²:k ÷yxÙgH2/ÂžQ,[@ûIß
€|Ë²Ï2,ôCÆ
~~~
{: .output}

> `if` is used for writing a conditional statement.
{: .notes}

If the FASTA contains non-standard characters, and we only want to count and 
filter using standard bases, we can use `gsub` function.

~~~
$ awk -vRS=">" -vORS="\n" -vFS="\n" -vOFS="\t" '
    NR>1 {a=$1; $1 = ""; print a, gsub(/[aAtTgGcCnN]/, $0)}
  ' files/fasta/example.fasta

$ awk -vRS=">" -vORS="" -vFS="\n" -vOFS="\n" '
    NR>1 {a=$0; $1 = ""; if(gsub(/[aAtTgGcCnN]/, $0)>=300) {print ">"a}}
  ' files/fasta/example.fasta
~~~
{: .language-bash}

~~~
Xanthomonas perforans strain GV872 contig_50    212
Xanthomonas perforans strain GV872 contig_51    431
Xanthomonas perforans strain GV872 contig_52    227
Xanthomonas perforans strain GV872 plasmid_contig_53    214
Xanthomonas perforans strain GV872 plasmid_contig_54    75

>Xanthomonas perforans strain GV872 contig_51
CCCTCCACAGGAAGCGGACCGGACAATGACGCTTGATTGCCGCGCAGGAAATGCCCCCGTTCCCGCAGGCATACACCCCA
TTTGCCGCGGACATACGCCGCACATTTCTCGGAGCCTGGACAAGTCGCCACAGCTTTGATGGAACGACGGAAAAAACGCC
CCGGAGTTGACCAGAACAACCGCCAGTACTGAAATCTGGACCTTAAGTTCCAGTTCCTGCAGCTCTGGACCTGGTTCATC
AGCGAGCTTCTATCGCTTCGCGGCGGCTTCGATCGCAGCGTGTGCGCGCCATCACTTGGTCGTGCGAAATCGTCGGCCGC
GTGTCGGCCACGGCGCGCTCGATCTTCGCCCGGAACTACCAGTCGTAGTCCTCGGCGGCTTCGGTCGACTCAAATTCGGA
AACAGGTGGTCGAGTGATACGCCCAGGGTTC
~~~
{: .output}

> Find the list of all IUPAC nucleotide and amino acid codes [here](https://en.wikipedia.org/wiki/International_Union_of_Pure_and_Applied_Chemistry#Amino_acid_and_nucleotide_base_codes).
{: .tips}

> ## Exercise title
> 
> Write a `awk` code to extract first 2 FASTA records from a multiFASTA file and
> write it to file 
> 
> <details markdown=1> <!--collapsible region-->
>   <summary></summary>
> 
> ~~~
> # 1
> $ awk -vRS=">" -vORS="" -vFS="\n" -vOFS="\n" '
>    NR>1 && NR <=3 {$1=$1; print ">"$0}
>  ' files/fasta/example.fasta
> ~~~
> {: .language-bash}
> 
> </details>
{: .challenge}

Using what you learned from optional lessons 1 and 2, 
you can now write a script to filter FASTA files and distribute it to others.
