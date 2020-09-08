---
title: "UNIX command line"
teaching: 30
exercises: 15
questions:
- "Key question (FIXME)"
objectives:
- "First learning objective. (FIXME)"
keypoints:
- "First key point. Brief Answer to questions. (FIXME)"
---

Command line and graphical user interaface are different ways of communicating with computer’s operating system. The shell is a program that provides the command line interface and allowing to control the computer using keyboard commands. For bioinformatics tools, limited software have graphical user interface and you will have to use shell. The shell is a powerful method of communicating with the computer that can help you to do your work more efficienty and understanding how to use shell will be transformative for you to apply in bioinformatics. It can be used to connect to remote and cloud computers.

## Terminal Command line

Once you login to HiperGator through SSH, you will start using a bash shell.

```
[<username>@login1 ~]$
```

The `$` prompt shows that the bash shell is ready to accept bash commands.
Before learning some basic commands, there are a few recommendations regarding UNIX systems.

- The ternminal syntax hates spaces between names. For long or complex names, use connectors such “_" or ‘-’ instead of spaces.
- Uppercase is different from lowercase. ‘R’ is not the same as ‘r’ in commands, paths and arguments.
- UNIX system uses `/` for path, unlike windows, which uses `\`.

Tip: We will be working on ‘blue’ storage in HiperGator for our workshop under a group name ‘general_workshop’. Each user has a profile in the ‘general_workshop’ folder. There is also a ‘share’ folder where all the datasets and information for this workshop are stored. Please remember to only copy requested files from shared folder to your user folder and run the analyses in the folder(directory) with your username only.

To go to your personal folder which is named as your username for the gatorlink account provided to you, enter the following command (we will talk about `cd` shortly). Do not forget to replace &lt;username&gt; with the username provided to you.

```sh
$ cd /blue/general_workshop/<username>
```

## Basic Commands

### Displaying current path/location

`pwd` displays your "path" (where you are located in the cluster).

```sh
$ pwd
/home/<username>
```

### Displaying files and folders in current location

`ls` command dipslays the files and folders in the current location. 
Adding argument `-l` to `ls` displays additional details. 

```sh
$ ls

$ ls -l

total 0
```

Tip: Your directory is currently empty. Run `ls` and `ls -l` 
later after we create some files and directories.
The command displays the properties of the files 
`r` means read only permission 
`w` means write permission 
`x` means read execution permission

### Copying files

`cp` is used for copying files. `cp -R` can be used for copying entire directories. `-R` stands for recursive.
Lets copy a demo file and a demo folder from share folder to you current working directory.

```sh
$ cp ../share/file1.txt ./file1.txt

$ cp file1.txt file2/txt

$ cp -R ../share/dir1 ./

$ ls
dir1     file1.txt     file2.txt
```

Tip: `..` stands for parent directory and `.` stands for current directory.

### Moving files

`mv` is used for moving files or directories. Unlike copying, moving deletes the original copy.

```sh
$ mv file2.txt newfile.txt

$ ls
dir1     file1.txt     newfile.txt
```

### Deleting files

`rm` can be used for deleting files (and directories too with `-r` recursive argument)

```sh
$ rm newfile.txt

$ ls
dir1     file1.txt
```

Tip: You can use <kbd>Tab</kbd> to autocomplete paths and filenames.

### Creating directories

`mkdir` creates a new directory called in current path.
Lets create a new directory called `newdir` and 
then use `ls` command to check if the directory was succesfully created.

```sh
$ mkdir newdir

$ ls
dir1     file1.txt     newdir
```

### Changing directories

`cd` changes your current path. 
Lets change out current path to the directory we just created. 
Use `pwd` to check the current path.

```sh
$ mkdir newdir

$ cd newdir

$ pwd
/blue/general_workshop/<username>/newdir
```

Tip: 
cd also accepts some special path shortcuts in linux. 
`cd ..` moves current path to parent directory, 
`cd ~` changes path to home directory and 
`cd /` changes path to root directory.

```sh
$ cd ..

$ pwd
/blue/general_workshop/<username>
```

### Removing directories

`rmdir` removes the specified directory from current path. 
Lets remove the `newdir` we created earlier and 
check if it is removed using `ls`.
```sh
$ ls
dir1    file1.txt   newdir

$ rmdir newdir

$ ls
dir1     file1.txt
```

### File content handling

There are a set of commands to read the contents of a file. 

`cat` reads the entire content of a file and returns to command line prompt. Let's see what `file1.txt` file contains.

```sh
$ cat file1.txt
CM008465.1      4979077 A       C       intergenic_region       T459_26891-T459_26892   gene26890-gene26891     4979077A>C              2       3       0       11
CM008458.1      97095206        A       G       intergenic_region       T459_11976-T459_11977   gene11975-gene11976     97095206A>G             11      6       12      8
...
...
CM008463.1      123496326       T       A       intergenic_region       T459_23389-T459_23390   gene23388-gene23389     123496326T>A            4       12      8       8
CM008457.1      21333612        G       T       intergenic_region       T459_07945-T459_07946   gene7944-gene7945       21333612G>T             31      1       24      3
```

`less` and `more` display small chunks of the file content at a time.
In `less` output, you can scroll up and down the file content using 
<kbd>&uparrow;</kbd> and <kbd>&downarrow;</kbd> keys.
In `more` output, you can scroll down the file using <kbd>enter</kbd> key.
You can return to command line prompt by pressing <kbd>q</kbd> key.
Lets view `file1.txt` again but using `more` or `less`.

```sh
$ less file1.txt

$ more file1.txt
```

`head` and `tail` can be used to read the start and end of the file respectively.
`-n` argument can be used to specify the the number of lines to read (default is 10 lines).
We can now only read the start or the end of `file1.txt`.

```sh
$ head file1.txt
CM008465.1      4979077 A       C       intergenic_region       T459_26891-T459_26892   gene26890-gene26891     4979077A>C              2       3       0       11
CM008458.1      97095206        A       G       intergenic_region       T459_11976-T459_11977   gene11975-gene11976     97095206A>G             11      6       12      8
CM008459.1      72668492        G       A       intergenic_region       T459_14464-T459_14465   gene14463-gene14464     72668492G>A             1       8       3       17
CM008465.1      57718962        G       A       downstream_gene_variant T459_27577      gene27576       *184G>A         14      9       10      5
CM008465.1      225524501       C       T       intergenic_region       T459_28204-T459_28205   gene28203-gene28204     225524501C>T            6       12      4       11
CM008464.1      76483552        T       A       intergenic_region       T459_25248-T459_25249   gene25247-gene25248     76483552T>A             15      12      8       9
CM008464.1      22690967        C       T       intergenic_region       T459_25013-T459_25014   gene25012-gene25013     22690967C>T             8       13      10      9
CM008463.1      130960788       A       G       intergenic_region       T459_23413-T459_23414   gene23412-gene23413     130960788A>G            1       13      3       12
CM008457.1      19966217        A       C       intergenic_region       T459_07919-T459_07920   gene7918-gene7919       19966217A>C             8       2       11      11
CM008466.1      61441376        G       T       intergenic_region       T459_30076-T459_30077   gene30075-gene30076     61441376G>T             5       7       2       14

$ tail -n 3 file1.txt
CM008463.1      16380489        T       A       intergenic_region       T459_22791-T459_22792   gene22790-gene22791     16380489T>A             14      2       20      3
CM008463.1      123496326       T       A       intergenic_region       T459_23389-T459_23390   gene23388-gene23389     123496326T>A            4       12      8       8
CM008457.1      21333612        G       T       intergenic_region       T459_07945-T459_07946   gene7944-gene7945       21333612G>T             31      1       24      3
```

Lines from middle of a file can be extracted usng `sed` as follows. To extract 5th to 7th lines:
```sh
$ sed -n 5,7p file1.txt
CM008465.1      225524501       C       T       intergenic_region       T459_28204-T459_28205   gene28203-gene28204     225524501C>T            6       12      4       11
CM008464.1      76483552        T       A       intergenic_region       T459_25248-T459_25249   gene25247-gene25248     76483552T>A             15      12      8       9
CM008464.1      22690967        C       T       intergenic_region       T459_25013-T459_25014   gene25012-gene25013     22690967C>T             8       13      10      9
```

### File length

`wc` can be used for reading length of the file. 
`wc` accepts argument 
`-l` for number of lines, 
`-c` for number of characters and 
`-w` for number of words.

```$
$ wc -l file1.txt
100 file1.txt

$ wc -c file1.txt
10546 file1.txt

$ wc -w file1.txt
1200 file1.txt
```

### Writing a file

`>` operator can be used to write the output into a file.
Lets write some files.

```sh
$ head -n 5 file1.txt > head.txt

$ tail -n 5 file1.txt > tail.txt

$ sed -n 56,60p file1.txt > middle.txt
```

You can even save the lists of files in current directory into a file.

```sh
$ ls > files.txt

$ cat files.txt
dir1
file1.txt
files.txt
head.txt
middle.txt
tail.txt
```

Tip: `>>` operator is used to append to a file.

### Concatenate/Join files

`cat` command is used for concatenation of multiple files. 
Lets concatenate the new files created in previous step. 
We can verify concatenation by checking the number of lines in concatenated file.

```sh
$ cat head.txt middle.txt tail.txt > concat.txt

$ wc -l concat.txt
15 concat.txt
```

### Print/output to screen

`echo` command can be used to display result to the screen.

```sh
$ echo "hi"
hi
```

### Variables

Variables can be assigned values with `=` operator. Do not use space around `=`.

```sh
$ a="Hello"

$ echo $a
Hello

$ b="World"

$ echo "$a $b!"
Hello World!
```

### Text manipulation

`cut` is used to extract fields of data from a string or a file. 
Some useful arguments that `cut` takes are listed below.
- `-c`: specify vertical character position to extract
- `-f`: specify columns to extract.
- `-d`: specify delimiter to break line into columns (default is Tab).

`file1.txt` contains part of the output of a real sequence analysis. 
The first part (CM00084xx) is the names of the chromosomes, which we need to extract.
Since chromosome name is first 8 characters, `cut -c` can extract chromosome name.   

```sh
$ cut -c 1-8 file1.txt
CM008465
CM008458
...
...
CM008463
CM008457
```

Looks like the columns in `files1.txt` are tab separated. 
Alternatively, we can just extract the first column using `cut`. 

```sh
$ cut -f 1 file1.txt
CM008465.1
CM008458.1
...
...
CM008463.1
CM008457.1
```

We don't want the `.1` part. In this case, 
we can get the chromosome name by separating lines 
into columns using `.` and taking the first column. 

```sh
$ cut -f 1 -d "." file1.txt
CM008465
CM008458
...
...
CM008463
CM008457
```

### Sorting

`sort` is used to sort lines in a file. 
By default, `sort` sorts in alphabetical order.
`-r` argument reverses the order.

```sh
$ sort file1.txt
CM008455.1      27261226        C       A       intergenic_region       T459_01015-T459_01016   gene1014-gene1015       27261226C>A             15      14      24      18
CM008455.1      289184259       AT      A       intron_variant  T459_03640      gene3639        159+2881delT            11      5       7       7
CM008455.1      289401193       G       A       intergenic_region       T459_03642-T459_03643   gene3641-gene3642       289401193G>A            9       6       7       8
...
...

$ sort -r file1.txt
CM008466.1      95058324        T       C       intergenic_region       T459_30200-T459_30201   gene30199-gene30200     95058324T>C             28      4       29      6
CM008466.1      71138968        G       A       intergenic_region       T459_30125-T459_30126   gene30124-gene30125     71138968G>A             6       8       3       9
CM008466.1      61600177        C       T       intergenic_region       T459_30076-T459_30077   gene30075-gene30076     61600177C>T             12      1       5       3
...
...
```

`-n` argument is used to sort in numerical ascending order, and is useful for numeric columns. 
For e.g. 2 comes before 15 in numerical order, but 2 comes after 15 in alphabetical order
as first character "2" comes after first character "1".

`-k` can be used to sort by a specific column.

```sh

$ sort -k2 file1.txt
108474056
109077809
...
...
97095206
9871059

$ sort -k2n file1.txt
3542278
3999373
...
...
289401193
294972840

$ sort -k2nr file1.txt
294972840
289401193
...
...
3999373
3542278
```

Tip: `-k2nr` is shortcut for `-k 2 -n -r`.

### Replacing text

`sed` command can be used for replacing text in a file. 
`sed 's/old/new/g' file1` replaces all instances of `old` 
with `new` in the file `file1`.

```sh
$ sed 's/CM0084/Chr_/g' file1.txt > prettyfile.txt

$ head -n 3 prettyfile.txt
Chr_65.1        4979077 A       C       intergenic_region       T459_26891-T459_26892   gene26890-gene26891     4979077A>C              2       3       0       11
Chr_58.1        97095206        A       G       intergenic_region       T459_11976-T459_11977   gene11975-gene11976     97095206A>G             11      6       12      8
Chr_59.1        72668492        G       A       intergenic_region       T459_14464-T459_14465   gene14463-gene14464     72668492G>A             1       8       3       17
```

### Search and Extract

`grep` command is used to find a string in a file and return the matching line.
Argument `-c` is used to return the number of matches.

```sh
$ grep "downstream_gene_variant" file1.txt
CM008465.1      57718962        G       A       𝐝𝐨𝐰𝐧𝐬𝐭𝐫𝐞𝐚𝐦_𝐠𝐞𝐧𝐞_𝐯𝐚𝐫𝐢𝐚𝐧𝐭 T459_27577      gene27576       *184G>A         14      9       10      5
CM008466.1      236800990       C       A       𝐝𝐨𝐰𝐧𝐬𝐭𝐫𝐞𝐚𝐦_𝐠𝐞𝐧𝐞_𝐯𝐚𝐫𝐢𝐚𝐧𝐭 T459_31014      gene31013       *1241C>A                11      9       4       15
CM008458.1      226074184       T       TA      𝐝𝐨𝐰𝐧𝐬𝐭𝐫𝐞𝐚𝐦_𝐠𝐞𝐧𝐞_𝐯𝐚𝐫𝐢𝐚𝐧𝐭 T459_12904      gene12903       *4901_*4902iT           5       9       7       14
CM008458.1      214028749       A       C       𝐝𝐨𝐰𝐧𝐬𝐭𝐫𝐞𝐚𝐦_𝐠𝐞𝐧𝐞_𝐯𝐚𝐫𝐢𝐚𝐧𝐭 T459_12654      gene12653       *2350T>G                5       5       9       8
CM008460.1      5909607 A       T               𝐝𝐨𝐰𝐧𝐬𝐭𝐫𝐞𝐚𝐦_𝐠𝐞𝐧𝐞_𝐯𝐚𝐫𝐢𝐚𝐧𝐭 T459_15880      gene15879       *2684A>T                16      6       20      5
CM008460.1      63781473        A       T       𝐝𝐨𝐰𝐧𝐬𝐭𝐫𝐞𝐚𝐦_𝐠𝐞𝐧𝐞_𝐯𝐚𝐫𝐢𝐚𝐧𝐭 T459_16347      gene16346       *2958A>T                1       20      1       13
CM008465.1      138360632       A       G       𝐝𝐨𝐰𝐧𝐬𝐭𝐫𝐞𝐚𝐦_𝐠𝐞𝐧𝐞_𝐯𝐚𝐫𝐢𝐚𝐧𝐭 T459_27891      gene27890       *3582T>C                6       15      6       16

$ grep -c "downstream_gene_variant" file1.txt
7
```

### Piping commands

`|` operators can be used for piping. 
Piping means that the output of first command serves as input of second command and so on. 
This eliminates need for saving intermediate results as a file.

Lets sort `file1` by chromosome (first column) and position (second column), 
and only display the first 5 columns in top of the sorted lines.

```sh
$ sort -k1,1 -k2,2n output.txt | cut -f1-5 | head -n5
CM008455.1      27261226        C       A       intergenic_region
CM008455.1      38680477        T       G       intergenic_region
CM008455.1      39264338        G       A       intergenic_region
CM008455.1      75147386        C       A       intergenic_region
CM008455.1      82264380        ACC     A       intergenic_region
```

## Exercise: Finding “alien genes” in the plant pathogen *Streptomyces scabies*.

Let’s use the commands in a simple real case. Copy a file in shared folder ‘aliens_in_scabies’ using following command

warning: make sure you are in /blue/general_workshop/&lt;username&gt; folder, if not use 'cd /blue/general_workshop/username' to go to your personal directory. 

$ cp ../share/strep/aliens_in_scabies ./

*Streptomyces scabies* is a plant pathogen that produces necrosis in potatoes. Most of the virulence factors are locate regions that have low GC content. Also, virulence is highly expressed during the interaction with roots. The file “aliens_in_scabies” contains a tabular table with more than 1000 genes (first column) with the level of change expression when comparing growth in rich medium vs. interaction with roots (second column). Besides, the gene sequences, GC% content is provided in the third column of the table.

Using UNIX commands 
1. Sort the table by expression levels (second column). 
2. Create a new file that contains the top 10 highest expressed genes. 
3. Using the new file, create a table to sort the table by gene GC content. 
4. Create a table with only the gene names but replace the name SCAB with SCABIES.

<details markdown=1>
  <summary> Click here for answer.</summary>

```sh
# 1
$ sort -k2n aliens_in_scabies

# 2
$ sort -k2nr aliens_in_scabies | head > newfile.txt

# 3
$ sort -k3n newfile.txt

# 4
$ cut -f1 aliens_in_scabies | sed 's/SCAB/SCABIES/g'
```

</details>

## An intro to loops
Many tasks are repetitive. It is not necessary to repeat the same command multiple times. We learn how to pipe and re-direct outputs of the commands.
let’s learn how to use loops to repeat a task

### Wildcards

There are certain symbols in UNIX that represent mutiple values, 
and are called “wildcards”. `*` is a universal wildcard that stands for anything.

```sh
$ ls *.fasta
...
...

$ cat *.fasta
...
...
```

Other wildcards are 

| Wildcard | Value |
|----------|-------|
| * | any character or characters |
| ? | any single character |
| [0-9] | any single number |
| [A-Z] | any single capital alphabets |
| [a-z] | aany single small alphabets |
| [!x] | not x |
| {x, y} | x or y |



### The ‘FOR’ loop

`For` loop iterates over a specific range or numbers or an array. For example.

```sh
$ for x in {0..3}
$ do
$   echo $x
$ done
0
1
2
3
```

The process between DO and DONE enclose task to be repeated.
The variable x that will take consecutive values from array in each iteration.

Tip: `echo` command is equivalent to print to screen.


### The ‘WHILE’ loop

`While` loop iterates as long as a condition is true. Lets use the example.

```sh
$ x=1
$ while(($x <5))
$ do 
$   echo $x
$   ((++x)) 
$ done
1
2
3
4
```

The process between DO and DONE enclose task to be repeated. 

Tips: 
More about Unix commands
https://swcarpentry.github.io/shell-novice/reference/
