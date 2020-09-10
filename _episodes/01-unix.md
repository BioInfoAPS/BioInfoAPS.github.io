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

~~~
[<username>@login1 ~]$
~~~
{: .terminal}

The `$` prompt shows that the bash shell is ready to accept bash commands.
Before learning some basic commands, there are a few recommendations regarding UNIX systems.

- The ternminal syntax hates spaces between names. 
For long or complex names, use connectors such “_" or ‘-’ instead of spaces.
- Uppercase is different from lowercase. 
‘R’ is not the same as ‘r’ in commands, paths and arguments.
- UNIX system uses `/` for path, unlike windows, which uses `\`.

We will be working on ‘blue’ storage in HiperGator 
for our workshop under a group name ‘general_workshop’. 
Each user has a directory in the ‘general_workshop’ folder. 
There is also a ‘share’ directory where all the datasets and 
information for this workshop are stored. 
Please remember to only copy requested files from 
shared folder to your user folder and run the analyses 
in the folder(directory) with your username only.

Your personal folder is named as your gatorlink username. 
Enter the following command to go to your work directory 
(we will talk about `cd` shortly). 
Do not forget to replace &lt;username&gt; with the username provided to you.

~~~
$ cd /blue/general_workshop/<username>
~~~
{: .language-bash}

> Do not copy `$` sign.
{: .caution}

## Basic Commands

### Displaying current path/location

`pwd` displays your "path" (where you are located in the cluster).

~~~
$ pwd
~~~
{: .language-bash}

~~~
/blue/general_workshop/<username>
~~~
{: .output}

### Displaying files and folders in current location

`ls` command dipslays the files and folders in the current location.

~~~
$ ls
~~~
{: .language-bash}

~~~

~~~
{: .output}

> Your directory is currently empty. 
> However, we can check somewhere else.
{: .notes}

~~~
$ ls /blue/general_workshop/
~~~
{: .language-bash}

~~~
anujsharma     Intro_slides.pptx     share     sujan.timilsina
...
...
~~~
{: .output}

Adding argument `-l` to `ls` displays additional details such as permissions, file owner, 
size, date modified etc. 

~~~
$ ls -l /blue/general_workshop/
~~~
{: .language-bash}

~~~
drwxr-sr-x 3 anujsharma      general_workshop    4096 Sep 10 03:47 anujsharma
-rw-r----- 1 jhuguet         general_workshop 1290989 Sep  8 15:11 Intro_slides.pptx
drwxr-sr-x 7 sujan.timilsina general_workshop    4096 Sep  8 15:27 share
drwxr-sr-x 2 sujan.timilsina general_workshop    4096 Sep  8 11:01 sujan.timilsina
~~~
{: .output}

> ## Permissions in linux
>
> Linux permissions look like this: `d  r w x  r w −  r − −`
>
> - 1st character represents special flag: file `−`, directory `d` or link `l`. 
> - The rest of the characters show permissions in set of three: `r` for read, `w` for write and `x` for execute. `−` means permission denied. 
>   - 2nd to 4th characters: permission for file owner
>   - 5th to 7th characters: permission for the group
>   - 8th to 10th characters: permission for others
{: .tips}

### Copying files

`cp` is used for copying files. `cp -R` can be used for copying entire directories. `-R` stands for recursive.
Lets copy a demo file and a demo folder from share folder to you current working directory.

~~~
$ cp ../share/file1.txt ./file1.txt

$ ls
~~~
{: .language-bash}

~~~
file1.txt
~~~
{: .output}

~~~
$ cp file1.txt file2.txt

$ ls
~~~
{: .language-bash}

~~~
file1.txt     file2.txt
~~~
{: .output}

~~~
$ cp -R ../share/demo ./

$ ls
~~~
{: .language-bash}

~~~
demo     file1.txt     file2.txt
~~~
{: .output}

> ## Common path symbols in linux
> Linux uses some symbols to represent commonly used paths.
> - `..` stands for parent directory.
> - `.` stands for current directory.
> - `/` at the beginning stands for root directory.
> - `~` stands for home directory.
{: .tips}

### Moving files

`mv` is used for moving files or directories. Unlike copying, moving deletes the original copy.

~~~
$ mv file2.txt newfile.txt

$ ls
~~~
{: .language-bash}

~~~
demo     file1.txt     newfile.txt
~~~
{: .output}

### Deleting files

`rm` can be used for deleting files (and directories too with `-r` recursive argument)

~~~
$ rm newfile.txt

$ ls
~~~
{: .language-bash}

~~~
demo     file1.txt
~~~
{: .output}

### Creating directories

`mkdir` creates a new directory called in current path.
Lets create a new directory called `newdir` and 
then use `ls` command to check if the directory was succesfully created.

~~~
$ mkdir newdir

$ ls
~~~
{: .language-bash}

~~~
demo     file1.txt     newdir
~~~
{: .output}

### Changing directories

`cd` changes your current path. 
Lets change out current path to the directory we just created. 
Use `pwd` to check the current path.

~~~
$ cd newdir

$ pwd
~~~
{: .language-bash}

~~~
/blue/general_workshop/<username>/newdir
~~~
{: .output}

~~~
$ cd ..

$ pwd
~~~
{: .language-bash}

~~~
/blue/general_workshop/<username>
~~~
{: .output}

### Removing directories

`rmdir` removes the specified directory from current path. 
Lets remove the `newdir` we created earlier and 
check if it is removed using `ls`.
~~~
$ rmdir newdir

$ ls
~~~
{: .language-bash}

~~~
demo     file1.txt
~~~
{: .output}

### File content handling

There are a set of commands to read the contents of a file. 

`cat` reads the entire content of a file and returns to command line prompt. Let's see what `file1.txt` file contains.

~~~
$ cat file1.txt
~~~
{: .language-bash}

~~~
CM008465.1      4979077 A       C       intergenic_region       gene26890-gene26891     4979077A>C
CM008458.1      97095206        A       G       intergenic_region       gene11975-gene11976     97095206A>G
CM008459.1      72668492        G       A       intergenic_region       gene14463-gene14464     72668492G>A
...
...
CM008463.1      16380489        T       A       intergenic_region       gene22790-gene22791     16380489T>A
CM008463.1      123496326       T       A       intergenic_region       gene23388-gene23389     123496326T>A
CM008457.1      21333612        G       T       intergenic_region       gene7944-gene7945       21333612G>T
~~~
{: .output}

`less` and `more` display small chunks of the file content at a time.
In `less` output, you can scroll up and down the file content using 
<kbd>&uparrow;</kbd> and <kbd>&downarrow;</kbd> keys.
In `more` output, you can scroll down the file using <kbd>enter</kbd> key.
You can return to command line prompt by pressing <kbd>q</kbd> key.
Lets view `file1.txt` again but using `more` or `less`.

~~~
$ less file1.txt

$ more file1.txt
~~~
{: .language-bash}

`head` and `tail` can be used to read the start and end of the file respectively.
`-n` argument can be used to specify the the number of lines to read (default is 10 lines).
We can now only read the start or the end of `file1.txt`.

~~~
$ head file1.txt
~~~
{: .language-bash}

~~~
CM008465.1      4979077 A       C       intergenic_region       gene26890-gene26891     4979077A>C
CM008458.1      97095206        A       G       intergenic_region       gene11975-gene11976     97095206A>G
CM008459.1      72668492        G       A       intergenic_region       gene14463-gene14464     72668492G>A
CM008465.1      57718962        G       A       downstream_gene_variant gene27576       *184G>A
CM008465.1      225524501       C       T       intergenic_region       gene28203-gene28204     225524501C>T
CM008464.1      76483552        T       A       intergenic_region       gene25247-gene25248     76483552T>A
CM008464.1      22690967        C       T       intergenic_region       gene25012-gene25013     22690967C>T
CM008463.1      130960788       A       G       intergenic_region       gene23412-gene23413     130960788A>G
CM008457.1      19966217        A       C       intergenic_region       gene7918-gene7919       19966217A>C
CM008466.1      61441376        G       T       intergenic_region       gene30075-gene30076     61441376G>T
~~~
{: .output}

~~~
$ tail -n 3 file1.txt
~~~
{: .language-bash}

~~~
CM008463.1      16380489        T       A       intergenic_region       gene22790-gene22791     16380489T>A
CM008463.1      123496326       T       A       intergenic_region       gene23388-gene23389     123496326T>A
CM008457.1      21333612        G       T       intergenic_region       gene7944-gene7945       21333612G>T
~~~
{: .output}

Lines from middle of a file can be extracted usng `sed -n` as follows. To extract 5th to 7th lines:
~~~
$ sed -n 5,7p file1.txt
~~~
{: .language-bash}

~~~
CM008465.1      225524501       C       T       intergenic_region       gene28203-gene28204     225524501C>T
CM008464.1      76483552        T       A       intergenic_region       gene25247-gene25248     76483552T>A
CM008464.1      22690967        C       T       intergenic_region       gene25012-gene25013     22690967C>T
~~~
{: .output}

> ## `sed`
> `sed` is a very powerful tool in in bash and 
> can be used to do wide range of text editing tasks.
> Here, `-n` argument directs `sed` to pick only the lines that 
> match the parameter `5,7p`.
{: .tips}

### File length

`wc` can be used for reading length of the file. 
`wc` accepts argument 
`-l` for number of lines, 
`-c` for number of characters and 
`-w` for number of words.

~~~
$ wc -l file1.txt
~~~
{: .language-bash}

~~~
100 file1.txt
~~~
{: .output}

~~~
$ wc -c file1.txt
~~~
{: .language-bash}

~~~
7502 file1.txt
~~~
{: .output}

~~~
$ wc -w file1.txt
~~~
{: .language-bash}

~~~
700 file1.txt
~~~
{: .output}

### Writing a file

`>` operator can be used to write the output into a file.
Lets write some files.

~~~
$ head -n 5 file1.txt > head.txt

$ tail -n 5 file1.txt > tail.txt

$ sed -n 56,60p file1.txt > middle.txt

$ ls
~~~
{: .language-bash}

~~~
demo     file1.txt     head.txt     middle.txt     tail.txt
~~~
{: .output}

You can even save the lists of files in current directory into a file.

~~~
$ ls > files.txt

$ cat files.txt
~~~
{: .language-bash}

~~~
demo
file1.txt
files.txt
head.txt
middle.txt
tail.txt
~~~
{: .output}

> ## `>` vs `>>`
> Writing to an existing file with `>` removes its existing contents. Use `>>` operator to append new contents to the end of exsiting file.
{: .tips}

### Concatenate/Join files

`cat` command is used for concatenation of multiple files. 
Lets concatenate the new files created in previous step. 
We can verify concatenation by checking the number of lines in concatenated file.

~~~
$ cat head.txt middle.txt tail.txt > concat.txt

$ wc -l concat.txt
~~~
{: .language-bash}

~~~
15 concat.txt
~~~
{: .output}

> `head.txt`. `middle.txt` and `tail.txt` each had 5 lines each. To confirm, you can check like this: `wc -l head.txt`.
{: .notes}

### Print/output to screen

`echo` command can be used to display result to the screen.

~~~
$ echo "hi"
~~~
{: .language-bash}

~~~
hi
~~~
{: .output}

### Variables

Variables can be assigned values with `=` operator. Do not use space around `=`.

~~~
$ a="Hello"

$ echo $a
~~~
{: .language-bash}

~~~
Hello
~~~
{: .output}

~~~
$ b="World"

$ echo "$a $b"
~~~
{: .language-bash}

~~~
Hello World
~~~
{: .output}

> ## `"` vs `'` (Double quotes vs single quotes)
> `"` and `'` mean different things in unix and should not be use interchangeably.
> Check yourself how the output of `echo '$a $b'` differs from that of `echo "$a $b"`.
{: .tips}

### Text manipulation

`cut` is used to extract fields of data from a string or a file.

`head.txt` contains part of the output of a real sequence analysis. 
The first part (CM00084xx) contains the names of the chromosomes, 
which we need to extract.

Since chromosome name is first 8 characters, 
`cut -c` to extract chromosome name. 
`-c` specifies vertical position of characters to extract  

~~~
$ cut -c 1-8 head.txt
~~~
{: .language-bash}

~~~
CM008465
CM008458
CM008459
CM008465
CM008465
~~~
{: .output}

The columns in `head.txt` are tab separated. You can verify with `cat head.txt`.
So, we can alternatively just extract the first **column** using `cut`.
`-f`: specify columns to extract.

~~~
$ cut -f 1 head.txt
~~~
{: .language-bash}

~~~
CM008465.1
CM008458.1
CM008459.1
CM008465.1
CM008465.1
~~~
{: .output}

We don't want the `.1` part. In this case, 
we can get the chromosome name by separating columns 
by `.` instead of space. `-d` arguemnt allows us to 
specify delimiter to break line into columns (default is Tab).

~~~
$ cut -f 1 -d "." head.txt
~~~
{: .language-bash}

~~~
CM008465
CM008458
CM008459
CM008465
CM008465
~~~
{: .output}

### Sorting

`sort` is used to sort lines in a file. 
By default, `sort` sorts in alphabetical order.

~~~
$ sort file1.txt
~~~
{: .language-bash}

~~~
CM0084𝟓𝟓.1      𝟐𝟕261226        C       A       intergenic_region       gene1014-gene1015       27261226C>A
CM0084𝟓𝟓.1      𝟐𝟖𝟗𝟏84259       AT      A       intron_variant  gene3639        159+2881delT
CM0084𝟓𝟓.1      𝟐𝟖𝟗𝟒01193       G       A       intergenic_region       gene3641-gene3642       289401193G>A
...
...
CM0084𝟲𝟲.1      𝟔1600177        C       T       intergenic_region       gene30075-gene30076     61600177C>T
CM0084𝟲𝟲.1      𝟕1138968        G       A       intergenic_region       gene30124-gene30125     71138968G>A
CM0084𝟲𝟲.1      𝟗5058324        T       C       intergenic_region       gene30199-gene30200     95058324T>C
~~~
{: .output}


`-r` argument reverses the order od sorting.

~~~
$ sort -r file1.txt
~~~
{: .language-bash}

~~~
CM0084𝟲𝟲.1      95058324        T       C       intergenic_region       gene30199-gene30200     95058324T>C
CM0084𝟲𝟲.1      71138968        G       A       intergenic_region       gene30124-gene30125     71138968G>A
CM0084𝟲𝟲.1      61600177        C       T       intergenic_region       gene30075-gene30076     61600177C>T
...
...
CM0084𝟓𝟓.1      289401193       G       A       intergenic_region       gene3641-gene3642       289401193G>A
CM0084𝟓𝟓.1      289184259       AT      A       intron_variant  gene3639        159+2881delT
CM0084𝟓𝟓.1      27261226        C       A       intergenic_region       gene1014-gene1015       27261226C>A
~~~
{: .output}

`-n` argument is used to sort in numerical ascending order.

> ## Numerical vs alphabetical order
> Sorting by numerical order is useful when dealing with numbers.
> In alphabetical order, 2 comes after 15 because the first character 
> "2" comes after first character "1".
> In numerical order, 2 comes before 15, since 2 is smaller than 15.
{: .tips}

`-k` can be used to sort by a specific column.

~~~
$ sort -k2 file1.txt
~~~
{: .language-bash}

~~~
CM008465.1      𝟭08474056       G       A       intergenic_region       gene27823-gene27824     108474056G>A
CM008463.1      𝟭09077809       C       T       intergenic_region       gene23339-gene23340     109077809C>T
...
...
CM008458.1      𝟵7095206        A       G       intergenic_region       gene11975-gene11976     97095206A>G
CM008461.1      𝟵871059 T       C       intergenic_region       gene18895-gene18896     9871059T>C
~~~
{: .output}

~~~
$ sort -k2n file1.txt
~~~
{: .language-bash}

~~~
CM008457.1      3542278 T       G       intergenic_region       gene7359-gene7360       3542278T>G
CM008465.1      3999373 G       A       upstream_gene_variant   gene26864       -4845C>T
...
...
CM008455.1      289401193       G       A       intergenic_region       gene3641-gene3642       289401193G>A
CM008455.1      294972840       C       T       intergenic_region       gene3691-gene3692       294972840C>T
~~~
{: .output}

~~~
$ sort -k2nr file1.txt
~~~
{: .language-bash}

~~~
CM008455.1      294972840       C       T       intergenic_region       gene3691-gene3692       294972840C>T
CM008455.1      289401193       G       A       intergenic_region       gene3641-gene3642       289401193G>A
...
...
CM008465.1      3999373 G       A       upstream_gene_variant   gene26864       -4845C>T
CM008457.1      3542278 T       G       intergenic_region       gene7359-gene7360       3542278T>G
~~~
{: .output}

> ## Argument shortcuts in bash
> Multiple argument, and sometimes expected values of arguments
> can be written together in bash.
> In above example, `-k2nr` is shortcut for `-k 2 -n -r`. 
> Guess what `sort -k2nr file1.txt` does.
{: .tips}

### Replacing text

`sed` command can be used for replacing text in a file. 
`sed 's/old/new/g'` replaces all instances of `old` with `new`.

~~~
$ sed 's/CM0084/Chr_/g' tail.txt > prettyfile.txt

$ cat prettyfile.txt
~~~
{: .language-bash}

~~~
Chr_66.1        28255843        A       C       intergenic_region       gene29731-gene29732     28255843A>C
Chr_65.1        236852266       C       T       intergenic_region       gene28293-gene28294     236852266C>T
Chr_63.1        16380489        T       A       intergenic_region       gene22790-gene22791     16380489T>A
Chr_63.1        123496326       T       A       intergenic_region       gene23388-gene23389     123496326T>A
Chr_57.1        21333612        G       T       intergenic_region       gene7944-gene7945       21333612G>T
~~~
{: .output}

### Search and Extract

`grep` command is used to find a string in a file and return the matching line.

~~~
$ grep "downstream_gene_variant" file1.txt
~~~
{: .language-bash}

~~~
CM008465.1      57718962        G       A       𝐝𝐨𝐰𝐧𝐬𝐭𝐫𝐞𝐚𝐦_𝐠𝐞𝐧𝐞_𝐯𝐚𝐫𝐢𝐚𝐧𝐭 gene27576       *184G>A 
CM008466.1      236800990       C       A       𝐝𝐨𝐰𝐧𝐬𝐭𝐫𝐞𝐚𝐦_𝐠𝐞𝐧𝐞_𝐯𝐚𝐫𝐢𝐚𝐧𝐭 gene31013       *1241C>A
CM008458.1      226074184       T       TA      𝐝𝐨𝐰𝐧𝐬𝐭𝐫𝐞𝐚𝐦_𝐠𝐞𝐧𝐞_𝐯𝐚𝐫𝐢𝐚𝐧𝐭 gene12903       *4901_*4902iT
CM008458.1      214028749       A       C       𝐝𝐨𝐰𝐧𝐬𝐭𝐫𝐞𝐚𝐦_𝐠𝐞𝐧𝐞_𝐯𝐚𝐫𝐢𝐚𝐧𝐭 gene12653       *2350T>G 
CM008460.1      5909607 A       T       𝐝𝐨𝐰𝐧𝐬𝐭𝐫𝐞𝐚𝐦_𝐠𝐞𝐧𝐞_𝐯𝐚𝐫𝐢𝐚𝐧𝐭 gene15879       *2684A>T
CM008460.1      63781473        A       T       𝐝𝐨𝐰𝐧𝐬𝐭𝐫𝐞𝐚𝐦_𝐠𝐞𝐧𝐞_𝐯𝐚𝐫𝐢𝐚𝐧𝐭 gene16346       *2958A>T
CM008465.1      138360632       A       G       𝐝𝐨𝐰𝐧𝐬𝐭𝐫𝐞𝐚𝐦_𝐠𝐞𝐧𝐞_𝐯𝐚𝐫𝐢𝐚𝐧𝐭 gene27890       *3582T>C 
~~~
{: .output}

Argument `-c` is used to return the number of matches.

~~~
$ grep -c "downstream_gene_variant" file1.txt
~~~
{: .language-bash}

~~~
7
~~~
{: .output}

### Piping commands

`|` operators can be used for piping. 
Piping means that the output of first command serves as input of second command and so on. 
This eliminates need for saving intermediate results as a file.

Lets sort `file1.txt` by chromosome (first column) and 
by position (second column) within a chromosome, 
and only display the first 5 columns in top of the sorted lines.

~~~
$ sort -k1,1 -k2,2n file1.txt | cut -f1-5 | head -n5
~~~
{: .language-bash}

~~~
CM008455.1      27261226        C       A       intergenic_region
CM008455.1      38680477        T       G       intergenic_region
CM008455.1      39264338        G       A       intergenic_region
CM008455.1      75147386        C       A       intergenic_region
CM008455.1      82264380        ACC     A       intergenic_region
~~~
{: .output}

> `,1` and `,2` specify the order/hierarchy of sorting.
{: .tips}

> ## Exercise: Finding “alien genes” in the plant pathogen *Streptomyces scabies*.
> 
> Let’s use the commands in a simple real case. Copy a file `aliens_in_scabies` inside `share` directory using following command
> 
> > Make sure you are in `/blue/general_workshop/<username>` directory, 
> > if not use `cd /blue/general_workshop/<username>` to go to your personal working directory. 
> {: .caution}
> 
> ~~~
> $ cp ../share/strep/aliens_in_scabies ./
> ~~~
> {: .language-bash}
> 
> *Streptomyces scabies* is a plant pathogen that produces necrosis in potatoes. 
> Most of the virulence factors are located in regions that have low GC content. 
> Also, virulence is highly expressed during the interaction with roots. 
> The file “aliens_in_scabies” contains a tabular data with more than 1000 genes 
> (first column) along with the level of change expression when growth in rich medium 
> vs. interaction with roots (second column). 
> GC content% is provided in the third column of the table.
> 
> Using UNIX commands 
> 1. Sort the table by expression levels (second column).  <input type="checkbox">
> 2. Create a new file that contains the top 10 highest expressed genes.  <input type="checkbox">
> 3. Using the new file, create a table to sort the table by gene GC content.  <input type="checkbox">
> 4. Create a table with only the gene names but replace the name SCAB with SCABIES. <input type="checkbox">
> 
> <details markdown=1>
>   <summary></summary>
> 
> ~~~
> # 1
> $ sort -k2n aliens_in_scabies
> 
> # 2
> $ sort -k2nr aliens_in_scabies | head > newfile.txt
> 
> # 3
> $ sort -k3n newfile.txt
> 
> # 4
> $ cut -f1 aliens_in_scabies | sed 's/SCAB/SCABIES/g'
> ~~~
> {: .language-bash}
> 
> </details>
{: .challenge}

## An intro to loops
Many tasks are repetitive. It is not necessary to repeat the same command multiple times. Instead, we can use wildcards and loops to repeat a task.

### Wildcards

There are certain symbols in UNIX that represent mutiple values, 
and are called “wildcards”. `*` is a universal wildcard that stands for anything.

~~~
$ ls *.txt
~~~
{: .language-bash}

~~~
concat.txt     files.txt     middle.txt         tail.txt
file1.txt      head.txt      prettyfile.txt
~~~
{: .output}

~~~
$ mv *.txt demo

$ ls
~~~
{: .language-bash}

~~~
demo
~~~
{: .output}

~~~
$ ls demo
~~~
{: .language-bash}

~~~
concat.txt     files.txt     middle.txt         tail.txt
file1.txt      head.txt      prettyfile.txt
~~~
{: .output}

Other wildcards are 

| Wildcard | Value |
|----------|-------|
| `*` | any character or characters |
| `?` | any single character |
| `[0-9]` | any single number |
| `[A-Z]` | any single capital alphabets |
| `[a-z]` | aany single small alphabets |
| `[!x]` | not x |
| `{x, y}` | x or y |



### The ‘FOR’ loop

`For` loop iterates over a specific range or numbers or an array. For example.

~~~
$ for x in {0..3}
> do
>   echo $x
> done
~~~
{: .language-bash}

~~~
0
1
2
3
~~~
{: .output}

> Be careful when copying code. 
> If you are copying the whole code block, 
> first paste it to a text editor, 
> remove staring `$` or `>` sign,
> copy edited command and only then, 
> paste it into terminal or powershell.
>
> Or click [here](/code/demo_for.txt){: target="_blank"} to open a separate page to copy the code above.
{: .caution}

The process between `do` and `done` enclose task to be repeated.
The variable `x` will take consecutive values from array `{0..3}` in each iteration.

> Try `echo {0..3}` to see what `{0..3}` represents.
{: .tips}

> Bash prompts
> `$` prompt specifies ready to start a commaond.
> `#` replaces `$` for root user.
> `>` specifies continuation of multiline command from previous line.
{: .tips}


### The ‘WHILE’ loop

`while` loop iterates as long as a condition is true. Lets use the example.

~~~
$ x=0
$ while(($x <4))
> do 
>   echo $x
>   ((++x)) 
> done
~~~
{: .language-bash}

~~~
0
1
2
3
~~~
{: .output}

As in `for` loop, the commands to be looped over are enclosed between `do` and `done`. 

> Click [here](/code/demo_while.txt){: target="_blank"} to open a separate page to copy the code above.
{: .caution}


> ## Learn more about unix commands
> [Software carpentary reference](https://swcarpentry.github.io/shell-novice/reference/)
{: .notes}
