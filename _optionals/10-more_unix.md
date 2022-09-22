---
title: "More UNIX commands"
objectives:
- To understand UNIX conditional statements.
- To learn about functions in bash.
- To get started with bash scripting.
- Parallelization in bash with GNU parallel.

---

The optional sections are targeted for more advanced users 
(psa: more advanced than regular sections, but still for beginners). 
Feel free to try these sections once you are done with regular lesson.

## Conditionals

A conditional expression consists of commands that are executed only if some
conditions are satisfied. 
In bash, the general syntax for conditional statement looks like:

~~~
$ if [ condition ]
  then
      commands if given condition is true ...
      ...
  elif [ another condition]
      commands if new condition is true but previous conditions are false ...
      ...
  else
      commands if all conditions is false ...
      ...
  fi
~~~
{: .terminal}

> `elif` and `else` sections are optional. 
> `elif` section can be repeated multiple times.
{: .notes}

Let's look at a simple example:

~~~
$ if [ 1<2 ]
  then
      echo "true"
  else
      echo "false"
  fi
~~~
{: .language-bash}

~~~
true
~~~
{: .output}

~~~
$ if [ 1>2 ]
  then
      echo "true"
  else
      echo "false"
  fi
~~~
{: .language-bash}

~~~
false
~~~
{: .output}

### Operators

Instead of commonly used operator symbols such as `=` or `>=`, 
it is preferable to use test operators for numeric comparisons, 
especially when working with variables that have not been declared as numbers.
- `-eq` for **equal to**, 
- `-gt` for **greater than**
- `-lt` for **greater than**, 
- `-ge` for **greater than or equal to**
- `-le` for **less than or equal to**

~~~
$ a=1

$ b=2

$ if [ $a -gt $b ]
  then
      echo "a is greater than b."
  elif [ $a -lt $b ]
      echo "a is lesser than b."
  else
      echo "a is equal to b."
  fi
~~~
{: .language-bash}

~~~
a is lesser than b.
~~~
{: .output}

Other useful operators are
- `!condition` for checking if the condition is false.
- `-z STRING` for checking if string (or variable) is empty.
- `-v VARIABLE` for checking if variable is defined.
- `-d /path/to/DIR` for checking if a directory exists.
- `-f /path/to/FILE` for checking if a file exists.
- `-s /path/to/FILE` for checking if a file exists and is not empty.

### Single line conditional statement

For shorter commands, conditional statements can be written as one-liner
using double brackets `[[...]]`.

The general syntax is:
~~~
$ [[ CONDITION ]] && COMMANDS IF TRUE || COMMANDS IF FALSE
~~~
{: .language-bash}

Let's look at some examples.

~~~
$ ls files/fasta
~~~
{: .language-bash}

~~~
example.fasta
~~~
{: .output}

~~~
$ [[ -f files/demo/example.txt ]] && echo "file exists".
~~~
{: .language-bash}

~~~
~~~
{: .output}

~~~
$ [[ -f files/demo/example.txt ]] && echo "file exists" || echo "file does not exist"
~~~
{: .language-bash}

~~~
file does not exist
~~~
{: .output}

~~~
$ [[ -f files/demo/example.fasta ]] && echo "file exists" || echo "file does not exist"
~~~
{: .language-bash}

~~~
file exists
~~~
{: .output}

~~~
$ [[ ! -f files/demo/example.fasta ]] && echo "file exists" || echo "file does not exist"
~~~
{: .language-bash}

~~~
file does not exist
~~~
{: .output}

## An intro to functions

Functions are scoped, reusable and often portable chunks of code.
Reusability is an important aspect to a function, 
as you can refer back to the same piece of code multiple times throughout
your pipeline once you define it once.

### Function syntax

In bash, functions are defined as follows:

~~~
$ function_name () 
  {
      function code ...
      ...
  }
~~~
{: .terminal}

> `function_name` is the name of the function, and it will be used to call the function later.
>
> `()` declares the variable `function_name` stores a function, 
> rather than string or another data.
>
> `{...}` delimits the scope of the function, i.e., it encloses the code that are
> part of the function.
{: .tips}

Let's write a simple function.

~~~
$ myfunc ()
  {
      echo "Hello, world!"
  }
~~~
{: .language-bash}

### Calling a function

As mentioned earlier, a function is called by its name. 
Let's call the `myfunc` function defined earlier.

~~~
$ myfunc
~~~
{: .language-bash}

~~~
Hello, world!
~~~
{: .output}

The code `echo ...` was executed when you called the function.

### Arguments

Another important component of functions are arguments. 
Arguments are values passed to the function
from interactive prompt or from parent script. 
Arguments enable the function to use different input data in different calls.

You can pass arguments to a function by simply 
supplying the values after the function name. 
The function stores the arguments passed to it in variables 
$1 (first argument), $2 (second argument), ... and so on.

~~~
$ myfunc2()
  {
      echo $1
      echo $2
  }

$ myfunc2 "Hello" "World"
~~~
{: .language-bash}

~~~
Hello
World
~~~
{: .output}

Let's call the function with new arguments.

~~~
$ myfunc2()
  {
      echo $1
      echo $2
  }

$ myfunc2 "This workshop" "is awesome."
~~~
{: .language-bash}

~~~
This workshop
is awesome.
~~~
{: .output}

> ## Return value
>
> Unlike many other programming languages,
> bash function do not return a data back,
> but rather return a status.
>
> Instead values can be saved to a global variable or a file.
> This part will not be discussed here.
>
> Alternatively, standard output, e.g., `echo`, can be captured as follows.
>
> ~~~
> $ echo $(myfunc2) | wc -c
>
> $ wc -c <<< $(myfunc2)
> ~~~
> {: .language-bash}
>
{: .notes}

> ## Exercise: Write your own script
> 
> 1. Write a function to find the maximum of two input numbers.
> 2. Run the function with 1 and 2.
> 3. Run the function with 4 and 3.
> 
> <details markdown=1>
>   <summary></summary>
> 
> ~~~
> # 1
> $ max_num() {
>     if [ $1 -gt $2 ]
>     then
>         echo $1
>     else
>         echo $2
>     fi
>   }
>
> #1 alternative one-liner
> $ max_num() {
>     [[ $1 -gt $2 ]] && echo $1 || echo $2
>   }
> ~~~
> {: .language-bash}
> 
> ~~~
> # 2
> $ max_num 1 2
> ~~~
> {: .language-bash}
> 
> ~~~
> 2
> ~~~
> {: .output}
> 
> ~~~
> # 3
> $ max_num 4 3
> ~~~
> {: .language-bash}
> 
> ~~~
> 4
> ~~~
> {: .output}
> </details>
{: .challenge}

## An intro to scripting

Scripts are another way of writing reusable portable code.
Scripts are files which contain code to be executed and 
are called interactively or from another script.

Scripts are generally considered more standalone than functions.
Scripts are preferrable to functions when the code is reused in
many other scripts.

### Writing a script

A typical UNIX script looks like this:

~~~
#!/bin/sh

your code here ...
...
~~~
{: .terminal}

> ## Script shell
> The first line of a script is a shebang, i.e., `#!` , 
> followed by the shell to run the script.
>
> A shell is an interpreter, which converts user command to machine language.
> - The Bourne shell `sh` is available in all modern UNIX systems 
> and is available from path `/bin/sh`.
> - The Bourne-Again shell `bash` is present in most modern UNIX systems 
> and is the default shell in many of them. 
> It is [slightly feature rich](https://www.gnu.org/software/bash/manual/html_node/Major-Differences-From-The-Bourne-Shell.html) than `sh` and 
> is available from path `/bin/bash`.
> - `zsh`, `fish`, `ksh`, `dash` etc. are other popular UNIX shells.
> - Many programming languages provide their shells. 
> For example, python scripts can be run from `/usr/bin/env python` or 
> `/usr/bin/env python3` etc. and R scripts from `/usr/bin/env Rscript`.
>
> If you are not sure about shell selection, 
> `bash` is the safest choice for UNIX scripting.
{: .tips}

Let's write a simple script. 
Use [`nano`](/02-slurm/#adding-information-to-slurm-script) to write a simple script as follows. 
Once you are done, press <kbd>Ctrl</kbd>+<kbd>x</kbd> to return to bash prompt.
Press <kbd>Y</kbd> and <kbd>Enter</kbd> to save the changes made to the file.

~~~
$ nano code.sh
~~~
{: .language-bash}

~~~
-----------------------------------------------------------------------------------------------
 GNU nano 3.3 beta 02                     File: code.sh
-----------------------------------------------------------------------------------------------
#!/bin/bash

echo "Hello, World!"




-----------------------------------------------------------------------------------------------
^G Get Help     ^O WriteOut     ^R Read File     ^Y Prev Page     ^K Cut Text       ^C Cur Pos
^X Exit         ^J Justify      ^W Where Is      ^V Next Page     ^U UnCut Text     ^T To Spell
-----------------------------------------------------------------------------------------------
~~~
{: .terminal}

### Making script executable

So far, your script is a regular text file and does not have **execution** permission.
You can make it executable by using `chmod +x`.

~~~
$ chmod +x code.sh
~~~
{: .language-bash}

### Calling the script

You can run the script by using its filepath (path+filename). 

~~~
$ ./code.sh
~~~
{: .language-bash}

~~~
Hello, World!
~~~
{: .output}

### Arguments

Just like functions, arguments to a script can be passed by simply 
supplying the values after the call. 

Similarly, argument values can be accessed from within script with 
$1, $2 and so on.

> ## Parsing the argument
>
> Argument parsing in a script allows
> - passing arguments without a fixed order
> - pass options without any arguments
>
> An example is `echo -h` for help, where `-h` is an option.
>
> We do not cover argument parsing here, but you can learn more about it 
> [here](https://www.shellscript.sh/tips/getopts/).
{: .notes}

> If you are publishing your own script for others to use, 
> it is considered best practice to adequately comment each step of the script.
{: .tips}

> ## Exercise: Write your own script
> 
> 1. Use what you learnt so far to write a script that does the following:  
>     a) Create a directory named in the first argument in current directory.  <input type="checkbox">  
>     b) Create a file that named in the second argument within directory from step 1.  <input type="checkbox">  
>     c) Write the content of third argument to the file from step 2.  <input type="checkbox">  
>     d) Print content of file created in step 2.  <input type="checkbox">  
> 
> 2. Next, run the script with following arguments:
>    - first argument: demo
>    - second argument: text
>    - third argument: "Success!"
> 
> <details markdown=1>
>   <summary></summary>
> 
> ~~~
> # 1
> nano script.sh
> ~~~
> {: .language-bash}
> 
> ~~~
> ------------------------------------------------------------------------------------------
>  GNU nano 3.3 beta 02                   File: script.sh
> ------------------------------------------------------------------------------------------
> #!/bin/bash
> 
> # 1a - creating a directory from first argument
> mkdir -p $1
> 
> # 1b - creating a file from second argument within the directory
> touch $1/$2
> 
> # 1c - writing third argument to the file created
> echo $3 > $1/$2
> 
> # 1d - reading the content of the file
> cat $1/$2
>
> ------------------------------------------------------------------------------------------
> ^G Get Help    ^O WriteOut    ^R Read File    ^Y Prev Page    ^K Cut Text      ^C Cur Pos
> ^X Exit        ^J Justify     ^W Where Is     ^V Next Page    ^U UnCut Text    ^T To Spell
> ------------------------------------------------------------------------------------------
> ~~~
> {: .terminal}
>
> ~~~
> # 2
> $ chmod +x $1/$2     # making script executable
> $ ./script.sh demo text "Success!"
> ~~~
> {: .language-bash}
> 
> ~~~
> Success!
> ~~~
> {: .output}
>
> > Using `-p` argument with `mkdir` eliminates errors related to creating 
> >  nested or pre-existing directories.
> {: .tips}
> 
> </details>
{: .challenge}


## Parallelization

In computer science, parallelization is the execution of multiple, 
often similar/repeated tasks concurrently.

Parallelization is important in bioinformatics because of:
- large dataset leading to lengthy computation time, e.g., 50 million reads, 30 Gb file
- similar processing for each unit of data, e.g., alignment of each read to reference genome using identical processes.
- availability of computing clusters with large computing resources.

Let's look at a demo example of parallelization using [GNU parallel](https://www.gnu.org/software/parallel/).

First let's define function that takes some time to complete.

~~~
$ myfunc3()
  {
    sleep 2   # suspend computing for 2 seconds
    echo $1
  }
~~~
{: .language-bash}

We can execute this n times function sequentially, i.e., not in parallel using a for loop.
~~~
$ for i in {1..10}; do myfunc3 $1; done
~~~
{: .language-bash}

We need some way to monitor the time taken. Let use some simple math for this, 
i.e., subtract start time from end time.

~~~
$ start=$SECONDS   # Start time
$ for i in {1..10}; do myfunc3 $i; done
$ end=$SECONDS     # End time
$ echo -e "---\nExecution time is $(echo "$end - $start" | bc -l) seconds."
~~~
{: .language-bash}

~~~
1
2
3
4
5
6
7
8
9
10
---
Execution time is 20 seconds.
~~~
{: .output}

Now let us parallelize.

## Parallelization with background tasks

The easiest way to parallelize is bash is by using background task. 
This is accomplished using `&` at the end of the command (within loop).

~~~
$ set +m           # Disable monitoring mode to prevent job status output
$ start=$SECONDS   # Start time
$ for i in {1..10}; do myfunc3 $i & done
$ wait             # Wait for all tasks to complete
$ end=$SECONDS     # End time
$ echo -e "---\nExecution time is $(echo "$end - $start" | bc -l) seconds."
~~~
{: .language-bash}

~~~
1
2
3
4
5
6
7
8
9
10
---
Execution time is 2 seconds.
~~~
{: .output}

> Your screen output might look slightly different. 
> The background task management in RedHat Linux (used by HiperGator)
> prints some more job id information on the screen. 
> Nonetheless, output of `myfunc3` function should be same as above.
{: .notes}

## GNU parallel

However, inbuilt parallelization methods in UNIX system are not
very user-intuitive and flexible, and sometime vary slightly between different
UNIX distros. GNU `parallel` is a great alternative for parallelization.
Among other benefits, it can set number of CPU cores per job based on number of jobs
and has options to specify minimum memory for a job,
which make sit a useful tool in bioinformatics.

In a personal computer,
parallel needs to be installed [source](https://ftp.gnu.org/gnu/parallel/). 
However, GNU `parallel` is available in HiperGator –
it only needs to be loaded.

~~~
$ ml parallel
~~~
{: .language-bash}

User-defined objects (variables and functions) are not accessible 
within parallel by default,
so the function `myfunc3` has to be made global (globally available) first.

~~~
$ export -f myfunc3
~~~
{: .language-bash}

Now, let's run 5 jobs at the same time.

~~~
$ start=$SECONDS   # Start time
$ parallel -j5 --ungroup --will-cite myfunc3 ::: {1..10}
$ end=$SECONDS     # End time
$ echo -e "---\nExecution time is $(echo "$end - $start" | bc -l) seconds."
~~~
{: .language-bash}

> ## Breakdown of GNU parallel
> - `parallel` is the command name.
> - `-j5` specifies 5 jobs at a time.
> - `--ungroup` specifies printing output as soon as one input is processed (optional).
> - `--will-cite` is added to suppress prompt to cite parallel (optional).
> - `myfunc3` is the function to execute.
> - `{1..10}` is an array of 10 arguments.
> 
> The function `myfunc3` will be run 10 time,
> each with a number between 1 to 10 as an argument.
{: .notes}

~~~
1
2
3
4
5
6
7
8
9
10
---
Execution time is 5 seconds.
~~~
{: .output}

Thus, with parallelization, we were able to reduce the execution time.


> ## Become a shellscript pro
> [Advanced Bash-Scripting Guide](https://tldp.org/LDP/abs/html/)
{: .notes}
