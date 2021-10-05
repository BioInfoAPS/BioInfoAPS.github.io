---
title: "More UNIX commands"
teaching: 0
exercises: 0
questions:
- 
objectives:
- 
keypoints:
- 

---

The optional sections are targeted for more advanced users 
(psa: more advanced than regular sections, but still for beginners). 
Feel free to try these sections once you are done with regular lesson.

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
{: .language-bash}

> `function_name` is the name of the function, that will be used to call it later.
>
> `()` decalres the variable `function_name` stores a function, 
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
(from interactive prompt or from ) parent script. 
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
> This part won't be discussed here.
>
> Alternatively, standard output, e.g. `echo`, can be captured as follows.
>
> ~~~
> $ echo $(myfunc2) | wc -c
>
> $ wc -c <<< $(myfunc2)
> ~~~
> {: .language-bash}
>
{: .notes}

## An intro to scripting

Scripts are another way of writing reusable portable code.
Scripts are files which contain code to be executed and 
are called interactively or from another script.

Scripts are generally considered more standalone than functions.
Scripts are preferrable to functions when the code is reused in
many other scripts.

### Writing a script

A typical unix script looks like this:

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
> - The Bourne shell `sh` is available in all modern unix systems 
> and is available from path `/bin/sh`.
> - The Bourne-Again shell `bash` is present in most modern unix systems 
> and is the default shell in many of them. 
> It is [slightly feature rich](https://www.gnu.org/software/bash/manual/html_node/Major-Differences-From-The-Bourne-Shell.html) than `sh` and 
> is available from path `/bin/bash`.
> - `zsh`, `fish`, `ksh`, `dash` etc. are other popular unix shells.
> - Many programming languages provide their shells. 
> Python scripts can be run from `/usr/bin/env python` or 
> `/usr/bin/env python3` etc. and R scripts from `/usr/bin/env Rscript`.
>
> If you are not sure about shell selection, 
> `bash` is the safest choice for unix scripting.
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

In computer science, parallelization is the execucation of multiple, 
often similar/repeated tasks concurrently.

Parallelization is important in bioinformatics because of:
- large dataset leading to lengthy computation time, e.g., 50 million reads, 30 Gb file
- similar processing for each unit of data, e.g., alignment of each read to reference genome using identical processes.
- availability of computing clusters with large computing resources.

Let's look at a demo example of parallelization useing [GNU parallel](https://www.gnu.org/software/parallel/).

First let's define function that takes some time to complete.

~~~
$ myfunc3()
  {
    sleep 2   # suspend computing for 2 seconds
    echo $1
  }
~~~
{: .language-bash}

We can execute this n times function sequenctially, i.e., not in parallel using a for loop.
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

Now let's parallelize.

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

However, inbuilt parallelization methods in unix system are not
very user-intuitive and flexible, and sometime vary slightly between different
unix distros. GNU `parallel` is a great alternative for parallelization.
Among other benefits, it can set number of cpu cores per job based on number of jobs
and has options to specify minimum memory for a job,
which make sit a useful tool in bioinformatics.

In a personal computer,
parallel needs to be installed [source](https://ftp.gnu.org/gnu/parallel/). 
However, GNU `parallel` is available in HiperGator –
it only needs to to be loaded.

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