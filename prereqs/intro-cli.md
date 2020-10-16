# Introduction to the Command Line Interface

Please take the opportunity to ensure that you can work through all of this material before the workshop starts. In this section, we will cover logging into a remote server and basic command line concepts, which you will need to be successful during the workshop.

While we will be available to help if you are still struggling with this material during the workshop itself, we find that participants who are comfortable with command line basics get more out of the workshop experience. If you require assistance with any of the material contained in this tutorial before the workshop begins, please let us know in the Slack channel.

## What is the command line?

The command line interface (CLI) is a text-based interface into which the user types commands in order to type commands and perform tasks. Most of the time when using the CLI, your mouse will not work.

<img src="figures/cli/what_is_cli.png" alt="terminal" width="800px"/>

## Syntax of a command

A command is made up of components separated by spaces:

* Prompt
* Command (or program)
* Parameters
* Arguments
* The separator used in issuing a command is space, number of spaces does not matter

<img src="figures/cli/cli_prompt.png" alt="prompt" width="800px"/>

A greater than sign (>) instead of a prompt means the shell is expecting more input. You can use Cntr-c to cancel the operation and return to a prompt.

## Logging into a remote server

For this we use the command Secure SHell ... SSH. Replace 'username' with your login name in the following and enter your password when prompted. Use the username / password combination you set up on the [Genome Center computing accounts portal](https://computing.genomecenter.ucdavis.edu). You will not see your password being typed.

    ssh username@tadpole.genomecenter.ucdavis.edu

for example my login is

    ssh hslyman@tadpole.genomecenter.ucdavis.edu

Once you're done working on the command line, you can exit. Anything that follows the character '#' is ignored.

    exit  # kills the current shell!

<img src="figures/cli/login.png" alt="logging in" width="800px"/>

Go ahead and log back into the server.

After opening, system messages are often displayed, followed by the "prompt".
A prompt is a short text message at the start of the command line and ends with '$' in bash shell, commands are typed after the prompt.
The prompt typically follows the form **username@server:location$**

## Directory Structure

Absolute path: always starts with ”/”

/share/workshop/epigenetics_workshop/

the folder (or file) epigenetics_workshop in the folder workshop in the folder share from root.

Relative path: always relative to our current location.

_a single dot (.) refers to the current directory_
_two dots (..) refers to the directory one level up_

<img src="figures/cli/linux_tree.png" alt="standard Linux directory tree" width="500px"/>

Usually, /home is where the user accounts reside, ie. user’s 'home' directory.
For example, for a user that has a username of "hslyman" (me): their home directory is /home/hslyman
It is the directory that a user is located after starting a new shell or logging into a remote server.

The tilde (~) is a short form of a user’s home directory.


## Command Line Basics

First some basics - how to look at your surroundings.

    pwd

present working directory ... where am I?

    ls .

list files here ... you should see nothing since your homes are empty

    ls /tmp/

list files somewhere else, like /tmp/


Because one of the first things that's good to know is *how to escape once you've started something you don't want*.

    sleep 1000  # wait for 1000 seconds!

Use Ctrl-c (shows as '^C' in screen) to exit (kill) a command. In some cases, a different key sequence is required (Ctrl-d).


### Options

Each command can act as a basic tool, or you can add 'options' or 'flags' that modify the default behavior of the tool. These flags come in the form of '-v' ... or, when it's a more descriptive word, two dashes: '\-\-verbose' ... that's a common (but not universal) one that tells a tool that you want it to give you output with more detail. Sometimes, options require specifying amounts or strings, like '-o results.txt' or '\-\-output results.txt' ... or '-n 4' or '\-\-numCPUs 4'. Let's try some, and see what the man page for the 'list files' command 'ls' is like.

    ls -R /

Lists directories and files *recursively*. How do I know which options do what?

    man ls

Within the man pages, navigate using the keyboard (spacebar, b, pgup, pgdn, g, G, /pattern, n, N, q). Look up and try the following parameters for the ls command. Remember, if you don't provide ls with arguments describing where to run the command, ls lists files in your current directory.

    ls -l
    ls -a
    ls -l -a
    ls -la  # option 'smushing' ... when no values need specifying
    ls -ltrha

And finally adding color:

    ls -ltrha --color  # single letter (smushed) vs word options (Linux)

**OR**

    ls -ltrhaG  # (MacOS)


Quick aside: what if I want to use same options repeatedly? And be lazy? You can create a shortcut to another command using 'alias'.

    alias ll='ls -lah'
    ll

## Getting Around

The filesystem you're working on is like the branching root system of a tree. The top level, right at the root of the tree, is called the 'root' directory, specified by '/' ... which is the divider for directory addresses, or 'paths'. We move around using the 'change directory' command, 'cd'. The command pwd return the present working directory.

    pwd
    cd  # no effect? that's because by itself it sends you home (to ~)
    pwd
    cd /  # go to root of tree's root system
    ls
    cd home  # go to where everyone's homes are
    pwd
    cd $USER  # $USER will be replaced with your username
    pwd
    cd /
    pwd
    cd ~  # a shortcut to home, from anywhere
    pwd
    cd .  # '.' always means this directory
    pwd
    cd ..  # '..' always means one directory up
    pwd

<img src="figures/cli/getting_around.png" alt="navigation using cd" width="800px"/>

You should also notice the location changes in your prompt.

## Absolute and Relative Paths

You can think of paths like addresses. You can tell your friend how to go to a particular store *from where they are currently* (a 'relative' path), or *from the main Interstate Highway that everyone uses* (in this case, the root of the filesystem, '/' ... this is an 'absolute' path). Both are valid. But absolute paths can't be confused, because they always start off from the same place. Relative paths, on the other hand, could be totally wrong for your friend *if you assume they're somewhere they're not*. With this in mind, let's try a few more:

    cd ~  # let's start at home

**relative** (start here, take two steps up, then down through share and workshop)

    cd ../../share/workshop/epigenetics_workshop
    pwd

**absolute** (start at root, take steps)

    cd /share/workshop/epigenetics_workshop
    pwd


Now, because it can be a real pain to type out, or remember these long paths, we need to discuss ...

## Tab Completion

Using tab-completion is a must on the command line. A single tab auto-completes file or directory names when there's only one name that could be completed correctly. If multiple files could satisfy the tab-completion, then nothing will happen after the first tab. In this case, press tab a second time to list all the possible completing names. Note that if you've already made a mistake that means that no files will ever be completed correctly from its current state, then tabs will do nothing.

The 'touch' command updates the timestamp on a file, here we use it to create three empty files.

    mkdir -p /share/workshop/epigenetics_workshop/$USER/tmp
    cd /share/workshop/epigenetics_workshop/$USER/tmp
    touch one seven september
    ls o

Pressing tab once, with no enter should complete to 'one.' Once the command is complete, pressing enter runs the command.

    ls s

Press tab. The command line completes up to 'se' since that's in common between seven and september. Press tab again and no enter, this second tab should cause listing of seven and september. Type 'v' then tab and no enter now it's unique to seven, and should complete to seven. Enter runs 'cat seven' command.

I can't overstate how useful tab completion is. You should get used to using it constantly. Watch experienced users type and they hit tab once or twice in between almost every character. You don't have to go that far, of course, but get used to constantly getting feedback from hitting tab and you will save yourself a huge amount of typing and trying to remember weird directory and filenames.

We don't really need these files, so let's clean up.

    cd ..
    rmdir tmp     # rmdir means remove directory, but this shouldn't work!
    rm tmp/one    # clear the directory first
    rm tmp/seven
    tmp/september
    rmdir tmp     # should work now that the directory is empty

#### CHALLENGE

After returning to your home directory (just enter 'cd' by itself), verify that the two following commands are equivalent (replacing 'username' with your actual username). Don't forget to tab complete!

    cd ../../home/username/
    cd ../../../../../../../home/username/

Why are these very different-looking commands equivalent?

## Create and Destroy

We already learned one command that will create a file, touch. Let's create a folder in /share/workshop/epigenetics_workshop for you to work in (you've done this already) and then another directory within that one called cli. We will use the environment variable $USER, that is your username. Wondering what the '-p' option does? Use the man page for mkdir to find out.

    mkdir -p /share/workshop/prereq_workshop/$USER/cli
    cd /share/workshop/prereq_workshop/$USER/cli

Inside the cli directory, make a new directory, "tmp."

    mkdir tmp
    ls -lh
    cd tmp/

echo text then redirect ('>') to a file.

    echo 'Hello, world!' > first.txt
    cat first.txt  # 'cat' means 'concatenate'

Why 'concatenate'? try this:

    cat first.txt first.txt first.txt > second.txt
    cat second.txt

OK, let's destroy what we just created:

    cd ../
    rmdir tmp  # 'rmdir' means 'remove directory', but this shouldn't work!
    rm tmp/first.txt
    rm tmp/second.txt  # clear directory first
    rmdir tmp  # should succeed now

So, 'mkdir' and 'rmdir' are used to create and destroy (empty) directories. 'rm' to remove files. To create a file can be as simple as using 'echo' and the '>' (redirection) character to put text into a file. Even simpler is the 'touch' command.

    touch newFile
    ls -ltra  # look at the time listed for the file you just created
    cat newFile  # it's empty!
    sleep 60  # go grab some coffee
    touch newFile
    ls -ltra  # same time?

So, 'touch' creates empty files, or updates the 'last modified' time. Note that the options on the 'ls' command you used here give you a Long listing, of All files, in Reverse Time order (l, a, r, t).

## Piping and Redirection

Pipes ('\|') allow commands to hand output to other commands, and redirection characters ('>' and '>>') allow you to put output into files.

    echo 'first' > test.txt
    cat test.txt
    echo 'second' > test.txt
    cat test.txt
    echo 'third' >> test.txt
    cat test.txt

The '>' character redirects output of a command that would normally go to the screen instead into a specified file. '>' replaces, '>>' appends.

Cut character one to three, from every line, from file 'test.txt'

    cut -c 1-3 test.txt

Same thing, piping output of one command into input of another

    cat test.txt | cut -c 1-3

Pipe cat to cut to sort (-r means reverse order sort), and then grep searches for pattern ('s') matches.

    cat test.txt | cut -c 1-3 | sort -r
    cat test.txt | cut -c 1-3 | sort -r | grep s

This is a great way to build up a set of operations while inspecting the output of each step in turn. We'll do more of this in a bit.

## History Repeats Itself

Linux remembers everything you've done (at least in the current shell session), which allows you to pull steps from your history, potentially modify them, and redo them. This can obviously save a lot of time and typing.

The 'head' and 'tail' commands view the first 10 (by default) lines of a file and last 10 lines of a file (type 'man head' or 'man tail' to consult their manuals).

    <up arrow>  # last command
    <up>  # next-to-last command
    <down>  # last command, again
    <down>  # current command, empty or otherwise
    history  # usually too much for one screen, so ...
    history | head
    history | tail
    history | tail -n 30
    history | less
    cat test.txt | cut -c 1-3 | sort -r | grep s > reallyImportantResult.txt
    rm reallyImportantResult.txt  # whoops! didn't mean to do that!
    history | tail
    !560  # re-executes 560th command (yours will have different numbers; choose the one that recreates your really important result!)

You can also search your history from the command line:

    <ctrl-r>fir  # should find most recent command containing 'fir' string: echo 'first' > test.txt
    <enter>  # to run command
    <ctrl-c>  # get out of recursive search
    <ctr-r>  # repeat <ctrl-r> to find successively older string matches

### CHALLENGE

What's the first command you executed today? How many times have you used the 'man' command today? Whatever that number is, it should be more! Just kidding. Sort of.

### CHALLENGE

How would you create a second text file - let's say 'test2.txt' - with the line that says 'third' *before* the line that says 'second'? Without directly editing the file with a text editor, of course ..

## Editing Yourself

Here are some more ways to make editing previous commands, or novel commands that you're building up, easier:

    <up><up>  # go to some previous command, just to have something to work on
    <ctrl-a>  # go to the beginning of the line
    <ctrl-e>  # go to the end of the line
    # now use left and right to move to a single word (surrounded by whitespace: spaces or tabs)
    <ctrl-k>  # delete from here to end of line
    <ctrl-w>  # delete from here to beginning of preceding word
    blah blah blah<ctrl-w><ctrl-w>  # leaves you with only one 'blah'

## Compression and Archives

As file sizes get large, you'll often see compressed files, or whole compressed folders. Note that **any good bioinformatics software** should be able to work with compressed file formats.

    gzip test.txt
    cat test.txt.gz

To uncompress a file

    gunzip -c test.txt.gz

The '-c' leaves the original file alone, but dumps expanded output to screen

    gunzip test.txt.gz  # now the file should change back to uncompressed test.txt

Tape archives, or .tar files, are one way to compress entire folders and all contained folders into one file. When they're further compressed they're called 'tarballs'. We can use wget (web get).

    wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/PhiX/Illumina/RTA/PhiX_Illumina_RTA.tar.gz

The .tar.gz and .tgz are *commonly used* extensions for compressed tar files, when gzip compression is used. The application tar is used to uncompress .tar files

    tar -xzvf PhiX_Illumina_RTA.tar.gz

Here -x = extract, -z = use gzip/gunzip, -v = verbose (show each file in archive), -f filename

Note that, unlike Windows, linux does not depend on file extensions to determine file behavior. So you could name a tarball 'fish.puppy' and the extract command above should work just fine. The only thing that should be different is that tab completion doesn't work within the 'tar' command if it doesn't see the 'correct' file extension.

## Forced Removal

When you're on the command line, there's no 'Recycle Bin'. Since we've expanded a whole directory tree, we need to be able to quickly remove a directory without clearing each subdirectory and using 'rmdir'.

    rmdir PhiX    # this should not work
    rm -rf PhiX

Here -r = recursively remove sub-directories, -f means *force*. We actually want to use those directories, so un-archive them again!

Obviously, be careful with 'rm -rf', there is no going back, if you delete something with rm, rmdir its gone!

## BASH Wildcard Characters

When we want to specify or operate on sets of files all at once, we can use 'wildcard' characters.

    ls ?hiX/Illumina

lists files in Illumina sub-directory of any directory ending in 'hiX'

    ls PhiX/Illumina/RTA/Sequence/*/*.fa

lists all .fa files a few directories down. So, '?' fills in for zero or one character, '\*' fills in for zero or more characters. find can be used to locate files of a server form.

    find . -name "*.f*"
    find . -name "*.f?"

How are these find commands different from the wildcard use in the previous ls commands?

#### Quick Note About the Quote(s)

The quote characters " and ' are different. In general, single quotes preserve the *literal* meaning of all characters between them. On the other hand, double quotes allow the shell to see what's between them and make substitutions when appropriate. For example:

    VRBL=someText
    echo '$VRBL'
    echo "$VRBL"

However, some commands try to be 'smarter' about this behavior, so it's a little hard to predict what will happen in all cases. It's safest to experiment first when planning a command that depends on quoting ... list filenames first, instead of changing them, etc. Finally, the 'backtick' characters \` (same key - unSHIFTED - as the tilde ~) causes the shell to interpret what's between them as a command, and return the result.

    echo `$VRBL`  # tries to execute a command with the name *someText*
    newVRBL=`echo $VRBL`
    echo $newVRBL

## Manipulation of a FASTA File

We just found the phiX-174 genome, so let's copy it to our current directory so we can play with it:

    cp ./PhiX/Illumina/RTA/Sequence/WholeGenomeFasta/genome.fa phix.fa

Similarly, we can also use the move command here, but then ./PhiX/Illumina/RTA/Sequence/WholeGenomeFasta/genome.fa will no longer be there:

    cp ./PhiX/Illumina/RTA/Sequence/WholeGenomeFasta/genome.fa  ./PhiX/Illumina/RTA/Sequence/WholeGenomeFasta/genome2.fa
    ls ./PhiX/Illumina/RTA/Sequence/WholeGenomeFasta/
    mv ./PhiX/Illumina/RTA/Sequence/WholeGenomeFasta/genome2.fa phix.fa
    ls ./PhiX/Illumina/RTA/Sequence/WholeGenomeFasta/

This functionality of mv is why it is used to rename files.

Note how we copied the 'genome.fa' file to a different name: 'phix.fa'

    wc -l phix.fa

counts the number of lines in the file using 'wc' (word count) and parameter '-l' (lines).

We can use the 'grep' command to search for matches to patterns. 'grep' comes from '**g**lobally search for a **r**egular **e**xpression and **p**rint'.

    grep -c '>' phix.fa

Only one FASTA sequence entry, since only one header line ('>gi\|somethingsomething...'). Using either single or double quotes around the > symbol is very important here, since '>' by itself has a meaning on the command line already: redirect.

    cat phix.fa

This may not be useful for anything larger than a virus! Let's look at start codon and 2 following:

    grep --color "ATG......" phix.fa

'.' characters are the single-character wildcards for grep

Use the '-o' option to **o**nly print the pattern matches, one per line

    grep -o "ATG......" phix.fa

Use the 'cut' command with '-c' to select characters 4-6, the second codon

    grep -o "ATG......" phix.fa | cut -c4-6

'sort' the second codon sequences (default order is same as ASCII table; see 'man ascii')

    grep -o "ATG......" phix.fa | cut -c4-6 | sort

Combine successive identical sequences, but count them ('-c' option)

    grep -o "ATG......" phix.fa | cut -c4-6 | sort | uniq -c

Finally, sort using reverse numeric order ('-rn')

    grep -o "ATG......" phix.fa | cut -c4-6 | sort | uniq -c | sort -rn

... which gives us the most common codons first

This may not be a particularly useful thing to do with a genomic FASTA file, but it illustrates the process by which one can build up a string of operations, using pipes, in order to ask quantitative questions about sequence content. More generally than that, this process allows one to ask questions about files and file contents and the operating system, and verify at each step that the process so far is working as expected. The command line is, in this sense, really a modular workflow management system.

### CHALLENGE

Many programs and data archives contain files named something like 'readme' or 'README' that contains important information for the user. How many of these files are there in the PhiX directory tree? How would you look at their contents?


## Symbolic Links

Since copying or even moving large files (like sequence data) around your filesystem may be impractical, we can use links to reference 'distant' files without duplicating the data in the files. Symbolic links are disposable pointers that refer to other files, but behave like the referenced files in commands.

**You should, by default, always use a symbolic (-s) link.**

Hard links point to the location of the data on the hard drive while symbolic links point to an secondary location of the data on the hard drive (or the original file itself)

![](./hard_vs_symbolic.png)

Draw back of hard links is that deleting the original file become much more difficult because all hard links need to be deleted first.

    ln -s PhiX/Illumina/RTA/Sequence/WholeGenomeFasta/genome.fa .
    ls -ltrhaF  # notice the symbolic link pointing at its target
    grep -c ">" genome.fa

## STDOUT & STDERR

Programs can write to two separate output streams, 'standard out' (STDOUT), and 'standard error' (STDERR). The former is generally for direct output of a program, while the latter is supposed to be used for reporting problems. I've seen some bioinformatics tools use STDERR to report summary statistics about the output, but this is probably bad practice. Default behavior in a lot of cases is to dump both STDOUT and STDERR to the screen, unless you specify otherwise. In order to nail down what goes where, and record it for posterity:

    wc -c genome.fa 1> chars.txt 2> any.err

The 1st output, STDOUT, goes to 'chars.txt'
The 2nd output, STDERR, goes to 'any.err'

    cat chars.txt

Contains the character count of the file genome.fa

    cat any.err

Empty since no errors occurred.

Saving STDOUT is pretty routine (you want your results, yes?), but remember that explicitly saving STDERR is important on a remote server, since you may not directly see the 'screen' when you're running jobs.

## Paste Command

The paste command is useful in creating tables.

    echo 'WT1' > b
    echo 'WT2' >> b
    echo 'control1' >> b
    echo 'control2' >> b
    cat b

Now we can number our four samples to conveniently refer to them in order

    for i in {1..4}; do echo $i >> a; done
    cat a
    paste a b > c
    cat c


## Running in the Background

Sometimes it's useful to continue working on the command line, even after you've executed a command that's going to take a while to finish. Normally this command would occupy the shell, and prevent you from typing in commands and receiving results. But we can 'put jobs in the background' so that they don't occupy your shell directly:

    sleep 1000000

Ctrl-Z to pause (stop) the command

    bg

To restart the last command in the background

'^Z' first suspends the sleep command. Then, 'bg' resumes running that command *in the background*, so that it doesn't occupy the terminal. The output of the 'bg' command tells you that you have one command running in the background. You could start more, suspend them, then resume them in the background, and query what background jobs are running or are suspended, not running:

    jobs

<div class="output">[1]+  Running                 sleep 1000000 &
</div>

We can also start a job in the background in one step, without having to suspend then resume it, using the '&' character at the end of the command:

    sleep 5000000 &

If we want to delete these jobs for any reason, we can kill them using the numbering that 'jobs' reveals:

    jobs

<div class="output">[1]-  Running                 sleep 1000000 &
[2]+  Running                 sleep 5000000 &
</div>

    kill %1
    jobs

<div class="output">[1]-  Terminated              sleep 1000000
[2]+  Running                 sleep 5000000 &
</div>

    kill %2
    jobs

<div class="output">[2]+  Terminated              sleep 5000000
</div>

Finally, the 'nohup' command (from 'no hangup'!) makes jobs extra resistant to lost connections or terminal problems. In other words, even jobs running in the background can be terminated if one's shell dies. 'nohup' separates the running job from the shell, so it'll keep running until it dies or is actively killed by you.

    nohup sleep 1000000 &

<div class="output">[1] 34993
username@c4-0:~/CLB$ nohup: ignoring input and appending output to nohup.out
</div>

    jobs

<div class="output">[1]+  Running                 nohup sleep 1000000 &
</div>

Output is dumped into the 'nohup.out' file unless specifically redirected in your command

    kill %1

If you used bg to send a process to the background you can use 'disown' to "nohup" the process

## Table of Processes (top)

The 'top' command prints a self-updating table of running processes and system stats. Use 'q' to exit top, 'z' to toggle better color contrast, 'M' to sort by memory use, 'P' to sort by processor use, and 'c' to toggle display of the full commands. Hit '1' to toggle display of all processors, and hit 'u' followed by typing in a username in order to only show processes (jobs) owned by that user.

<img src="figures/cli/top.png" alt="table of processes" width="800px"/>

## Shell Scripts, File Permissions

Often it's useful to define a whole string of commands to run on some input, so that (1) you can be sure you're running the same commands on all data, and (2) so you don't have to type the same commands in over and over! Let's use the 'nano' text editor program that's pretty reliably installed on most linux systems.

    nano test.sh

<img src="figures/cli/nano.png" alt="nano text editor" width="800px"/>

nano now occupies the whole screen; see commands at the bottom
Type/paste in the following (note that '#!' is an interpreted command to the shell, not a comment):

<div class="script">
#!/bin/bash
grep -o . $1 | \
    sort | \
    uniq -c | \
    sort -rn -k1,1
</div>

Use control-X to exit. You will be prompted to save (type 'y' to save) and select a name for the document (to keep the original name, press enter). Follow the instruction at the bottom of the screen.

Note that '$1' means 'the value of the 1st argument to the shell script.' In other words, the text that follows the shell script name when we run it (see below).

Though there are ways to run the commands in test.sh right now, it's generally useful to give yourself (and others) 'execute' permissions for test.sh, really making it a shell script. Note the characters in the first (left-most) field of the file listing:

    ls -lh test.sh

<div class="output">-rw-rw-r-- 1 msettles biocore 79 Aug 19 15:05 test.sh
</div>

The first '-' becomes a 'd' if the 'file' is actually a directory. The next three characters represent **r**ead, **w**rite, and e**x**ecute permissions for the file owner (you), followed by three characters for users in the owner's group, followed by three characters for all other users. Run the 'chmod' command to change permissions for the 'test.sh' file, adding execute permissions ('+x') for the user (you) and your group ('ug'):

    chmod ug+x test.sh
    ls -lh test.sh

<div class="output">-rwxr-xr-- 1 msettles biocore 79 Aug 19 15:05 test.sh
</div>

The first 10 characters of the output represent the file and permissions.

The first character is the file type, the next three sets of three represent the file permissions for the user, group, and everyone respectively.
- r = read
- w = write
- x = execute

OK! Let's run this script, feeding it the phiX genome. When we put the genome file 1st after the name of the script, this filename becomes variable '1', which the script can access by specifying '$1'. We have to provide a relative reference to the script './' because its not our our "PATH".

    ./test.sh genome.fa

<div class="output">msettles@tadpole:/share/workshop/msettles/cli$ ./test.sh genome.fa
    1686 T
    1292 A
    1253 G
    1155 C
       1 x
       1 p
       1 i
       1 h
       1 >
</div>

The script's grep command splits out every character in the file on a separate line, then sorts them so it can count the occurrences of every unique character and show the most frequent characters first ... a quick and dirty way to get at GC content.

### A CLI cheat sheet

<object data="https://files.fosswire.com/2007/08/fwunixref.pdf" type="application/pdf" width="700px" height="700px">
<embed src="https://docs.conda.io/projects/conda/en/4.6.0/_downloads/52a95608c49671267e40c689e0bc00ca/conda-cheatsheet.pdf">
      <p>This browser does not support PDFs. Please download the PDF to view it: <a href="http://yoursite.com/the.pdf">Download PDF</a>.</p>
    </embed>
</object>
