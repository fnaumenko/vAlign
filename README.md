# vAlign
**V**erify **Align**ment is a fast verifier of reads forming the aligned DNA sequence, 
which is recalled from an artificial [FastQ](https://en.wikipedia.org/wiki/FASTQ_format) sequence. 
It compares the original and actual coordinates of each read and prints statistics of right and wrong mappings.

To do this each read in an initial artificial sequence should keep its location as a part of its name. 
Such template sequence can be generated by [**isChIP**](https://github.com/fnaumenko/isChIP) software. 
Both single end and paired end reads are possible. 
If initial reads do not answer this requirement, the treatment will be rejected.

It runs on the command line under Linux and Windows.

## Installation
### Executable file

**Linux**<br>
Go to the desire directory and type commands:<br>
```wget -O vAlign.gz https://github.com/fnaumenko/vAlign/releases/download/1.0/vAlign-Linux-x64.gz```<br>
```gzip -d vAlign.gz```<br>
```chmod +x vAlign```

**Windows**<br>
Download archive from [here](https://github.com/fnaumenko/vAlign/releases/download/1.0/vAlign-Windows-x64.zip) 
and unzip by any archiver, for instance [WinRar](https://www.win-rar.com/download.html?&L=0).

### Compiling in Linux
Required libraries:<br>
g++<br>
zlib (optionally)

Go to the desired directory and type commands:<br>
```wget -O vAlign.zip https://github.com/fnaumenko/vAlign/archive/1.0.zip```<br>
```unzip vAlign.zip```<br>
```cd vAlign-1.0```<br>
```make```

If **zlib** is not installed on your system, a message will be displayed from the linker.<br>
In that case you can compile the program without the ability to work with .gz files. 
To do this, open *makefile* in any text editor, uncomment last macro in the second line, comment third line, save *makefile*, and try again ```make```.<br>
To be sure about **zlib** on your system, type ```whereis zlib```.

## Usage
```vAlign [options] -g|--gen <name> sequence```

### Help
```
Input:
  -g|--gen <name>       reference genome library or single nucleotide sequence. Required
  -c|--chr <chars>      treat specified chromosome only
  --min-scr <int>       score threshold for treated reads
  --char-case <OFF|ON>  recognize uppercase and lowercase characters in template and test as different [OFF]
Output:
  -i|--info <NM|CNT|STAT>       print information about file:
                        NM - name only, CNT - number of reads, STAT - statistics [STAT]
  -w|--warn             print each read ambiguity, if they exist
  -o|--out              duplicate standard output to vAlign_out.txt file
Other:
  -t|--time             print run time
  -v|--version          print program's version and exit
  -h|--help             print usage information and exit 
```

## Details

### Input
Aligned DNA sequence in [BED](https://www.ensembl.org/info/website/upload/bed.html) format.

Compressed files in gzip format (.gz) are acceptable.

### Output
**vAlign** outputs number of exactly matched reads, and number of wrong placed reads with 0, 1, 2, … N mismatches, where N is length of read.

### Options
Enumerable option values are case insensitive.

```-g|--gen <name>```<br>
Reference genome library or single nucleotide sequence.<br>
Genome library is a directory containing nucleotide sequences for each chromosome in [FASTA](https://en.wikipedia.org/wiki/FASTA_format) format.<br>
If ```name``` is a .fa[.gz] file, **vALign** accepts the corresponding chromosome as the only treated.<br>
Otherwise first the program searches for .fa files in the directory ```name```. If there are no such files in this directory, **vALign** searches for .fa.gz files.<br>
If chromosome is specified by option ```–c|--chr```, the program searches for the corresponding .fa[.gz] file.

One can obtain a genome library in UCSC ftp://hgdownload.soe.ucsc.edu/goldenPath/ or in Ensemble ftp://ftp.ensembl.org/pub/release-73/fasta storage. 
In the second case please copy genomic sequences with the same masked type only, f.e. unmasked (‘dna'), since program does not recognise mask’s types.<br>
This option is required.

```-c|--chr <chars>```<br>
Treat specified chromosome only. Samples of option’s value: 1, 20, X.<br>
Reduces run time by 1.5-20 times depending on how far the chromosome is placed in an alignment.<br>
Default: all.

```--min-scr <int>```<br>
Score threshold for treated reads. Reads with the score equal or less then specified will be ignored.<br>
Default: all reads are accepted.

```--char-case <OFF|ON>```<br>
Recognize uppercase and lowercase characters in template and test as different.<br>
Default: ```OFF```.

```-i|--info <NM|CNT|STAT>```<br>
Output information about number of items (features/reads/intervals).<br>
```NM```:  brief output. Prints file names without any additional information.<br>
```CNT```:  prints file names and number of all and accepted items, if they are different.<br>
```STAT```: prints item ambiguities statistics, if they exist.<br>
For the alignments, such ambiguities can be duplicated reads or reads with different length (size).<br>
Thus, not all reads present in the file can be accepted.<br>
In some circumstances you need to be aware of these issues. 
There are two methods used to identify them: detailed and summary.<br>
The ```STAT``` value provides the summary method. 
It forces to display number of all recognized certain type ambiguities, and appropriate treatment.<br>
In particular, this is a simple way to know the number of duplicated reads.<br>
The detailed method is managed by option ```-w|--warn```.

```-w|--warn```<br>
Output ambiguity for each read, if they exist.<br>
If it is specified, information about type of ambiguity, number of line where it occurs, and resulting treatment would be printed each time it happens.<br>
Duplicated reads are not printed when using this method as they are considered normal for the alignment, 
but they are reported in the summary method, see ```-i|--info``` option.

```-o|--out```<br>
Duplicate standard output to **vAling_out.txt** file. 
It is an analogue of the *tee* Linux command and is rather useful by calling **vAling** under Windows.

##
If you face to bugs, incorrect English, or have commentary/suggestions, please do not hesitate to write me on fedor.naumenko@gmail.com
