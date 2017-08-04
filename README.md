# vAlign
<b>v</b>erify <b>Align</b>ment is a fast verifier of reads of aligned DNA sequence, recalled from initial artificial FastQ sequence. It compares the original and actual coordinates of each read and prints statistics of right and wrong mappings.
To do this each read in an initial artificial sequence should keep its location as a part of its name. Such template sequence can be generated by isChIP package. Both single-end and paired-end reads are possible. If initial reads do not answer this requirement, the treatment will be rejected.

## Usage
vAlign [options] -g|--gen <name> sequence

## Help
```
Input:
  -g|--gen <name>       genome size file, or genome library, or single nucleotide sequence<br>
  -c|--chr <chars>      treat stated chromosome only (all)<br>
  --min-scr <int>       score threshold for treated reads (lack)<br>
  --char-case <OFF|ON>  recognize uppercase and lowercase characters in template and test
                        as different [OFF]<br>
Output:
  --alarm               output features ambiguities, if they exist
  --stat                output features ambiguities statistics, if they exist
  -o|--out              duplicate standard output to vAlign_out.txt file
Other:
  -t|--time             output run time
  -v|--version          print program's version and quit
  -h|--help             print usage information and quit
  ```

## Details

### Input
Aligned DNA sequence in BED format.
Zipped files (.gz) are accepted too.

### Output
vAlign outputs number of exactly matched reads, and number of wrong placed reads with 0, 1, 2, … N mismatches, where N is length of read.
The output can be duplicated into a file (see ```-o|--out``` option).

### Options description
```-g|--gen <file>```<br>
Genome size file, or genome library, or single nucleotide sequence. 
Genome library is a directory contained nucleotide sequences for each chromosome in FASTA format.
The difference between genome size file and genome library/file is that in the last case all the undefined regions in reference genome (gaps), will be excluded from calculation. 
Undefined regions are regions with only ambiguous reference characters ‘N’ in them.
The minimal length of accounting gaps is managed by ```--gap-len``` option.
For example, chromosome 1 from mm9 library contains 14 regions, separated by gaps with length more then 400 bps, and 10 regions, separated by gaps with length more then 1000.
Indicating genome library has the same effect as ```-f|--fbed``` option, where ‘template’ is a set of defined regions.
You can obtain genome library in UCSC or in ensemble storage. In the second case please copy genomic sequences with the same masked type only, f.e. unmasked (‘dna'), since program does not recognise mask’s types.
Zipped .fa files can be mixed with unzipped.
The single pointed FASTA file has the same effect as ```-c|--chr``` option.
This option is required.

```-c|--chr <chars><br>```<br>
Treat stated chromosome only. Samples of option’s value: 1, 20, X.
Reduces run time on 1.5-20 times depends of how far this chromosome is placed in an alignment.<br>
Default: all.

```--min-scr <int><br>```<br>
Score threshold for treated reads. Reads with the score equal or less then stated will be ignored.<br>
Default: all reads are accepted.

```--char-case <OFF|ON><br>```<br>
Recognize uppercase and lowercase characters in template and test as different.<br>
Default: OFF.

```--alarm```<br>
Output ambiguities, if they exist.
For alignments such ambiguities can be duplicated reads or reads with different length (size).
In some circumstances you need to be aware about these issues. There are two ways to know about them: detailed and summary.
This option provides a detailed way. If it is pointed, information about type of ambiguity, number of line where it occurs, and resulting treatment would be printed each time when it happens.
Summary way is managed by option ```–-stat```.
Duplicated reads are not printed this way since they are considered as normal for alignment, but they are reported in summary way as well.

```--stat```<br>
Output features ambiguities statistics, if they exist. Prints number of all recognized certain type ambiguities, and appropriate treatment.
In particular, this is a simple way to know the number of duplicated reads in alignments.
For more details about ambiguities see --alarm option.

```-o|--out```<br>
Duplicate standard output to bioCC_out.txt file. It is analogue of tee Linux command and rather is useful by calling under Windows.

