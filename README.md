# LongTarget: Predicting long noncoding RNAs' epigenetic target genes

## Contents

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
    + [Compilation](#compilation)
    + [Running](#running)
    + [Help information](#help-information)
    + [Time consumption](#time-consumption)
- [Demo](#demo)
    + [Inputs and their formats](#inputs-and-their-formats)
    + [Results](#results)
    + [Example datasets](#example-datasets)
 - [Instructions for use](#instructions-for-use)
    + [How to run LongTarget on your data](#how-to-run-longtarget-on-your-data)
    + [Bug reports](#bug-reports)
- [Other codes](#other-codes)
- [License](./LICENSE)
- [Citation](#citation)

# Overview
Since the pioneering genome-wide discovery of mouse lncRNAs in the FANTOM consortium, experimental studies have identified abundant lncRNAs in humans, mice, and other mammals. Many lncRNAs can bind to both DNA sequences and DNA- and histone-modifying enzymes, thus recruiting these enzymes to specific genomic sites to epigenetically regulate the expression of genes at these sites. Genomic imprinting is a specific kind of epigeneticregulation. 

LongTarget was developed to predict one or many lncRNA’s DNA binding motifs and binding sites in one or many genome regions based on all known ncRNA/DNA base pairing rules (Hoogsteen and reverse Hoogsteen base pairing rules). LongTarget consists of a few C/C++ programs, and is distributed under the AGPLv3 license. It can be used as a standalone program, and we have also integrated it into the lncRNA database LongMan at the website http://lncRNA.smu.edu.cn, making it a web service. Our use of LongTarget indicates that it can satisfactorily predict lncRNAs’ DNA binding motifs and binding sites, and multiple pipelines have been developed to seamless bridge database search and lncRNA/DNA binding prediction.

# Repo Contents
- [longtarget.cpp](./longtarget.cpp): The main program for generating the executable "LongTarget".
- [rules.h](./rules.h): Base-pairing rules and codes for handling these rules.
- [sim.h](./sim.h):   The SIM program for local alignment.
- [stats.h](./stats.h): Michael Farrar's code (with SSE2) for local alignment.
- [H19.fa](./H19.fa):  A sample lncRNA sequence.  
- [testDNA.fa](./testDNA.fa): A sample DNA sequence. 

# System Requirements
- OS: Linux, we compile and run the LongTarget under CentOS 6.0. 
- System software:	g++.
- RAM: 16G or above, depending on the number of lncRNAs and length of genome region.
- CPU: 4 cores or above, depending on the number of lncRNAs and length of genome region.
- To use our web service, Google Chrome and Mozilla Firefox are recommended, because functions were tested under these browsers. 

# Installation Guide
## Compilation
Typically, this command will generate an executable LongTarget program: 

```
g++ longtarget.cpp -O -msse2 -o LongTarget.
```

## Running 
A simple case is:

```
./LongTarget -f1 testDNA.fa -f2 H19.fa -r 0
```

A more complex case is:

```
./LongTarget -f1 testDNA.fa -f2 H19.fa -r 0 -O /home/test/example -c 6000 -i 70 -S 1.0 -ni 25 -na 1000 -pc 1 -pt -500 -ds 10 -lg 60
```

In this command, output path is /home/test/example, cut sequence's length is 6000, identity is 70%, stability is 1.0, ntMin is 25 nt, ntMax is 1000 nt, penaltyC is 1, penaltyT is -500, distance between TFOs is 10, min length of triplexes is 60 (for more details about parameters, visit the website http://lncRNA.smu.edu.cn).

## Help information
Here is a brief explanation of the command line arguments:

```
Options   Parameters      Functions
f1   DNA sequence file  A string, indicate the DNA sequence file name.
f2   RNA sequence file  A string, indicate the RNA sequence file name.
r    rules              An integer, indicate base-pairing rules, "0" indicates all rules. 
O    Output path        A string, indicate the directory into which the results are outputted.
c    Cutlength          An integer, indicate the length of each segment, the default value is 5000.
i    identity           An integer, indicate the criterion of alignment output, the default value is 60.
S    stability          A floating point, indicate the criterion of base-pairing, the default value is 1.0.
ni   ntmin              An integer, indicate the min length of triplexes, the default value is 20.
na   ntmax              An integer, indicate the max length of triplexes, the default is 100000 but is rarely used.
pc   penaltyC           An integer, indicate penalty, the default value is 0.
pt   penaltyT           An integer, indicate penalty, the default value is -1000.
ds   c_dd               An integer, indicate the distance between TFOs, the default value is 15.
lg   c_length           An integer, indicate the min length of triplexes, the default value is 50.
```

## Time consumption
This depends on the number and length of lncRNAs and the length of genome regions. The expected running time for the H19 demo should be no more than ten minutes even on a normal desktop computer. 

# Demo
## Inputs and their formats
H19.fa and testDNA.fa in the directory give a demo example. To obtain more details, go to our website http://lncRNA.smu.edu.cn and/or check files in the "examples" subdirectory.

The H19.fa indicates that the lncRNA sequence file should have a title line in the format ">species_lncRNA" without any space within letters, and the lncRNA sequence should be in a new line.

The testDNA.fa indicates that the DNA sequence file should have a title line in the format ">species|chr|start-end" without any space between letters, and the DNA sequence should be in a new line.

## Results
The results include three files whose filenames ending with: (1)*TFOsorted, (2)*TFOclass1, 
(3)*TFOclass2. The TFOsorted file contains the details of all triplexes, the TFOclass1 file contains the TTS distribution of TFO1 in the genome region, and the TFOclass2 file contains the TTS distribution of TFO2. 

## Example datasets
An example dataset giving detailed results of examples is given in the subdirectory "examples".

# Instructions for use
## How to run LongTarget on your data
To run LongTarget using the web service, both lncRNAs and genome sequences are available in the database LongMan. To run it as a standalone program, you should obtain lncRNAs and genome sequences from websites such as https://www.gencodegenes.org/ , http://genome.ucsc.edu/ , http://www.noncode.org/ , and http://asia.ensembl.org/index.html . 
  
## Bug reports
Please send comments and bug reports to: zhuhao@smu.edu.cn.

# Other codes
The [DatabaseScripts](./DatabaseScripts) subdirectory contains the scripts written in Python and Perl for building LongMan database. The web application of LongMan is available at http://lncrna.smu.edu.cn/show/info .

# License
The program is distributed under the AGPLv3 license.

# Citation
