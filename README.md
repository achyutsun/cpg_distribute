# cpg_distribute

cpg_distribute as written by Tim Cutts and Gos Micklem
    
This directory contains the source code of program written by Tim Cutts
and Gos Micklem (in 2004) that detects CpG islands in nucleic acid
sequences. The program is referenced in the Ensembl database, but
is not available elsewhere to my knowledge.

The code provided in this commit was provided to my by Gos Micklem.

Usage:
```
./cpg -h

Usage: ./cpg [options] [file]

        Options:

                -h      This message

                -a      .ace format output

                -l nn   length threshold
                -g nn   %GC threshold
                -o nn   observed/expected threshold

        Stringent CpG island: -l 1000 -g 50 -o 0.6


* Sensible values to run are:
./cpg -l 400 -g50 -o0.6 14b7.seq



"acedb output":
./cpg -l 400 14b7.seq
EM:HS14B7       1950    2383    143     CpG:32  62.4    0.71
seqname         start   end     score  #_of_CpGs  %GC   observed/ expected ratio

Default output:
./cpg -al 400 14b7.seq
Feature Predicted_CpG_island 1950 2383 0.71 "Predicted CpG island: 62.4 %GC,   o/e=0.71,   #CpGs= 32"
```
