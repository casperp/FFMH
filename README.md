# FFMH
FFMH stands for Find FIMO motif hits.
This tool is made to quickly find genes after an FIMO motif hit.
___

## Usage
This script sshould work on most linux and mac osx systems with python 3 installed.
The only file's needed to work are the fimo.txt with the hits and the full genbank file form the organism

### syntax 
Python3 FFMH_scan.py {txt file from FIMO} {gb(ff) from the organism} {number of bp before the gene} {type of output}

As example:

```python3 FFMH_scan.py fimo.txt genome.gb 300 0```

This command will give feedback in the terminal how far the process is.

## Output
The defaut output is  an extensive overview per hit. all the info from the gen and all the motif info.
the second output is a small one line output with:
 >gene_loc	product	id	start_motif	end_motif	streng	motif
 
| output  | command  |
|:--------: |:--------:| 
| extensive      | **0** | 
| small      | **1** |

____________
Build by:
Casper Prins