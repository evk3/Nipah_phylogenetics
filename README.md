# Nipah_phylogenetics
Collection of scripts used for "Inference of Nipah virus Evolution, 1999-2015" 


## Enrichment Probes

How does one make enrichment probes?  It's surprisingly simple.  The user generates a fasta file containing sequences of interest.  A multiline or single line fasta can be used as input.  These sequences can represent only a limited subset of the available sequences, or perhaps representative sequences (maybe for each viral clade?).  Some limited data suggests that an enrichment probe set works better if there are many diverse sequences in your fasta.

The [enrichment probe script](/python_multi_seq_fasta_test_hybrid_oligos_V2.py) will slice a fasta sequence to generate enrichment probe sequences.  The probe size is controled by "j" and the spacing between the probes is controled by "i".  Currently the script will generate probes of 80 bp long with a spacing of 200 bp between the probes.  If you want to change the probe size (100 or 120 bp length?) just change the value of j.  If you want the probes to have no gaps between them, then change the value of i to j+1.

Also, this was written in Python 2.  If you want to use it with Python 3, you probably just need to change the "print" syntax.

Script usage:
```
python python_multi_seq_fasta_test_hybrid_oligos_V2.py > enrichment_probes.txt
```

## Enrichment probes used in "Inference of Nipah virus Evolution, 1999-2015"

File containing probe sequences can be found [here](/Nipah_oligos.xls)
Probes can be purchase from [Twist Biosciences](https://twistbioscience.com).  Please note that probes need to be 5' biotinylated.

My last conversation with the folks at Twist Biosciences indicated that they could make approximately 6,500 probes at a total probe concentration of 30ug.

## Assembly script

File for our in-house assembly pipeline can be found [here](/mapping_scripts).  The script that does the heavy lifting for genome assembly is [here](/mapping_scripts/map_and_dedup_Bangladesh_NEBNextV1.sh).  The other scrips are wrapped scripts for input to CDC's High performance computing cluster (uses Sun Grid Engine).

Usage:
```
sh array_script_NEBNext_V1.sh
```

## Read Strandedness

This script was originally from from Jason [Ladner](https://github.com/jtladner/Scripts/blob/master/strand_specific_coverage/strandspec_covplot/strandspec_covplot_v1.1.py), but I made some edits and used this [verison](/strandspec_covplot_v1.2.py).  

Usage:
Please see Jason' [readme](https://github.com/jtladner/Scripts/tree/master/strand_specific_coverage/strandspec_covplot).  The usage is the same.

## Chimeric Reads

This script is from Jason Ladner.  Check out his github [page](https://github.com/jtladner/Scripts/tree/master/chimeric_reads) for original [script](https://github.com/jtladner/Scripts/blob/master/chimeric_reads/chimeric_reads_v3.6.2.py).

## Visualizing Geographic Spread from a BEAST Continuous Diffusion Analysis

Hold on to your hats, folks.  This script has a lot of moving parts.  The original geographic spread script was made by Gytis Dudas and can be found here.  Kudos to him for this excellent piece of work!  This script was used to make the really cool Ebola virus spread video found here.

######How does my modified curionia script work?
Well, let's break it down into several steps:
1. 
