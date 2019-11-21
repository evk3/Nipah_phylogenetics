# Nipah_phylogenetics
Collection of scripts used for "Inference of Nipah virus Evolution, 1999-2015" 


## Enrichment Probes

How does one make enrichment probes?  It's surprisingly simple.  The user generates a fasta file containing sequences of interest.  A multiline or single line fasta can be used as input.  These sequences can represent only a limited subset of the available sequences, or perhaps representative sequences (maybe for each viral clade?).  Some limited data suggests that an enrichment probe set works better if there are many diverse sequences in your fasta.

The [enrichment probe script](/python_multi_seq_fasta_test_hybrid_oligos_V2.py) will slice a fasta sequence to generate enrichment probe sequences.  The probe size is controled by "j" and the spacing between the probes is controled by "i".  Currently the script will generate probes of 80 bp long with a spacing of 200 bp between the probes.  If you want to change the probe size (100 or 120 bp length?) just change the value of j.  If you want the probes to have no gaps between them, then change the value of i to j+1.

Also, this was written in Python 2.  If you want to use it with Python 3, you probably just need to change the "print" syntax.

Script usage:
'''
python python_multi_seq_fasta_test_hybrid_oligos_V2.py > enrichment_probes.txt
'''

## Enrichment probes used in "Inference of Nipah virus Evolution, 1999-2015"

File containing probe sequences can be found [here](/Nipah_oligos.xls)
Probes can be purchase from [Twist Biosciences](https://twistbioscience.com).  Please note that probes need to be 5' biotinylated.

My last conversation with the folks at Twist Biosciences indicated that they could make approximately 6,500 probes at a total probe concentration of 30ug.

##
