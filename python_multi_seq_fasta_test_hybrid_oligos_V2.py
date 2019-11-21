#!/apps/x86_64/python/2.7.3/bin/python
# File created on March 14, 2017 by Shannon Whitmer, evk3@cdc.gov

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
from Bio.SeqUtils import MeltingTemp as mt
import numpy as np


for myseq in SeqIO.parse("./hantavirus_M_1200_3700_segments_seq_1st_filter_oligos.fasta", "fasta"):
	#print(seq_record.id)
	#print 'seq %s is %i bases long\toligo_sequence\tPercent_GC\tTm_mean' % (myseq.id[0:59], len(myseq) )

	i=0
	j=80
	while j < len(myseq):
		a = [ mt.Tm_NN(myseq.seq[i:j], strict=False), mt.Tm_GC(myseq.seq[i:j], strict=False) ]
		Tm_mean = ( sum( a ) / 2 )
	
	
		print '%s_M_%i_%i\t%s\t%0.2f\t%0.2f' % (myseq.id[0:29], i, j, myseq.seq[i:j], GC(myseq.seq[i:j]), Tm_mean)
		#print 'GC content is: %i' % (GC(myseq.seq[i:j]))
		#print('%0.2f' % mt.Tm_NN(myseq.seq[i:j]) )
		#print('%0.2f' % mt.Tm_GC(myseq.seq[i:j]) )
		i+= 200
		j = i + 80

