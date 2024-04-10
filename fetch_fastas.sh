for i in P99999 P00004 P0C0X8 P00091 Q93VA3
do
  wget https://www.uniprot.org/uniprot/$i.fasta
done

cat P99999.fasta P00004.fasta P0C0X8.fasta P00091.fasta Q93VA3.fasta > cytc_aln.clw
