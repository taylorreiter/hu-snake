# downloaded pfam profiles from pfam. 
# build hmm
hmmbuild PF00204_full_gyrb.hmm PF00204_full_gyrb.txt

# search plass amino acid sequences for matches
hmmsearch PF00204_full_gyrb.hmm plass-gyra.fas > plass-gyra.out
hmmsearch globins4.hmm uniprotsprot.fasta > globins4.out

hmmbuild PF00521_full_gyra.hmm PF00521_full_gyra.txt
hmmsearch PF00521_full_gyra.hmm plass-gyra.fas > plass-gyra.out


# do on all bin prots, and on all sb1 plass prots
ln -s ../../../inputs/hu-genomes-plass100/sb1/all_sb1.nostop.c100.fas .
hmmpress PF00521_full_gyra.hmm

hmmscan --tblout plass-gyra-seq-hmmscan.out --domtblout plass-gyra-dom-hmmscan.out --pfamtblout plass-gyra-pfam-hmmscan.out PF00521_full_gyra.hmm all_sb1.nostop.c100.fas

ln -s ~/github/hu-snake/outputs/hu-bins/prokka/all-bins.faa .

hmmscan --tblout bin-gyra-seq-hmmscan.out --domtblout bin-gyra-dom-hmmscan.out --pfamtblout bin-gyra-pfam-hmmscan.out PF00521_full_gyra.hmm all-bins.faa

hmmpress PF00204_full_gyrb.hmm

hmmscan --tblout bin-gyrb-seq-hmmscan.out --domtblout bin-gyrb-dom-hmmscan.out --pfamtblout bin-gyrb-pfam-hmmscan.out PF00204_full_gyrb.hmm all-bins.faa

hmmscan --tblout plass-gyrb-seq-hmmscan.out --domtblout plass-gyrb-dom-hmmscan.out --pfamtblout plass-gyrb-pfam-hmmscan.out PF00204_full_gyrb.hmm all_sb1.nostop.c100.fas
