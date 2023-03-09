1. FLT3's protein tyrosine kinase domain (610-943) was downloaded directly from NCBI as FASTA.

2. Deltablast:
deltablast -db swissprot -query flt3_610_943.fa -out flt3_610_943.deltablast.2000.out -outfmt 7 -remote -max_target_seqs 2000 

3. Getting a list of hits (filtering out those with >95% identity to the query):
grep -v ^# flt3_610_943.deltablast.2000.out | awk '$3<95{print $2}' | awk -F\| '{print $NF}' > flt3_610_943.deltablast.2000.hits

4.  Downloading the hits from swissprot
swissprot_batch_retrieve.pl flt3_610_943.deltablast.2000.hits  > flt3_610_943.deltablast.2000.hits.fa

5.  Removing similar (95%) seqs:
skipredundant -sequences  flt3_610_943.deltablast.2000.hits.fa  -mode 1 -outseq flt3_610_943.deltablast.2000.hits.filt.fa -auto

6. Preparing the input file:
cp flt3_610_943.fa flt3_610_943.fa_and_deltablast_hits_filt_2000.fa
cat flt3_610_943.deltablast.2000.hits.filt.fa  >>  flt3_610_943.fa_and_deltablast_hits_filt_2000.fa

7.  Running mafft:
mafft --clustalout flt3_610_943.fa_and_deltablast_hits_filt_2000.fa > flt3_610_943_and_deltablast_hits_filt_2000.aln

8. Testing the variations:
~/scripts/MolEvol/get_seq_pos_variation.pl --alnfile flt3_610_943_and_deltablast_hits_filt_2000.aln --shift 610 --mut FLT3_resistance_mutations

To gain some more info into closely related seqs (whole protein):
9. deltablast -db swissprot -query NP_004110.2.fa -out flt3.dblast.50.out -outfmt 7 -remote -max_target_seqs 50

10. grep -v ^# flt3.dblast.50.out | awk '$3<95{print $2}' | awk -F\| '{print $NF}' > flt3.dblast.50.hits

11. swissprot_batch_retrieve.pl flt3.dblast.50.hits  > flt3.dblast.50.hits.fa

12. skipredundant -sequences flt3.dblast.50.hits.fa  -mode 1 -outseq flt3.dblast.50.hits.filt.fa -auto

13. cp NP_004110.2.fa flt3_and_dblast.50.hits.filt.fa ; cat flt3.dblast.50.hits.filt.fa >> flt3_and_dblast.50.hits.filt.fa

14. mafft --auto --clustalout flt3_and_dblast.50.hits.filt.fa >> flt3_and_dblast.50.hits.filt.aln

15. Making a blast DB for hits:
makeblastdb -in flt3_610_943.deltablast.2000.hits.filt.fa -parse_seqids -dbtype prot

16. Generating a PSSM:
deltablast -db flt3_610_943.deltablast.2000.hits.filt.fa -query flt3_610_943.fa -out flt3_db.out -outfmt 7 -out_pssm flt3_2000.pssm -out_ascii_pssm flt3_2000_pssm.txt

17. Making the PSSM according to the real seq:
perl -ne 'if (/^\s+(\d+)(\s+\S+.*)/){printf "%5d%s\n",$1+609,$2}else{print}' flt3_2000_pssm.txt > flt3_610_943.pssm.txt

18. Extracting residues from the PSSM (e.g., residue 842):
awk '$1=="A"||$1==842' flt3_610_943.pssm.txt



