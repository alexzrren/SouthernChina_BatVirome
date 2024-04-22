cat viralfamilies2abbr.tsv | while read -r vf vfabbr
do
echo $vf $vfabbr
#ls -lh 1.MSA/$vf.afa
#ls -lh 2.HMM_align/$vf.afa.hmm
#ls -lh 2.HMM_align/${vfabbr}_extracted_sequences.fasta
echo "hmmalign --trim --mapali 1.MSA/$vf.afa --outformat Afa -o 3.aligned_extracted_sequences_trim_merged/aligned_${vfabbr}_extracted_sequences_trim_merged.fasta 2.HMM_align/${vf}.afa.hmm 2.HMM_align/${vfabbr}_extracted_sequences.fasta" > haln_$vfabbr.sh
qsub -clear -cwd -q st.q -binding linear:1 -P P20Z10200N0206 -l vf=16g,p=1 haln_$vfabbr.sh
done
