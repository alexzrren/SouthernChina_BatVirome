import Bio.SeqIO as SeqIO
import numpy as np
import pandas as pd
from tqdm import tqdm
import sys

inpath, outpath = sys.argv[1:]

aln = SeqIO.to_dict(list(SeqIO.parse(inpath, format='fasta')))
aln

# %%
aln_cleaned = { seqid: ''.join([ char for char in list(seqrecord.seq) if (char == '-') or (char.isupper()) ]) for seqid, seqrecord in aln.items() }

# %%
def calc_aai(seq1, seq2, gap_rate=0.8):
    if len(seq1) != len(seq2):
        raise Exception('Sequence not aligned properly!')
    else:
        seq1_np = np.array(list(seq1))
        seq2_np = np.array(list(seq2))

        valid_positions = ~((seq1_np == '-') & (seq2_np == '-') & (seq1_np == seq2_np))
        match_count = np.sum(seq1_np[valid_positions] == seq2_np[valid_positions])
        gap_count = np.sum((seq1_np[valid_positions] == '-') | (seq2_np[valid_positions] == '-'))
        aligned_length = np.sum(valid_positions)
        
        try:
            pident = match_count / (aligned_length) * 100
        except ZeroDivisionError:
            pident = 'NA'
        
        return match_count, (aligned_length - match_count), gap_count, aligned_length, len(seq1), pident

# %%
result = []
for i in tqdm(range(len(aln_cleaned))):
    for j in range(i+1, len(aln_cleaned)):
        seqid_i = list(aln_cleaned.keys())[i]
        seqid_j = list(aln_cleaned.keys())[j]
        seq_i = aln_cleaned[seqid_i]
        seq_j = aln_cleaned[seqid_j]
        match, mismatch, gap, alnlen, msalen, pident = calc_aai(seq_i, seq_j, gap_rate=0.8)
        result.append([seqid_i, seqid_j, match, mismatch, gap, alnlen, msalen, pident])

# %%
resultdf = pd.DataFrame(result, columns=['seq1', 'seq2', 'match', 'mismatch', 'gap', 'alnlen', 'msalen', 'pident'])
resultdf.to_csv(outpath, index=False)


