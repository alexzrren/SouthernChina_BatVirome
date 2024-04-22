import Bio.SeqIO as SeqIO
import sys
import numpy as np

def calc_aai(seq1, seq2, gap_rate=0.8):
    if len(seq1) != len(seq2):
        raise Exception('Sequence not aligned properly!')
    
    length = len(seq1)
    
    # Convert sequences to NumPy arrays for vectorized operations
    seq1_array = np.array(list(seq1))
    seq2_array = np.array(list(seq2))
    
    # Identify positions where neither seq1 nor seq2 has a gap
    non_gap_positions = (seq1_array != '-') & (seq2_array != '-')
    
    # Count matches at non-gap positions
    match = np.sum(seq1_array[non_gap_positions] == seq2_array[non_gap_positions])
    
    # Count gaps
    gap = length - np.count_nonzero(non_gap_positions)
    
    try:
        pident = match / (length - gap) * 100
    except ZeroDivisionError:
        pident = 'NA'
    
    return match, (length - gap - match), gap, length, pident
    
    
if __name__=='__main__':
    assert len(sys.argv) == 3, "Invalid length of input parameters! 2 is required but %d is given" % (len(sys.argv) - 1)
    inputpath, outputpath = sys.argv[1:]
    vsname_abbr = inputpath.split('/')[-1].split('.')[0]
    
    seqdict = SeqIO.to_dict(SeqIO.parse(inputpath, format='fasta'))

    seqids = list(seqdict.keys())
    
    result_list = []
    with open(outputpath, 'w') as fdout:
        fdout.write('sn1,sn2,mismatch,length,vsname_abbr\n')
        for i in range(len(seqdict)):
            seqid_i = seqids[i]
            for j in range(i+1, len(seqdict)):
                seqid_j = seqids[j]
                seq_i = seqdict[seqid_i].seq
                seq_j = seqdict[seqid_j].seq
                match, mismatch, gap, length, pident = calc_aai(seq_i, seq_j, gap_rate=0.8)
                fdout.write(f"{seqid_i},{seqid_j},{mismatch},{length},{vsname_abbr}\n")
