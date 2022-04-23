import random
from Bio import SeqIO
import Bio
from collections import defaultdict
import pandas as pd
import numpy as np
from tqdm import tqdm

def read_seq(directory, fname):
    '''
    Read reference TRB genes
    '''
    # remove gaps from the sequences
    seq_list = [str(seq.seq).replace(".","") for seq in SeqIO.parse(directory+"/"+fname,"fasta")]
    random.shuffle(seq_list)
    return seq_list

def sample_gene(V, J, D):
    '''
    Sample one sequence from gene repository
    '''
    v = random.choice(V)
    j = random.choice(J)
    d = random.choice(D)
    return [v,j,d]

def mutate(seq, mut):
    '''
    Mutate mut number of loci in sequences
    '''
    for m in range(mut):
        idx = random.choice(range(3))
        loc = random.choice(range(len(seq[idx])))
        a = list(seq[idx])
        char = random.choice([c for c in ["a","t","c","g"] if c != seq[idx][loc]])
        act = random.choice(range(3)) # 0 - sub, 1 - ins, 2 - del
        if act == 0:
            a[loc] = char
        if act == 1:
            a.insert(loc, char)
        if act == 2:
            a = a[:loc] + a[loc+1:]
        seq[idx] = "".join(a)
    return seq

def random_partition(total,num_part):
    '''
    Randomly split 100
    '''
    parts = sorted(random.sample(range(total),k=num_part-1))+[total]
    sizes = [b-a for a,b in zip([0]+parts[:-1], parts)]
    return sizes

if __name__ == "__main__":
    base_dir = "/home/yutongq/gtedemedtest/gtedemedtest/data/"
    seq_dir = base_dir + "TCR_ref_seq/"
    target_dir = base_dir +"TCR_simulated/"

    TRB_fnames = ["TRBV.fasta", "TRBJ.fasta", "TRBD.fasta"]

    TRBV, TRBJ, TRBD = [read_seq(seq_dir, f) for f in TRB_fnames]

    # get number of ref seqs for each gene group
    TRBV_sizes = list(range(5, len(TRBV), len(TRBV)//5))
    TRBJ_sizes = list(range(5, len(TRBJ), len(TRBJ)//5))+[len(TRBJ)-1]

    all_seqs = defaultdict(set)

    # Generate all sequences in the gene group
    for stop in range(5):
        repo_TRBV = random.sample(TRBV, k=TRBV_sizes[stop])
        repo_TRBJ = random.sample(TRBJ, k=TRBJ_sizes[stop])
        repo_TRBD = TRBD
        for num in range(1,100):
            for mut in [1,3,5,8,10]:
                while True:
                    seq = sample_gene(repo_TRBV, repo_TRBJ, repo_TRBD)
                    seq = "".join(mutate(seq, mut))
                    if seq not in all_seqs[stop]:
                        all_seqs[stop].add(seq)
                        break

    # Write generated sequences to file 
    for i, stop in enumerate(all_seqs):
        seq_records = []
        for idx,s in enumerate(all_seqs[stop]):
            record = SeqIO.SeqRecord(Bio.Seq.Seq(s), "simulated_%i" % (idx+1), "", "")
            seq_records.append(record)
        SeqIO.write(seq_records, target_dir+"/all_seqs_"+str(i)+".fasta", "fasta")
    
    all_seqs = [read_seq(target_dir, "all_seqs_"+str(i)+".fasta") for i in range(5)]
    
    # Select and construct sequence sets
    all_seq_sets = []
    for i in range(5):
        for num_part in range(2,12,1):
            part_sizes = random_partition(100, num_part)
            seq_set = []
            for size in part_sizes:
                group_id = i % 5
                seq_id= random.sample(range(len(all_seqs[group_id])),k=1)[0]
                seq = all_seqs[group_id][seq_id]
                seq_set.append([group_id, seq_id, size, seq])
            all_seq_sets.append(seq_set)
    
    # Write constructed sequence set to file
    for i,seq_set in enumerate(all_seq_sets):
        seq_records = []
        for idx, seq in enumerate(seq_set):
            record = SeqIO.SeqRecord(Bio.Seq.Seq(seq[3]), "simulated_seq_%i_%i_%i" %(seq[0], seq[1], seq[2]), "", "")
            seq_records.append(record)
        SeqIO.write(seq_records, target_dir+"/constructed_set_"+str(i)+".fasta","fasta")


