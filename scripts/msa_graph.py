from collections import defaultdict
from Bio import SeqIO
import numpy as np

def read_msa(fname):
    records = SeqIO.parse(fname,"fasta")
    seqid, seqs = np.array([(s.id, str(s.seq).upper()) for s in records if s.id.split("_")[-1]!="0"]).T
    return seqid, seqs

def read_msa_graph(fname):
    with open(fname, "rb") as ifile:
        return pickle.load(ifile)

class MSA_Graph:
    
    def __init__(self, seqid, seqs):
        
        self.nodes = dict({0:"$",1:"#"})                  # key: u1, value = label
        self.adj_list = defaultdict(dict)    # key: u1, value = [dict((u2, weight), (u3, weight))]
        self.edge_to_idx = dict()             # key: (u,v), value = idx
        self.idx_to_edge = dict()             # key: idx, value = (u,v)
        self.paths = [[0] for i in range(len(seqs))]
        self.seqs = seqs
        self.seqid = seqid
        
        self.construct_paths()
        self.add_edges()
        self.add_edge(1, 0, sum(self.adj_list[0].values()))
            
    def add_edge(self, u,v,weight):
        edge_idx = len(self.edge_to_idx)
        if v not in self.adj_list[u]:
            self.adj_list[u][v] = weight
            self.edge_to_idx[(u,v)] = edge_idx
            self.idx_to_edge[edge_idx] = (u,v, weight)
        else:
            self.adj_list[u][v] += weight        
            curr_edge = self.idx_to_edge[self.edge_to_idx[(u,v)]]
            edge_idx = self.edge_to_idx[(u,v)]
            curr_edge = (curr_edge[0], curr_edge[1], curr_edge[2] + weight)
            self.idx_to_edge[edge_idx] = curr_edge
            
    def add_edges(self):
        for idx,path in enumerate(self.paths):
            weight = int(self.seqid[idx].split("_")[-1])
            for i in range(len(path)-1):
                u = path[i]
                v = path[i+1]
                self.add_edge(u,v,weight) 
    
    def construct_paths(self):
        for i in range(len(self.seqs[0])):
            curr_nodes = dict()
            for j in range(len(self.seqs)):
                c = self.seqs[j][i]
                idx = 0
                if c != "-":
                    if c not in curr_nodes:
                        curr_nodes[c] = len(self.nodes)
                        self.nodes[len(self.nodes)] = c
                    idx = curr_nodes[c]
                    self.paths[j].append(idx)
        for p in self.paths:
            p.append(1)