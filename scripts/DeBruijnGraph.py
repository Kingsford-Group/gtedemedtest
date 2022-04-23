from collections import defaultdict
import pickle
import sys, os
from Bio import SeqIO

def read_fasta(fname, get_seq = False):
    set_1 = [("_".join(seq.id.split("_")[:-1]), seq.id.split("_")[-1]) for seq in SeqIO.parse(fname,"fasta")]
    seqs = []
    if get_seq:
        seqs = [str(seq.seq) for seq in SeqIO.parse(fname,"fasta")]
        return set_1, seqs
    return set_1

def read_dbg_graph(fname):
    with open(fname, "rb") as ifile:
        return pickle.load(ifile)

class DeBruijnGraph(object):

    def __init__(self, seqid, seqs, k):

        self.nodes = dict({0:"$",1:"#"})                  # key: u1, value = label
        self.adj_list = defaultdict(dict)     # key: u1, value: [dict((u2, weight), (u3, weight))]
        self.edge_to_idx = dict()             # key: (u,v), value: idx
        self.idx_to_edge = dict()             # key: idx, value: (u,v, weight)
        self.paths = [[0] for i in range(len(seqs))]
        self.seqs = seqs
        self.seqid = seqid

        self.construct_kmers(seqid, seqs,k)
        self.add_edges([1,0], sum(self.adj_list[0].values()))


    def add_edges(self, node_list, weight):
        for i in range(len(node_list)-1):
            edge_idx = len(self.idx_to_edge)
            u = node_list[i]
            v = node_list[i+1]
            edge = (u,v,weight)

            if v not in self.adj_list[u]:
                self.adj_list[u][v] = weight
                self.edge_to_idx[(u,v)] = edge_idx
                self.idx_to_edge[edge_idx] = edge
            else:
                self.adj_list[u][v] += weight
                edge_idx = self.edge_to_idx[(u,v)]
                self.idx_to_edge[edge_idx] = (u,v,self.adj_list[u][v])


    def construct_kmers(self, seqid, seq, k):
        self.kmers = dict()  # map string to [[node_idx], source]

        pathid = 0
        for idx, s in zip(seqid, seq):
            prev_node = 0
            weight = int(idx[1])
            for i in range(len(s) - k + 1):
                source = (i == 0)
                substr = s[i:i+k]
                node_idx = []

                if substr in self.kmers:
                    node_idx = self.kmers[substr][0]
                    # add nodes if kmer becomes a source
                    if self.kmers[substr][1] == False and source == True:
                        self.kmers[substr][1] = True
                        node_idx = [len(self.nodes)+ j for j in range(k-1)] + node_idx
                        self.kmers[substr][0] = node_idx
                        for j in range(k-1):
                            self.nodes[node_idx[j]] = substr[j]

                # kmer to set if new kmer
                else:
                    if source:
                        node_idx = [len(self.nodes) + j for j in range(k)]
                        for j in range(k):
                            self.nodes[node_idx[j]] = substr[j]
                    else:
                        node_idx = [len(self.nodes)]
                        self.nodes[len(self.nodes)] = substr[-1]
                    self.kmers[substr] = [node_idx, source]

                # add edge
                node_list = [prev_node] + node_idx
                if not source:
                    node_list = [prev_node, node_idx[-1]]
                self.add_edges(node_list, weight)
                for n in node_list[1:]:
                    self.paths[pathid].append(n)
                
                prev_node = node_idx[-1]

            self.add_edges([prev_node, 1], weight)
            self.paths[pathid].append(1)
            pathid += 1


def test():
    seq = ["ATCATC", "TCAT"]
    seqid = [["1", "3"], ["2", "2"]]

    graph = DeBruijnGraph(seqid, seq, 3)
    print(graph.kmers)
    print(graph.nodes)
    print(graph.adj_list)
    print(graph.edge_to_idx)
    print(graph.idx_to_edge)
    print(graph.paths)

if __name__ == "__main__":

    # test()
    
    if len(sys.argv) < 4:
        print("USAGE: python3.6 %s <fasta_dir> <k-mer size> <graph output dir>" % sys.argv[0])
        exit()

    fasta_dir = sys.argv[1]
    k = int(sys.argv[2])
    target_dir = sys.argv[3]

    for f in os.listdir(fasta_dir):
        if "fasta" in f:
            prefix=f.split(".")[0]
            print(f, prefix)
            seq_id, seqs = read_fasta(fasta_dir + "/" + f, True)
            graph = DeBruijnGraph(seq_id, seqs, k)
            print(len(graph.nodes), len(graph.idx_to_edge))
            with open(target_dir+"/"+prefix+"."+str(k)+".dbg.graph", "wb") as handle:
                pickle.dump(graph, handle)
             
