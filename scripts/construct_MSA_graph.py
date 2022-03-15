from msa_graph import *
import sys
import pickle

msa_fname = sys.argv[1]
prefix = msa_fname.split(".")[0]

seqid, seqs = read_msa(msa_fname)
graph = MSA_Graph(seqid, seqs)

print(len(graph.nodes))
print(prefix+".msa.graph")
with open(prefix+".msa.graph", "wb") as handle:
    pickle.dump(graph, handle)
