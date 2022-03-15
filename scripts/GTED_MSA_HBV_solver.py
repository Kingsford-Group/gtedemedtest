####
# Main function for running FGTED on MSA graphs built from HBV sequences
####

from Bio import SeqIO
import os
import Bio
import numpy as np
import random
import pickle
from collections import defaultdict
from GTED_LP_solver_2 import *
from msa_graph import *
import sys


if __name__ == "__main__":
    # seq_dir = "simulated/1018_fewer_parts/corrected/"

    if len(sys.argv) < 4:
        print("!!!! If you do not know how to deal with the output, try running pipeline_hbv.sh to run everything")
        print("!!!! USAGE: python3.6 GTED_MSA_HBV_solver.py <pairs.txt> <msa_graph_dir> <gted_msa_dir>")
        print("pairs.txt --- txt file that contains comma-delimited filename for each pair of msa graphs")
        print("msa_graph_dir --- Graph directory that contains all msa files")
        print("gted_msa_dir --- output directory for all the .lp and .sol files from gurobi")
        exit()

    pairs_file = sys.argv[1]
    msa_graph_dir = sys.argv[2]
    gted_msa_dir = sys.argv[3]

    gp.setParam("Threads", 12)
    gp.setParam("LogToConsole", 1)
    gp.setParam("OutputFlag", 1)
    gp.setParam("TimeLimit",12000)

    pairs = [l.rstrip("\n").split(",") for l in open(pairs_file)]
    for p in pairs:
        name1 = p[0]
        name2 = p[1]

        prefix1 = name1.split(".")[0]
        prefix2 = name2.split(".")[0]
        fname = "GTED_MSA_%s_%s" % (prefix1,prefix2)
        print(fname)

        graph1 = read_msa_graph(msa_graph_dir+name1)
        graph2 = read_msa_graph(msa_graph_dir+name2)

        aln_graph = Alignment_Graph(graph1, graph2)

        write_LP_exact(gted_msa_dir + fname + ".lp", aln_graph, graph1, graph2)

        model = gp.read(gted_msa_dir + fname + ".lp")

        model._objfname = gted_msa_dir + fname + ".obj"
        model._solfname = gted_msa_dir + fname + ".sol"

        model.optimize(callback_lp)

        if model.status == GRB.INF_OR_UNBD:
            # Turn presolve off to determine whether model is infeasible or unbounded
            model.setParam(GRB.Param.Presolve, 0)
            model.optimize()

        if model.status == GRB.OPTIMAL:
            with open(model._objfname,"w") as obj_file:
                obj_file.write(str(model.objVal))
                model.write(gted_msa_dir + fname + ".sol")
