####
# Main function for running FGTED on MSA graphs built from TCR sequences
####

from Bio import SeqIO
import Bio
import numpy as np
import random
import pickle
from collections import defaultdict
from GTED_LP_solver_2 import *
import sys

if __name__ == "__main__":

    if len(sys.argv) < 4:
        print("!!!! If you do not know how to deal with the output, try running pipeline_tcr.sh to run everything")
        print("!!!! USAGE: python3.6 GTED_MSA_TCR_solver.py <num1> <num2> <msa_graph_dir> <gted_msa_dir>")
        print("num1 --- starting index of file prefix (inclusive)")
        print("num2 --- ending index of file prefix (exclusive)")
        print("msa_graph_dir --- Graph directory that contains all msa files")
        print("gted_msa_dir --- output directory for all the .lp and .sol files from gurobi")
        exit()

    msa_graph_dir = "MSA/graphs/"
    gted_msa_dir = "GTED_MSA_LP_solver_file/1030/"

    gp.setParam("Threads", 12)
    gp.setParam("LogToConsole", 1)
    gp.setParam("OutputFlag", 1)
    gp.setParam("TimeLimit",6100)

    num1 = int(sys.argv[1])
    num2 = int(sys.argv[2])

    for i in range(num1, num2):
        for j in range(i+1, num2):
            fname = "GTED_MSA_%i_%i" % (i,j)
            print(fname)

            graph1 = read_msa_graph(msa_graph_dir+"constructed_set_%i.msa.graph" % i)
            graph2 = read_msa_graph(msa_graph_dir+"constructed_set_%i.msa.graph" % j)

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
