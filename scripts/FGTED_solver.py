import os
import Bio
import numpy as np
import pickle
# from GTED_LP_solver_2 import *
from FGTED_utils import *
from msa_graph import *
from DeBruijnGraph import *
import time
import sys

def print_help(flag=0):
    # general help
    if flag==0:
        print("!!!! Run this script to output FGTED for pairs of graphs specified")
        print("!!!! USAGE: python3.6 FGTED_solver.py <pairs.txt> <msa_graph_dir> <outprefix>")
        print("pairs.txt --- txt file that contains comma-delimited filename for each pair of msa graphs")
        print("graph_dir --- Graph directory that contains all graoh files")
        print("outprefix -- directory/prefix.csv will store pairwise FGTED")
    # graph type error
    if flag==1:
        print("!!! Check the graph type/extension: it should either be *.dbg.graph or *.msa.graph !!!")

    exit()
    
def read_graph(fname):
    try:
        with open(fname, "rb") as ifile:
            return pickle.load(ifile)
    except:
        print_help(1)
        exit()

if __name__ == "__main__":

    if len(sys.argv) < 4:
        print_help()
    
    pairs_file = sys.argv[1]
    graph_dir = sys.argv[2]
    outprefix = sys.argv[3]

    ofilename = outprefix + ".csv"
    
    pairs = [l.rstrip("\n").split(",") for l in open(pairs_file)]

    gp.setParam("LogToConsole", 0)
    gp.setParam("OutputFlag", 0)
    gp.setParam("Threads", 12)
    gp.setParam("TimeLimit",12000)

    with open(ofilename, "w") as ofile:
        ofile.write(",".join(["idx1", "idx2", "FGTED", "time_alnGraph", "time_FGTED"])+"\n")

    for p in pairs:
        name1 = p[0]
        name2 = p[1]

        prefix1 = ".".join(name1.split(".")[:-1])
        prefix2 = ".".join(name2.split(".")[:-1])

        lp_fname = "FGTED_%s_%s.lp" % (prefix1, prefix2)

        print("Computing FGTED between")
        print("---", prefix1,prefix2)

        graph1 = read_graph(graph_dir+name1)
        graph2 = read_graph(graph_dir+name2)
        
        print("Constructing Alignment graph...")
        start = time.time()
        aln_graph = Alignment_Graph(graph1, graph2)
        end = time.time()
        t1 = end - start

        print("Computing FGTED...")
        start = time.time()
        objValue = FGTED_lp(aln_graph,graph1, graph2,lp_fname)
        end = time.time()
        t2 = end - start

        if objValue == -1:
            print("!!! ERROR: LP model did not return a solution !!!")
            exit()
        
        with open(ofilename, "a") as ofile:
            ofile.write(",".join([prefix1, prefix2,str(objValue/100),str(t1),str(t2)])+"\n")
