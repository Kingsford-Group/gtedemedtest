from datetime import datetime
import random
from GTED_MSA_solver import *
import pickle
import gurobipy as gp
from gurobipy import GRB
import Levenshtein
from tqdm import trange
import seaborn as sns
import matplotlib.pyplot as plt
import os,sys

sys.setrecursionlimit(5500)

def sample_path(u,adj_list, path,width,sink):
    path.append(u)
    if u == sink:
#         total_paths.append(path)
        for i in range(len(path)-1):
            adj_list[path[i]][path[i+1]] -= width
        return path,width
    
    for i in range(len(path) - 1):
        n1 = path[i]
        n2 = path[i+1]
        adj_list[n1][n2] -= width
    
    neighbors = [e for e in adj_list[u] if adj_list[u][e] > 0]
    neighbors = sorted(neighbors, key=lambda x: adj_list[u][x])
    
    v = 0
    if len(neighbors) > 0:
        v = neighbors[-1]
        if random.choice([0,1]) == 0:
            v = random.choice(neighbors)
        flow = adj_list[u][v]
    
        for i in range(len(path) - 1):
            n1 = path[i]
            n2 = path[i+1]
            adj_list[n1][n2] += width

        width = min(width, adj_list[u][v])    

    try:
        return sample_path(v,adj_list, path, width, sink)
    except BaseException as err:
        print(err)
        return []

def sample_all_paths(source, adj_list, sink):
    total_paths = []
    times = 0
    while True:
        if sum(adj_list[source].values()) > 0:
            path = sample_path(source, adj_list, [], 1000, sink)
            if len(path) > 0:
                total_paths.append(path)
        if times > 1000 or sum([n[1] for n in total_paths]) == 100:
            return total_paths
        times += 1
    return total_paths

def read_msa_graph(fname):
    with open(fname, "rb") as ifile:
        return pickle.load(ifile)

def get_msa_adj_list(graph):
    adj_list = defaultdict(dict)
    incoming = defaultdict(list)
    outgoing = defaultdict(list)
    
    for u in graph.adj_list:
        for v in graph.adj_list[u]:
            incoming[v].append(u)
            outgoing[u].append(v)
            adj_list[u][v] = graph.adj_list[u][v]
    return adj_list, incoming, outgoing
'''
Writes LP to file specified
fname -- file that stores LP
dist_mat -- edit distance between s1 from set 1, s2 from set 2
set_1, set_2 -- list of tuples(sname, weight)
'''
def write_EMD_LP(fname, set_1, set_2):
    out_file = open(fname, 'w')
    
    # objective
    out_file.write("Minimize\n")
    out_str=""
    for i,s1 in enumerate(set_1[1]):
        for j,s2 in enumerate(set_2[1]):
            dist = Levenshtein.distance(s1, s2)
            out_str += " + " + str(dist) + " x" + str(i) + "_" + str(j)
    out_file.write(out_str[3:]+"\n")
    
    # constraints
    out_file.write("Subject To\n")
    
    # for all s1
    for i,s1 in enumerate(set_1[0]):
        out_str = "v"+str(i)+": "
        curr_len = len(out_str)
        for j,s2 in enumerate(set_2[0]):
            out_str += " + " + " x" + str(i) + "_" + str(j)
        out_str += " = " +  str(s1)
        out_file.write(out_str[:curr_len] + out_str[curr_len+3:]+"\n")
    
    # for all s2
    for j,s2 in enumerate(set_2[0]):
        out_str = "w"+str(j)+": "
        curr_len = len(out_str)
        for i,s1 in enumerate(set_1[0]):
            out_str += " + " + " x" + str(i) + "_" + str(j)
        out_str += " = " + str(s2)
        out_file.write(out_str[:curr_len] + out_str[curr_len+3:]+"\n")
    
    out_file.write("End")
    out_file.close

def sample_diameter(graph_path):
    #graph1 = read_msa_graph("MSA/graphs/constructed_set_%i.msa.graph" % graph_idx)
    graph1 = read_msa_graph(graph_path)
    source = [n for n in graph1.nodes if graph1.nodes[n] == "$"][0]
    sink = [n for n in graph1.nodes if graph1.nodes[n] == "#"][0]
    max_dist = 0

    for counter in range(50):
        now = datetime.now()
        curr_time = now.strftime("%H:%M:%S")
        print(graph_path, counter, curr_time, flush=True)
        
        adj_list, incoming, outgoing = get_msa_adj_list(graph1)
        total_paths1 = sample_all_paths(source, adj_list,sink)
        total_satisfied = sum([n[1] for n in total_paths1])

        adj_list, incoming, outgoing = get_msa_adj_list(graph1)
        total_paths2 = sample_all_paths(source, adj_list,sink)
        total_satisfied = sum([n[1] for n in total_paths2])
        strings1 = [[n[1] for n in total_paths1], [path_to_string(graph1,total_paths1[i][0]) for i in range(len(total_paths1))]]
        strings2 = [[n[1] for n in total_paths2], [path_to_string(graph1,total_paths2[i][0]) for i in range(len(total_paths2))]]
        
        graph_idx = graph_path.split("/")[-1].split(".")[0]
        fname = output_emd_lp+"curr_emed_lp_%s.lp"%graph_idx

        write_EMD_LP(fname, strings1, strings2)

        gp.setParam("LogToConsole", 0)
        gp.setParam("OutputFlag", 0)
        model = gp.read(fname)
        model.optimize()

        dist = model.objVal
        max_dist = max(max_dist, dist)

    return max_dist

def path_to_string(graph,path):
    return "".join([graph.nodes[i] for i in path])


if __name__ == "__main__":
    
    if len(sys.argv) < :
        print("!!! USAGE: python3.6 diameter_flow_decomp.py <all_graphs.txt> <lp_output> <graph_dir>")
        print("all_graphs.txt --- txt file that contains all graph names, one per line")
        print("lp_output --- output directory for .lp files from gurobi")
        print("graph_dir --- directory that contains all the graphs; also output directory for the final sampled_diameter.txt file")
        exit()

    all_graphs = sys.argv[1]
    output_emd_lp = sys.argv[2]
    graph_dir = sys.argv[3]

    files = [l.rstrip("\n") for l in open(all_graphs)]
    print(files)
    suff = all_graphs.split(".")[0].split("_")[-1]
    sampled_diameter_file = output_emd_lp+"sampled_diameters_%s.txt" % suff
    print(sampled_diameter_file)

    with open(sampled_diameter_file,"w") as ofile:
        ofile.write("Graph idx\tDiameter\n")

        #for graph_path in [f for f in os.listdir(graph_dir) if "msa.graph" in f and "pan" not in f and "hylobates" not in f] :
        for graph_path in files:
            diameter = sample_diameter(graph_dir + graph_path)
            print(diameter, flush=True)
            ofile.write(graph_path.split(".")[0] + "\t" + str(diameter) + "\n")
            ofile.flush()
