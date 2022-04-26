import time
import random
import pickle
from DeBruijnGraph import *
from msa_graph import *
import gurobipy as gp
from gurobipy import GRB
from collections import Counter
import Levenshtein
from tqdm import trange
import os, sys

def EMED_lp(set_1, set_2, ofilename):
    '''
    Objective: min cost(e) w(e)
    Each e is a pair of string from set 1 and set 2 respectively
    (min cost max flow)

    Input: set_1, set_2 --- [[weights], [strings]]
    '''

    ofile = open(ofilename, "w")

    try:
        
        # create a new model
        model = gp.Model("EMED")

        # create variable sets
        X = [(i,j) for i in range(len(set_1[1])) for j in range(len(set_2[1]))]
        x = model.addVars(X, vtype=GRB.CONTINUOUS, name="x")

        # find pair-wise distance
        dist = dict([((i,j),Levenshtein.distance(set_1[1][i], set_2[1][j])) for (i,j) in X])



        # set objective
        model.setObjective(sum([dist[(i,j)]* x[(i,j)] for (i,j) in X], GRB.MINIMIZE))

        # constraints
        for i,w1 in enumerate(set_1[0]):
            model.addConstr(sum([x[(i,j)] for j in range(len(set_2[0]))]) == w1)
        for j,w2 in enumerate(set_2[0]):
            model.addConstr(sum([x[(i,j)] for i in range(len(set_1[0]))]) == w2)

        model.optimize()

        if model.status == GRB.OPTIMAL:
            return model.objVal

        if model.status == GRB.INFEASIBLE:
            print("!!! INFEASIBLE MODEL?")
            model.write(ofilename)
            return -1

        if model.status == GRB.TIME_LIMIT:
            print("!!! TIME LIMIT EXCEEDS")
            model.write(ofilename)
            return -1
    
    except gp.GurobiError as e:
        sys.exit('!!! GurobiError: Error code ' + str(e.errno) + ': ' + str(e))   
        
def process_path(path, graph):
    '''
    Convert a path into corresponding string
    path -- a list of node ids
    graph -- the original input graph
    '''
    #path = [str(graph.nodes[path[i][1]]) for i in range(len(path)-1,-1,-1)]
    # get the correct ordering
    path = [path[i][1] for i in range(len(path)-1,-1,-1)]
    strings = [str(graph.nodes[i]) for i in path]

    #print(path)
    string_sets = Counter()
    current_string = ""
    for char in strings:
       if char == "#":
           string_sets[current_string[1:]] += 1
           current_string = ""
       else:
           current_string += char
    if current_string != "":
        string_sets[current_string[1:]] += 1
    return string_sets


def Eulerian_path(G, source):
    '''
    Find an Eulerian cycle in the graph starting from source provided
    '''

    # get degree and shuffle edges in the adjacency list
    degree = defaultdict(int)
    edges = defaultdict(list)
    for u in G.adj_list:
        for v in G.adj_list[u]:
            degree[u] += G.adj_list[u][v]
            for i in range(G.adj_list[u][v]):
                edges[u].append((u,v,i))
        random.shuffle(edges[u])

    vertex_stack = [(source, None)]
    last_vertex = None
    last_key = None

    output_stack = []
    while vertex_stack:
        current_vertex, current_key = vertex_stack[-1]
        if degree[current_vertex] == 0:
            if last_vertex is not None:
                yield (last_vertex, current_vertex, last_key)
            last_vertex, last_key = current_vertex, current_key
            vertex_stack.pop()
        else:
            triple = edges[current_vertex][-1]
            edges[current_vertex].pop()
            _, next_vertex, next_key = triple
            vertex_stack.append((next_vertex, next_key))
            degree[current_vertex] -= 1

def sample_eulerian_cycles(graph, number=25):
    all_sets = []
    for i in range(number):
        path = list(Eulerian_path(graph, 1))
        string_sets = process_path(path, graph)

        strings = [s for s in string_sets]
        weights = [string_sets[s] for s in strings]
        all_sets.append([weights, strings])
    return all_sets


if __name__ == "__main__":

    if len(sys.argv) < 4:
        print("!!! USAGE: python3.6 diameter_flow_decomp.py <all_graphs.txt> <lp_output> <graph_dir>")
        print("all_graphs.txt --- txt file that contains all graph names, one per line")
        print("lp_output --- output directory for .lp files from gurobi if error occurs")
        print("graph_dir --- directory that contains all the graphs; also output directory for the final sampled_diameter.txt file")
        exit()

    gp.setParam("LogToConsole", 0)
    gp.setParam("OutputFlag", 0)
    gp.setParam("Threads", 12)
    gp.setParam("TimeLimit",12000)

    # load graph
    graph_file = sys.argv[1]
    lp_output = sys.argv[2]
    graph_dir = sys.argv[3]

    graph_names = [graph_dir +"/" + f.rstrip("\n") for f in open(graph_file).readlines()]
    graph_names = sorted(graph_names)
    graphs = [pickle.load(open(name,'rb')) for name in graph_names ]

    graph_dict = dict(zip(graph_names, graphs))

    ofilename = os.path.basename(graph_file).split(".")[0] + "_diameters.csv"
    with open(ofilename, "w") as ofile:
        ofile.write("idx,diameter\n")

    all_dists = []
    for name in graph_names:
        print(os.path.basename(name))
        graph = graph_dict[name]

        max_dist = 0
        # for each graph, find two decompositions for 25 iterations
        for _ in trange(25):
            s1,s2 = sample_eulerian_cycles(graph, 2)
            lp_name = lp_output + "/" + ".".join(os.path.basename(name).split(".")[:-1])+".lp"

            # compute EMED between two decompositions
            dist = EMED_lp(s1,s2,lp_name)
            if dist == -1:
                print("!!! ERROR: LP model did not return a solution !!!")
                exit()
            max_dist = max(max_dist, dist)
        print("Diameter:", max_dist/100)

        # output to file
        with open(ofilename, "a") as ofile:
            ofile.write(os.path.basename(name)+","+str(max_dist/100)+"\n")
    