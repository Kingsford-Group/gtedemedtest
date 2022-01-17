
###
# Functions necessary for building alignment graph and computing FGTED
###

import numpy as np
from collections import defaultdict
from Bio import SeqIO
import gurobipy as gp
from gurobipy import GRB
from collections import Counter
import pandas as pd
import sys

rev_dict = {"A":"T", "C":"G","G":"C","T":"A"}
gap_penalty = 1
mismatch_penalty = 1
match_reward = 0
large_penalty = 100000

def rev_comp(string):
    new_string = ""
    for c in string[::-1]:
        new_string += rev_dict[c]
    return new_string


def check_source_sink(adj_set):
    '''
    Check if weights to and from sink and source are the same
    '''
    weight_source = 0
    weight_sink = 0
    for u1 in adj_set:
        for u2 in adj_set[u1]:
            if u2 == 2:
                weight_sink += adj_set[u1][u2]
            if u1 == 1:
                weight_source += adj_set[u1][u2]
    assert(weight_source == weight_sink)
    return weight_source


def read_fasta(fname, get_seq = False):
    set_1 = [("_".join(seq.id.split("_")[:-1]), seq.id.split("_")[-1]) for seq in SeqIO.parse(fname,"fasta")]
    seqs = []
    if get_seq:
        seqs = [str(seq.seq) for seq in SeqIO.parse(fname,"fasta")]
        return set_1, seqs
    return set_1

def read_gfa(fname, weights, graph_type="rlz"):
    '''
    Read gfa into node set and adjacency lists. 
    Add node's reverse complement if reverse edge exists
    Parameters: gfa file name
    Return: node_set, adj_set -- dictionary: [(u1, u2)]: weight
    '''
    node_set = dict()
    node_weight_dict = defaultdict(int)
    adj_set = defaultdict(dict)
    reverse_adj_set = defaultdict(dict)
    rev_comp_dict = dict()
    path_list = []
    node_counter = 1
    for line in open(fname):
        l = line.rstrip("\n").split("\t")
        if l[0] == "S":
            node_set[int(l[1])] = l[2]
            node_counter += 1
        
        # read edge from path for rlzgraphs
        if graph_type == "rlz":
            if l[0] == "P":
                path = l[2].split(",")
                path_id = int(l[1])
                
                new_path = [1]
                
                w = int(weights[path_id][1])
                
                curr_node = int(path[0][:-1])
                next_node = 0
                if path[0][-1] == "-":
                    if curr_node not in rev_comp_dict:
                        node_set[node_counter] = rev_comp(node_set[curr_node])
                        rev_comp_dict[curr_node] = node_counter
                        curr_node = node_counter
                        node_counter += 1
                    else:
                        curr_node = rev_comp_dict[curr_node]
                
                new_path.append(curr_node)
                
                # source node
                if curr_node not in adj_set[1]:
                    adj_set[1][curr_node] = w
                else:
                    adj_set[1][curr_node] += w

                if 1 not in reverse_adj_set[curr_node]:
                    reverse_adj_set[curr_node][1] = w
                else:
                    reverse_adj_set[curr_node][1] += w

                    
                for n in path[1:]:
                    next_node = int(n[:-1])
                    if n[-1] == "-":
                        if next_node not in rev_comp_dict:
                            node_set[node_counter] = rev_comp(node_set[next_node])
                            rev_comp_dict[next_node] = node_counter
                            next_node = node_counter
                            node_counter += 1
                        else:
                            next_node = rev_comp_dict[next_node]
                    
                    new_path.append(next_node)
                    
                    if next_node not in adj_set[curr_node]:
                        adj_set[curr_node][next_node] = w
                    else:
                        adj_set[curr_node][next_node] += w
                        
                    if curr_node not in reverse_adj_set[next_node]:
                        reverse_adj_set[next_node][curr_node] = w
                    else:
                        reverse_adj_set[next_node][curr_node] += w    
                    
                    curr_node = next_node
        
                # sink node
                if 2 not in adj_set[next_node]:
                    adj_set[next_node][2] = w
                else:
                    adj_set[next_node][2] += w
#                 path_list.append(path)
                new_path.append(2)
                path_list.append(new_path)
          
    return node_set, adj_set, reverse_adj_set, path_list


class GTED_Graph(object):
    
    def __init__(self, node_set, adj_set, reverse_adj_set, total_weight):
        self.nodes = dict()                  # key: u1, value = label
        self.adj_list = defaultdict(dict)    # key: u1, value = [dict((u2, weight), (u3, weight))]
        self.edge_to_idx = dict()             # key: (u,v), value = idx
        self.idx_to_edge = dict()             # key: idx, value = (u,v)
        
        self.new_to_old = dict()
        self.old_to_new = dict()
        
        self.source = -1
        self.sink = -1
        
        node_counter = 1
        for u in adj_set:
            u_node_list, node_counter = self.split_node(u, node_set, reverse_adj_set, node_counter)
            
            # find source
            if u_node_list[0][1] == "$":
                self.source = u_node_list[0]
                
            for v in adj_set[u]:
                
                v_node_list, node_counter = self.split_node(v, node_set, reverse_adj_set, node_counter)
                
                # find sink
                if self.sink == -1 and v_node_list[0][1] == "#":
                    self.sink = v_node_list[0]

                weight = adj_set[u][v]
                self.create_edges([u_node_list[-1], v_node_list[0]], weight)
        
        # add edge from sink to source
        self.create_edges([self.sink, self.source], total_weight)
                
    def create_edges(self, node_list, weight):
        '''
        Given a list of nodes, create edge between each ordered pair with weight specified.
        Adds the created edge to edge_to_id and id_to_edge
        Return: None
        '''
        edge_idx = len(self.idx_to_edge)
        for i in range(len(node_list)-1):
            u = node_list[i][0]
            v = node_list[i+1][0]
            edge = (u, v, weight)
            self.adj_list[u][v] = weight
            self.edge_to_idx[(u,v)] = edge_idx
            self.idx_to_edge[edge_idx] = edge
            edge_idx += 1
            
    def split_node(self, node, node_set, reverse_adj_list, node_counter):
        '''
        Given a node, split into single-character nodes
        return: 
        node_list -- list of [nodeid, label] pairs
        '''
        
        edge_weight = sum([reverse_adj_list[node][i] for i in reverse_adj_list[node]])
        
        node_list = []
        if node not in self.old_to_new:        
            label = node_set[node]
            node_list = []
            for i in range(len(label)):
                node_list.append((node_counter, label[i]))
                self.nodes[node_counter] = label[i]
                node_counter += 1
            self.old_to_new[node] = node_list
            self.create_edges(node_list, edge_weight)
        else:
            node_list = self.old_to_new[node]
        
        return node_list, node_counter
    
def matched(label_1, label_2):
    '''
    Determine match score
    1) if characters: return match or mismatch
    2) free edges from both sink to both source, otherwise forbidden
    '''
    if label_1 in ["$", "#"] or label_2 in ["$","#"]:
        if label_1 == label_2:
            return 0
        else:
            return large_penalty
    if label_1 == label_2:
        return 0
    else:
        return 1
    
def find_gap(u1,u2,v1,v2):

    '''
    Impose large gap penalty if any node involves a sink
    In other words, only allow edges to sink if it is a match edge
    '''
    if (u1 == "#") or (u2 == "#") or (v1 == "#") or (v2 == "#"):
        return large_penalty
    else:
        return gap_penalty

class Alignment_Graph(object):
    
    def __init__(self, graph_1, graph_2):
        
        self.nodes = dict()                   # key: (u1, u2) value: idx
        self.project_1 = defaultdict(set)    # maps edge from graph 1 to its projection on alignment graph
        self.project_2 = defaultdict(set)    # maps edge from graph 1 to its projection on alignment graph
        self.edge_to_idx = dict()         # key: (n1,n2,cost,weight) value: idx
        self.idx_to_edge = dict()         # key: idx value:(n1,n2,cost, weight)
        
        self.incoming = defaultdict(set)    # key: node_id1, value: [node_id2] for edge (node2, node1)
        self.outgoing = defaultdict(set)    # key: node_id1, value: [node_id2] for edge (node1, node2)
        
        node_counter = 0
        # create aligned nodes
        for u1 in graph_1.nodes:
            for u2 in graph_2.nodes:
                self.nodes[(u1,u2)] = node_counter
                node_counter += 1

        edge_counter = 0
        # create aligned edges
        for u1 in graph_1.nodes:
            for v1 in graph_1.adj_list[u1]:
                for u2 in graph_2.nodes:
                    for v2 in graph_2.adj_list[u2]:
                        
                        # obtain edge in original graphs
                        edge1 = graph_1.edge_to_idx[(u1,v1)]
                        edge2 = graph_2.edge_to_idx[(u2,v2)]
                        weight1 = graph_1.adj_list[u1][v1]
                        weight2 = graph_2.adj_list[u2][v2]
                        
                        
                        # horizontal edge
                        # free edge for the last step to reach sink
                        g = find_gap(graph_1.nodes[u1], graph_2.nodes[u2], graph_1.nodes[v1], graph_2.nodes[u2])
                        h_edge = (self.nodes[(u1,u2)], self.nodes[(v1,u2)], g, weight1)

                        # vertical edge
                        g = find_gap(graph_1.nodes[u1], graph_2.nodes[u2], graph_1.nodes[u1], graph_2.nodes[v2])
                        v_edge = (self.nodes[(u1,u2)], self.nodes[(u1,v2)], g, weight2)

                        # diagonal edge
                        cost = matched(graph_1.nodes[v1], graph_2.nodes[v2])
                        d_edge = (self.nodes[(u1,u2)], self.nodes[(v1,v2)], cost, min(weight1, weight2))

#                         if self.nodes[(u1,u2)] == 0:
#                             print(u1,v1,u2,v2)
#                             print(h_edge, v_edge, d_edge)
                        edge_counter = self.create_edges([h_edge,v_edge,d_edge], edge_counter, edge1, edge2)
                        
    def create_edges(self, edges, edge_counter, e1, e2):
        '''
        Updates the data structures with newly created AG edges
        Parameters: edges -- [h_edge, v_edge, d_edge] ~ must follow the same order!
                    edge_counter -- keep tracks of the number of edges
                    e1, e2 -- edges in the original edges
        '''
        for i,edge in enumerate(edges):
            
            # don't add edge if edge has large penalty
            if edge[2] == large_penalty:
                continue

            # check if edge already exists
            if edge not in self.edge_to_idx:
                self.idx_to_edge[edge_counter] = edge
                self.edge_to_idx[edge] = edge_counter
                self.incoming[edge[1]].add(edge)
                self.outgoing[edge[0]].add(edge)
                edge_counter += 1

            edge_idx = self.edge_to_idx[edge]

            # update projection
            if i == 0: # horizontal
                self.project_1[e1].add(edge_idx)
            if i == 1: # vertical
                self.project_2[e2].add(edge_idx)
            if i == 2: # diagonal
                self.project_1[e1].add(edge_idx)
                self.project_2[e2].add(edge_idx)

        return edge_counter

def write_LP_exact(ofile_name, aln_graph, graph_1, graph_2):
    '''
    Objective: min cost(e) f(e)
    s.t. flow symmetry and flow coverage constraints
    '''
    ofile = open(ofile_name, "w")
    ofile.write("Minimize\n")
    
    # objective
    out_str = ""
    for edge_idx in aln_graph.idx_to_edge:
        edge = aln_graph.idx_to_edge[edge_idx]
        cost = edge[2]
        out_str += " + " + str(cost) + " x" + str(edge_idx)
    ofile.write(out_str[3:]+"\n")
    
    # constraints
    ofile.write("Subject To\n")
    
    # flow symmetry (v)
    constraint_idx = 0
    for node in aln_graph.nodes.values():
        out_str = ""
        # only add constraint if the node is not isolated
        if node in aln_graph.incoming:
            for prec in aln_graph.incoming[node]:
                edge_id = aln_graph.edge_to_idx[prec]
                out_str += " + x" + str(edge_id)
            for succ in aln_graph.outgoing[node]:
                edge_id = aln_graph.edge_to_idx[succ]
                out_str += " - x" + str(edge_id)
            out_str += " = 0\n"
            ofile.write("v" + str(constraint_idx) + ": " + out_str[3:])
            constraint_idx += 1
        
    # flow capacity (w)
    constraint_idx = 0
    out_str = ""
    for edge_id in aln_graph.idx_to_edge:
        out_str = "w" + str(constraint_idx) + ": "
        weight = aln_graph.idx_to_edge[edge_id][3]
        out_str += "x" + str(edge_id) + " <= " + str(weight) + "\n"
        
        out_str2 = "t" + str(constraint_idx) + ": " + "x" + str(edge_id) + " >= 0\n"
        ofile.write(out_str)
        ofile.write(out_str2)
        constraint_idx += 1
        
    # flow coverage (u)
    constraint_idx = 0
    for e1 in aln_graph.project_1:
        cons_str = "u" + str(constraint_idx) + ": "
        total_weight = graph_1.idx_to_edge[e1][2]
        out_str = ""
        for aln_edge_idx in aln_graph.project_1[e1]:
            out_str += " + x" + str(aln_edge_idx)
        out_str += " = " + str(total_weight) + "\n"
        ofile.write(cons_str + out_str[3:])
        constraint_idx += 1
        
    for e2 in aln_graph.project_2:
        cons_str = "u" + str(constraint_idx) + ": "
        total_weight = graph_2.idx_to_edge[e2][2]
        out_str = ""
        for aln_edge_idx in aln_graph.project_2[e2]:
            out_str += " + x" + str(aln_edge_idx)
        out_str += " = " + str(total_weight) + "\n"
        ofile.write(cons_str + out_str[3:])
        constraint_idx += 1
    
    ofile.write("End")
    ofile.close()    
    

def callback_lp(model, where):
    if where == GRB.Callback.BARRIER:
        time = model.cbGet(GRB.Callback.RUNTIME)
        iter_cnt = model.cbGet(GRB.Callback.BARRIER_ITRCNT)
        if iter_cnt > 50 or time > 6000:
            f = open(model._objfname,"w")
            f.write(str(model.cbGet(GRB.Callback.BARRIER_PRIMOBJ)))
            model.write(model._solfname)
        if iter_cnt > 51 or time > 6060:
            print("Here, barrier")
            model.terminate()

if __name__ == "__main__":
    test_dir = "test_dir/"
    seq_dir = "simulated/1018_fewer_parts/corrected/"
    rlz_dir = "rlz_graphs/1018/"
    gted_rlz_dir = "GTED_RLZ_LP_solver_file/1027/"

    # seq_dir = "test_dir/test_loop/"
    # rlz_dir = "test_dir/test_loop/"
    # gted_rlz_dir = "test_dir/test_loop/"

    num1 = int(sys.argv[1])
    num2 = int(sys.argv[2])

    gp.setParam("Threads", 12)
    gp.setParam("LogToConsole", 1)
    gp.setParam("OutputFlag", 1)
    gp.setParam("TimeLimit",6100)

    for i in range(num1, num2):
        print("Constructing GTED graph for %i" % i)
        weights1, seqs1 = read_fasta(seq_dir+"constructed_set_"+str(i)+".fasta", True)
        node_set1, adj_set1, reverse_adj_set1, path_list1 = read_gfa(rlz_dir + "constructed_set_"+str(i)+".gfa", weights1)
        total_weight1 = check_source_sink(adj_set1)
        graph1 = GTED_Graph(node_set1, adj_set1, reverse_adj_set1, total_weight1)

        for j in range(i+1, num2):
            fname = "GTED_RLZ_%i_%i" % (i,j)
            print(fname)

            print("Constructing GTED graph for %i" % j)

            weights2, seqs2 = read_fasta(seq_dir+"constructed_set_"+str(j)+".fasta", True)
            node_set2, adj_set2, reverse_adj_set2, path_list2 = read_gfa(rlz_dir + "constructed_set_"+str(j)+".gfa", weights2)
            total_weight2 = check_source_sink(adj_set2)

            graph2 = GTED_Graph(node_set2, adj_set2, reverse_adj_set2, total_weight2)

            print("Constructing Alignment Graph.")

            aln_graph = Alignment_Graph(graph1, graph2)


            print("Writing LP to file.")
            write_LP_exact(gted_rlz_dir + fname + ".lp", aln_graph, graph1, graph2)

            model = gp.read(gted_rlz_dir + fname + ".lp")
            model._objfname = gted_rlz_dir+ fname + ".obj"
            model._solfname = gted_rlz_dir + fname + ".sol"

            model.optimize(callback_lp) 

            if model.status == GRB.INF_OR_UNBD:
                # Turn presolve off to determine whether model is infeasible or unbounded
                model.setParam(GRB.Param.Presolve, 0)
                model.optimize()

            if model.status == GRB.OPTIMAL:
                with open(model._objfname,"w") as obj_file:
                    obj_file.write(str(model.objVal))
                model.write(gted_rlz_dir + fname + ".sol")

