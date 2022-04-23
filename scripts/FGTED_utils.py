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
        
        self.nodes = dict()                  # key: (u1, u2) value: idx
        self.project_1 = defaultdict(set)    # maps edge from graph 1 to its projection on alignment graph
        self.project_2 = defaultdict(set)    # maps edge from graph 1 to its projection on alignment graph
        self.edge_to_idx = dict()            # key: (n1,n2,cost,weight) value: idx
        self.idx_to_edge = dict()            # key: idx value:(n1,n2,cost, weight)
        
        self.incoming = defaultdict(set)     # key: node_id1, value: [node_id2] for edge (node2, node1)
        self.outgoing = defaultdict(set)     # key: node_id1, value: [node_id2] for edge (node1, node2)
        
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

def incoming_edges(aln_graph, node):
    '''
    helper function for FGTED_lp
    outputs a list of edge id for edges coming into node
    '''
    return [aln_graph.edge_to_idx[n] for n in aln_graph.incoming[node]]

def outgoing_edges(aln_graph, node):
    '''
    helper function for FGTED_lp
    outputs a list of edge id for edges coming out of node
    '''
    return [aln_graph.edge_to_idx[n] for n in aln_graph.outgoing[node]]

def FGTED_lp(aln_graph, graph_1, graph_2, ofilename):
    '''
    Objective: min cost(e) f(e)
    s.t. flow symmetry and flow coverage constraints
    '''
    try:
        # create variables set
        E = [idx for idx in aln_graph.idx_to_edge]  # each variable corresponds to an edge in alignment graph

        # create a new model
        model = gp.Model("FGTED")
        
        # creat variables
        x = model.addVars(E, vtype=GRB.CONTINUOUS, name="x")

        # set objective
        model.setObjective(sum([aln_graph.idx_to_edge[idx][2] * x[idx] for idx in E]), GRB.MINIMIZE)

        # flow conservation/symmetry
        for node in aln_graph.nodes.values():
            model.addConstr(sum([x[i] for i in incoming_edges(aln_graph, node)]) - sum([x[i] for i in outgoing_edges(aln_graph, node)]) == 0)

        # flow capacity
        for edge_id in aln_graph.idx_to_edge:
            model.addConstr(x[edge_id] <= aln_graph.idx_to_edge[edge_id][3])

        # flow coverage
        for e1 in aln_graph.project_1:
            total_weight = graph_1.idx_to_edge[e1][2]
            model.addConstr(sum([x[i] for i in aln_graph.project_1[e1]]) == total_weight)
        for e2 in aln_graph.project_2:
            total_weight = graph_2.idx_to_edge[e2][2]
            model.addConstr(sum([x[i] for i in aln_graph.project_2[e2]]) == total_weight) 

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


