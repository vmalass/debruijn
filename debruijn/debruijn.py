#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics
import matplotlib.pyplot as plt
matplotlib.use("Agg")

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="k-mer size (default 22)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    with open(fastq_file, "r") as fasta:
        for i in fasta:
            yield next(fasta).strip()
            next(fasta)
            next(fasta)


def cut_kmer(read, kmer_size):
    for i in range(len(read)):
        if i + kmer_size <= len(read):
            yield read[i:i + kmer_size]


def build_kmer_dict(fastq_file, kmer_size):
    kmer_dict = {}
    for i in read_fastq(fastq_file):
        for j in cut_kmer(i, kmer_size):
            if j in kmer_dict:
                kmer_dict[j] += 1
            else:
                kmer_dict[j] = kmer_dict.get(j, 0) + 1
    return kmer_dict


def build_graph(kmer_dict):
    graph = nx.DiGraph()
    for key, value in kmer_dict.items():
        graph.add_edge(key[:-1], key[1:], weight = value)
    return graph


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    for i in path_list:
        if delete_entry_node and delete_sink_node is not True:
            graph.remove_nodes_from(i[:-1])
        elif delete_entry_node is not True and delete_sink_node:
            graph.remove_nodes_from(i[1:])
        elif delete_entry_node is not True and delete_sink_node is not True:
            graph.remove_nodes_from(i[1:-1])
        else:
            graph.remove_nodes_from(i)
        return graph 


def std(data):
    return statistics.stdev(data)


def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    list_path = []
    for first_path in range(len(path_list)):
        for second_path in range((first_path + 1), len(path_list)):
            std_length = std([path_length[first_path], path_length[second_path]])
            std_weight = std([weight_avg_list[first_path], weight_avg_list[second_path]])
            if std_length == 0 and std_weight == 0:
                rdm = randint([first_path, second_path])
                list_path.append(path_list[rdm])
            elif std_length != 0 and std_weight == 0:
                if path_length[first_path] > path_length[second_path]:
                    list_path.append(path_list[second_path])
                else:
                    list_path.append(path_list[first_path])
            else:
                if weight_avg_list[first_path] > weight_avg_list[second_path]:
                    list_path.append(path_list[second_path])
                else:
                    list_path.append(path_list[first_path])
    remove_paths(graph, list_path, delete_entry_node, delete_sink_node)
    return graph


def path_average_weight(graph, path):
    return statistics.mean([d["weight"] for (_, _, d) in graph.subgraph(path).edges(data=True)])


def solve_bubble(graph, ancestor_node, descendant_node):
    all_paths = list(nx.all_simple_paths(graph,ancestor_node, descendant_node))
    weights = []
    lengths = []
    for path in all_paths:
        weights.append(path_average_weight(graph, path))
        lengths.append(len(path))
    return select_best_path(graph, all_paths, lengths, weights)


def simplify_bubbles(graph):
    bubble = False
    for node in graph.nodes:
        list_predecessor = list(graph.predecessors(node))
        if len(list_predecessor) > 1:
            for i in range(len(list_predecessor) - 1):
                for j in range(i + 1, len(list_predecessor)):
                    ancestor_node = nx.lowest_common_ancestor(graph, list_predecessor[i], list_predecessor[j])
                    if ancestor_node:
                        bubble = True
                        break
            if bubble:
                graph = simplify_bubbles(solve_bubble(graph, ancestor_node, node))
                break
    return graph


def solve_entry_tips(graph, starting_nodes):
    lst_path = []
    wei_path = []
    len_path = []
    for node in starting_nodes:
        for descendant in nx.descendants(graph, node):  # descendants(G, source) Returns all nodes reachable from source in G.
            predecessor = list(graph.predecessors(descendant))
            if len(predecessor) > 1:
                for path in nx.all_simple_paths(graph, node, descendant):
                    lst_path.append(path)
                    len_path.append(len(path))
                    wei_path.append(path_average_weight(graph, path))
    graph = select_best_path(graph, lst_path, len_path, wei_path, True, False)
    return(graph) 


def solve_out_tips(graph, ending_nodes):
    lst_path = []
    wei_path = []
    len_path = []
    for node in ending_nodes:
        for ancestor in nx.ancestors(graph, node):  # ancestors(G, source) Returns all nodes having a path to source in G.
            successor = list(graph.successors(ancestor))
            if len(successor) > 1:
                for path in nx.all_simple_paths(graph, node, ancestor):
                    lst_path.append(path)
                    len_path.append(len(path))
                    wei_path.append(path_average_weight(graph, path))
    graph = select_best_path(graph, lst_path, len_path, wei_path, False, True)
    return(graph)


def get_starting_nodes(graph):
    starting = []
    for node in graph.nodes:
        predecessor = list(graph.predecessors(node))
        if not predecessor:
            starting.append(node)
    return starting


def get_sink_nodes(graph):
    sink = []
    for node in graph.nodes:
        successor = list(graph.successors(node))
        if not successor:
            sink.append(node)
    return sink


def get_contigs(graph, starting, sink):
    contigs = []
    for i in starting:
        for j in sink:
            for path in nx.all_simple_paths(graph, i, j):
                counter = path[0]
                for k in range(1, len(path)):
                    counter = counter + path[k][-1]
                contigs.append((counter, len(counter)))
    return(contigs)


def save_contigs(contigs, output_file):
    with open(output_file, "w") as filout:
        for i, j in enumerate(contigs):
            filout.write(">contig_"+str(i)+" len="+str(j[1])+"\n")
            filout.write(fill(j[0])+"\n")


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def draw_graph(graph, graphimg_file):
    """Draw the graph
    """                                    
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


def save_graph(graph, graph_file):
    """Save the graph with pickle
    """
    with open(graph_file, "wt") as save:
            pickle.dump(graph, save)


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    #Lecture du fichier et construction du graphe:
    dic_kmer = build_kmer_dict(args.fastq_file, args.kmer_size)
    graph = build_graph(dic_kmer)
    starting_nodes = get_starting_nodes(graph)
    sink_nodes = get_sink_nodes(graph)
    #Résolution des bulles:
    simplify_bubbles = simplify_bubbles(graph)
    #Resolution des points d'entrée et de sortie
    entry_tips = solve_entry_tips(simplify_bubbles, starting_nodes)
    out_tips = solve_out_tips(entry_tips, sink_nodes)
    #Ecriture des contigs:
    contigs = get_contigs(out_tips, starting_nodes, sink_nodes)
    save_contigs(contigs, args.output_file)

    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit 
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #     draw_graph(graph, args.graphimg_file)
    # Save the graph in file
    # if args.graph_file:
    #     save_graph(graph, args.graph_file)


if __name__ == '__main__':
    main()
