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

# from random import randint
# from operator import itemgetter
import statistics
import argparse
import os
import sys
import random
import networkx as nx
# import matplotlib
# import matplotlib.pyplot as plt
random.seed(9001)

__author__ = "Debbah Nagi"
__copyright__ = "Universite de Paris"
__credits__ = ["Debbah Nagi"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Debbah Nagi"
__email__ = "debbah.nagi@gmail.fr"
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
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()

def read_fastq(fastq):
    """
    Read a fastq file and return an iterator containing each sequence
    on this fastq file.
    Parameters
    ----------
    fastq : filepath
        Filepath of the fastq file.

    Returns
    -------
    iterator
        Iterator containing sequences from the fastq file.
    """

    # bases = ['A','T','C','G']
    # dict_seq = {}
    # Opening and reading the file
    with open(fastq, "r") as f_read:
        # Version with yield that i never used before
        for line in f_read:
            yield next(f_read).strip('\n')
            next(f_read)
            next(f_read)

def cut_kmer(sequence, k_mer):
    """
    Take a sequence and cut it according a size given into k_mer.
    Parameters
    ----------
    sequence : Fasta.
        Nucleotidic sequence.
    k_mer : Integer
        Size of the k_mer chosen.

    Returns
    -------
    Iterator.
        Iterator containing all k_mer possible from the sequnece.

    """
    for i in range(0, len(sequence)-k_mer + 1):
        yield sequence[i:i+k_mer]

def build_kmer_dict(fastq, k_mer):
    """
    Take a fastq file, read it, and cut each sequence according to
    a kmer size chosen. Return a dictionnary of k_mer and their count.
    Parameters
    ----------
    fastq : filepath
        filepath of the fastq file.
    k_mer : integer
        size of the kmer chosen.

    Returns
    -------
    dict
        dict of kmer produced according the k_mer size, and containing the
        count of every kmer.

    """
    km_dict = {}
    for seq in read_fastq(fastq):
        kmerseq = cut_kmer(seq, k_mer)
        for kmer in kmerseq:
            if kmer not in km_dict.keys():
                km_dict[kmer] = 1
            else:
                km_dict[kmer] += 1
    return km_dict

def build_graph(kmer_dict):
    """
    BUilding a graph with NetworkX from a k_mer dictionnary

    Parameters
    ----------
    kmer_dict : dict
        Dictionnary of k_mer containing their sequence and count.

    Returns
    -------
    grap : diGraph
        Network diGraph.

    """
    # We create the graph
    grap = nx.DiGraph()
    keys = kmer_dict.keys()
    # We create our suffix and prefix from the key of the dictionnary
    for key in keys:
        prefix = key[0:len(key)-1]
        suffix = key[1:len(key)]
        # We add nodes and edge of our graph according the suffix and prefix
        grap.add_node(prefix)
        grap.add_node(suffix)
        grap.add_edge(prefix, suffix, weight=kmer_dict[key])
    return grap

def get_starting_nodes(grap):
    """
    Taking out starting nodes from a networkX graph

    Parameters
    ----------
    grap : networkX graph
        Graph obtained from the NetworkX module.

    Returns
    -------
    starting_nd : list
        list of starting nodes.

    """
    starting_nd = []
    for node in grap.nodes():
        pred_list = list(grap.predecessors(node))
        if len(pred_list) == 0:
            starting_nd.append(node)
    return starting_nd

def get_sink_nodes(grap):
    """
    Taking out sinking nodes from a networkX graph

    Parameters
    ----------
    grap : networkX graph
        Graph obtained from the NetworkX module.

    Returns
    -------
    sink_nd : list
        list of sinking nodes.

    """
    sink_nd = []
    for node in grap.nodes():
        pred_list = list(grap.successors(node))
        if len(pred_list) == 0:
            sink_nd.append(node)
    return sink_nd

def get_contigs(grap, starting_nd, sink_nd):
    """
    Getting all conting between starting nodes and sinkings nodes.

    Parameters
    ----------
    g : networkX graph
        Graph obtained from the NetworkX module.
    starting_nd : list
        list of starting nodes.
    sink_nd : list
        list of sinking nodes.

    Returns
    -------
    contig_list : list
        list containing contigs.

    """
    contig_list = []
    for starting in starting_nd:
        for sinking in sink_nd:
            for path in nx.all_simple_paths(grap, starting, sinking):
                contig = ""
                for i in range(0, len(path)):
                    if not contig:
                        contig += path[i]
                    else:
                        contig += path[i][-1]
                contig_list.append((contig, len(contig)))
    return contig_list

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def save_contigs(contig_list, output_name):
    """
    Saving all contings into a fasta file format.

    Parameters
    ----------
    contig_list : list
        list containing contigs.
    output_name : string
        name and path of the output, by default "contigs.fasta".

    Returns
    -------
    None.

    """
    with open(output_name, "w") as fil_out:
        counter = 0
        for tpl in contig_list:
            fil_out.write(">contig_{} len={}\n".format(counter, tpl[1]))
            fil_out.write(fill(tpl[0]) + "\n")
            counter += 1

def std(val_list):
    """
    Perform the calcul of the standard deviation from a provided list.

    Parameters
    ----------
    val_list : list
        List with values.

    Returns
    -------
    float
        standard deviation calculated from a list.

    """
    return statistics.stdev(val_list)

def path_average_weight(grap, path):
    """
    Perfom the calcul of the average weight from a path.

    Parameters
    ----------
    grap : networkX graph
        Graph obtained from the NetworkX module.
    path : path
        path obtained from differents nodes of the graph.

    Returns
    -------
    avg_weight : float
        average weight of the path.

    """
    weight = 0
    for i in range(len(path)-1):
        weight += grap.edges[path[i], path[i+1]]["weight"]
    avg_weight = weight/(len(path)-1)
    return avg_weight

def remove_paths(grap, paths, delete_entry_node, delete_sink_node):
    """
    Remove a path from a graph with/without his entry and sinking nodes.

    Parameters
    ----------
    grap : networkX graph
        Graph obtained from the NetworkX module.
    paths : list
        list of path obtained from differents nodes of the graph.
    delete_entry_node : Boolean
        Boolean value of "keeping or not" the entry node.
    delete_sink_node : Boolean
        Boolean value of "keeping or not" the sinking node.

    Returns
    -------
    grap : networkX graph
        Graph obtained from the NetworkX module with a removed path.

    """
    for i in range(0, len(paths)):
        if delete_entry_node:
            grap.remove_node(paths[i][0])
        if delete_sink_node:
            grap.remove_node(paths[i][-1])
        grap.remove_nodes_from(paths[i][1:-1])
    return grap

def select_best_path(grap, paths, len_paths, weight_paths,
                     delete_entry_node=False,
                     delete_sink_node=False):
    """
    Removing unwanted paths from a graph, while/without removing entry and
    sinking nodes.

    Parameters
    ----------
    grap : networkX graph
        Graph obtained from the NetworkX module.
    paths : list
        list of path obtained from differents nodes of the graph.
    len_paths : list
        list of length of paths.
    weight_paths : list
        list of weigh of paths.
    delete_entry_node : Boolean, optionnal
        Boolean value of "keeping or not" the entry node (by default, it is False).
    delete_sink_node : Boolean, optionnal
        Boolean value of "keeping or not" the sinking node (by default, it is False).

    Returns
    -------
    grap : networkX graph
        Graph obtained from the NetworkX module with unwanted paths removed.

    """

    wm_paths = []
    wm_weight = []
    wm_length = []
    lm_paths = []
    lm_weight = []
    lm_length = []
    undes_paths = []
    for i in range(0, len(paths)):
        if weight_paths[i] == max(weight_paths):
            wm_paths.append(paths[i])
            wm_weight.append(weight_paths[i])
            wm_length.append(len_paths[i])

    for i in range(0, len(wm_paths)):
        if wm_length[i] == max(wm_length):
            lm_length.append(wm_length[i])
            lm_paths.append(wm_paths[i])
            lm_weight.append(wm_weight[i])

    for path in paths:
        if path not in lm_paths:
            undes_paths.append(path)
    grap = remove_paths(grap, undes_paths, delete_entry_node, delete_sink_node)

    return grap

def solve_bubble(grap, ancestor, descendant):
    """
    CHoosing the best path between two nodes containing differents paths
    (bubble).

    Parameters
    ----------
    grap : networkX graph
        Graph obtained from the NetworkX module.
    ancestor : node object
        node from a networkX graph, here the starting one of the bubble, or
        ancestor.
    descendant : node object
                node from a networkX graph, here the ending one of the bubble, or
        descendant.

    Returns
    -------
    grap : networkX graph
        Graph obtained from the NetworkX module, with the best path chosen
        from the bubble.

    """

    paths = []
    path_l = []
    path_w = []
    simple_paths = nx.all_simple_paths(grap, ancestor, descendant)

    for path in simple_paths:
        path_l.append(len(path))
        path_w.append(path_average_weight(grap, path))
        paths.append(path)

    grap = select_best_path(grap, paths, path_l, path_w)

    return grap

def simplify_bubbles(grap):
    """
    Simplifying the bubbles of a given graph

    Parameters
    ----------
    grap : networkX graph
        Graph obtained from the NetworkX module.

    Returns
    -------
    grap : networkX graph
        Graph obtained from the NetworkX module with bubbles simplified.

    """

    bad_nd = []
    for des_nd in grap.nodes:
        pred_list = list(grap.predecessors(des_nd))
        leng = len(pred_list)
        if leng > 1:
            anc_nd = nx.lowest_common_ancestor(grap, pred_list[0], pred_list[1])
            bad_nd.append([anc_nd, des_nd])
    for anc_des in bad_nd:
        grap = solve_bubble(grap, anc_des[0], anc_des[1])
    return grap

def solve_entry_tips(grap, starting_nd):
    """
    Removing tips in starting not interesting, for keeping only pertinent ones.

    Parameters
    ----------
    grap : networkX graph
        Graph obtained from the NetworkX module.
    starting_nd : list
        list of starting nodes.
    sink_nd : list
        list of sinking nodes.

    Returns
    -------
    grap : networkX graph
        Graph obtained from the NetworkX module without entry tips uninteresting.

    """
    ancestors = []
    paths = []
    path_l = []
    path_w = []


    for node in starting_nd:
        for des in nx.descendants(grap, node):
          # with while it's tricky so we'll go with a for loop
          # while len(g.pred[des]) >= 2:
          #     if  des not in ancestors:
          #          ancestors.append(des)
            n_predecessor = grap.pred[des]
            if len(n_predecessor) >= 2 and des not in ancestors:
                ancestors.append(des)
        for anc in ancestors:
            for path in nx.all_simple_paths(grap, node, anc):
                path_w.append(path_average_weight(grap, path))
                path_l.append(len(path))
                paths.append(path)

        grap = select_best_path(grap, paths, path_l, path_w,
                                delete_entry_node=True,
                                delete_sink_node=False)

    return grap

def solve_out_tips(grap, sink_nd):
    """
    Removing tips in sinking not interesting, for keeping only pertinent ones.

    Parameters
    ----------
    grap : networkX graph
        Graph obtained from the NetworkX module.
    sink_nd : list
        list of sinking nodes.

    Returns
    -------
    grap : networkX graph
        Graph obtained from the NetworkX module without out tips uninteresting.

    """

    descendants = []
    paths = []
    path_l = []
    path_w = []

    for node in sink_nd:
        for nd_nxt in nx.ancestors(grap, node):
            n_successor = grap.succ[nd_nxt]
            if nd_nxt not in descendants and len(n_successor) >= 2:
                descendants.append(nd_nxt)
        for des in descendants:
            for path in nx.all_simple_paths(grap, des, node):
                path_w.append(path_average_weight(grap, path))
                path_l.append(len(path))
                paths.append(path)

        grap = select_best_path(grap, paths, path_l, path_w,
                                delete_entry_node=False,
                                delete_sink_node=True)

    return grap


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    return args

if __name__ == '__main__':
    # Getting arguments
    ARGS = main()
    # GEtting sequence from fastq
    SEQ = read_fastq(ARGS.fastq_file)
    # Getting the k_mer dictionnary
    KMER_DICT = build_kmer_dict(ARGS.fastq_file, ARGS.kmer_size)
    # Building the graph
    G = build_graph(KMER_DICT)
    # Getting the starting nodes and sinking nodes
    STARTING_ND = get_starting_nodes(G)
    SINK_ND = get_sink_nodes(G)
    # Taking out bubbles
    G = simplify_bubbles(G)
    # Taking out bad entries and exits nodes
    G = solve_entry_tips(G, STARTING_ND)
    G = solve_out_tips(G, SINK_ND)
    # Saving new starting nodes and sinking nodes
    SINK_ND = get_starting_nodes(G)
    STARTING_ND = get_sink_nodes(G)
    # Saving contigs
    CONTIGS = get_contigs(G, STARTING_ND, SINK_ND)
    save_contigs(CONTIGS, ARGS.output_file)

    # Blast result with eva71.fna and contigs.fasta (kmer_size = 21):
    # E-Value = 0, Perc.ident = 100.0%
    # Pylint test = 9.51/10