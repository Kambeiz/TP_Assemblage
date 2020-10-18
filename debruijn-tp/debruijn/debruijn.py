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
import matplotlib.pyplot as plt
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics

__author__ = "Debbah Nagi"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Debbah Nagi"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
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
    bases = ['A','T','C','G']
    dict_seq = {}
    # Opening and reading the file
    with open(fastq, "r") as f_read:
        ''' 
        Version with dict and list
        for line in f_read:
            # Lines starting with @ are identifiant of a sequence
            if line.startswith("@"):
                seq_id = line[1:]
            # A line with a sequence will contain only bases    
            elif all(base in line for base in bases):
            # We create with that id and the sequence a key and id into the sequence
                dict_seq[seq_id] = line
        return dict_seq
'''
        # Version with yield that i never used before
        for line in f_read:
            yield next(f_read).strip('\n')
            next(f_read)
            next(f_read)
    
def cut_kmer(sequence, k_mer):
    ''' Dict/list version
    # We create a list for storing our k_mer
    k_mer_list=[]
    # The total k_mer poissible is len(sequence) minus the k_mer size plus one
    for i in range(0, len(sequence)-k_mer+1):
        k_mer_list.append(sequence[i:i+k_mer])
    return k_mer_list'''
    for i in range(0,len(sequence) - k_mer + 1 ):
        yield(sequence[i:i+k_mer])

def build_kmer_dict(fastq, k_mer):
    ''' Dict version
    k_mer_dict = {}
    dict_seq = read_fastq(fastq)
    for value in dict_seq.values():
        kmerseq = cut_kmer(value, k_mer)
        for kmer in kmerseq :
            if kmer not in k_mer_dict.keys():
                k_mer_dict[kmer] = 1
            elif kmer in k_mer_dict.keys():
                k_mer_dict[kmer] += 1
    return k_mer_dict
    '''
    # Other version
    kmer_dict ={}
    for seq in read_fastq(fastq):
        kmerseq = cut_kmer(seq, k_mer)
        for kmer in kmerseq :
            if kmer not in kmer_dict.keys():
                kmer_dict[kmer] = 1
            else:
                kmer_dict[kmer] += 1
    return kmer_dict    

'''
def build_graph(fastq, k_mer):
    plt.figure(figsize=(20,40))
    gen_seq = read_fastq(fastq) 
    g = nx.DiGraph()
    kmer_dict = build_kmer_dict(fastq, k_mer)
    for seq in gen_seq:
        kmerseq = list(cut_kmer(seq,k_mer))
        for i in range(0, len(kmerseq)-1):
            weight = kmer_dict[kmerseq[i]]
            g.add_edge(kmerseq[i], kmerseq[i+1], weight=weight)
    nx.draw(g, with_labels=True)
    plt.show()
    '''
def build_graph(kmer_dict):
    # BUILD
    g = nx.DiGraph()
    keys = kmer_dict.keys()
    for key in keys:
        prefix = key[0:len(key)-1]
        suffix = key[1:len(key)] 
        g.add_node(prefix)
        g.add_node(suffix)
        g.add_edge(prefix, suffix, weight=kmer_dict[key])
    return g

def get_starting_nodes(g):
    starting_nd = []
    for node in g.nodes():
        pred_list = list(g.predecessors(node))
        if len(pred_list) == 0:
            starting_nd.append(node)
    return starting_nd

def get_sink_nodes(g):
    sink_nd= []
    for node in g.nodes():
        pred_list = list(g.successors(node))
        if len(pred_list) == 0:
            sink_nd.append(node)
    return sink_nd

def get_contigs(g, starting_nd, sink_nd):
    contig_list =[]
    for starting in starting_nd:
        for sinking in sink_nd:
            for path in nx.all_simple_paths(g, starting, sinking):
                # print("Voici une liste de chemin: ", path)
                contig = ""
                for i in range(0, len(path)) :
                    if not contig:
                        contig += path[i]
                    else:
                        contig += path[i][-1]
                contig_list.append( (contig, len(contig)) )
    return contig_list

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def save_contigs(contig_list, output_name):
    with open(output_name, "w") as fil_out:
        counter = 0
        for tpl in contig_list:
            fil_out.write(">contig_{} len={}\n".format(counter, tpl[1]))
            fil_out.write(fill(tpl[0]) + "\n")
            counter +=1

def std(val_list):
    return statistics.stdev(val_list)

def path_average_weight(g, path):
    weight=0
    for i in range(len(path)-1):
        weight += g.edges[path[i], path[i+1]]["weight"]    
    avg_weight = weight/(len(path)-1)
    return avg_weight

def remove_paths(g, paths, delete_entry_node, delete_sink_node):
    for i in range(0, len(paths)):
        if delete_entry_node==True:
            g.remove_node(paths[i][0])
        if delete_sink_node==True:
            g.remove_node(paths[i][-1])
        g.remove_nodes_from(paths[i][1:-1])
    return g

def select_best_path(g, paths, len_paths,weight_paths, delete_entry_node=False, 
                     delete_sink_node=False):

    wm_paths = []
    wm_weight =[]
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

    for i in range(0,len(wm_paths)):
        if wm_length[i] == max(wm_length):
            lm_length.append(wm_length[i])
            lm_paths.append(wm_paths[i])
            lm_weight.append(wm_weight[i])

    for path in paths:
        if path not in lm_paths:
            undes_paths.append(path)
    g = remove_paths(g, undes_paths, delete_entry_node, delete_sink_node)
    
    return g

    
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
    args =  main()
    fastq = args.fastq_file
    kmer_size = args.kmer_size
    test = build_kmer_dict(fastq,kmer_size)
    g = build_graph(test)
    t = get_contigs(g, get_starting_nodes(g), get_sink_nodes(g))
