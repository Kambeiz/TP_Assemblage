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

def read_fastq(fastq_file):
    """Reads a Fastq file.
    Arguments
    ---------
    fastq_file: Path to the Fastq file
    Yield
    -----
    Generator yielding Fasta reads
    """
    with open(fastq_file, "r") as f_in:
        for _ in f_in:
            yield next(f_in).strip()
            next(f_in)
            next(f_in)


def cut_kmer(fasta_seq, k):
    """Reads a Fastq file.
    Arguments
    ---------
    fasta_seq: Fasta sequence
    Yield
    -----
    Generator yielding Fasta sequence k-mers
    """
    for i in range(len(fasta_seq) - k+1):
        yield fasta_seq[i:i+k]


def build_kmer_dict(fastq_file, k):
    """Reads a Fastq file.
    Arguments
    ---------
    fasta_seq: Fasta sequence
    k: k-mer size
    Returns
    -------
    kmer_dict: k-mer dictionary
    """
    kmer_dict = {}
    for fasta_seq in read_fastq(fastq_file):
        for kmer in cut_kmer(fasta_seq, k):
            if kmer not in kmer_dict:
                kmer_dict[kmer] = 0
            kmer_dict[kmer] += 1
    return kmer_dict


def build_graph(kmer_dict):
    """Build a NetworkX graph from a k-mer dictionary.
    Arguments
    ---------
    kmer_dict: k-mer dictionary
    Returns
    -------
    graph: NetworkX graph
    """
    graph = nx.DiGraph()
    for kmer in kmer_dict:
        graph.add_edge(kmer[:-1], kmer[1:], weight=kmer_dict[kmer])
    return graph


def get_starting_nodes(graph):
    """Get the starting nodes from the graph.
    Arguments
    ---------
    graph: NetworkX graph
    Returns
    -------
    start_nodes: Starting nodes list
    """
    start_nodes = []
    for node in graph.nodes:
        if len(list(graph.predecessors(node))) == 0:
            start_nodes.append(node)
    return start_nodes


def get_sink_nodes(graph):
    """Get the sinking nodes from the graph.
    Arguments
    ---------
    graph: NetworkX graph
    Returns
    -------
    sink_nodes: Sinking nodes list
    """
    sink_nodes = []
    for node in graph.nodes:
        if len(list(graph.successors(node))) == 0:
            sink_nodes.append(node)
    return sink_nodes


def get_contigs(graph, start_nodes, sink_nodes):
    """Get all the contigs between the staring nodes and the sinking nodes.
    Arguments
    ---------
    graph: NetworkX graph
    start_nodes: Start nodes list
    sink_nodes: Sink nodes list
    Returns
    -------
    contigs_list: Contig list
    """
    contigs_list = []

    for start in start_nodes:
        for sink in sink_nodes:
            for path in nx.all_simple_paths(graph, source=start, target=sink):
                contig = ""
                for kmer in path:
                    if not contig:
                        contig += kmer
                    else:
                        contig += kmer[-1]
                contigs_list.append((contig, len(contig)))

    return contigs_list


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def save_contigs(contig_list, contig_filename):
    """Save all the contigs of a list into a Fasta format.
    Arguments
    ---------
    contig_list: Contig list
    contig_filename: Contif file name
    Returns
    -------
    Fasta file in the working directory
    """
    with open(contig_filename, "w") as f_out:
        for i, contig in enumerate(contig_list):
            f_out.write(f">contig_{i} len={contig[1]}\n")
            f_out.write(fill(contig[0]) + "\n")


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
    save_contigs(t, "test_y.fasta")