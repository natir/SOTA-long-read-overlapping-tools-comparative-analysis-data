#!/usr/bin/env python3

from sys import intern

import pickle
import itertools

import os
import multiprocessing
import begin

from enum import Enum
from types import SimpleNamespace

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.style.use("dark_background")
import venn


def __compute_overhang(strand, len_a, beg_a, end_a, len_b, beg_b, end_b):
    if strand == "+":
        return min(beg_a, beg_b) + min(len_a - end_a, len_b - end_b)
    return min(beg_a, len_b - end_b) + min(beg_b, len_a - end_a)


def good_overlap(strand, len_a, beg_a, end_a, len_b, beg_b, end_b):
    overhang = __compute_overhang(strand, len_a, beg_a, end_a, len_b, beg_b, end_b)
    maplen = max(end_a - beg_a, end_b - beg_b)
    
    # print("overhang", overhang)
    # print("maplen", maplen)
    # print("param", strand, len_a, beg_a, end_a, len_a - end_a, len_b, beg_b, end_b, len_b - end_b)
    if overhang > maplen*0.8 or overhang > 1000:
    # Strange overlap
        return False
    elif strand == "+" and beg_a <= beg_b and len_a - end_a <= len_b - end_b:
    # B containe A
        return False 
    elif strand == "-" and beg_a <= len_b - end_b and len_a - end_a <= beg_b:
    # B containe A
        return False 
    elif strand == "+" and beg_a >= beg_b and len_a - end_a >= len_b - end_b:
        # A containe B
        return False
    elif strand == "-" and beg_a >= len_b - end_b and len_a - end_a >= beg_b:
        # A containe B
        return False

    return True

def parse_mhap_gene(file_path):
    with open(file_path) as fh_in:
        for line in fh_in:
            line = line.strip().split()
            if good_overlap("-+"[int(line[4] == line[8])], int(line[7]),
                            int(line[5]), int(line[6]), int(line[11]),
                            int(line[9]), int(line[10])):
                if int(line[0]) < int(line[1]):
                    yield (int(line[0]), int(line[1]))
                else:
                    yield (int(line[1]), int(line[0]))

def parse_mhap(file_path):
    return set(parse_mhap_gene(file_path))

def parse_paf_gene(file_path):
    with open(file_path) as fh_in:
        for line in fh_in:
            line = line.split("\t")
            if good_overlap(line[4], int(line[1]), int(line[2]), int(line[3]),
                            int(line[6]),int(line[7]), int(line[8])):
                if int(line[0]) < int(line[5]):
                    yield (int(line[0]), int(line[5]))
                else:
                    yield (int(line[5]), int(line[0]))

def parse_paf(file_path):
    return set(parse_paf_gene(file_path))

def parse_hisea_gene(file_path):
    with open(file_path) as fh_in:
        for line in fh_in:
            line = line.strip().split()
            if good_overlap("+-"[int(line[8])], int(line[7]), int(line[5]), int(line[6]), int(line[11]), int(line[9]), int(line[10])):
                if int(line[0]) < int(line[1]):
                    yield (int(line[0]), int(line[1]))
                else:
                    yield (int(line[1]), int(line[0]))

def parse_hisea(file_path):
    return set(parse_hisea_gene(file_path))


def parse(file_path):
    if os.path.splitext(file_path)[1] == ".mhap":
        return parse_mhap(file_path)
    elif os.path.splitext(file_path)[1] == ".paf":
        return parse_paf(file_path)
    elif os.path.splitext(file_path)[1] == ".hisea":
        return parse_hisea(file_path)
    else:
        raise NameError("File extension isn't support "+file_path)

def read_overlap(file_path):
    if os.path.isfile(file_path+".pickle"):
        return file_path, pickle.load(open(file_path+".pickle", "rb"))
    else:
        obj = parse(file_path)
        pickle.dump(obj, open(file_path+".pickle", "wb"))
        return file_path, obj


def rename_name(label):
    if label.startswith("hisea"):
        return "hisea"
    elif label.startswith("mhap1.6"):
        return "mhap1.6"
    elif label.startswith("mhap2.1"):
        return "mhap2.1"
    elif label.startswith("graphmap"):
        return "graphmap"
    elif label.startswith("minimap2"):
        return "minimap2"
    elif label.startswith("minimap"):
        return "minimap"
    else:
        return label

@begin.start
def main(out_prefix, *file_list):

    select = "1101"
    n = 10
    
    if len(file_list) == 2:
        venn_generator = venn.venn2
    elif len(file_list) == 3:
        venn_generator = venn.venn3
    elif len(file_list) == 4:
        venn_generator = venn.venn4
    elif len(file_list) == 5:
        venn_generator = venn.venn5
    elif len(file_list) == 6:
        venn_generator = venn.venn6
    else:
        print("Error we just support ensembl between 2 and 6 value")
        return -1

    pool = multiprocessing.Pool(len(file_list))
    result = pool.map_async(read_overlap, file_list)
    file2set = {k: v for k, v in result.get()}

    labels, s = venn.get_labels(list(file2set.values()), select=select, n=n)
   
    for (key_a, val_a), (key_b, val_b) in itertools.combinations(file2set.items(), 2):
        print(key_a, key_b, len(val_a & val_b)/len(val_a | val_b))

    fig, ax = venn_generator(labels, names=[rename_name(name) for name in file2set.keys()])
    fig.savefig(out_prefix)

