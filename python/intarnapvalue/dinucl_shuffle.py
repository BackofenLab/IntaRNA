#!/usr/bin/env python3

# dinucl_shuffle.py
# Original code by P. Clote, Oct 2003
# Modified for python3 by Fabio Gutmann <https://github.com/fabio-gut>, 2019
# NOTE: One cannot use function "count(s,word)" to count the number
# of occurrences of dinucleotide word in string s, since the built-in
# function counts only non-overlapping words, presumably in a left to
# right fashion.


import sys
import random


def compute_count_and_lists(s):
    # WARNING: Use of function count(s,'UU') returns 1 on word UUU
    # since it apparently counts only nonoverlapping words UU
    # For this reason, we work with the indices.

    # Initialize lists and mono- and dinucleotide dictionaries
    nucl_dict = {
      'A': [],
      'G': [],
      'C': [],
      'U': []
    }  # List is a dictionary of lists
    nucl_list = ["A", "C", "G", "U"]
    s = s.upper().replace("T", "U")
    nucl_cnt = {}  # empty dictionary
    dinucl_cnt = {}  # empty dictionary
    for x in nucl_list:
        nucl_cnt[x] = 0
        dinucl_cnt[x] = {}
        for y in nucl_list:
            dinucl_cnt[x][y] = 0

    # Compute count and lists
    nucl_cnt[s[0]] = 1
    nucl_total = 1
    dinucl_total = 0
    for i in range(len(s)-1):
        x = s[i]
        y = s[i+1]
        nucl_dict[x].append(y)
        nucl_cnt[y] += 1
        nucl_total += 1
        dinucl_cnt[x][y] += 1
        dinucl_total += 1
    assert (nucl_total == len(s))
    assert (dinucl_total == len(s)-1)
    return nucl_cnt, dinucl_cnt, nucl_dict
 
 
def choose_edge(x, dinucl_cnt):
    num_in_list = 0
    for y in ['A', 'C', 'G', 'U']:
        num_in_list += dinucl_cnt[x][y]
    z = random.random()
    denom = dinucl_cnt[x]['A'] + dinucl_cnt[x]['C'] + dinucl_cnt[x]['G'] + dinucl_cnt[x]['U']
    numerator = dinucl_cnt[x]['A']
    if z < float(numerator)/float(denom):
        dinucl_cnt[x]['A'] -= 1
        return 'A'
    numerator += dinucl_cnt[x]['C']
    if z < float(numerator)/float(denom):
        dinucl_cnt[x]['C'] -= 1
        return 'C'
    numerator += dinucl_cnt[x]['G']
    if z < float(numerator)/float(denom):
        dinucl_cnt[x]['G'] -= 1
        return 'G'
    dinucl_cnt[x]['U'] -= 1
    return 'U'


def connected_to_last(edge_list, nucl_list, last_ch):
    d = {}
    for x in nucl_list:
        d[x] = 0
    for edge in edge_list:
        a = edge[0]
        b = edge[1]
        if b == last_ch:
            d[a] = 1
    for i in range(2):
        for edge in edge_list:
            a = edge[0]
            b = edge[1]
            if d[b] == 1:
                d[a] = 1
    for x in nucl_list:
        if x != last_ch and d[x] == 0:
            return 0
    return 1
 

def eulerian(s):
    nucl_cnt, dinucl_cnt, nucl_dict = compute_count_and_lists(s)
    # compute nucleotides appearing in s
    nucl_list = []
    for x in ["A", "C", "G", "U"]:
        if x in s:
            nucl_list.append(x)
    # compute numInList[x] = number of dinucleotides beginning with x
    num_in_list = {}
    for x in nucl_list:
        num_in_list[x] = 0
        for y in nucl_list:
            num_in_list[x] += dinucl_cnt[x][y]
    # create dinucleotide shuffle L
    last_ch = s[-1]
    edge_list = []
    for x in nucl_list:
        if x != last_ch:
            edge_list.append([x, choose_edge(x, dinucl_cnt)])
    ok = connected_to_last(edge_list, nucl_list, last_ch)
    return ok, edge_list, nucl_list, last_ch


def shuffle_edge_list(l):
    n = len(l)
    barrier = n
    for i in range(n-1):
        z = int(random.random() * barrier)
        tmp = l[z]
        l[z] = l[barrier - 1]
        l[barrier - 1] = tmp
        barrier -= 1
    return l


def dinucl_shuffle(s):
    ok = 0
    edge_list, nucl_list = [], []  # initialize arrays
    while not ok:
        ok, edge_list, nucl_list, last_ch = eulerian(s)
    nucl_cnt, dinucl_cnt, nucl_dict = compute_count_and_lists(s)

    # remove last edges from each vertex list, shuffle, then add back
    # the removed edges at end of vertex lists.
    for x, y in edge_list:
        nucl_dict[x].remove(y)
    for x in nucl_list:
        shuffle_edge_list(nucl_dict[x])
    for x, y in edge_list:
        nucl_dict[x].append(y)

    # construct the eulerian path
    l = [s[0]]
    prev_ch = s[0]
    for i in range(len(s)-2):
        ch = nucl_dict[prev_ch][0]
        l.append(ch)
        del nucl_dict[prev_ch][0]
        prev_ch = ch
    l.append(s[-1])
    t = ''.join(l)
    return t
 

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: ./dinucl_shuffle.py <RNA sequence>")
        sys.exit(1)
    seq = sys.argv[1]
    for _ in range(4000):
        print(dinucl_shuffle(seq))
