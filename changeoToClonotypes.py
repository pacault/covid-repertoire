#!/usr/bin/env python

import argparse
import os, sys
from natsort import natsorted, ns
from collections import Counter


def define_clone_dict():
    return {
        "most_fq_junction" : (),
        "freq" : -1,
        "v_call" : {},
        "j_call" : {},
        "seq_ids" : set(),
        "junctions": {}
        }


# def get_most_frequent_junction(junctions_list):
#     """ return a tupe with the most frequent element and his count"""
#     c = Counter(junctions_list)
#     most_common = c.most_common(1)
#     if len(most_common) > 0:
#         return most_common[0]
#     else:
#         return None


def get_all_allele(field, keep_allele=True):
    # Homsap IGHV3-25*04 ORF,Homsap IGHV3-72*01 F,Homsap IGHV3-73*01 F,Homsap IGHV3-73*02 F
    call = []
    for c in field.split(" "):
        if c.startswith("IGH"):
            if not keep_allele:
                call.append(c.split("*")[0])
            else:
                call.append(c)
            break
    return call


def get_most_frequent_item(dict):
    """ return a tupe with the most frequent element and his count"""
    occ = -1
    most_freq_item = ""
    for item, count in dict.items():
        if count > occ:
            most_freq_item = item
            occ = count
    return (most_freq_item, occ)


def main():         
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, help="Output of the changeo DefineClones.py command", required=True)
    parser.add_argument("-o", "--output", type=str, help="output file", required=True)
    args = parser.parse_args()

    clones = {}
    total_seq = 0
    with open(args.input, 'r') as f, open(args.output, 'w') as o:
        for i, line in enumerate(f):
            if i > 0:
                total_seq += 1
                cols = line.rstrip().split("\t")
                # print(cols)
                if cols[-2]:
                    clone_id = int(cols[-2])
                    if not clone_id in clones:
                        clones[clone_id] = define_clone_dict()
                    clones[clone_id]["seq_ids"].add(cols[0])
                    
                    # junctions
                    if not cols[10] in clones[clone_id]["junctions"]:
                        clones[clone_id]["junctions"][cols[10]] = 1
                    else:
                        clones[clone_id]["junctions"][cols[10]] += 1
                    
                    #v_call
                    v_calls = get_all_allele(cols[4], False)
                    for call in v_calls:
                        if not call in clones[clone_id]["v_call"]:
                            clones[clone_id]["v_call"][call] = 1
                        else:
                            clones[clone_id]["v_call"][call] += 1
                    
                    #j_call
                    v_calls = get_all_allele(cols[6])
                    for call in v_calls:
                        if not call in clones[clone_id]["j_call"]:
                            clones[clone_id]["j_call"][call] = 1
                        else:
                            clones[clone_id]["j_call"][call] += 1


        header = ["clone_id", "most_freq_junction", "count", "freq", "v_call", "j_call", "junctions"]
        o.write("\t".join(header) + "\n")
        i = 0
        for clone_id, clone in sorted(clones.items(), key=lambda x: len(x[1]["seq_ids"]), reverse=True): 
            if clone_id != 0:
                clone["most_fq_junction"] = get_most_frequent_item(clone["junctions"])
                clone["freq"] = len(clone["seq_ids"]) / total_seq
                junctions_str = ""
                d = ""
                for junction, count in clone["junctions"].items():
                    junctions_str += d + junction + "," + str(count)
                    d = ";"
                val = [clone_id,
                    clone["most_fq_junction"][0],
                    len(clone["seq_ids"]),
                    clone["freq"],
                    get_most_frequent_item(clone["v_call"])[0],
                    get_most_frequent_item(clone["j_call"])[0],
                    junctions_str]
                o.write("\t".join(map(str, val)) + "\n")


main()