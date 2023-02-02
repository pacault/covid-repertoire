#!/usr/bin/env python

import argparse
import os
import shutil
import math
import collections
import statistics
import gzip
from natsort import natsorted, ns

MAX_DIFF = 1
MAX_DIST = 1
MIN_READ = 5
MIN_PATIENT_NB = 4
HASH_NAMES = {}
CONTROLES_ADN = []
PATIENTS_ADN = []
PATIENTS_ARN = []
CONTROLES_ARN = []

patients_corr = {
    "04R" : ("Patient25R", "4d4d4dff"),
    "05O" : ("Patient26O", "5da4daff"),
    "06R" : ("Patient27R", "faa53bff"),
    "09R" : ("Patient28R", "60bd67ff"),
    "10R" : ("Patient12R", "f27cb1ff"),
    "11O" : ("Patient220", "b1912eff"),
    "12O" : ("Patient29O", "b277b2ff"),
    "14O" : ("Patient30O", "dfcf40ff"),
    "15_VAL_FIL" : ("Control01", "5a5a5aff"),
    "16R" : ("Patient31R", "f15955ff"),
    "18_LAU_MAR" : ("Control02", "5a5a5aff"),
    "18O" : ("Patient24O", "386e97ff"),
    "19_PHA_CHE" : ("Control03", "5a5a5aff"),
    "20O" : ("Patient32O", "bc6b09ff"),
    "24_GOR_GEO" : ("Control04", "5a5a5aff"),
    "24R" : ("Patient33R", "298130ff"),
    "25R" : ("Patient20R", "d81f72ff"),
    "27_SUC_GUY" : ("Control05", "5a5a5aff"),
    "29_MAR_MIC" : ("Control06", "5a5a5aff"),
    "29R" : ("Patient34R", "7b6930ff"),
    "31_DUM_CLA" : ("Control07", "5a5a5aff"),
    "32R" : ("Patient17R", "be37beff"),
    "33_BOR_COL" : ("Control08", "5a5a5aff"),
    "33R" : ("Patient04R", "b9ab33ff"),
    "35R" : ("Patient10R", "b91f1bff"),
    "36O" : ("Patient140", "ea8c89ff"),
    "37R" : ("Patient01R", "a3c2daff"),
    "37_ROC_ANT" : ("Control09", "5a5a5aff"),
    "38_COI_GIL" : ("Control10", "5a5a5aff"),
    "39_DAC_JOS" : ("Control11", "b186ffff"),
    "42_MAS_MON" : ("Control12", "c77cffff"),
    "AB_53" : ("Patient01R", "a3c2daff"),
    "BAY_FRE" : ("Control13", "e76bf3ff"),
    "BEL_AZI" : ("Control14", "f265e8ff"),
    "BG_10" : ("Patient02", "fa62dbff"),
    "BG_11" : ("Patient03", "ff61ccff"),
    "BLO_HEN" : ("Control15", "ff62bcff"),
    "BZ_42" : ("Patient04R", "b9ab33ff"),
    "CG_12" : ("Patient05", "ff6a98ff"),
    "CHA_FRE" : ("Control16", "ff6c91ff"),
    "CHE_ALA" : ("Control17", "5a5a5aff"),
    "CM_06" : ("Patient06", "ff65aaff"),
    "DC_02" : ("Patient07", "d973fcff"),
    "DJ_04" : ("Patient08", "a3a500ff"),
    "DOU_DAN" : ("Control18", "5a5a5aff"),
    "FENNICHE" : ("Patient09", "00c091ff"),
    "FIL_JEA" : ("Control19", "5a5a5aff"),
    "GM_45" : ("Patient10R", "b91f1bff"),
    "GORER" : ("Patient11", "7cae00ff"),
    "HH_64" : ("Patient12R", "f27cb1ff"),
    "JC_07" : ("Patient13", "39b600ff"),
    "KM_74" : ("Patient140", "ea8c89ff"),
    "KUEBO" : ("Patient15", "00c1a3ff"),
    "LAC_CAM" : ("Control20", "5a5a5aff"),
    "LR_026" : ("Patient16", "00bb4eff"),
    "MA_66" : ("Patient17R", "be37beff"),
    "MM_08" : ("Patient18", "00bfc4ff"),
    "PE_80" : ("Patient19", "00bae0ff"),
    "PH_73" : ("Patient20R", "d81f72ff"),
    "RA_03" : ("Patient21", "7099ffff"),
    "RIG_REN" : ("Control21", "5a5a5aff"),
    "SJ_50" : ("Patient220", "b1912eff"),
    "TES_ALA" : ("Control22", "5a5a5aff"),
    "TOSUN" : ("Patient23", "9590ffff"),
    "VS_03_VM" : ("Control23", "5a5a5aff"),
    "VS_09_GJ" : ("Control24", "5a5a5aff"),
    "VS_13FR" : ("Control25", "5a5a5aff"),
    "VS_30" : ("Control26", "5a5a5aff"),
    "YM_65" : ("Patient24O", "b5a628ff"),
    "BM_67" : ("Patient35", "d17b10ff"),
    "CD_65" : ("Patient36", "8898a5ff"),
    "DA_49" : ("Patient37", "b24946ff"),
    "LT_70" : ("Patient38", "82c088ff"),
    "MJ_71" : ("Patient39", "7c201dff"),

}

def hamming_distance(s1, s2):
    #Return the Hamming distance between equal-length sequences
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2)) / len(s1)


def parse_sample_file(sample_file):
    with open(sample_file, 'r') as f:
        for i, line in enumerate(f):
            group, name, id, day, imgt_dir = line.rstrip().split("\t")
            patient = {
                "name": name,
                "id": id,
                "day": day,
                "imgt_dir": imgt_dir
            }
            if group == "CONTROLES_ADN":
                CONTROLES_ADN.append(patient)
            elif group == "CONTROLES_ARN":
                CONTROLES_ARN.append(patient)
            elif group == "PATIENTS_ARN":
                PATIENTS_ARN.append(patient)
            elif group == "PATIENTS_ADN":
                PATIENTS_ADN.append(patient)


def get_patients_clolor(pat):
    if pat in patients_corr:
        return "#" + patients_corr[pat][1][:-2].upper()
    else:
        return "#D4D4D4"

def get_samples():
    samples = {}
    
    for sample in PATIENTS_ARN:
        print(sample["name"])
        # print(sample["imgt_dir"])
        if not sample["name"] in samples:
            samples[sample["name"]] = {}
            samples[sample["name"]]["name"] = sample["name"]
            samples[sample["name"]]["id"] = sample["id"]
        if not "ARN" in samples[sample["name"]]:
            samples[sample["name"]]["ARN"] = {}
        samples[sample["name"]]["ARN"][sample["day"]] = {}
        samples[sample["name"]]["ARN"][sample["day"]]["imgt_dir"] = sample["imgt_dir"]
        # print(sample["name"])
        for file in os.listdir(samples[sample["name"]]["ARN"][sample["day"]]["imgt_dir"]):
            # print(file)
            if file.endswith("_clonotypes_changeo.tsv.gz"):
                samples[sample["name"]]["ARN"][sample["day"]]["clonofile"] = os.path.join(samples[sample["name"]]["ARN"][sample["day"]]["imgt_dir"], file)
            elif file.endswith("clones_id.tsv.gz"):
                samples[sample["name"]]["ARN"][sample["day"]]["clonesfile"] = os.path.join(samples[sample["name"]]["ARN"][sample["day"]]["imgt_dir"], file)
            elif file == "8_V-REGION-nt-mutation-statistics.txt.gz":
                samples[sample["name"]]["ARN"][sample["day"]]["mutationsfile"] = os.path.join(samples[sample["name"]]["ARN"][sample["day"]]["imgt_dir"], file)
        samples[sample["name"]]["ARN"][sample["day"]]["clonotypes"] = load_clonotypes(samples[sample["name"]]["ARN"][sample["day"]]["clonofile"])
        samples[sample["name"]]["ARN"][sample["day"]]["mutation_rates"], samples[sample["name"]]["ARN"][sample["day"]]["read_names"] = get_all_mutation_rate(samples[sample["name"]]["ARN"][sample["day"]]["mutationsfile"], samples[sample["name"]]["ARN"][sample["day"]]["clonesfile"])

    for sample in PATIENTS_ADN:
        print(sample["name"])
        # print(sample["imgt_dir"])
        if not sample["name"] in samples:
            samples[sample["name"]] = {}
            samples[sample["name"]]["name"] = sample["name"]
            samples[sample["name"]]["id"] = sample["id"]
        if not "ADN" in samples[sample["name"]]:
            samples[sample["name"]]["ADN"] = {}
        samples[sample["name"]]["ADN"][sample["day"]] = {}
        samples[sample["name"]]["ADN"][sample["day"]]["imgt_dir"] = sample["imgt_dir"]
        for file in os.listdir(samples[sample["name"]]["ADN"][sample["day"]]["imgt_dir"]):
            if file.endswith("_clonotypes2.tsv.gz") or file.endswith("_clonotypes_changeo.tsv.gz"):
                samples[sample["name"]]["ADN"][sample["day"]]["clonofile"] = os.path.join(samples[sample["name"]]["ADN"][sample["day"]]["imgt_dir"], file)
            elif file.endswith("clones_id.tsv.gz"):
                samples[sample["name"]]["ADN"][sample["day"]]["clonesfile"] = os.path.join(samples[sample["name"]]["ADN"][sample["day"]]["imgt_dir"], file)
            elif file == "8_V-REGION-nt-mutation-statistics.txt.gz":
                samples[sample["name"]]["ADN"][sample["day"]]["mutationsfile"] = os.path.join(samples[sample["name"]]["ADN"][sample["day"]]["imgt_dir"], file)
        samples[sample["name"]]["ADN"][sample["day"]]["clonotypes"] = load_clonotypes(samples[sample["name"]]["ADN"][sample["day"]]["clonofile"], False)
        samples[sample["name"]]["ADN"][sample["day"]]["mutation_rates"], samples[sample["name"]]["ADN"][sample["day"]]["read_names"] = get_all_mutation_rate(samples[sample["name"]]["ADN"][sample["day"]]["mutationsfile"], samples[sample["name"]]["ADN"][sample["day"]]["clonesfile"])

    return samples


def get_controls():
    controls = {}
    

    for sample in CONTROLES_ARN:
        print(sample["name"])
        # print(sample["imgt_dir"])
        if not sample["name"] in controls:
            controls[sample["name"]] = {}
            controls[sample["name"]]["name"] = sample["name"]
            controls[sample["name"]]["id"] = sample["id"]
        if not "ARN" in controls[sample["name"]]:
            controls[sample["name"]]["ARN"] = {}
        controls[sample["name"]]["ARN"]["imgt_dir"] = sample["imgt_dir"]
        for file in os.listdir(controls[sample["name"]]["ARN"]["imgt_dir"]):
            if file.endswith("_clonotypes_changeo.tsv.gz"):
                controls[sample["name"]]["ARN"]["clonofile"] = os.path.join(controls[sample["name"]]["ARN"]["imgt_dir"], file)
            elif file.endswith("clones_id.tsv.gz"):
                controls[sample["name"]]["ARN"]["clonesfile"] = os.path.join(controls[sample["name"]]["ARN"]["imgt_dir"], file)
            elif file == "8_V-REGION-nt-mutation-statistics.txt.gz":
                controls[sample["name"]]["ARN"]["mutationsfile"] = os.path.join(controls[sample["name"]]["ARN"]["imgt_dir"], file)
        # print(sample["name"])
        controls[sample["name"]]["ARN"]["clonotypes"] = load_clonotypes(controls[sample["name"]]["ARN"]["clonofile"])
        controls[sample["name"]]["ARN"]["mutation_rates"], controls[sample["name"]]["ARN"]["read_names"] = get_all_mutation_rate(controls[sample["name"]]["ARN"]["mutationsfile"], controls[sample["name"]]["ARN"]["clonesfile"])
        
    for sample in CONTROLES_ADN:
        print(sample["name"])
        # print(sample["imgt_dir"])
        if not sample["name"] in controls:
            controls[sample["name"]] = {}
            controls[sample["name"]]["name"] = sample["name"]
            controls[sample["name"]]["id"] = sample["id"]
        if not "ADN" in controls[sample["name"]]:
            controls[sample["name"]]["ADN"] = {}
        controls[sample["name"]]["ADN"]["imgt_dir"] = sample["imgt_dir"]
        for file in os.listdir(controls[sample["name"]]["ADN"]["imgt_dir"]):
            if file.endswith("_clonotypes2.tsv.gz") or file.endswith("_clonotypes_changeo.tsv.gz"):
                controls[sample["name"]]["ADN"]["clonofile"] = os.path.join(controls[sample["name"]]["ADN"]["imgt_dir"], file)
            elif file.endswith("clones_id.tsv.gz"):
                controls[sample["name"]]["ADN"]["clonesfile"] = os.path.join(controls[sample["name"]]["ADN"]["imgt_dir"], file)
            elif file == "8_V-REGION-nt-mutation-statistics.txt.gz":
                controls[sample["name"]]["ADN"]["mutationsfile"] = os.path.join(controls[sample["name"]]["ADN"]["imgt_dir"], file)
        controls[sample["name"]]["ADN"]["clonotypes"] = load_clonotypes(controls[sample["name"]]["ADN"]["clonofile"], False)
        controls[sample["name"]]["ADN"]["mutation_rates"], controls[sample["name"]]["ADN"]["read_names"] = get_all_mutation_rate(controls[sample["name"]]["ADN"]["mutationsfile"], controls[sample["name"]]["ADN"]["clonesfile"])
    
    return controls


def load_clonotypes(file, rna=True):
    clonos = {}
    with gzip.open(file, 'rb') as f:
        for i, line in enumerate(f):
            if i > 0:
                cols = line.decode().rstrip().split("\t")
                if cols[0] != "" and cols[0] != "1" and (rna or int(cols[2]) >= MIN_READ):
                    call_hash = cols[4] + "_" + cols[5].split("*")[0]
                    cdr3 = cols[1][1:-1]
                    cdr3_length = len(cdr3)
                    if not call_hash in clonos:
                        clonos[call_hash] = {}
                    if not cdr3_length in clonos[call_hash]:
                        clonos[call_hash][cdr3_length] = []

                    clono = {
                            "id":int(cols[0]),
                            "cdr3":cols[1],
                            "count":int(cols[2]),
                            "freq":float(cols[3]),
                            "v":cols[4],
                            "j":cols[5],
                            "sequences": cols[6],
                            "hash": cols[4] + "_" + cols[5] + "_" + cols[1]                        
                        }
                    if len(cols) == 9:
                        clono["isotypes"] = cols[7]
                        clono["subclasses"] = cols[8]
                    clonos[call_hash][cdr3_length].append(clono)

    return clonos


def get_clonotypes_roi(samples):
    shared_clonotypes = {}
    shared_idx = 1
    combinaison_found = {}

    for samplename in samples.keys():
        if "ADN" in samples[samplename]:
            days = natsorted(samples[samplename]["ADN"].keys())
            
            if len(days) > 1:
                for i in range(len(days)):
                    for call_hash in samples[samplename]["ADN"][days[i]]["clonotypes"].keys():
                        for cdr3_length in samples[samplename]["ADN"][days[i]]["clonotypes"][call_hash]:
                            for clonotype in samples[samplename]["ADN"][days[i]]["clonotypes"][call_hash][cdr3_length]:
                                cdr3_list = clonotype["sequences"].split(";")
                                cdr3_seq = set()
                                for s in cdr3_list:
                                    cdr3_seq.add(s.split(",")[0][1:-1])
                                results = []
                                for j in range(i+1, len(days)):
                                    if i != j:
                                        if call_hash in samples[samplename]["ADN"][days[j]]["clonotypes"].keys():
                                            if cdr3_length in samples[samplename]["ADN"][days[j]]["clonotypes"][call_hash].keys():
                                                for clonotype_target in samples[samplename]["ADN"][days[j]]["clonotypes"][call_hash][cdr3_length]:
                                                    if hamming_distance(clonotype["cdr3"][1:-1], clonotype_target["cdr3"][1:-1]) <= MAX_DIST:
                                                        results.append(clonotype_target)
                                
                                
                                if len(results) > 0:
                                    if not call_hash in shared_clonotypes:
                                        shared_clonotypes[call_hash] = {}

                                    if not cdr3_length in shared_clonotypes[call_hash]:
                                        shared_clonotypes[call_hash][cdr3_length] = {}

                                    shared_clonotypes[call_hash][cdr3_length][shared_idx] = {}
                                    shared_clonotypes[call_hash][cdr3_length][shared_idx]["sequences"] = set()
                                    shared_clonotypes[call_hash][cdr3_length][shared_idx]["clonos"] = []
                                    shared_clonotypes[call_hash][cdr3_length][shared_idx]["clonos"].append(clonotype)
                                    most_freq_cdr3 = ""
                                    most_freq_cdr3_count = 0
                                    for cdr3 in clonotype["sequences"].split(";"):
                                        seq, count = cdr3.split(",")
                                        if int(count) > most_freq_cdr3_count:
                                            most_freq_cdr3_count = int(count)
                                            most_freq_cdr3 = seq
                                        shared_clonotypes[call_hash][cdr3_length][shared_idx]["sequences"].add(seq[1:-1])
                                    for c in results:
                                        shared_clonotypes[call_hash][cdr3_length][shared_idx]["clonos"].append(c)
                                        for cdr3 in c["sequences"].split(";"):
                                            seq, count = cdr3.split(",")
                                            if int(count) > most_freq_cdr3_count:
                                                most_freq_cdr3_count = int(count)
                                                most_freq_cdr3 = seq
                                            shared_clonotypes[call_hash][cdr3_length][shared_idx]["sequences"].add(seq[1:-1])
                                    shared_clonotypes[call_hash][cdr3_length][shared_idx]["most_freq_cdr3"] = most_freq_cdr3

                                    shared_idx += 1

    return shared_clonotypes



def get_all_mutation_rate(mutation_file, clones_file):
    reads_names = {}
    clones_id_seq = {}
    with gzip.open(clones_file, 'rt') as f:
        for i , line in enumerate(f):
            if i > 0:
                cols = line.split("\t")
                clone_id = cols[-2]
                seqname = cols[0]
                # if not clone_id in clones_id_seq:
                clones_id_seq[seqname] = clone_id

    mutations_rates = {}
    with gzip.open(mutation_file, 'rt') as f:
        for i, line in enumerate(f):
            if i > 0:
                cols = line.split("\t")
                if cols[1] in clones_id_seq:
                    clone_id = clones_id_seq[cols[1]]
                    
                    if not clone_id in reads_names:
                        reads_names[clone_id] = []
                    reads_names[clone_id].append(cols[1])
                    if not clone_id in mutations_rates:
                        mutations_rates[clone_id] = {
                            "mut_total" : [],
                            "mut_silent" : [],
                            "mut_nonsilent" : []
                        }
                        if len(cols) >= 5 and int(cols[5].split()[0]) > 0:
                            v_nucl = int(cols[5].split()[0])
                            v_mut_total = int(cols[7].split()[0]) / v_nucl * 100
                            v_mut_silent = int(cols[8].split()[0]) / v_nucl * 100
                            v_mut_nonsilent = int(cols[9].split()[0]) / v_nucl * 100
                            mutations_rates[clone_id]["mut_total"].append(v_mut_total)
                            mutations_rates[clone_id]["mut_silent"].append(v_mut_silent)
                            mutations_rates[clone_id]["mut_nonsilent"].append(v_mut_nonsilent)
    
    for clone_id in mutations_rates.keys():
        if len(mutations_rates[clone_id]["mut_total"]) > 0:
            mutations_rates[clone_id]["mut_total"] = statistics.mean(mutations_rates[clone_id]["mut_total"])
        if len(mutations_rates[clone_id]["mut_silent"]) > 0:
            mutations_rates[clone_id]["mut_silent"] = statistics.mean(mutations_rates[clone_id]["mut_silent"])
        if len(mutations_rates[clone_id]["mut_nonsilent"]) > 0:
            mutations_rates[clone_id]["mut_nonsilent"] = statistics.mean(mutations_rates[clone_id]["mut_nonsilent"])

    return mutations_rates, reads_names


def get_mutation_rate(sample, clone_ids_list):

    mut_total_list = []
    mut_silent_list = []
    mut_nonsilent_list = []

    if not "clonotypes" in sample:
        for day in clone_ids_list.keys():
            clones_ids_file = sample[day]["clonesfile"]
            mutation_file = sample[day]["mutationsfile"]
            reads_names = set()

            with gzip.open(clones_ids_file, 'rb') as f:
                for i, line in enumerate(f):
                    if i > 0:
                        cols = line.decode().strip().split("\t")
                        if int(cols[-2]) in clone_ids_list[day]:
                            reads_names.add(cols[0])

            # print(reads_names)
            with gzip.open(mutation_file, 'rb') as f:
                for i, line in enumerate(f):
                    if i > 0:
                        cols = line.decode().split("\t")

                        if cols[1] in reads_names:
                            if len(cols) >= 5 and int(cols[5].split()[0]) > 0:
                                v_nucl = int(cols[5].split()[0])
                                v_mut_total = int(cols[7].split()[0]) / v_nucl * 100
                                v_mut_silent = int(cols[8].split()[0]) / v_nucl * 100
                                v_mut_nonsilent = int(cols[9].split()[0]) / v_nucl * 100
                                mut_total_list.append(v_mut_total)
                                mut_silent_list.append(v_mut_silent)
                                mut_nonsilent_list.append(v_mut_nonsilent)
            

    else:
        clones_ids_file = sample["clonesfile"]
        mutation_file = sample["mutationsfile"]
        reads_names = set()
        with gzip.open(clones_ids_file, 'rb') as f:
            for i, line in enumerate(f):
                if i > 0:
                    cols = line.decode().strip().split("\t")
                    if int(cols[-2]) in clone_ids_list:
                        reads_names.add(cols[0])

        with gzip.open(mutation_file, 'rb') as f:
            for i, line in enumerate(f):
                if i > 0:
                    cols = line.decode().split("\t")

                    if cols[1] in reads_names:
                        if len(cols) >= 5 and int(cols[5].split()[0]) > 0:
                            v_nucl = int(cols[5].split()[0])
                            v_mut_total = int(cols[7].split()[0]) / v_nucl * 100
                            v_mut_silent = int(cols[8].split()[0]) / v_nucl * 100
                            v_mut_nonsilent = int(cols[9].split()[0]) / v_nucl * 100
                            mut_total_list.append(v_mut_total)
                            mut_silent_list.append(v_mut_silent)
                            mut_nonsilent_list.append(v_mut_nonsilent)
                            
    mut_total = 0
    mut_silent = 0
    mut_nonsilent = 0
    if len(mut_total_list) > 0:
        mut_total = statistics.mean(mut_total_list)
    if len(mut_silent_list) > 0:
        mut_silent = statistics.mean(mut_silent_list)
    if len(mut_nonsilent_list) > 0:
        mut_nonsilent = statistics.mean(mut_nonsilent_list)
    return mut_total, mut_silent, mut_nonsilent


def search_clonotypes(clono_queries, samples, controls):

    clonotypes_partages = {}

    nowrite = False
    table = {}

    groups_name = ["J0", "J7", "J14", "M4"]
    clono_traited = set()
    for call_hash in clono_queries.keys():
        for cdr3_length in clono_queries[call_hash].keys():
            for queryid, query in clono_queries[call_hash][cdr3_length].items():                

                cdr3_seq = query["sequences"]
                clonos = query["clonos"]
                most_freq_cdr3 = query["most_freq_cdr3"]
                outfile = most_freq_cdr3 + "_" + call_hash + "_" + str(queryid) + ".tsv"
                results = []
                sample_share_clono = set()
                for samplename in samples.keys():
                    if "ADN" in samples[samplename]:
                        for day in groups_name:
                            if not day in samples[samplename]["ADN"]:
                                results.append([samplename, samples[samplename]["id"], "patient_ADN", day, "nan", "nan", 0])
                            else:
                                result_found = False
                                if call_hash in samples[samplename]["ADN"][day]["clonotypes"] and cdr3_length in samples[samplename]["ADN"][day]["clonotypes"][call_hash]:
                                    for clonotype in samples[samplename]["ADN"][day]["clonotypes"][call_hash][cdr3_length]:
                                        if hamming_distance(most_freq_cdr3[1:-1], clonotype["cdr3"][1:-1]) <= MAX_DIST:
                                            result_found = True
                                            sample_share_clono.add(samplename)
                                            results.append([samplename, samples[samplename]["id"], "patient_ADN", day, clonotype["freq"], clonotype["count"], clonotype])
                                if not result_found:
                                    results.append([samplename, samples[samplename]["id"], "patient_ADN", day, 0, 0, 0])
                
                if len(sample_share_clono) >= MIN_PATIENT_NB:
                    sharedclono_hash = most_freq_cdr3 + "_" + call_hash
                    
                    if not sharedclono_hash in clono_traited:
                        clono_traited.add(sharedclono_hash)

            
                        # find in RNA samples
                        for samplename in samples.keys():
                            if "ARN" in samples[samplename]:
                                for day in groups_name:
                                    if not day in samples[samplename]["ARN"]:
                                        results.append([samplename, samples[samplename]["id"], "patient_ARN", day, "nan", "nan", 0])
                                    else:
                                        result_found = False
                                        if call_hash in samples[samplename]["ARN"][day]["clonotypes"] and cdr3_length in samples[samplename]["ARN"][day]["clonotypes"][call_hash]:
                                            for clonotype in samples[samplename]["ARN"][day]["clonotypes"][call_hash][cdr3_length]:
                                                        
                                                if hamming_distance(most_freq_cdr3[1:-1], clonotype["cdr3"][1:-1]) <= MAX_DIST:
                                                    result_found = True
                                                    results.append([samplename, samples[samplename]["id"], "patient_ARN", day, clonotype["freq"], clonotype["count"], clonotype])
                                        if not result_found:
                                            results.append([samplename, samples[samplename]["id"], "patient_ARN", day, 0, 0, 0])
                                    

                        # find in DNA AND RNA controls
                        for controlname in controls.keys():
                            if "ADN" in controls[controlname]:
                                result_found = False
                                if call_hash in controls[controlname]["ADN"]["clonotypes"] and cdr3_length in controls[controlname]["ADN"]["clonotypes"][call_hash]:
                                    for clonotype in controls[controlname]["ADN"]["clonotypes"][call_hash][cdr3_length]:
                                        if hamming_distance(most_freq_cdr3[1:-1], clonotype["cdr3"][1:-1]) <= MAX_DIST:
                                            result_found = True
                                            results.append([controlname, controls[controlname]["id"], "control_ADN", "", clonotype["freq"], clonotype["count"], clonotype])
                                if not result_found:
                                    results.append([controlname, controls[controlname]["id"], "control_ADN", "", 0, 0, 0])
                            
                            if "ARN" in controls[controlname]:
                                result_found = False
                                if call_hash in controls[controlname]["ARN"]["clonotypes"] and cdr3_length in controls[controlname]["ARN"]["clonotypes"][call_hash]:
                                    for clonotype in controls[controlname]["ARN"]["clonotypes"][call_hash][cdr3_length]:
                                        if hamming_distance(most_freq_cdr3[1:-1], clonotype["cdr3"][1:-1]) <= MAX_DIST:
                                            result_found = True
                                            results.append([controlname, controls[controlname]["id"], "control_ARN", "", clonotype["freq"], clonotype["count"], clonotype])
                                if not result_found:
                                    results.append([controlname, controls[controlname]["id"], "control_ARN", "", 0, 0, 0])
                        

                        for i, r in enumerate(results):
                            pn = r[0]
                            pid = r[1]
                            ptype = r[2]
                            pday = r[3]
                            if pday == "":
                                pday = "Ctrl"
                            cfreq = r[4]
                            ccount = r[5]
                            clonotype = r[6]
                            if clonotype != 0:
                                if not pn in clonotypes_partages:
                                    clonotypes_partages[pn] = {
                                        "id": pid,
                                        "type": ptype,
                                        "days": {}
                                    }
                                if not pday in clonotypes_partages[pn]["days"]:
                                    clonotypes_partages[pn]["days"][pday] = {}

                                clonotypes_partages[pn]["days"][pday][clonotype["hash"]] = set()
                                for j, r2 in enumerate(results):
                                    if i != j:
                                        pn2 = r2[0]
                                        pid2 = r2[1]
                                        ptype2 = r2[2]
                                        pday2 = r2[3]
                                        ptag = pn2 + "_" + pid2 + "_" + ptype2 + "_" + pday2
                                        clonotypes_partages[pn]["days"][pday][clonotype["hash"]].add(ptag)

                        outfile = sharedclono_hash + ".tsv"
                        with open(outfile, 'w') as o:
                            o.write("sample\tsample_anonyme\tcategory\tgroup\tvalue\tcolor\n")
                            for r in results:
                                o.write(r[0] + "\t" + r[1] + "\t" + r[2] + "\t" + r[3] + "\t" + str(r[4]) + "\t" + get_patients_clolor(r[0]) + "\n")

                        isooutfile = sharedclono_hash + "iso.tsv"
                        with open(isooutfile, 'w') as isout:
                            isout.write("sample\tsample_anonyme\tcategory\tday\tgroup\tvalue\n")
                            for r in results:
                                if r[2] == "patient_ARN":
                                    isotypes = {
                                        "IgM":0,
                                        "IgA":0,
                                        "IgG":0,
                                        "Undetermined":0
                                    }
                                    subclasses = {
                                        "IgG1":0,
                                        "IgG2":0,
                                        "IgG3":0,
                                        "IgG4":0,
                                        "pseudoIgG":0,
                                        "Unknow":0
                                    }
                                    if r[-1] != 0:
                                        isotypeString = r[-1]["isotypes"]
                                        for token in isotypeString.split(";"):
                                            key, value = token.split(":")
                                            isotypes[key] += int(value)
                                        subclassesString = r[-1]["subclasses"]
                                        for token in subclassesString.split(";"):
                                            key, value = token.split(":")
                                            subclasses[key] += int(value)
                                    for key, count in isotypes.items():
                                        isout.write(r[0] + "\t" + r[1] + "\tisotype\t" + r[3] + "\t" + key + "\t" + str(count) + "\n")
                                    for key, count in subclasses.items():
                                        isout.write(r[0] + "\t" + r[1] + "\tsubclasse\t" + r[3] + "\t" + key + "\t" + str(count) + "\n")

                        control_dna_traited = set()
                        control_rna_traited = set()
                        for clono in results:
                            if clono[2].startswith("control") and clono[-1] != 0:
                                if clono[2] == "control_ADN":
                                    control_dna_traited.add(clono[0])
                                elif clono[2] == "control_ARN":
                                    control_rna_traited.add(clono[0])
                            
                            
                                    
                        for samplename in samples.keys():
                            if "ADN" in samples[samplename]:
                                sample_tag = samplename + "_ADN"
                                
                                if not sample_tag in table:
                                    table[sample_tag] = {}
                                    table[sample_tag]["id"] = samples[samplename]["id"]
                                    table[sample_tag]["type"] = "ADN"
                                    table[sample_tag]["clonos"] = {}
                                
                                # init clono columns
                                table[sample_tag]["clonos"][sharedclono_hash] = {}
                                table[sample_tag]["clonos"][sharedclono_hash]["freq"] = ""
                                table[sample_tag]["clonos"][sharedclono_hash]["days"] = []
                                table[sample_tag]["clonos"][sharedclono_hash]["isotypes"] = ""
                                table[sample_tag]["clonos"][sharedclono_hash]["subclasses"] = ""
                                table[sample_tag]["clonos"][sharedclono_hash]["v"] = ""
                                table[sample_tag]["clonos"][sharedclono_hash]["j"] = ""
                                table[sample_tag]["clonos"][sharedclono_hash]["mutation_rate_total"] = ""
                                table[sample_tag]["clonos"][sharedclono_hash]["mutation_rate_silent"] = ""
                                table[sample_tag]["clonos"][sharedclono_hash]["mutation_rate_nonsilent"] = ""
                                table[sample_tag]["clonos"][sharedclono_hash]["nb_DNA_ctrl"] = 0
                                table[sample_tag]["clonos"][sharedclono_hash]["nb_RNA_ctrl"] = 0

                                days = {"J0":"", "J7":"", "J14":"", "M4":""}
                                freq = []
                                reads = []
                                table[sample_tag]["clonos"][sharedclono_hash]["nb_DNA_ctrl"] = len(control_dna_traited)
                                table[sample_tag]["clonos"][sharedclono_hash]["nb_RNA_ctrl"] = len(control_rna_traited)

                                clono_in_patient = {}

                                for clono in results:
                                    if samplename == clono[0] and clono[2] == "patient_ADN":
                                        if clono[-1] != 0:
                                            if not clono[3] in clono_in_patient:
                                                clono_in_patient[clono[3]] = set()
                                            clono_in_patient[clono[3]].add(clono[-1]["id"])

                                            days[clono[3]] = "OK"
                                            table[sample_tag]["clonos"][sharedclono_hash]["v"] = clono[-1]["v"]
                                            table[sample_tag]["clonos"][sharedclono_hash]["j"] = clono[-1]["j"]
                                            freq.append(clono[4])
                                            reads.append(clono[5])

                                        elif clono[4] == "nan":
                                            days[clono[3]] = "NaN"
                                        else:
                                            days[clono[3]] = "ABS"

                                if len(freq) > 0:
                                    table[sample_tag]["clonos"][sharedclono_hash]["freq"] = statistics.mean(freq)
                                
                                if len(freq) > 0:
                                    table[sample_tag]["clonos"][sharedclono_hash]["reads"] = statistics.mean(reads)

                                for d in ["J0", "J7", "J14", "M4"]:
                                    table[sample_tag]["clonos"][sharedclono_hash]["days"].append(days[d])

                

                                if len(clono_in_patient) > 0:
                                    mut_rate_total, mut_rate_silent, mut_rate_nonsilent = get_mutation_rate(samples[samplename]["ADN"], clono_in_patient)
                                    table[sample_tag]["clonos"][sharedclono_hash]["mutation_rate_total"] = mut_rate_total
                                    table[sample_tag]["clonos"][sharedclono_hash]["mutation_rate_silent"] = mut_rate_silent
                                    table[sample_tag]["clonos"][sharedclono_hash]["mutation_rate_nonsilent"] = mut_rate_nonsilent
                


                            if "ARN" in samples[samplename]:
                                sample_tag = samplename + "_ARN"
                                
                                if not sample_tag in table:
                                    table[sample_tag] = {}
                                    table[sample_tag]["id"] = samples[samplename]["id"]
                                    table[sample_tag]["type"] = "ARN"
                                    table[sample_tag]["clonos"] = {}
                                
                                # init clono columns
                                table[sample_tag]["clonos"][sharedclono_hash] = {}
                                table[sample_tag]["clonos"][sharedclono_hash]["freq"] = ""
                                table[sample_tag]["clonos"][sharedclono_hash]["days"] = []
                                table[sample_tag]["clonos"][sharedclono_hash]["isotypes"] = ""
                                table[sample_tag]["clonos"][sharedclono_hash]["subclasses"] = ""
                                table[sample_tag]["clonos"][sharedclono_hash]["v"] = ""
                                table[sample_tag]["clonos"][sharedclono_hash]["j"] = ""
                                table[sample_tag]["clonos"][sharedclono_hash]["mutation_rate_total"] = ""
                                table[sample_tag]["clonos"][sharedclono_hash]["mutation_rate_silent"] = ""
                                table[sample_tag]["clonos"][sharedclono_hash]["mutation_rate_nonsilent"] = ""
                                table[sample_tag]["clonos"][sharedclono_hash]["nb_DNA_ctrl"] = 0
                                table[sample_tag]["clonos"][sharedclono_hash]["nb_RNA_ctrl"] = 0

                                days = {"J0":"", "J7":"", "J14":"", "M4":""}
                                freq = []
                                reads = []
                                isotypes = {"IgA":0, "IgG":0, "IgM":0}
                                subclasses = {"IgG1":0, "IgG2":0, "IgG3":0,"IgG4":0,"pseudoIgG":0}

                                
                                table[sample_tag]["clonos"][sharedclono_hash]["nb_DNA_ctrl"] = len(control_dna_traited)
                                table[sample_tag]["clonos"][sharedclono_hash]["nb_RNA_ctrl"] = len(control_rna_traited)
                                
                                clono_in_patient = {}

                                for clono in results:
                                    if samplename == clono[0] and clono[2] == "patient_ARN":
                                        if clono[-1] != 0:
                                            
                                            if not clono[3] in clono_in_patient:
                                                clono_in_patient[clono[3]] = set()
                                            clono_in_patient[clono[3]].add(clono[-1]["id"])
                                            
                                            days[clono[3]] = "OK"
                                            table[sample_tag]["clonos"][sharedclono_hash]["v"] = clono[-1]["v"]
                                            table[sample_tag]["clonos"][sharedclono_hash]["j"] = clono[-1]["j"]
                                            freq.append(clono[4])
                                            reads.append(clono[5])
                                            if "isotypes" in clono[-1]:
                                                for isotoken in clono[-1]["isotypes"].split(";"):
                                                    iso, nb = isotoken.split(":")
                                                    if iso in isotypes:
                                                        isotypes[iso] += int(nb)
                                            if "subclasses" in clono[-1]:
                                                for subctoken in clono[-1]["subclasses"].split(";"):
                                                    subc, nb = subctoken.split(":")
                                                    if subc in subclasses:
                                                        subclasses[subc] += int(nb)
                                        elif clono[4] == "nan":
                                            days[clono[3]] = "NaN"
                                        else:
                                            days[clono[3]] = "ABS"

                                if len(freq) > 0:
                                    table[sample_tag]["clonos"][sharedclono_hash]["freq"] = statistics.mean(freq)
                                    table[sample_tag]["clonos"][sharedclono_hash]["reads"] = statistics.mean(reads)
                                    d = ""
                                    for iso in ["IgM", "IgA", "IgG"]:
                                        table[sample_tag]["clonos"][sharedclono_hash]["isotypes"] += d + iso + ":" + str(isotypes[iso])
                                        d = ";"
                                    d = ""
                                    for subc in ["IgG1", "IgG2", "IgG3", "IgG4", "pseudoIgG"]:
                                        table[sample_tag]["clonos"][sharedclono_hash]["subclasses"] += d + subc + ":" + str(subclasses[subc])
                                        d = ";"

                                for d in ["J0", "J7", "J14", "M4"]:
                                    table[sample_tag]["clonos"][sharedclono_hash]["days"].append(days[d])

                                if len(clono_in_patient) > 0:
                                    mut_rate_total, mut_rate_silent, mut_rate_nonsilent = get_mutation_rate(samples[samplename]["ARN"], clono_in_patient)
                                    table[sample_tag]["clonos"][sharedclono_hash]["mutation_rate_total"] = mut_rate_total
                                    table[sample_tag]["clonos"][sharedclono_hash]["mutation_rate_silent"] = mut_rate_silent
                                    table[sample_tag]["clonos"][sharedclono_hash]["mutation_rate_nonsilent"] = mut_rate_nonsilent


                        for controlname in controls.keys():
                            if "ADN" in controls[controlname]:
                                control_tag = controlname + "_ADN"
                                
                                if not control_tag in table:
                                    table[control_tag] = {}
                                    table[control_tag]["id"] = controls[controlname]["id"]
                                    table[control_tag]["type"] = "ADN"
                                    table[control_tag]["clonos"] = {}
                                
                                # init clono columns
                                table[control_tag]["clonos"][sharedclono_hash] = {}
                                table[control_tag]["clonos"][sharedclono_hash]["freq"] = ""
                                table[control_tag]["clonos"][sharedclono_hash]["days"] = []
                                table[control_tag]["clonos"][sharedclono_hash]["isotypes"] = ""
                                table[control_tag]["clonos"][sharedclono_hash]["subclasses"] = ""
                                table[control_tag]["clonos"][sharedclono_hash]["v"] = ""
                                table[control_tag]["clonos"][sharedclono_hash]["j"] = ""
                                table[control_tag]["clonos"][sharedclono_hash]["mutation_rate_total"] = ""
                                table[control_tag]["clonos"][sharedclono_hash]["mutation_rate_silent"] = ""
                                table[control_tag]["clonos"][sharedclono_hash]["mutation_rate_nonsilent"] = ""
                                table[control_tag]["clonos"][sharedclono_hash]["nb_DNA_ctrl"] = 0
                                table[control_tag]["clonos"][sharedclono_hash]["nb_RNA_ctrl"] = 0
            
                                table[control_tag]["clonos"][sharedclono_hash]["nb_DNA_ctrl"] = len(control_dna_traited)
                                table[control_tag]["clonos"][sharedclono_hash]["nb_RNA_ctrl"] = len(control_rna_traited)
                                
                                freq = []
                                reads = []
                                clono_in_patient = set()
                                for clono in results:
                                    if controlname == clono[0] and clono[2] == "control_ADN":
                                        if clono[-1] != 0:
                                            clono_in_patient.add(clono[-1]["id"])

                                            table[control_tag]["clonos"][sharedclono_hash]["v"] = clono[-1]["v"]
                                            table[control_tag]["clonos"][sharedclono_hash]["j"] = clono[-1]["j"]
                                            freq.append(clono[4])
                                            reads.append(clono[5])
                                
                                if len(freq) > 0:
                                    table[control_tag]["clonos"][sharedclono_hash]["freq"] = statistics.mean(freq)
                                    table[control_tag]["clonos"][sharedclono_hash]["reads"] = statistics.mean(reads)

                                for d in ["J0", "J7", "J14", "M4"]:
                                    table[control_tag]["clonos"][sharedclono_hash]["days"].append("")

                                if len(clono_in_patient) > 0:
                                    mut_rate_total, mut_rate_silent, mut_rate_nonsilent = get_mutation_rate(controls[controlname]["ADN"], clono_in_patient)
                                    table[control_tag]["clonos"][sharedclono_hash]["mutation_rate_total"] = mut_rate_total
                                    table[control_tag]["clonos"][sharedclono_hash]["mutation_rate_silent"] = mut_rate_silent
                                    table[control_tag]["clonos"][sharedclono_hash]["mutation_rate_nonsilent"] = mut_rate_nonsilent


                            if "ARN" in controls[controlname]:
                                control_tag = controlname + "_ARN"
                                
                                if not control_tag in table:
                                    table[control_tag] = {}
                                    table[control_tag]["id"] = controls[controlname]["id"]
                                    table[control_tag]["type"] = "ARN"
                                    table[control_tag]["clonos"] = {}
                                
                                # init clono columns
                                table[control_tag]["clonos"][sharedclono_hash] = {}
                                table[control_tag]["clonos"][sharedclono_hash]["freq"] = ""
                                table[control_tag]["clonos"][sharedclono_hash]["days"] = []
                                table[control_tag]["clonos"][sharedclono_hash]["isotypes"] = ""
                                table[control_tag]["clonos"][sharedclono_hash]["subclasses"] = ""
                                table[control_tag]["clonos"][sharedclono_hash]["v"] = ""
                                table[control_tag]["clonos"][sharedclono_hash]["j"] = ""
                                table[control_tag]["clonos"][sharedclono_hash]["mutation_rate_total"] = ""
                                table[control_tag]["clonos"][sharedclono_hash]["mutation_rate_silent"] = ""
                                table[control_tag]["clonos"][sharedclono_hash]["mutation_rate_nonsilent"] = ""
                                table[control_tag]["clonos"][sharedclono_hash]["nb_DNA_ctrl"] = 0
                                table[control_tag]["clonos"][sharedclono_hash]["nb_RNA_ctrl"] = 0
            
                                table[control_tag]["clonos"][sharedclono_hash]["nb_DNA_ctrl"] = len(control_dna_traited)
                                table[control_tag]["clonos"][sharedclono_hash]["nb_RNA_ctrl"] = len(control_rna_traited)
                                
                                freq = []
                                reads = []
                                isotypes = {"IgA":0, "IgG":0, "IgM":0}
                                subclasses = {"IgG1":0, "IgG2":0, "IgG3":0,"IgG4":0,"pseudoIgG":0}
                                
                                clono_in_patient = set()
                                for clono in results:
                                    if controlname == clono[0]:
                                        if clono[-1] != 0 and clono[2] == "control_ARN":
                                            clono_in_patient.add(clono[-1]["id"])

                                            table[control_tag]["clonos"][sharedclono_hash]["v"] = clono[-1]["v"]
                                            table[control_tag]["clonos"][sharedclono_hash]["j"] = clono[-1]["j"]
                                            freq.append(clono[4])
                                            reads.append(clono[5])
                                            if "isotypes" in clono[-1]:
                                                for isotoken in clono[-1]["isotypes"].split(";"):
                                                    iso, nb = isotoken.split(":")
                                                    if iso in isotypes:
                                                        isotypes[iso] += int(nb)
                                            if "subclasses" in clono[-1]:
                                                for subctoken in clono[-1]["subclasses"].split(";"):
                                                    subc, nb = subctoken.split(":")
                                                    if subc in subclasses:
                                                        subclasses[subc] += int(nb)
                                
                                if len(freq) > 0:
                                    table[control_tag]["clonos"][sharedclono_hash]["freq"] = statistics.mean(freq)
                                    table[control_tag]["clonos"][sharedclono_hash]["reads"] = statistics.mean(reads)

                                    d = ""
                                    for iso in ["IgM", "IgA", "IgG"]:
                                        table[control_tag]["clonos"][sharedclono_hash]["isotypes"] += iso + ":" + str(isotypes[iso]) + d
                                        d = ";"
                                    d = ""
                                    for subc in ["IgG1", "IgG2", "IgG3", "IgG4", "pseudoIgG"]:
                                        table[control_tag]["clonos"][sharedclono_hash]["subclasses"] += subc + ":" + str(subclasses[subc]) + d
                                        d = ";"
                        
                                for d in ["J0", "J7", "J14", "M4"]:
                                    table[control_tag]["clonos"][sharedclono_hash]["days"].append("")
                                
                                if len(clono_in_patient) > 0:
                                    mut_rate_total, mut_rate_silent, mut_rate_nonsilent = get_mutation_rate(controls[controlname]["ARN"], clono_in_patient)
                                    table[control_tag]["clonos"][sharedclono_hash]["mutation_rate_total"] = mut_rate_total
                                    table[control_tag]["clonos"][sharedclono_hash]["mutation_rate_silent"] = mut_rate_silent
                                    table[control_tag]["clonos"][sharedclono_hash]["mutation_rate_nonsilent"] = mut_rate_nonsilent

    header1 = ["", "", ""]
    header2 = ["Patient", "id", "type"]

    if nowrite:
        return

    with open("table.tsv", 'w') as tableout:
        for i, pn in enumerate(table.keys()):
            cols = []
            cols.append(pn)
            cols.append(table[pn]["id"])
            cols.append(table[pn]["type"])

            for clonohash in table[pn]["clonos"].keys():
                    header1.append(clonohash)
                    for j in range(13):
                        header1.append("")
                    header2.append("J0")
                    header2.append("J7")
                    header2.append("J14")
                    header2.append("M4")
                    header2.append("Freq")
                    header2.append("Reads")
                    header2.append("V")
                    header2.append("J")
                    header2.append("Isotypes")
                    header2.append("Subclasses")
                    header2.append("Mutation rate total")
                    header2.append("Mutation rate silent")
                    header2.append("Mutation rate non silent")
                    header2.append("Nb control ADN")
                    header2.append("Nb control ARN")

                    for d in table[pn]["clonos"][clonohash]["days"]:
                        cols.append(d)
                    
                    cols.append(table[pn]["clonos"][clonohash]["freq"])

                    if 'reads' in table[pn]["clonos"][clonohash]:
                        cols.append(table[pn]["clonos"][clonohash]["reads"])
                    else:
                        cols.append("")
                    cols.append(table[pn]["clonos"][clonohash]["v"])
                    cols.append(table[pn]["clonos"][clonohash]["j"])
                    cols.append(table[pn]["clonos"][clonohash]["isotypes"])
                    cols.append(table[pn]["clonos"][clonohash]["subclasses"])
                    cols.append(table[pn]["clonos"][clonohash]["mutation_rate_total"])
                    cols.append(table[pn]["clonos"][clonohash]["mutation_rate_silent"])
                    cols.append(table[pn]["clonos"][clonohash]["mutation_rate_nonsilent"])
                    cols.append(table[pn]["clonos"][clonohash]["nb_DNA_ctrl"])
                    cols.append(table[pn]["clonos"][clonohash]["nb_RNA_ctrl"])

            if i == 0:
                tableout.write("\t".join(header1) + "\n")
                tableout.write("\t".join(header2) + "\n")
            tableout.write("\t".join(map(str, cols)) + "\n")
            
    # write clonotypes
    with open("clonotypes_totaux.tsv", 'w') as o:
        header = ["#samplename", "barcode", "group", "type", "day", "V", "J", "Junction", "cdr3_length", "isotypes", "subclasse", "reads", "freq", "cdr3_list", "total_clonos_count", "total_reads_count", "mutation_rate_total (%)", "mutation_rate_silent (%)", "mutation_rate_nonsilent (%)", "shared", "shared_list"]
        o.write("\t".join(header) + "\n")

        for samplename in samples.keys():
            for stype in ["ADN", "ARN"]:
                if stype in samples[samplename]:
                    for day in samples[samplename][stype].keys():
                        clonotypes = samples[samplename][stype][day]["clonotypes"]
                        clonotypes_count = 0
                        reads_count = 0
                        for call_hash in clonotypes.keys():
                            for cdr3_length in clonotypes[call_hash].keys():
                                clonotypes_count += len(clonotypes[call_hash][cdr3_length])
                                for clono in clonotypes[call_hash][cdr3_length]:
                                    reads_count += clono["count"]
                        for call_hash in clonotypes.keys():
                            for cdr3_length in clonotypes[call_hash].keys():
                                
                                for clono in clonotypes[call_hash][cdr3_length]:
                                    mut_total = samples[samplename][stype][day]["mutation_rates"][str(clono["id"])]["mut_total"]
                                    mut_silent = samples[samplename][stype][day]["mutation_rates"][str(clono["id"])]["mut_silent"]
                                    mut_nonsilent = samples[samplename][stype][day]["mutation_rates"][str(clono["id"])]["mut_nonsilent"]
                                    reads_names = ";".join(samples[samplename][stype][day]["read_names"][str(clono["id"])])
                                    clono_to_write = []
                                    clono_to_write.append(samplename)
                                    clono_to_write.append(samples[samplename]["id"])
                                    clono_to_write.append("Patient")
                                    clono_to_write.append(stype)
                                    clono_to_write.append(day)
                                    clono_to_write.append(clono["v"])
                                    clono_to_write.append(clono["j"].split("*")[0])
                                    clono_to_write.append(clono["cdr3"])
                                    clono_to_write.append(len(clono["cdr3"])-2)
                                    if "isotypes" in clono:
                                        clono_to_write.append(clono["isotypes"])
                                    else:
                                        clono_to_write.append("na")
                                    if "subclasses" in clono:
                                        clono_to_write.append(clono["subclasses"])
                                    else:
                                        clono_to_write.append("na")
                                    clono_to_write.append(clono["count"])
                                    clono_to_write.append(clono["freq"])
                                    clono_to_write.append(clono["sequences"])
                                    clono_to_write.append(clonotypes_count)
                                    clono_to_write.append(reads_count)
                                    clono_to_write.append(round(mut_total, 5))
                                    clono_to_write.append(round(mut_silent, 5))
                                    clono_to_write.append(round(mut_nonsilent, 5))
                                    shared_list = ""
                                    if samplename in clonotypes_partages:
                                        if day in clonotypes_partages[samplename]["days"]:
                                            if clono["hash"] in clonotypes_partages[samplename]["days"][day]:
                                                shared_list = ";".join(clonotypes_partages[samplename]["days"][day][clono["hash"]])
                                    if len(shared_list) > 0:
                                        clono_to_write.append("1")
                                        clono_to_write.append(shared_list)
                                    else:
                                        clono_to_write.append("0")
                                        clono_to_write.append("na")
                                    clono_to_write.append(reads_names)
                    
                                    o.write("\t".join(map(str, clono_to_write)) + "\n")

        for samplename in controls.keys():
            for stype in ["ADN", "ARN"]:
                if stype in controls[samplename]:
                    clonotypes = controls[samplename][stype]["clonotypes"]
                    clonotypes_count = 0
                    reads_count = 0
                    for call_hash in clonotypes.keys():
                        for cdr3_length in clonotypes[call_hash].keys():
                            clonotypes_count += len(clonotypes[call_hash][cdr3_length])
                            for clono in clonotypes[call_hash][cdr3_length]:
                                reads_count += clono["count"]
                    for call_hash in clonotypes.keys():
                        for cdr3_length in clonotypes[call_hash].keys():
                            
                            for clono in clonotypes[call_hash][cdr3_length]:
                                mut_total = controls[samplename][stype]["mutation_rates"][str(clono["id"])]["mut_total"]
                                mut_silent = controls[samplename][stype]["mutation_rates"][str(clono["id"])]["mut_silent"]
                                mut_nonsilent = controls[samplename][stype]["mutation_rates"][str(clono["id"])]["mut_nonsilent"]
                                reads_names = ";".join(controls[samplename][stype]["read_names"][str(clono["id"])])
                                clono_to_write = []
                                clono_to_write.append(samplename)
                                clono_to_write.append(controls[samplename]["id"])
                                clono_to_write.append("Control")
                                clono_to_write.append(stype)
                                clono_to_write.append("na")
                                clono_to_write.append(clono["v"])
                                clono_to_write.append(clono["j"].split("*")[0])
                                clono_to_write.append(clono["cdr3"])
                                clono_to_write.append(len(clono["cdr3"])-2)
                                if "isotypes" in clono:
                                    clono_to_write.append(clono["isotypes"])
                                else:
                                    clono_to_write.append("na")
                                if "subclasses" in clono:
                                    clono_to_write.append(clono["subclasses"])
                                else:
                                    clono_to_write.append("na")
                                clono_to_write.append(clono["count"])
                                clono_to_write.append(clono["freq"])
                                clono_to_write.append(clono["sequences"])
                                clono_to_write.append(clonotypes_count)
                                clono_to_write.append(reads_count)
                                clono_to_write.append(round(mut_total, 5))
                                clono_to_write.append(round(mut_silent, 5))
                                clono_to_write.append(round(mut_nonsilent, 5))
                                shared_list = ""
                                
                                if samplename in clonotypes_partages:
                                    if "Ctrl" in clonotypes_partages[samplename]["days"]:
                                        if clono["hash"] in clonotypes_partages[samplename]["days"]["Ctrl"]:
                                            shared_list = ";".join(clonotypes_partages[samplename]["days"]["Ctrl"][clono["hash"]])
                                if len(shared_list) > 0:
                                    clono_to_write.append("1")
                                    clono_to_write.append(shared_list)
                                else:
                                    clono_to_write.append("0")
                                    clono_to_write.append("na")
                                clono_to_write.append(reads_names)
                                o.write("\t".join(map(str, clono_to_write)) + "\n")


def main():                 
    global MAX_DIFF
    global MAX_DIST
    global MIN_PATIENT_NB
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--samples", type=str, help="Samples file. TSV file with columns : group\tname\tid\tday\timgt_dir.\n\nExample : Patient_RNA\tAA01\tpatient01\tJ7\t/path/to//patient01/imgtdir\nControl_DNA\tcc01\tControl01\tNA\t/path/to/control01/imgt_dir/", required=True)
    parser.add_argument("-d", "--max-diff", type=int, help="Minimum mismatch for 10 aa allowed between 2 clonotypes. Default : 1", default=1)
    parser.add_argument("-n", "--min-patient-number", type=int, help="Minimum patient number that share a clonotype", default=4)

    args = parser.parse_args()

    MAX_DIFF = args.max_diff
    MAX_DIST = MAX_DIFF / 10
    MIN_PATIENT_NB = args.min_patient_number

    parse_sample_file(args.samples)
    samples = get_samples()
    controls = get_controls()
    

    clono_queries = get_clonotypes_roi(samples)
    search_clonotypes(clono_queries, samples, controls)

main()