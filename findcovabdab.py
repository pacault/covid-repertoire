#!/usr/bin/env python
# coding: utf8`

import argparse
import sys, os
import gzip
import csv


MIN_READ = 5
CONTROLES_ADN = []
PATIENTS_ADN = []
PATIENTS_ARN = []
CONTROLES_ARN = []


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


def get_samples():
    samples = {}
    for sample in PATIENTS_ADN:
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
            if file.endswith("_clonotypes2.tsv.gz"):
                samples[sample["name"]]["ADN"][sample["day"]]["clonofile"] = os.path.join(samples[sample["name"]]["ADN"][sample["day"]]["imgt_dir"], file)
        samples[sample["name"]]["ADN"][sample["day"]]["clonotypes"] = load_clonotypes(samples[sample["name"]]["ADN"][sample["day"]]["clonofile"], False)


    for sample in PATIENTS_ARN:
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
            if file.endswith("_clonotypes2.tsv.gz") or file.endswith("_clonotypes_changeo.tsv.gz"):
                samples[sample["name"]]["ARN"][sample["day"]]["clonofile"] = os.path.join(samples[sample["name"]]["ARN"][sample["day"]]["imgt_dir"], file)
        samples[sample["name"]]["ARN"][sample["day"]]["clonotypes"] = load_clonotypes(samples[sample["name"]]["ARN"][sample["day"]]["clonofile"])

    return samples


def get_controls():
    controls = {}
    for sample in CONTROLES_ADN:
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
        controls[sample["name"]]["ADN"]["clonotypes"] = load_clonotypes(controls[sample["name"]]["ADN"]["clonofile"], False)

    for sample in CONTROLES_ARN_MIG1:
        # print(sample["imgt_dir"])
        if not sample["name"] in controls:
            controls[sample["name"]] = {}
            controls[sample["name"]]["name"] = sample["name"]
            controls[sample["name"]]["id"] = sample["id"]
        if not "ARN" in controls[sample["name"]]:
            controls[sample["name"]]["ARN"] = {}
        controls[sample["name"]]["ARN"]["imgt_dir"] = sample["imgt_dir"]
        for file in os.listdir(controls[sample["name"]]["ARN"]["imgt_dir"]):
            if file.endswith("_clonotypes2.tsv.gz") or file.endswith("_clonotypes_changeo.tsv.gz"):
                controls[sample["name"]]["ARN"]["clonofile"] = os.path.join(controls[sample["name"]]["ARN"]["imgt_dir"], file)
        controls[sample["name"]]["ARN"]["clonotypes"] = load_clonotypes(controls[sample["name"]]["ARN"]["clonofile"])

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
                            "cdr3":cdr3,
                            "count":int(cols[2]),
                            "freq":float(cols[3]),
                            "v":cols[4],
                            "j":cols[5],
                            "sequences": cols[6]                            
                        }
                    if len(cols) == 9:
                        clono["isotypes"] = cols[7]
                        clono["subclasses"] = cols[8]
                    clonos[call_hash][cdr3_length].append(clono)

    return clonos



def main():
    parser = argparse.ArgumentParser("Recherche des clonotypes communs entre les repertoires covid et covabdab. (Même V, même J CDR3 10%)")
    parser.add_argument("-s", "--samples", type=str, help="Samples file. TSV file with columns : group\tname\tid\tday\timgt_dir.\n\nExample : Patient_RNA\tAA01\tpatient01\tJ7\t/path/to//patient01/imgtdir\nControl_DNA\tcc01\tControl01\tNA\t/path/to/control01/imgt_dir/", required=True)
    parser.add_argument("-i", "--covabdabfile", type=str, help="IMGT High VQuest output dir", required=True)
    args = parser.parse_args()

    parse_sample_file(args.samples)
    
    samples = get_samples()
    controls = get_controls()

    groups_name = ["J0", "J7", "J14", "M4"]
    types = ["ADN", "ARN"]
    
    with open(args.covabdabfile, 'r') as f:
        csv_reader = csv.reader(f, delimiter=',')
        for i, row in enumerate(csv_reader):
            if i > 0:
                if row[10].strip() != "" and row[11].strip() != "" and row[14].strip() != "":
                    if "Human" in row[10] and "Human" in row[11]:
                        v = row[10].split()[0]
                        j = row[11].split()[0]
                        cdr3 = row[14]
                        ref = row[18]
                        name = row[0]
                        neutrilisingVs = "Neutralising vs " + row[4] if row[4] != "" else ""
                        notNeutrilisingVs = "Not neutralising vs " + row[5] if row[5] != "" else ""
                        hashcall = v + "_" + j

                        complete_hash = hashcall + "_" + cdr3
                        results = []

                        for samplename in samples.keys():
                            for sample_type in types:
                                if sample_type in samples[samplename]:
                                    for day in samples[samplename][sample_type].keys():
                                        clonotypes = samples[samplename][sample_type][day]["clonotypes"]
                                        if hashcall in clonotypes:
                                            if len(cdr3) in clonotypes[hashcall]:
                                                for clono in clonotypes[hashcall][len(cdr3)]:
                                                    if hamming_distance(cdr3, clono["cdr3"]) <= 0.1:
                                                        results.append(
                                                            [samplename, sample_type, day, clono]
                                                        )
                        
                        if len(results) > 0:
                            print("#" + name + "\t"+ neutrilisingVs + "\t" + notNeutrilisingVs + "\t" + v + "\t" + j + "\t" + cdr3 + "\t" + ref)
                            for r in results:
                                isotypes = ""
                                subclasses = ""
                                if "isotypes" in r[3]:
                                    isotypes = r[3]["isotypes"]
                                if "subclasses" in r[3]:
                                    subclasses = r[3]["subclasses"]
                                print("\t" + r[0] + "\t" + r[1] + "\t" + r[2] + "\t" + r[3]["cdr3"] + "\t" + str(r[3]["count"]) + " ("+ str(round(r[3]["freq"]*100, 2)) + "%)\t" + isotypes + "\t" + subclasses)
                        


if __name__ == '__main__':  main()
