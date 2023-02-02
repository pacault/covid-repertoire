#!/usr/bin/env python

import argparse
import os, sys
import subprocess
import tempfile
import docker
from argparse import RawTextHelpFormatter
from pathlib import Path


DISTANCE_MODEL = ""

def runCommand(cmd, silent=False):
    if not silent:
        print(cmd)
    child = subprocess.Popen(cmd , shell=True, stdout=subprocess.PIPE, executable='/bin/bash')
    output, err = child.communicate()
    output = output.decode('utf-8').strip()
    return output, err


def run_job (cmd_line, errorMess, silent = False):
    """
        Run an program.

        :param frameinfo: From the inspect module. Allow to print the line number of the script in a error message.
        :param cmd_line: The command line to run.
        :type cmd_line: str
        :errorMess: A message error.
        :type errorM: str
        :param silent: Print or not the commande line in the standard output.
        :type silent: bool
        :return: void
    """
    if not silent:
        print (cmd_line)
    try:
        tmp = tempfile.NamedTemporaryFile().name
        error = open(tmp, 'w')
        proc = subprocess.Popen( args=cmd_line, shell=True, stderr=error)
        returncode = proc.wait()
        error.close()
        error = open( tmp, 'rb' )
        stderr = ''
        buffsize = 1048576
        try:
            while True:
                stderr += str(error.read( buffsize ))
                if not stderr or len( stderr ) % buffsize != 0:
                    break
        except OverflowError:
            pass
        error.close()
        os.remove(tmp)
        if returncode != 0:
            raise Exception(stderr)
    except Exception as e:
        print (e)
        raise Exception(errorMess)


def get_samples_from_config(config_file):
    samples = {}
    with open(config_file, 'r') as f:
        for i, line in enumerate(f):
            if not line.startswith('#'):
                if line.strip():
                    cols = line.rstrip().split('\t')
                    sample_name = cols[0]
                    
                    if len(cols) != 3 and len(cols) != 2:
                        raise ValueError(F"Bad configuration file, line {i+1}, {line.rstrip()}")
                    if sample_name in samples:
                        raise ValueError(F"Duplicate sample name : {sample_name}")
                    imgt_dir = cols[1]
                    if not os.path.exists(imgt_dir):
                        raise ValueError(F"The imgt output dir {imgt_dir} doesn't exists")
                    if not os.path.isdir(imgt_dir):
                        raise ValueError(F"The imgt output dir {imgt_dir} isn't a valid directory")

                    if len(cols) == 3:
                        fasta_file = cols[2]
                    else:
                        ff = ""
                        for file in os.listdir(imgt_dir):
                            if file.endswith(".fasta"):
                                if ff != "":
                                    raise ValueError(F"Multiple fasta files found in {imgt_dir} directory")
                                else:
                                    ff = os.path.join(imgt_dir, file)
                        fasta_file = ff

                    if not os.path.exists(fasta_file):
                        print(line)
                        raise ValueError(F"The fasta file {fasta_file} doesn't exists")
                        
                    samples[sample_name] = {"imgt" : imgt_dir.replace(" ", "\\ "), "fasta" : fasta_file.replace(" ", "\\ ")}

    return samples


def worker(sample_name, sample_files):
    """ pipeline worker"""
    
    # MakeDB
    print(F"### Parsing IMGT output of {sample_name}")
    outfile_makedb = os.path.join(sample_files['imgt'], F"{sample_name}_db-pass.tsv")
    cmd = F"MakeDb.py imgt -i {sample_files['imgt']} -s {sample_files['fasta']} -o {outfile_makedb} --extended"
    if not os.path.exists(outfile_makedb):
        run_job(cmd, F"error when parsing {sample_name}")

    # ParseDb
    print("### Filter out unproductives sequences")
    outfile_parsedb_name = F"{sample_name}_db-pass_prod.tsv"
    outfile_parsedb = os.path.join(sample_files['imgt'], outfile_parsedb_name)
    cmd = F"ParseDb.py select -d {outfile_makedb} -o {outfile_parsedb} -f productive -u T"
    if not os.path.exists(outfile_parsedb):
        run_job(cmd, F"error when filter out non productive sequences for {sample_name}")
    # if os.path.exists(outfile_makedb.replace("\\ ", " ")):
    #     os.remove(outfile_makedb.replace("\\ ", " "))
    
    # compute Clustering Threshold
    print("### Define distance threshold for clustering")
    cmd = f"./computeClusteringThreashold.R -f {outfile_parsedb} -o {os.path.abspath(sample_files['imgt'])}"
    if not os.path.exists(os.path.join(sample_files['imgt'].replace("\\ ", " "), "threshold.txt")):
        run_job(cmd, F"error when computeClusteringThreashold for {sample_name}")
    with open(os.path.join(sample_files['imgt'].replace("\\ ", " "), "threshold.txt"), 'r') as threshold_file:
        threshold = float(threshold_file.read())

    # Assign clones
    print("### Clustering clones")
    outfile_clones = os.path.join(sample_files['imgt'], F"{sample_name}_clones_id.tsv")
    cmd = F"DefineClones.py -d {outfile_parsedb} -o {outfile_clones} --act first --model {DISTANCE_MODEL} --norm len --dist {threshold}"
    if not os.path.exists(outfile_clones):
        run_job(cmd, F"error when compute distance threshold for {sample_name}")

    # # Set clonotype
    print("### Compute clonotypes")
    outfile_clonotypes = os.path.join(sample_files['imgt'], F"{sample_name}_clonotypes_changeo.tsv")

    # get fasta isotypes
    for file in os.listdir(sample_files['imgt']):
        # if file.endswith("_isotypes.fasta"):
        if file.endswith(".fasta"):
            fasta_isotypes = os.path.join(sample_files['imgt'], file)

    cmd = F"./changeoToClonotypes.py -i {outfile_clones} -o {outfile_clonotypes}"
    run_job(cmd, F"error when computing clonotype of sample : {sample_name}")
    if os.path.exists(outfile_parsedb.replace("\\ ", " ")):
        os.remove(outfile_parsedb.replace("\\ ", " "))


def main():
    global DISTANCE_MODEL
    parser = argparse.ArgumentParser(description='Clusterize sequences from IMGT output and define clonotypes', formatter_class=RawTextHelpFormatter)
    parser.add_argument("-c", "--config_file", type=str, help="Config file with one line by sample. Example :\n#sample_name\timgt\tfatsa\nsample1\tpath/to/imgt-oudir1\tpath/to/sequence1.fasta\nsample2\tpath/to/imgt-oudir2\tpath/to/sequence2.fasta", required=True)
    parser.add_argument("-m", "--model", type=str, help="Model to use for clusterize sequences (aa : aminoacid distance, ham : nulceotides distance). Default : aa", default="aa")
    args = parser.parse_args()

    DISTANCE_MODEL = args.model

    samples = get_samples_from_config(args.config_file)
    
    for sample_name, files in samples.items():
        print(F"##### {sample_name}")
        worker(sample_name, files)
    
main()