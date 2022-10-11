import subprocess
from collections import defaultdict
from igem_tool import compute, parse_msa
from Bio import SeqIO, Seq
from os.path import exists
import os
import pandas as pd
from datetime import datetime
import numpy as np

def parse_output(data):
    for key, value in data.items():
        value = value.replace("NA", "0.0")
        data[key] = [np.float32(i) for i in value.split(" ")]
    data = pd.DataFrame(data).T
    z = data.to_numpy().tolist()
    x = data.columns
    y = data.index
    #sns.heatmap(data)
    return list(x), list(y), z

## run ssh from python and then run the tool from there
def run(type_inp, fasta_sequence, pos, matrix, preserve, database, number_of_mutant=None, mandatory_mutation=None):
    now = str(datetime.now()).split(".")[1]
    print(database)
    if database == "pdb70":
        database = "./uniclust/pdb70"
    elif database == "pfama":
        database = "./uniclust/pfam"
    elif database == "scop70":
        database = "./uniclust/scop70_1.75"
    else:
        database = "./uniclust/UniRef30_2022_02"
    print(database)
    bashCommand = "mkdir " + now
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    out, error = process.communicate()
    if exists(fasta_sequence):
        gene = fasta_sequence.split(".")[0]
    else:
        temp = fasta_sequence.replace(" ", "")
        temp = temp.split("\n")
        temp = "".join(temp[1:]).replace("\r", "")
        protein = Seq.Seq(temp)
        if type_inp == "nucleotide":
            protein = protein.translate()
            p = str(protein).find("*")
            protein = protein[:p]
        gene = "temp"
        with open("./" + now + "/" + now + "temp_protein.fasta", "w") as output_handle :
           temp = SeqIO.SeqRecord(protein, id="WT", name="", description="")
           SeqIO.write(temp, output_handle, "fasta")
        fasta_sequence = "./" + now + "/" + now + "temp_protein.fasta"
    print(protein)
    print("--> Running PSIBlast")
    try:
        bashCommand = "singularity exec /home/igem/iGEM_software/hh-suite_latest.sif hhblits -e 1e-10 -i 1 -p 40 -b 1 -B 20000 -i ./" + now + "/" + now + "temp_protein.fasta -o ./" + now + "/try-hhlist.txt -oa3m ./" + now + "/try-hhlist.a3m -d " + database + " -cpu 20"
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        out, error = process.communicate()
        bashCommand = "singularity exec /home/igem/iGEM_software/hh-suite_latest.sif  reformat.pl a3m fas ./" + now + "/try-hhlist.a3m ./" + now + "/try-hhlist.fas"
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        out, error = process.communicate()
        parse_msa(now)
    except FileNotFoundError:
        print("No psiblast")
    print("--> Running our tool")
    print(pos)
    if pos != "0":
        output = compute(type_inp, fasta_sequence, pos, matrix, preserve, now, number_of_mutant=number_of_mutant, mandatory_mutation=mandatory_mutation)
    else:
        output = defaultdict(list)
        output["nc"] = {}
        output["prot"] = {}
    print("--> Running gemme")
    try:
        path="/home/igem/iGEM_software/" + now + "/"
        os.chdir(path)
        if pos != "0":
            bashCommand = "singularity exec /home/igem/iGEM_software/gemme_gemme.sif python2.7 /opt/GEMME/gemme.py " + now + ".fasta -r input -f " + now + ".fasta -m " + "./" + now + "for_gemme.txt"
            process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
            out, error = process.communicate()
            with open("WT_normPred_evolCombi.txt", "r") as f:
                lines = f.readlines()
                lines = [i.split(" ") for i in lines[1:]]
            output["mutant_predictions"] = dict(sorted({i[0][1 :-1] : float(i[1][:-1]) for i in lines}.items(), key=lambda item: item[1]))
        else:
            bashCommand = "singularity exec /home/igem/iGEM_software/gemme_gemme.sif python2.7 /opt/GEMME/gemme.py " + now + ".fasta -r input -f " + now + ".fasta"
            process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
            out, error = process.communicate()
            with open("WT_normPred_evolCombi.txt", "r") as f:
                lines = f.readlines()
                lines = [i.split(" ") for i in lines[1 :]]
            x, y, z = parse_output({i[0][1 :-1] : " ".join(i[1 :-1]) for i in lines})
            output["mutant_predictions"] = {"x": x, "y": y, "z": z}
    except FileNotFoundError:
        if pos == "0":
            output = defaultdict(list)
        print("No gemme")
        bashCommand = "singularity exec /home/igem/iGEM_software/gemme_gemme.sif python2.7 /opt/GEMME/gemme.py " + now + ".fasta -r input -f " + now + ".fasta"
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        out, error = process.communicate()
        with open("WT_normPred_evolCombi.txt", "r") as f :
            lines = f.readlines()
            lines = [i.split(' ') for i in lines[1:]]
        x, y, z = parse_output({i[0][1 :-1] : " ".join(i[1 :-1]) for i in lines})
        output["mutant_predictions"] = {"x" : x, "y" : y, "z" : z}
    print("--> We are done")
    bashCommand = "rm -r " + now
    process = subprocess.Popen(bashCommand, stdout=subprocess.PIPE, shell=True)
    out, error = process.communicate()
    return output
