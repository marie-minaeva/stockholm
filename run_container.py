import subprocess
from collections import defaultdict
from igem_tool import compute, parse_msa
from Bio import SeqIO, Seq
from os.path import exists
import os
import numpy as np
import pandas as pd
from datetime import datetime
import argparse


def parse_output(data):
    """
    Parsing of GEMME screening mode outputs
    Parameters
    ----------
    data : dict
        GEMME output file
        Ex.: {"a": "2.0 4.5 3.6 8.2 -1.0 NA",
                "b": "0.0 3.0 NA NA 1.0 5.0"}

    Returns
    -------
    x : list
        List of amino acid positions
    y : list
        List of amino acids
    z : object
        List of list where z[i, j] = score of having amino acid y[j] at the position x[i]
    """
    # Parsing the string and converting to float
    for key, value in data.items():
        value = value.replace("NA", "0.0")
        data[key] = [np.float32(i) for i in value.split(" ")]
    data = pd.DataFrame(data).T
    z = data.to_numpy().tolist()
    x = data.columns
    y = data.index
    return list(x), list(y), z


def run(type_inp, fasta_sequence, pos, matrix, preserve, database, mandatory_mutation=None, number_of_mutant=None):
    """
    Runs the whole pipeline
    Parameters
    ----------
    type_inp : str
        Type of input. Either "nucleotide" or "protein"
    fasta_sequence : str
        Input sequence
    pos : str
        Comma separated list of positions to be mutated
    matrix : str
        Substitution matrix to be used
    preserve : bool
        Whether the closest (True) or the furthest (False) amino acid substituent is used
    database : str
        Database to be used for creating multiple sequence alignment (MSA)
    mandatory_mutation : str
        Comma separated list of positions ALWAYS to be mutated
    number_of_mutant : str
        Maximal number of mutations per sequence

    Returns
    -------
    output : collections.defaultdict(<class 'list')
        Contains all desired outputs
        Keys:
            - 'nc' : defaultdict
            - 'prot' : defaultdict
            - 'mutant_predictions' : dict

        Ex.: {'nc': defaultdict(<class 'str'>, {}), 'prot': defaultdict(<class 'str'>, {'M1L': 'LAWFALYLLSLLWATASTQT'}),
        'mutant_predictions': {'V28I,I57V,M97L': -1.2911568830377, 'V28I,I57V': -0.867803185985955,
        'I57V,M97L': -0.623305422209993}}
    """
    # Creating WD
    now = str(datetime.now()).split(".")[1]
    wd = os.getcwd()
    changeWD = "mkdir " + now
    process = subprocess.Popen(changeWD.split(), stdout=subprocess.PIPE)
    _, _ = process.communicate()

    # Parsing inputs
    if database == "pdb70":
        database = "./uniclust/pdb70"
    elif database == "pfama":
        database = "./uniclust/pfam"
    elif database == "scop70":
        database = "./uniclust/scop70_1.75"
    else:
        database = "./uniclust/UniRef30_2022_02"
    if not exists(fasta_sequence):
        temp = fasta_sequence.replace(" ", "")
        temp = temp.split("\n")
        temp = "".join(temp[1:]).replace("\r", "")
        protein = Seq.Seq(temp)
        if type_inp == "nucleotide":
            protein = protein.translate()
            p = str(protein).find("*")
            protein = protein[:p]
        with open("./" + now + "/" + now + "temp_protein.fasta", "w") as output_handle:
            temp = SeqIO.SeqRecord(protein, id="WT", name="", description="")
            SeqIO.write(temp, output_handle, "fasta")
        fasta_sequence = "./" + now + "/" + now + "temp_protein.fasta"

    # Constructing and parsing MSA
    print("--> Running HHBlits")
    try:
        # Running HHblits
        runHHBlist = \
            "singularity exec /home/igem/iGEM_software/hh-suite_latest.sif hhblits -e 1e-10 -i 1 -p 40 -b 1 -B 20000 -i ./" \
            + now + "/" + now + "temp_protein.fasta -o ./" + now + "/try-hhlist.txt -oa3m ./" + now + \
            "/try-hhlist.a3m -d " + database + " -cpu 20"
        process = subprocess.Popen(runHHBlist.split(), stdout=subprocess.PIPE)
        _, _ = process.communicate()
        # Converting MSA to FASTA format
        parseMSA = "singularity exec /home/igem/iGEM_software/hh-suite_latest.sif  reformat.pl a3m fas ./" + now + \
                   "/try-hhlist.a3m ./" + now + "/try-hhlist.fas"
        process = subprocess.Popen(parseMSA.split(), stdout=subprocess.PIPE)
        _, _ = process.communicate()
        # Ungapping MSA
        parse_msa(now)
    except FileNotFoundError:
        print("No HHBlist")

    # Running Mutant Generation
    print("--> Running Mutant Generator tool")

    # Prediction mode
    if pos != "0":
        output = compute(type_inp, fasta_sequence, pos, matrix, preserve, now, number_of_mutant=number_of_mutant,
                         mandatory_mutation=mandatory_mutation)
    # Screening mode
    else:
        output = defaultdict(list)
        output["nc"] = []
        output["prot"] = []

    # Running evolutionary score calculations
    print("--> Running GEMME")
    runGEMMEscreen = "singularity exec /home/igem/iGEM_software/gemme_gemme.sif python2.7 /opt/GEMME/gemme.py " \
                     + now + ".fasta -r input -f " + now + ".fasta"
    # Changing wd
    path = "/home/igem/iGEM_software/" + now + "/"
    os.chdir(path)
    try:
        # Prediction mode
        if pos != "0":
            # Running calculations
            runGEMMEmut = "singularity exec /home/igem/iGEM_software/gemme_gemme.sif python2.7 /opt/GEMME/gemme.py " + \
                          now + ".fasta -r input -f " + now + ".fasta -m " + "./" + now + "for_gemme.txt"
            process = subprocess.Popen(runGEMMEmut.split(), stdout=subprocess.PIPE)
            _, _ = process.communicate()
            # Parsing outputs
            with open("WT_normPred_evolCombi.txt", "r") as f:
                lines = f.readlines()
                lines = [i.split(" ") for i in lines[1:]]
            output["mutant_predictions"] = \
                dict(sorted({i[0][1:-1]: float(i[1][:-1]) for i in lines}.items(), key=lambda item: item[1]))
        # Screening mode
        else:
            # Running calculations
            process = subprocess.Popen(runGEMMEscreen.split(), stdout=subprocess.PIPE)
            _, _ = process.communicate()
            # Parsing outputs
            with open("WT_normPred_evolCombi.txt", "r") as f:
                lines = f.readlines()
                lines = [i.split(" ") for i in lines[1:]]
            x, y, z = parse_output({i[0][1:-1]: " ".join(i[1:-1]) for i in lines})
            output["mutant_predictions"] = {"x": x, "y": y, "z": z}
    # If too small MSA enable screening mode
    except FileNotFoundError:
        # Running calculations
        process = subprocess.Popen(runGEMMEscreen.split(), stdout=subprocess.PIPE)
        _, _ = process.communicate()
        # Parsing outputs
        with open("WT_normPred_evolCombi.txt", "r") as f:
            lines = f.readlines()
            lines = [i.split(' ') for i in lines[1:]]
        x, y, z = parse_output({i[0][1:-1]: " ".join(i[1:-1]) for i in lines})
        output["mutant_predictions"] = {"x": x, "y": y, "z": z}

    print("--> We are done")
    # Deleting wd
    os.chdir(wd)
    cleanTemps = "rm -r " + now
    process = subprocess.Popen(cleanTemps, stdout=subprocess.PIPE, shell=True)
    _, _ = process.communicate()
    return output

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process input arguments.')
    parser.add_argument('type_inp', type=str,
                        help='Type of input. Either "nucleotide" or "protein"')
    parser.add_argument('fasta', type=str,
                        help='Input sequence')
    parser.add_argument('pos', type=str,
                        help='Comma separated list of positions to be mutated')
    parser.add_argument('matrix', type=str, default="Blosum62",
                        help='Substitution matrix to be used. \n Options: "Blosum62", "Blosum80", "Blosum45", \n\t '
                             '"Blosum50", "Blosum90", "Pam30", "Pam90", "Pam250"')
    parser.add_argument('preserve', type=str, default="True",
                        help='Whether the closest (True) or the furthest (False) amino acid substituent is used')
    parser.add_argument('db', metavar="database", type=str, default="uniclust",
                        help='Database to be used for creating multiple sequence alignment (MSA)')
    parser.add_argument('--mand', metavar="mandatory_mutation", type=str, nargs="?",
                        help='Comma separated list of positions ALWAYS to be mutated')
    parser.add_argument('--num', metavar="number_of_mutant", type=str, nargs="?",
                        help='Maximal number of mutations per sequence')
    args = parser.parse_args()
    print(args)
    print(run(type_inp=args.type_inp, fasta_sequence=args.fasta, pos=args.pos, matrix=args.matrix, preserve=args.preserve, database=args.db, mandatory_mutation=args.mand, number_of_mutant=args.num))