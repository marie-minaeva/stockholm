from Bio import SeqIO, Seq
from Bio.SubsMat import MatrixInfo
import numpy as np
from collections import defaultdict
from itertools import chain, combinations
from os.path import exists

# Codons
my_dict = {
    'F': ['TTT', 'TTC'], 'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
    'I': ['ATT', 'ATC', 'ATA'], 'M': ['ATG'], 'V': ['GTT', 'GTC', 'GTA', 'GTG'],
    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
    'P': ['CCT', 'CCC', 'CCA', 'CCG'], 'T': ['ACT', 'ACC', 'ACA', 'ACG'],
    'A': ['GCT', 'GCC', 'GCA', 'GCG'], 'Y': ['TAT', 'TAC'], 'H': ['CAT', 'CAC'], 'Q': ['CAA', 'CAG'],
    'N': ['AAT', 'AAC'], 'K': ['AAA', 'AAG'], 'D': ['GAT', 'GAC'],
    'E': ['GAA', 'GAG'], 'C': ['TGT', 'TGC'], 'W': ['TGG'], 'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'G': ['GGT', 'GGC', 'GGA', 'GGG'],
    '*': ["TAA", "TAG", "TGA"], "B": ['AAT', 'AAC', 'GAT', 'GAC'], 'Z': ['CAA', 'CAG', 'GAA', 'GAG']
}
# Substitution matrices
matrices = {"Blosum62": MatrixInfo.blosum62, "Blosum80": MatrixInfo.blosum80, "Blosum45": MatrixInfo.blosum45,
            "Blosum50": MatrixInfo.blosum50, "Blosum90": MatrixInfo.blosum90,
            "Pam30": MatrixInfo.pam30, "Pam90": MatrixInfo.pam90, "Pam250": MatrixInfo.pam250}


def get_matrix(matrix):
    """
    Converts substitution matrix to a record
    Parameters
    ----------
    matrix : str
        Selected option of substitution matrix

    Returns
    -------
    blosum_mat : collections.defaultdict
        Substitution matrix in the records format
    """
    blosum = matrices[matrix]
    aas = list(set([x[0] for x in blosum.keys()]))
    blosum_mat = defaultdict()
    for a in aas:
        blosum_mat[a] = defaultdict(float)
        for b in aas:
            try:
                blosum_mat[a][b] = blosum[(a, b)]
            except KeyError:
                blosum_mat[a][b] = blosum[(b, a)]
    return blosum_mat


def parse_msa(wd):
    """
    Parses MSA hhblist outputs (initial .a3m was converted to .fas (FASTA) format beforehand) to ungap the MSA with
    respect to the query  sequence
    Parameters
    ----------
    wd : str
        Path to the working directory

    Returns
    ----------
    Resulting MSA is stored in a file ./wd/wd.fasta

    Files
    ----------
    ./wd/try-hhlist.fas :
        Input MSA
    ./wd/wd.fasta :
        Output MSA
    """
    seqs = defaultdict(str)
    query = False
    # Reading the MSA file
    for seq_record in SeqIO.parse(wd + "/" + "try-hhlist.fas", "fasta"):
        seqs[seq_record.name] = seq_record.seq
        if not query:
            query = seq_record.name
    # Getting gaps positions in the query sequence to be deleted
    pos = [i for i in range(len(seqs[query])) if seqs[query].startswith('-', i)]
    pos = pos[::-1]
    # Deleting positions
    for key, value in seqs.items():
        temp = value
        for p in pos:
            temp = value[:p] + temp[p + 1:]
        seqs[key] = temp
    # Writing resulting MSA to file
    with open(wd + "/" + wd + ".fasta", "w") as output_handle:
        for i in seqs.items():
            writer = SeqIO.FastaIO.FastaWriter(output_handle, wrap=50)
            writer.write_record(SeqIO.SeqRecord(Seq.Seq(i[1]), id=i[0], name=i[0], description=""))


def compute(type_inp, fasta_sequence, pos, matrix, preserve, wd, mandatory_mutation=None, number_of_mutant=None):
    """
    Generates mutant sequences
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
    wd : str
        Path to the working directory
    mandatory_mutation : str
        Comma separated list of positions ALWAYS to be mutated
    number_of_mutant : str
        Maximal number of mutations per sequence

    Returns
    -------
    output : collections.defaultdict(<class 'list')
        Dictionary that contains a dictionary of protein mutants and a dictionary of nucleotide mutants where applicable
        Ex.: {'nc': defaultdict(<class 'str'>, {}), 'prot': defaultdict(<class 'str'>, {'M1L': 'LAWFALYLLSLLWATASTQT'})}
    """
    # Arguements parsing
    # Positions to mutate
    pos = [int(x) - 1 for x in pos.split(",")]
    # Type of input
    if type_inp == "protein":
        is_nc = False
    else:
        is_nc = True

    # Reading input fasta
    if exists(fasta_sequence):
        for seq_record in SeqIO.parse(fasta_sequence, "fasta"):
            protein = seq_record.seq
    else:
        temp = fasta_sequence.split("\n")
        temp = "".join(temp[1:]).replace("\r", "")
        protein = Seq.Seq(temp)
        with open("./" + wd + "/" + wd + "temp_protein.fasta", "w") as output_handle:
            temp = SeqIO.SeqRecord(protein, id="WT", name="", description="")
            SeqIO.write(temp, output_handle, "fasta")
    seq = protein
    # Preproccessing of a substitution matrix
    blosum_mat = get_matrix(matrix)
    # Creating combinations of positions to mutate
    if number_of_mutant:
        to_mutate = chain.from_iterable(combinations(pos, r) for r in range(int(number_of_mutant) + 1))
    else:
        to_mutate = chain.from_iterable(combinations(pos, r) for r in range(len(pos) + 1))
    if mandatory_mutation:
        mand = [int(i) - 1 for i in mandatory_mutation.split(",")]
        to_mutate = [i for i in to_mutate if len(set(i).intersection(set(mand))) == len(mand)]

    # Initialization
    mut_prots = defaultdict(str)
    mut_nuc = defaultdict(str)

    # Mutant sequence creation
    for mut in list(to_mutate)[1:]:
        new_seq = protein
        name = ""
        name_nuc = ""
        new_seq_nuc = seq
        for p in mut:
            aa = protein[p]
            # Selecting a proper substituent
            if preserve == "True":
                sub = list(blosum_mat[aa].keys())[(-np.array(list(blosum_mat[aa].values()))).argsort()[1]]
            else:
                sub = list(blosum_mat[aa].keys())[(np.array(list(blosum_mat[aa].values()))).argsort()[1]]

            # Creating an amino acid mutant sequence
            new_seq = list(new_seq)
            new_seq[p] = sub
            new_seq = "".join(new_seq)
            name += aa + str(p + 1) + sub + ","

            # Creating a nucleotide mutant sequence
            flag = 3
            if is_nc:
                for s in my_dict[sub]:
                    prev = str(new_seq_nuc[p * 3:p * 3 + 3])
                    new_seq_nuc = list(new_seq_nuc)
                    ham = not (prev == s)
                    if ham <= flag:
                        new_seq_nuc[p * 3:p * 3 + 3] = s
                        new_seq_nuc = "".join(new_seq_nuc)
                        flag = ham
                        temp = prev + str(p) + s
                name_nuc += temp
        mut_prots[name[:-1]] = new_seq
        if is_nc:
            mut_nuc[name_nuc] = new_seq_nuc

    # Combining output dicctt
    output = defaultdict(defaultdict)
    output["nc"] = mut_nuc
    output["prot"] = mut_prots

    # Store mutations to pass to GEMME
    with open("./" + wd + "/" + wd + "for_gemme.txt", "w") as f:
        f.write("\n".join(output["prot"].keys()))
    return output