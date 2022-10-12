<img src="https://static.igem.wiki/teams/4214/wiki/home/soft-logo.png" alt="ProMutor" width="800"/>

![Code Coverage](https://img.shields.io/badge/Python-3.10.6-green)



# ProMutor 
**Pro**tein **M**utant genera**tor** : a platform for generating advantageous point mutants based on explicit modelling of evolutionary history.


## Description
ProMutor is a web-based tool that predicts the effect of suggested mutations based on explicit modelling of the evolutionary history of natural sequences. Given an input sequence, it either generates a complete landscape of protein mutations or predicts an epistatic effect of mutations of interest. It also incorporates mutant sequence generation step as well as provides direct access to ColabFold notebook for predicting protein structure.

### Main steps of pipeline
<img src="https://static.igem.wiki/teams/4214/wiki/software/promutor-scheme.png" alt="Pipeline scheme" width="800">

Depending on the selected mode, steps of the pipeline vary. However, both start with the generation of a multiple sequence alignment (MSA).

**Screening mode:**

1. MSA generation.
2. Mutation effect landsccape generation (using GEMME tool).

**Mutant generation mode:**

1. MSA generation.
2. Generation of all possible mutants containing desired mutation sites.
3. Prediction of their effect on protein structure based on evolutionary conservation (GEMME).

As an additional suggested step to procceed with, we propose using AlphaFold 2 or mutant protein structure prediction. We provide direct access to its notebook version ColabFold. 

## Installation
There are several options of running promutor remotely either using [ProMutor](https://promutor.com/) webpage (does not require any installation) or via downloading `promutor.html` web-form to submit request to the remote compute from a local device.

### Installations required to run ProMutor locally 
ProMutor is written in python. In addition, it utilises `singularity` for running `hhblits`. 

To run ProMutor locally the following tools are needed:
- `python` v3.9 or higher
- `singularity`

To use ProMutor, you have to download the repository and pull required containers from the DockerHub as well as install required python packages:
```
git clone https://gitlab.igem.org/2022/software-tools/stockholm.git
cd stockholm
singularity pull docker://elodielaine/gemme:gemme
singularity pull docker://soedinglab/hh-suite:latest
pip3 install -r requirements.txt
```
### Database Installation

For constructing an MSA, `hhlits` tool utilises protein sequence databases. Here, we provide commands to install Pfam-A database, however other databases could be easily installed by providing an appropriate link.
```
mkdir ./uniclust
cd uniclust/
wget  https://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/pfamA_35.0.tar.gz
tar xvfz pfamA_35.0.tar.gz
```
**List of available databases**
- [Uniclust30](https://uniclust.mmseqs.com) [[pub]](https://doi.org/10.1093/nar/gkw1081)
- [BFD](https://bfd.mmseqs.com) [[pub]](https://doi.org/10.1038/s41592-019-0437-4)
- [Pfam/SCOP/PDB70/dbCAN](http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/)

In addition, one can build its own database. Please check `hh-suite` [user guide](https://github.com/soedinglab/hh-suite/wiki#building-customized-databases) or more information.

## Usage

### Quick start

#### Running remotely using web-form
1. Go to [ProMutor](https://promutor.com/) webpage
2. Fill in the form acccording to tutorials provided 
    a. [Screening mode](https://2022.igem.wiki/stockholm/software#tag3_1)
    b. [Mutant Generation mode](https://2022.igem.wiki/stockholm/software#tag3_2)

#### Running remotely using local submition form
1. Download `promutor.html` file
2. Open file in any browser
3. Fill in the form acccording to tutorials provided 
    a. [Screening mode](https://2022.igem.wiki/stockholm/software#tag3_1)
    b. [Mutant Generation mode](https://2022.igem.wiki/stockholm/software#tag3_2)

#### Running locally
1. Clone the repository
```
cd stockholm
```
2. Run tool locally in Mutant generation mode:
```
python3 run_container.py protein example/P27352.fasta 1,2,3 Blosum62 True pfama
```
3. Run tool locally in Screening mode:
```
python3 run_container.py protein example/P27352.fasta 0 Blosum62 True pfama
```
### Command line arguments

A command to run ProMutor can be used as follows:
```
python3 run_container.py <type_inp> <fasta> <pos> <matrix> <preserve> <database> --mand <string> --num <string>
```

Arguments | Description | Valid value
:---   | :---    | :--
type_inp&nbsp;*string*  | Type of input. | *nucleotide* or *protein*    
fasta&nbsp;*string*  | Input sequence. | -- 
pos&nbsp;*string*     | Comma separated list of positions to be mutated. | -- 
matrix&nbsp;*string* |  Substitution matrix to be used. |  *Blosum45*, *Blosum50*, *Blosum62*, *Blosum80*, *Blosum90*, *Pam30*, *Pam90* or *Pam250* 
preserve&nbsp;*string*        | Whether the closest (True) or the furthest (False) amino acid substituent is used. | *True* or *False* 
database&nbsp;*string*          | Database to be used for creating multiple sequence alignment (MSA). | *uniclust*, *pdb70*, *scop70*, *pfama* 
&#x2011;&#x2011;mand&nbsp;*string*    | Comma separated list of positions ALWAYS to be mutated. | -- 
&#x2011;&#x2011;num&nbsp;*string*    |Maximal number of mutations per sequence. | --

**Example:**

Example input `fasta` file is provided in `example/` directory. This simple command will run ProMutor in the Mutant generation mode.
```
python3 run_container.py protein example/P27352.fasta 1,2,3 Blosum62 True pfama
```

## Contributing
We are open to any contributions and colaborations. When contributing to this repository, please first discuss the change you wish to make via issue, email, or any other method with the owners of this repository before making a change. 

## License

[License](LICENSE.txt)

## Authors and acknowledgment

### Acknowledgment

We would like to express our greatest gratitude to our principal technical supervisor Frédéric Mauron who has been of the greatest help on the technical parts of the software creation.
His availability, knowledge and faith in our project allowed us to bring it to life. More specifically, Mr Mauron guided us throughout the webform-server connection steps. Last but not least, we are remarkably thankful for the server he provided us with for free.

### Authors

[Raphaël Mauron](https://www.linkedin.com/in/raphael-mauron/) & [Mariia Minaeva](https://www.linkedin.com/in/mariia-minaeva-019615193), iGEM Stockholm Team 2022
