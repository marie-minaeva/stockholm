<img src="https://static.igem.wiki/teams/4214/wiki/home/soft-logo.png" alt="B-LORE logo" width="150"/>

![Code Coverage](https://img.shields.io/badge/Python-3.10.6-green)

# ProMutor 
**Pro**tein **M**utant genera**tor** : a platform for generating advantageous point mutants based on explicit modelling of evolutionary history.


## Description
ProMutor is a we-based tool that predicts the effect of suggested mutations based on explicit modelling of the evolutionary history of natural sequences. Given an input sequence, it either generates a complete landscape of protein mutations or predicts an epistatic effect of mutations of interest. It also incorporates mutant sequence generation step as well as provedes direct access to ColabFold notebook for predicting protein structure.

### Main steps of pipeline
Depending on the selected mode, steps o the pipeline vary. However, both start with the generation o a multiple sequence alignment (MSA).

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

To run ProMutor locally following tools are needed:
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

For constructing an MSA, `hhlits` tool utilises protein sequence databases. Here we provide commands to install Pfam-A database, however other databases could be easily installed by providing an appropriate link.
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
2. Fill in te form acccording to tutorials provided 
    a. [Screening mode](https://2022.igem.wiki/stockholm/software#tag3_1)
    b. [Mutant Generation mode](https://2022.igem.wiki/stockholm/software#tag3_2)

#### Running remotely using local submition form
1. Download `promutor.html` file
2. Open file in any browser
3. Fill in te form acccording to tutorials provided 
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

An command to run ProMutor can be used as follows:
ython3 cif2fasta.py -i <all_cifs> -o pdb100.fas -c <num_cores> -p pdb_filter.dat
```
python3 run_container.py <> example/P27352.fasta 0 Blosum62 True pfama
```
Option | Description | Priority | Default value
:---   | :---        |:---      | :--
&#x2011;&#x2011;gen&nbsp;*filename(s)*  | Input genotype file(s), all loci should have separate genotype files and specified here (wildcards allowed) | Required    | --
&#x2011;&#x2011;sample&nbsp;*filename*  | Input sample file | Required | --
&#x2011;&#x2011;pheno&nbsp;*string*     | Name of the phenotype as it appears in the header of the sample file| Optional | `pheno`
&#x2011;&#x2011;regoptiom               | If specified, the variance of the regularizer will be optimized, otherwise it will be N(0, σ<sup>2</sup>) where σ is specified by `--reg` | Optional | --
&#x2011;&#x2011;reg&nbsp;*float*        | Value of the standard deviation (σ) of the regularizer | Optional | 0.01
&#x2011;&#x2011;pca&nbsp;*int*          | Number of principal components of the genotype to be included as covariates | Optional | 0
&#x2011;&#x2011;cov&nbsp;*string(s)*    | Name of covariate(s) as they appears in the header of the sample file, multiple covariates can be specified as space-separated strings | Optional | None
&#x2011;&#x2011;out&nbsp;*directory*    | Name of the output directory where summary statistics will be created | Optional | directory of the genotype files
&#x2011;&#x2011;prefix&nbsp;*string* | Prefix for the summary statistics files | Optional | `_summary`

## Contributing
We are open to contributions from anyone! Please let us know via an issue if you find a problem with our design or would like to request a feature.
If you find some improvements or changes that work for you, consider opening a pull request with them!

State if you are open to contributions and what your requirements are for accepting them.

For people who want to make changes to your project, it's helpful to have some documentation on how to get started.
Perhaps there is a script that they should run or some environment variables that they need to set. Make these steps
explicit. These instructions could also be useful to your future self.

You can also document commands to lint the code or run tests. These steps help to ensure high code quality and reduce
the likelihood that the changes inadvertently break something. Having instructions for running tests is especially
helpful if it requires external setup, such as starting a Selenium server for testing in a browser.

## Authors and acknowledgment
Show your appreciation to those who have contributed to the project.
