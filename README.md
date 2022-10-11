<img src="https://static.igem.wiki/teams/4214/wiki/home/soft-logo.png" alt="B-LORE logo" width="150"/>

# ProMutor

**Pro**tein **M**utant genera**tor** 

A platform for generating advantageous point mutants based on explicit modelling of evolutionary history.


## Description
Let people know what your project can do specifically. Provide context and add a link to any reference visitors might
be unfamiliar with (for example your team wiki). A list of Features or a Background subsection can also be added here.
If there are alternatives to your project, this is a good place to list differentiating factors.

## Installation
Within a particular ecosystem, there may be a common way of installing things, such as using Yarn, NuGet, or Homebrew.
However, consider the possibility that whoever is reading your README is a novice and would like more guidance. Listing
specific steps helps remove ambiguity and gets people to using your project as quickly as possible. If it only runs in a
specific context like a particular programming language version or operating system or has dependencies that have to be
installed manually, also add a Requirements subsection.
```
git clone https://gitlab.igem.org/2022/software-tools/stockholm.git
cd stockholm
singularity pull docker://elodielaine/gemme:gemme
singularity pull docker://soedinglab/hh-suite:latest
pip3 install -r requirements.txt
```
### Database Installation
```
mkdir ./uniclust
cd uniclust/
wget  https://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/pfamA_35.0.tar.gz
tar xvfz pfamA_35.0.tar.gz
```
## Usage
Use examples liberally, and show the expected output if you can. It's helpful to have inline the smallest example of
usage that you can demonstrate, while providing links to more sophisticated examples if they are too long to reasonably
include in the README.
```
python3 run_container.py protein example/P27352.fasta 1,2,3 Blosum62 True pfama
```
## Contributing
State if you are open to contributions and what your requirements are for accepting them.

For people who want to make changes to your project, it's helpful to have some documentation on how to get started.
Perhaps there is a script that they should run or some environment variables that they need to set. Make these steps
explicit. These instructions could also be useful to your future self.

You can also document commands to lint the code or run tests. These steps help to ensure high code quality and reduce
the likelihood that the changes inadvertently break something. Having instructions for running tests is especially
helpful if it requires external setup, such as starting a Selenium server for testing in a browser.

## Authors and acknowledgment
Show your appreciation to those who have contributed to the project.
