#!/bin/bash
cd $HOME

git clone https://github.com/asimoesgomes/bsa_sim.git

cd bsa_sim

module load R/4.2.0
module load gsl/2.7

R -e "install.packages('renv')"
R -e "renv::init()"
