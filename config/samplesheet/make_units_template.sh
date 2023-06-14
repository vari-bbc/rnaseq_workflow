#!/bin/bash

set -e
set -u
set -o pipefail

module load bbc2/R/alt/R-4.2.1-setR_LIBS_USER
Rscript make_units_template.R 
