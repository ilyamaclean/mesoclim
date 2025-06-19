#!/bin/bash
#SBATCH --partition=short-serial
#SBATCH --job-name=run_mesoclim
#SBATCH -e /gws/nopw/j04/uknetzero/mesoclim/script_err/%j.err
#SBATCH -o /gws/nopw/j04/uknetzero/mesoclim/script_out/%j.out
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=25000
module add jasr
Rscript /gws/nopw/j04/uknetzero/mesoclim/mesoclim_scripts/jrun_newmesoclim_addtrees_10_1.R /gws/nopw/j04/uknetzero/mesoclim/mesoclim_inputs/parcels/focal_parcels.shp /gws/nopw/j04/uknetzero/mesoclim/mesoclim_outputs/focal_area 01 2011 2020
