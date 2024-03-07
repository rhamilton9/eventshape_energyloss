#!/bin/sh

#  Shell script for running all code
#  Created by Ryan Hamilton.

root -l -q 'runGlauber_local.c(5000, true)'
root -l -q compile_glauber.c > ../logs/glauber.log
root -l -q eventshape_summary.c > ../logs/glauber.log
