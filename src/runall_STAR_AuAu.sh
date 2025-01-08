#!/bin/sh

#  Shell script for running all code
#  for STAR AuAu at 0.20 TeV
#  Created by Ryan Hamilton.

nevent="4900"

root -l -q "runGlauber_local.c($nevent, true, 0.00, 3.31)"
root -l -q "runGlauber_local.c($nevent, true, 3.31, 4.68)"
root -l -q "runGlauber_local.c($nevent, true, 4.68, 6.61)"
root -l -q "runGlauber_local.c($nevent, true, 6.61, 8.10)"
root -l -q "runGlauber_local.c($nevent, true, 8.10, 9.35)"
root -l -q "runGlauber_local.c($nevent, true, 9.35, 10.50)"
root -l -q "runGlauber_local.c($nevent, true, 10.50, 11.50)"
root -l -q "runGlauber_local.c($nevent, true, 11.50, 12.40)"
root -l -q "runGlauber_local.c($nevent, true, 12.40, 13.20)"
