#!/usr/bin/python

""" A simple script for generating WHAM input files """

import os, sys
import re
import numpy as np

#----------------------#
# Function Definitions #
#----------------------#

def extractColVar(cv1, wham_inp_file):

  k1 = 4 #kcal/mol
  start_jobstp = 0
  end_jobstp   = 1

  base_name = 'lp_3hy5_refmac1_R234W-ions-wats'
  for jobstep in range(start_jobstp, end_jobstp):
    if (jobstep == start_jobstp):
      os.system('awk \' FNR > 1000 { if ($1 != "#") print $2; }\' /home/daniel/Downloads/trajectory_data/CRALBP/sampling/R234_lp/%s-%s.colvars.traj > ./timeseries/cv-%s-ts' % (base_name, cv1, cv1))
    else:
      os.system('awk \' { if ($1 != "#") print $2; }\' /home/daniel/Downloads/trajectory_data/CRALBP/sampling/R234_lp/%s-%s.colvars.traj >> ./timeseries/cv-%s-ts' % (base_name, cv1, cv1))


  wham_inp_file.write('./timeseries/cv-%s-ts\n' % cv1)
  wham_inp_file.write('%s   %03.1f\n' % (cv1, k1))

#--------------#
# Main Routine #
#--------------#

wham_inp = open('wham.inp', 'w')
wham_inp.write("""R234W_pmf.dat
R234W_rho.dat
R234W_bia.dat
R234W_fff.dat
-11.8 11.1  0.2  0.0
21 100000 100
0.001
""")

for cv1 in range(-10, 11):
    extractColVar(cv1, wham_inp)
