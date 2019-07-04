#! /usr/bin/env python
import os, sys

masses = ['127']
masses.extend([str(i) for i in range(150,1001,25)])

os.system("./compile.py")
for mass in masses:
    cardname = "datacard_SMS-TChiHH_mGluino-"+mass+"_mLSP-1_bfH-100_35p9ifb.txt"
    print "Processing", mass
    cmd = "./run/hig/scan_point.exe --cards -f cards/"+cardname + " >> limits_cards.txt"
    os.system(cmd)
