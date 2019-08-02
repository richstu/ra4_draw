#! /usr/bin/env python

import os, glob

card_dir = "reco_gen_avg"
outdir = "batch_t2tt"

cards = [ i.split("/")[1] for i in glob.glob(card_dir+"/datacard*txt")]


for card in cards:
    tmp_dir = outdir+"/"+card.strip("datacard_").strip(".txt")
    os.mkdir(tmp_dir)
    cmd = "./run/ra4/scan_point.exe -d "+card+' -i '+card_dir+' -o '+tmp_dir+' >> '+tmp_dir+'/limits.txt\n'   
    os.system(cmd)

print "Done"
