#!/usr/bin/env python
import os, sys, subprocess
import glob
import json
import string
import time

model = "T5tttt"
do_t2tt = True
mass_pts_job = 15

# dir used to collect gluino - LSP masses to run on 
tag = 'SMS-'+model+'_mGluino-XXX_mLSP-YYY_0_nom'
if (do_t2tt):
  tag = tag.replace('_nom','_t2tt')
example_dir = "/net/cms2/cms2r0/babymaker/babies/2019_07_16/"+model+"/skim_sys_abcd/*root"
release = 'CMSSW_8_1_0'

base_dir = "batch_"+model+"/" 
if do_t2tt:
  base_dir = "batch_"+model+"_t2tt/" 
if not os.path.exists(base_dir): 
  os.system("mkdir -p "+base_dir)

card_dir = os.path.join(base_dir,'cards')
if not os.path.exists(card_dir): 
  os.system("mkdir -p "+card_dir)

run_dir = os.path.join(base_dir,'run')
if not os.path.exists(run_dir): 
  os.system("mkdir -p "+run_dir)

wdir_tmpl = os.path.join(base_dir, 'scan_point_mGluino-XXX_mLSP-YYY')

# collect all mass pairs
allfiles = glob.glob(example_dir)
mass_pairs = []
for file in allfiles:
  tmp = file.split("mGluino-")[1]
  mglu = tmp.split("_mLSP-")[0]
  mlsp = tmp.split("_mLSP-")[1].split("_Tune")[0]
  mass_pairs.append([mglu, mlsp])
print 'Found a total of ',len(mass_pairs), 'mass pairs.'

os.system("JobSetup.csh")

ijob = 0
mpts_cards = []
mpts_limits = []
for mass in mass_pairs:
  # find which points do not have datacard and which do not have limit calculated
  impt_dir = wdir_tmpl.replace('XXX',mass[0]).replace('YYY',mass[1])
  if not os.path.exists(impt_dir): 
    os.system("mkdir -p "+impt_dir)
  card_path = os.path.join(card_dir,'datacard_'+tag.replace('XXX',mass[0]).replace('YYY',mass[1])+'.txt')
  if not os.path.exists(card_path):
    mpts_cards.append(mass)
  lim_path = os.path.join(impt_dir,'limits.txt')
  if not os.path.exists(lim_path):
    mpts_limits.append(mass)

  # submit job
  if (len(mpts_cards)>0 and len(mpts_cards)%mass_pts_job==0) or (len(mpts_limits)>0 and len(mpts_limits)%mass_pts_job==0) or mass==mass_pairs[-1]: 
    exename = os.path.join(run_dir,"find_limit_sig_"+str(ijob)+".sh")
    fexe = open(exename,"w")
    os.system("chmod u+x "+exename)
    fexe.write("#!/bin/bash\n\n")
    fexe.write(". /cvmfs/cms.cern.ch/cmsset_default.sh \n")
    fexe.write("cd ~/code/"+release+"/src/ \n")
    fexe.write("eval `scramv1 runtime -sh` \n")
    fexe.write("cd ~/code/ra4_draw/ ; \n\n")

    mpts_opt = ','.join([(impt[0]+'_'+impt[1]) for impt in mpts_cards])
    cmd = "./run/ra4/write_datacards.exe -u -y 0 -p "+mpts_opt+' -m '+model+ ' -o '+card_dir
    if (do_t2tt): 
      cmd += " -x t2tt"
    fexe.write(cmd+'\n')
    for impt in mpts_limits:
      impt_dir = wdir_tmpl.replace('XXX',impt[0]).replace('YYY',impt[1])
      card = 'datacard_'+tag.replace('XXX',impt[0]).replace('YYY',impt[1])+'.txt'
      fexe.write("./run/ra4/scan_point.exe -d "+card+' -i '+card_dir+' -o '+impt_dir+' >> '+impt_dir+'/limits.txt\n')
    fexe.close()
    cmd = "JobSubmit.csh ./run/wrapper.sh "+release+" ./"+exename
    os.system(cmd)
    ijob +=1
    mpts_cards = []
    mpts_limits = []

print 'Submitted', ijob, 'jobs.'
