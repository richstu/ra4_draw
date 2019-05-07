#!/usr/bin/env python

###### Script to send 58 jobs to the batch finding the limit for each signal mass point
import os, sys, subprocess
import pprint, operator
import glob
import json
import string
import time

# Setting folders
model = "T1tttt"
# dir used to collect gluino - LSP masses to run on 
example_dir = "/net/cms2/cms2r0/babymaker/babies/2018_12_17/T1tttt/skim_sys_abcd/*root"
release = 'CMSSW_8_1_0'
card = 'cards/datacard_SMS-T1tttt_mGluino-XXX_mLSP-YYY_0_nom.txt'

only_missed = False

runfolder = "batch_"+model+"/" 
if not os.path.exists(runfolder):
  os.system("mkdir -p "+runfolder)

allfiles = glob.glob(example_dir)

mass_pairs = []
if only_missed:
  for idir in glob.glob("scan_point/scan_point*"):
    tmp_ = idir.split("mGluino-")[1]
    mglu_ = tmp_.split("_mLSP-")[0]
    mlsp_ = tmp_.split("_mLSP-")[1].split("_")[0]
    mlsp_ = mlsp_[0:-6]
    if (len(glob.glob("scan_point/scan_point*mGluino-"+mglu_+"_mLSP-"+mlsp_+"*/limit*.txt"))==0):
      mass_pairs.append([mglu_, mlsp_])

  print "Found ", len(mass_pairs),"missing mass pairs:"
  mass_pairs.sort(key=operator.itemgetter(0,1))
  print mass_pairs
else:
  for file in allfiles:
    tmp = file.split("mGluino-")[1]
    mglu = tmp.split("_mLSP-")[0]
    mlsp = tmp.split("_mLSP-")[1].split("_Tune")[0]
    mass_pairs.append([mglu, mlsp])
  print 'Found',len(mass_pairs), 'mass pairs.'

os.system("JobSetup.csh")

len(mass_pairs)/njobs
if files_job*njobs<len(mass_pairs): files_job +=1
print 'Submitting ',njobs, 'jobs with ',files_job,' mass points per job.'

for ijob in range(njobs):
  exename = runfolder+"/find_limit_sig_"+str(ijob)+".sh"
  fexe = open(exename,"w")
  os.system("chmod u+x "+exename)
  fexe.write("#!/bin/bash\n\n")
  fexe.write(". /cvmfs/cms.cern.ch/cmsset_default.sh \n")
  fexe.write("cd ~/code/"+release+"/src/ \n")
  fexe.write("eval `scramv1 runtime -sh` \n")
  fexe.write("cd ~/code/ra4_draw/ ; \n\n")

  last = (ijob+1)*files_job
  if (last>len(mass_pairs)): last = len(mass_pairs)
  mpts = ','.join([(mass_pairs[i][0]+'_'+mass_pairs[i][1]) for i in range(ijob*files_job, last)])
  if (not only_missed):
    fexe.write("./run/ra4/write_datacards.exe -u -y 0 -p "+mpts+' \n')
  for i in range(ijob*files_job, last):
    fexe.write("./run/ra4/scan_point.exe -f "+card.replace('XXX',mass_pairs[i][0]).replace('YYY',mass_pairs[i][1])+' >> txt/limits_'+model+'_'+str(ijob)+'.txt\n')

  fexe.close()
  cmd = "JobSubmit.csh ./run/wrapper.sh "+release+" ./"+exename
  #print cmd
  os.system(cmd)

print "Done."
