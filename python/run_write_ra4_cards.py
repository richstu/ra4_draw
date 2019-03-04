#! /usr/bin/env python

import os, sys

lumi = 135
tag = "mjbins"

infolder  = "/net/cms29/cms29r0/babymaker/babies/2017_02_22_grooming/T1tttt/renormed/"
outfolder = "/net/cms29/cms29r0/babymaker/datacards/2018_09_06/T1tttt/"+tag+"_"+'{:.1f}'.format(lumi).\
replace('.','p')+'/'
runfolder = outfolder+"run/"
if not os.path.exists(runfolder):
  os.system("mkdir -p "+runfolder)

#input datasets                                                                                        
inputfiles = [i for i in os.listdir(infolder) if ".root" in i]

os.system("JobSetup.csh")
nfiles = 10
njobs = len(inputfiles)/nfiles
if (len(inputfiles)> njobs*nfiles): njobs +=1

for ibatch in range(njobs+1):
  exename = runfolder+"/run_cards_"+str(ibatch)+".sh"
  fexe = open(exename,"w")
  os.system("chmod u+x "+exename)
  fexe.write("#!/bin/bash\n\n")
  fexe.write("./run/ra4/write_datacards.exe -b "+str(ibatch)+' -n '+str(nfiles)+' -o '+outfolder+' -t \
'+tag+' -l '+str(lumi)+'\n')
  fexe.close()
  cmd = "JobSubmit.csh ./run/core/wrapper.sh CMSSW_8_0_25 "+exename
  if (ibatch==0): print cmd
  os.system(cmd)

print "\nSubmitted "+str(len(inputfiles))+" files in "+str(njobs)+" jobs. Output goes to "+outfolder+"\
\n"
sys.exit(0)