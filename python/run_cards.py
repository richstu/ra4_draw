#! /usr/bin/env python
import os, sys
import glob
import argparse

def parseArguments():
  parser = argparse.ArgumentParser(description="Runs all datacards in a folder that have .txt extension", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument("input_directory", default="cards", help="Input directory which holds card files.")
  parser.add_argument("-o", "--outputFile", default="card.limits", metavar="LIMIT_FILE", help="Output limit file to be ran by plot_limit")
  parser.add_argument('-b', '--runBatch', action='store_true', help='Run in batch')
  parser.add_argument('-w', '--cmsswWrapperPath', default='./run/cmsswWrapperAndLogging.sh', metavar="CMSSW_WRAPPER_FILE", help='CMSSW wrapper file for batch')
  args = parser.parse_args()
  # Check arguments
  if not os.path.isdir(args.input_directory):
    print ('[Error] There is no directory called '+args.input_directory)
    sys.exit()
  if len(glob.glob(args.input_directory+'/*.txt')) == 0:
    print ('[Error] There is no txt files in '+args.input_directory)
    sys.exit()
  if os.path.exists(args.outputFile):
    print ('[Error] There is already a file called '+args.outputFile)
    sys.exit()
  if args.runBatch and not os.path.isfile(args.cmsswWrapperPath):
    print ('[Error] Cmssw wrapper file does not exist => '+args.cmsswWrapperPath)
    sys.exit()
  return args

#def makeCmsswWrapperScript(cmsswWrapperPath, args):
#  cmssw_base = os.environ['CMSSW_BASE']
#  ra4drawDirectory = os.environ['RA4_DRAW_BASE']
#  cmsswWrapper = open(cmsswWrapperPath,'w')
#  cmsswWrapper.write('#!/bin/bash\n')
#  cmsswWrapper.write(". /cvmfs/cms.cern.ch/cmsset_default.sh \n")
#  cmsswWrapper.write('cd '+cmssw_base+'/src\n')
#  cmsswWrapper.write("eval `scramv1 runtime -sh` \n")
#  cmsswWrapper.write("echo Setup "+cmssw_base+'\n')
#  cmsswWrapper.write('cd '+ra4drawDirectory+'\n')
#  #cmsswWrapper.write('eval $1\n')
#  cmsswWrapper.write('./run/hig/scan_point.exe --cards -f $1 >> '+args.outputFile+'\n')
#  cmsswWrapper.write('echo Complete\n')
#  os.system("chmod u+x "+cmsswWrapperPath)

def makeCommands(dataCardPaths):
  commands = []
  for datacardPath in dataCardPaths:
    cmd = "./run/hig/scan_point.exe --cards -f "+datacardPath + " >> " + args.outputFile
    #os.system(cmd)
    commands.append(cmd)
  return commands

def runCommands(commands):
  for command in commands:
    print ("Running " + command)
    os.system(command)

def runCommandsInBatch(cmsswWrapperPath, command):
  for command in commands:
    if '>>' in command:
      runCommand, outFile = command.split('>>')
    bsubCommand = 'JobSubmit.csh -nomail '+cmsswWrapperPath+' ' + os.environ['SRT_CMSSW_VERSION_SCRAMRTDEL'] + ' ' + runCommand + ' '+outFile
    #bsubCommand = 'JobSubmit.py '+cmsswWrapperPath+' '+cardFile.rstrip()
    print (bsubCommand)
    os.system(bsubCommand)

if __name__ == "__main__":
  args = parseArguments()
  dataCardPaths = glob.glob(args.input_directory+'/*.txt')

  exitCode = os.system("./compile.py")
  if (exitCode): sys.exit()

  commands = makeCommands(dataCardPaths)
  #print(commands)
  #sys.exit()
  
  if (not args.runBatch): runCommands(commands)
  else:
    runCommandsInBatch(args.cmsswWrapperPath, commands)

  print('Run ./run/hig/plot_limit.exe -f '+args.outputFile)
