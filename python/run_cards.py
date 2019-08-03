#! /usr/bin/env python
import os, sys
import glob
import argparse

def parseArguments():
  parser = argparse.ArgumentParser(description="Runs all datacards in a folder that have .txt extension", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument("-i", "--inputDirectory", default="cards", metavar="INPUT_DIR", help="Directory which holds card files.")
  parser.add_argument("-o", "--outputFile", default="cards.limits", metavar="LIMIT_FILE", help="Output limit file to be ran by plot_limit")
  args = parser.parse_args()
  # Check arguments
  if not os.path.isdir(args.inputDirectory):
    print ('[Error] There is no directory called '+args.inputDirectory)
    sys.exit()
  if len(glob.glob(args.inputDirectory+'/*.txt')) == 0:
    print ('[Error] There is no txt files in '+args.inputDirectory)
    sys.exit()
  if os.path.exists(args.outputFile):
    print ('[Error] There is already a file called '+args.outputFile)
    sys.exit()
  return args

if __name__ == "__main__":
  args = parseArguments()

  dataCardPaths = glob.glob(args.inputDirectory+'/*.txt')

  exitCode = os.system("./compile.py")
  if (exitCode): sys.exit()

  for datacardPath in dataCardPaths:
    print ("Processes " + datacardPath)
    cmd = "./run/hig/scan_point.exe --cards -f "+datacardPath + " >> " + args.outputFile
    os.system(cmd)

  print('Run ./run/hig/plot_limit.exe -n -f '+args.outputFile)
