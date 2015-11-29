#!/usr/bin/env python

#
#Script for submitting toys in which new datasets are generated before being fit
#

print '=== Welcome to startToys.py ==='
print 'This script is used to generate and fit toy datasets, which is the procedure for most systematic checks (pretty much anything not involving the CLEO strong phase parameters).'
print 'If you want to rerun the fit to data but change some input parameters use startMultipleData.py instead.'
print 'For more information about which script to use, see ../README_sytsematics'
print

import os, time, sys, pwd, random, datetime, shutil

uname = pwd.getpwuid(os.getuid()).pw_name
todaysdate = datetime.date.today()

#Base directory

basedir = "/data/lhcb/users/" + uname + "/B02DKstar_Kshh_toys"
if not os.path.exists(basedir):
    print 'Creating', basedir
    os.mkdir(basedir)
    os.chmod(basedir, 0750)

#Toy description
    
desc = raw_input("Enter a description of the toys you are running (today's date will automatically be included): ")

if desc == '':
    print 'You entered a blank string so I am using "toys"'
    desc = 'toys'

toyDir = basedir + '/' + str(todaysdate) + "_" + desc

if os.path.exists(toyDir):
    print
    print '*** WARNING: directory', toyDir, 'already exists! ***'
    decision = raw_input('*** Input Y to overwrite, anything else to abort: ').lower()
    if decision != 'y':
        sys.exit()
    else:
        shutil.rmtree(toyDir)
        print 'Directory has been removed and will be recreated'
print 'Creating', toyDir
os.mkdir(toyDir)
print

#Copy code to new directory

dirsToCopy = ['Settings', 'Inputs', 'ToyReading', 'ToyReading_MINOS', 'src']
dirsToMake = ['bin', 'obj', 'figs', 'figs/fits', 'figs/residuals', 'output']
              
print 'Setting up code'
os.mkdir(toyDir + '/B02DKstar_Kshh_MassFit_copy')
for dtm in dirsToMake:
    print ' creating subdirectory', dtm
    os.mkdir(toyDir + '/B02DKstar_Kshh_MassFit_copy/' + dtm)
for dtc in dirsToCopy:
    print ' copying subdirectory', dtc
    shutil.copytree('../' + dtc, toyDir + '/B02DKstar_Kshh_MassFit_copy/' + dtc)
    #shutil.rmtree(toyDir + '/B02DKstar_Kshh_MassFit_copy/' + dtc + '/.svn')
print ' copying Makefile'    
shutil.copy('../Makefile', toyDir + '/B02DKstar_Kshh_MassFit_copy')
print

#remove useless .svn directory
os.system('find '+toyDir+' -name ".svn" ')
os.system('find '+toyDir+' -name ".svn" -exec rm -Rf "{}" \;')

#manual .svn removal in subdirectories
#shutil.rmtree(toyDir + '/B02DKstar_Kshh_MassFit_copy/Settings/PDFShapes/.svn')
#shutil.rmtree(toyDir + '/B02DKstar_Kshh_MassFit_copy/Settings/PDFShapes/Fit/.svn')
#shutil.rmtree(toyDir + '/B02DKstar_Kshh_MassFit_copy/Settings/PDFShapes/Gen/.svn')

#copy over binary or recompile?
cprun = raw_input('Input Y to recompile, anything else to use precompiled ../bin/run: ').lower()
if cprun != 'y':
    print ' copying ../bin/run'
    shutil.copy('../bin/run', toyDir + '/B02DKstar_Kshh_MassFit_copy/bin')
else:
    print ' compiling code'
    initDir = os.getcwd();
    os.chdir(toyDir + '/B02DKstar_Kshh_MassFit_copy')
    os.system('make clean && make')
    os.chdir(initDir)
print

#Job details

seedInput = int(raw_input("Enter random seed: ")) #could do seedInput = random.randint(0, 10000), but no real point
numBatchJobs = int(raw_input("Enter number of batch jobs: "))
numMCStudyToysPerBatchJob = int(raw_input("Enter number of toy studies per batch job (recommended: 2--10): "))
print

#Function to create submission scripts

def makeSubmitScript(toyDir, seed):
    print 'Making submission script for seed', seed
    initDir = os.getcwd();
    dirname = toyDir
    os.chdir(dirname)
    saved_stdout = sys.stdout
    sys.stdout = open('BatchSubmit.sh', 'w')
    print '#!/bin/bash'
    print '#PBS -j oe'
    print '#PBS -l ncpus=1'
    print '#PBS -l cput=1:59:59'
    print '. /data/lhcb/sw/scripts/lbsetup-cvmfs-osagnostic.sh'
    print '. SetupProject.sh Gaudi v25r0 ROOT\n '
    print 'cd ' + dirname + '/../B02DKstar_Kshh_MassFit_copy/'
    print './bin/run ../toySeed_' + seed + '/GeneralSettings.txt > ../toySeed_' + seed + '/runOutput.txt'
    sys.stdout = saved_stdout
    os.chdir(initDir)
    
#Start toys
    
jobCounter = 1
for toySeed in range(seedInput, numBatchJobs * numMCStudyToysPerBatchJob + seedInput, numMCStudyToysPerBatchJob):
    # Make toy directory
    toyDirName = toyDir + "/toySeed_" + str(toySeed) + "/"
    print 'Creating subdirectory', toyDirName
    os.mkdir(toyDirName)

    # Create copy of fit code parameter file there and set the seed
    shutil.copy('../Settings/GeneralSettings.txt', toyDirName + 'GeneralSettings.txt')
    f_fitParFile = open(toyDirName + 'GeneralSettings.txt') #original file
    f_newFitParFile = open(toyDirName + 'GeneralSettings_new.txt', 'w') #modified file
    for line in f_fitParFile.readlines():
        if(line.find('doFit') == 0):
            f_newFitParFile.write('doFit true\n')
        elif(line.find('drawProjections') == 0):
            f_newFitParFile.write('drawProjections false\n')
        elif(line.find('startSeed') == 0):
            f_newFitParFile.write('startSeed ' + str(toySeed) + "\n")
        elif(line.find('nToys') == 0):
            f_newFitParFile.write('nToys ' + str(numMCStudyToysPerBatchJob) + '\n')
        elif(line.find('genToys') == 0):
            f_newFitParFile.write('genToys true\n')
        elif(line.find('readToys') == 0):
            f_newFitParFile.write('readToys false\n')
        elif(line.find('toyLocation') == 0):
            f_newFitParFile.write('toyLocation ../toySeed_' + str(toySeed) + '/\n')
        elif(line.find('UNBLIND') == 0):
            f_newFitParFile.write('UNBLIND true\n')
        else:
            f_newFitParFile.write(line)
    f_fitParFile.close();
    f_newFitParFile.close();
    shutil.move(toyDirName + 'GeneralSettings_new.txt', toyDirName + 'GeneralSettings.txt')

    # Make batch submit script
    makeSubmitScript(toyDirName, str(toySeed))

    # Submit script from that directory
    initDir = os.getcwd();
    os.chdir(toyDirName)
    print 'Submitting job ' + str(jobCounter) + " with seed " + str(toySeed)
    time.sleep(0.1)
    os.system('qsub -N Toy_' + str(toySeed) + ' BatchSubmit.sh')
    os.chdir(initDir)
    jobCounter += 1
    print

