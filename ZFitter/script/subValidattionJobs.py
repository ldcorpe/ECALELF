#!/usr/bin/env python

from optparse import OptionParser
parser = OptionParser()
#parser.add_option("-b","--bkgfilename",help="Data and background workspace file")
#parser.add_option("-s","--sigfilename",help="Signal file (can be binned or parametric) or left blank")
#parser.add_option("-c","--cats",type="int",help="Number of categories to run")
#parser.add_option("-l","--catLabels",default="mk_default",help="Category labels (comma separated) default will use Category %cat")
#parser.add_option("-S","--sqrts",type='int',default=8,help="Sqrt(S) COM energy for finding strings etc")
#parser.add_option("--isMultiPdf",default=False,action="store_true",help="Use for multipdf workspaces")
#parser.add_option("--doBands",default=False,action="store_true",help="Use to draw bands")
#parser.add_option("--useBinnedData",default=False,action="store_true",help="Use binned data")
#parser.add_option("--makeCrossCheckProfPlots",default=False,action="store_true",help="Make some cross check plots - this is very slow!!")
#parser.add_option("--massStep",type="float",default=0.5,help="Mass step for calculating bands. Use a large number like 5 for quick running")
#parser.add_option("--nllTolerance",type="float",default=0.05,help="Tolerance for nll calc in %")
#parser.add_option("--blind",default=False,action="store_true",help="Blind the mass spectrum in the range [110,150]")
#parser.add_option("--runLocal",default=False,action="store_true",help="Run locally")

parser.add_option("-a","--auto",default="",help="Submit jobs from list (%default)")
parser.add_option("-b","--baseDir",default="test",help="Out directory for plots default: %default")
parser.add_option("-c","--validationConfig",default="data/validation/test_75x.dat",help="Out directory for plots default: %default")
parser.add_option("--floatTails",default="fixTails",action="store_true",help="fixTails or floatTails")
parser.add_option("-q","--queue",default="cmscaf1nd")
parser.add_option("-m","--mass",default="invMass_SC")
parser.add_option("-g","--gitTag",default="lcorpe:tag")
parser.add_option("-v","--verbose",default=False,action="store_true",help="Print more output")
parser.add_option("--check",default=False,action="store_true",help="check jobs and resubmit")
(options,args) = parser.parse_args()

import os
import time
fail=0
done=0
run=0
if not options.check:
  os.system('mkdir -p jobs')
  counter=0;
  os.system('rm jobs/*')
  if options.auto :
    with open('%s'%options.auto,'r') as f:
        for line in f:
                counter += 1
                if line[0]=="#" :
    							continue
                options.validationConfig=line.split(" ")[0]
                options.baseDir=line.split(" ")[1]
                options.floatTails=line.split(" ")[2]
                options.gitTag=line.split(" ")[3]
                options.mass=line.split(" ")[4]
							
    					
    						
                #print "job ",options.validationConfig, options.baseDir, options.floatTails
                os.system('mkdir -p %s'%options.baseDir)
                f = open('jobs/sub%d.sh'%(counter),'w')
                f.write('#!/bin/bash\n')
                f.write('cd %s\n'%os.getcwd())
                f.write('eval `scramv1 runtime -sh`\n')
                #execLine = '$CMSSW_BASE/src/Calibration/ZFitter/script/validation.sh -f %s --invMass_var invMass_SC --baseDir %s --validation \"%s\"'%(options.validationConfig,options.baseDir)
                execLine = '$CMSSW_BASE/src/Calibration/ZFitter/script/validation.sh  -f %s --invMass_var %s --baseDir %s --slides --validation'%(options.validationConfig,options.mass,options.baseDir)
                print execLine
                if options.floatTails=="floatTails" :
                	execLine += ' --floatTails'
                execLine += '\n'
                f.write('if ( %s ) then\n'%execLine);
                f.write('touch jobs/sub%d.sh.done\n'%counter);
                f.write('else\n');
                f.write('touch jobs/sub%d.sh.fail\n'%counter);
                f.write('fi\n');
                f.close()
                 
                os.system('chmod +x %s'%f.name)
                os.system('bsub -q %s -o %s.log %s'%(options.queue,os.path.abspath(f.name),os.path.abspath(f.name)))
                time.sleep(1) #wait for a few seconds so the jobs don't cause each other to fail by accessing stuff in tmp
else:
  print "[INFO] checking jobs, since --check was ", options.check
  counter=0;
  if options.auto :
    with open('%s'%options.auto,'r') as f:
        for line in f:
                counter += 1
                if line[0]=="#" :
    							continue
                options.validationConfig=line.split(" ")[0]
                options.baseDir=line.split(" ")[1]
                options.floatTails=line.split(" ")[2]
                options.gitTag=line.split(" ")[3]
                options.mass=line.split(" ")[4]
    					
                if os.path.exists('jobs/sub%d.sh.fail'%counter) or (not os.path.exists('jobs/sub%d.sh'%counter)) :
                  fail += 1
                  os.system('rm jobs/sub%d.sh*'%counter)
                  print "[RESUBMIT]resubmit job ", counter
                  os.system('mkdir -p %s'%options.baseDir)
                  f = open('jobs/sub%d.sh'%(counter),'w')
                  f.write('#!/bin/bash\n')
                  f.write('cd %s\n'%os.getcwd())
                  f.write('eval `scramv1 runtime -sh`\n')
                  execLine = '$CMSSW_BASE/src/Calibration/ZFitter/script/validation.sh  -f %s --invMass_var %s --baseDir %s --slides --validation'%(options.validationConfig,options.mass,options.baseDir)
                  print execLine
                  if options.floatTails=="floatTails" :
                	  execLine += ' --floatTails'
                  execLine += '\n'
                  f.write('if ( %s ) then\n'%execLine);
                  f.write('touch jobs/sub%d.sh.done\n'%counter);
                  f.write('else\n');
                  f.write('touch jobs/sub%d.sh.fail\n'%counter);
                  f.write('fi\n');
                  f.close()
                 
                  os.system('chmod +x %s'%f.name)
                  os.system('bsub -q %s -o %s.log %s'%(options.queue,os.path.abspath(f.name),os.path.abspath(f.name)))
                  time.sleep(1) #wait for a few seconds so the jobs don't cause each other to fail by accessing stuff in tmp
                elif (os.path.exists('jobs/sub%d.sh.done'%counter)):
									print "[SUCESSS] Job ", counter," finished OK!"
									done += 1
                else:
									print "[RUNNING] Job ", counter," still running!"
									run += 1
  print "SUMMARY"
  print "======="
  print "done ", done
  print "fail ", fail
  print "run ", run
  
  
exit



