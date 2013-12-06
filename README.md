ECALELF
=======

Code for electron calibration in CMS
Testing


======
Download instructions.

cp ~shervin/public/setup-newMaster.sh .
./setup.sh CMSSW_5_3_7_patch6
cd CMSSW_5_3_7_patch6/src/
cmsenv

NB: I'm not able to say to scram not to compile a particular directory, then the Calibration/ZFitter directory is not skipped by scram and the compilation fails. you should compile by hand all the packages except Calibration/ZFitter
In order to compile Calibration/ZFitter do

cd Calibration/ZFitter
make


======
Instructions for developments:
fork the repository in GIT to your own area
git@github.com:ECALELFS/ECALELF.git

Add your git repository to the list of remotes to be used like:
git remote add -f myfork git@github.com:hengne/ECALELF.git
git checkout newMaster

Create a new branch for your development (use a meaningful name)
git branch myNewBranch
Push it to your git repository
git branch push myNewBranch myfork

Then start to develop, remember to do commits as much as possible


======
Code documentation:
Once downloaded the code, in Calibration/ you can run the command

doxygen fulldoc

to have the code documentation produced by doxygen opening the doc/doxygen/fulldoc/html/index.html with your browser 

