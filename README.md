Use the code in this repository to export hit information 
to be used in the GPU-based reconstruction of trax.

This code works with the CMSSW releases 5_3_X

cd CMSSW_5_3_X/src/traxExport/Export
cmsDriver.py Configuration/Generator/python/TTbar_8TeV_cfi.py -n 10 -s GEN,SIM --conditions auto:mc --datatier GEN-SIM-RAW --eventcontent RAWSIM
cmsDriver.py RECO -s DIGI,L1,DIGI2RAW,RAW2DIGI,RECO --conditions auto:mc --number 10 --eventcontent FEVTDEBUG --mc --filein file:./TTbar_8TeV_cfi_py_GEN_SIM.root 
cmsRun python/exportHits.py



