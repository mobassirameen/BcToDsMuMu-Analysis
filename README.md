
## Intructions to run the examples.
```
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc700
cmsrel CMSSW_10_2_18
cd CMSSW_10_2_18/src/
cmsenv
voms-proxy-init -voms cms -valid 192:00
git clone -b newbcdstomumucode https://:@gitlab.cern.ch:8443/ckar/bctodsmumu-analysis.git

rm -rf bctodsmumu-analysis/JPsiKsPAT

scram b -j8

cd bctodsmumu-analysis/BcToDsMuMuPAT/test
open the BcToDsMuMuRootupler.py and change the line 
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(2000))
to 
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(number of event you want to run))


cmsRun BcToDsMuMuRootupler.py

```



