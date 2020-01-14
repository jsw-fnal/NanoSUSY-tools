# NanoSUSY-tools
Postprocessing script for Stop 0L analysis

See the [readme](python/processors/Condor/README.md) in "NanoSUSYTools/python/processors/Condor" for specific instructions for condor submission.

### Set up CMSSW

```tcsh
source /cvmfs/cms.cern.ch/cmsset_default.csh
setenv SCRAM_ARCH slc6_amd64_gcc700
cmsrel CMSSW_10_2_9
cd CMSSW_10_2_9/src/
cmsenv
```

### Set up NanoSusyTools framework
```tcsh
cd $CMSSW_BASE/src
cmsenv
git clone https://github.com/cms-tau-pog/TauIDSFs TauPOG/TauIDSFs
git clone -b Stop0l git@github.com:susy2015/nanoAOD-tools.git PhysicsTools/NanoAODTools
# For condor submission check the specific tag checkout instructions in [readme](python/processors/Condor/README.md)
git clone -b dev_v5 git@github.com:susy2015/NanoSUSY-tools.git PhysicsTools/NanoSUSYTools
git clone -b Stop0l_NanoAOD_production_V3.1 git@github.com:susy2015/TopTagger.git
scram b
cd $CMSSW_BASE/src/TopTagger/TopTagger/test
./configure
make
cmsenv
cd $CMSSW_BASE/src/PhysicsTools/NanoSUSYTools/python/processors
getTaggerCfg.sh -n -t DeepResolved_DeepCSV_GR_nanoAOD_2016_v1.0.3
getTaggerCfg.sh -n -t DeepResolved_DeepCSV_GR_nanoAOD_2017_v1.0.3
getTaggerCfg.sh -n -t DeepResolved_DeepCSV_GR_nanoAOD_2018_v1.0.3
```

### Run local MC test
```tcsh
python Stop0l_postproc.py -i file:[input file] -s [MC sample name (from sampleSet cfg file)] -e [year]
```

### Run local Data test
```tcsh
python Stop0l_postproc.py -i file:[input file] -d [data period] -e [year]
```


## To Do:
* PDF uncertainty module (Done)
    * weights stored in NanoAOD accordingly already, need code to extract the envelope
* lepton SF module follow [SUSY](https://twiki.cern.ch/twiki/bin/viewauth/CMS/SUSLeptonSF#Scale_Factors_for_SUSY_IDs) (Done with latest SF)
    * Need to update based upon the existing example code from NanoAOD-tools
* PU reweighting  (Done, to be tested)
    * Example code existed, but need recompute the pileup distribution from data
* btag SF update (instead of the btag weight stored during production) (done)
    * Need follow up with which method to be apply (iterative fit for DeepCSV?)
* JEC uncertainty update (done)
    * Need to understand the existing tool in NanoAOD-Tools
* DeepAK8/DeepResolved SF
    * Are they available yet?
* update Jet ID : [JetID Twiki](https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID13TeVRun2018) (Done with temp fix )
    * Do we need it in post-processing, or wait for next production?
* L1EcalPrefiring [twiki]* (https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1ECALPrefiringWeightRecipe#Call_the_producer_in_your_config)
    * Rumor is it will be included in next NanoAOD
    * If not, we will make a module for it
* Various systematics 
* Trigger path and efficiency:
    * Once Hui's study is finalized, we will store the bit and efficiency + systematic
