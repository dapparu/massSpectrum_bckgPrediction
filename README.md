# massSpectrum_bckgPrediction

This repositery is dedicated to HSCP analysis, for the background estimate of the mass spectrum. 

## Setup working area

```bash
export SCRAM_ARCH=slc7_amd64_gcc700
cmsrel CMSSW_10_6_27
cd CMSSW_10_6_27/src/
cmsenv
```

For the following step you should have a ssh key associated to your GitHub account.
For more information, see [connecting-to-github-with-ssh-key](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent).

```bash
git clone -b master git@github.com:dapparu/massSpectrum_bckgPrediction.git massSpectrum_bckgPrediction 
``` 

## Input file path

The path of input file is given here: 
-  Unblinded production: ```/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/crab_Analysis_SingleMuon_Run2017_CodeVUnB_v1_v1.root``` and ```/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/crab_Analysis_SingleMuon_Run2018_CodeVUnB_v1_v1.root```. You will also find files for the different eras. 
-  Blinded production: ```/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/crab_Analysis_SingleMuon_Run2017_CodeV73p3_v4.root``` and ```/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/crab_Analysis_SingleMuon_Run2018_CodeV73p3_v4.root```

## Run the background estimate code 

This part concerns the run of the background estimate method. 

```bash
python launchBckgEstimateOnAll.py
```

You can change on what you want to run directly in ```launchBckgEstimateOnAll.py``` and ```step2_backgroundPrediction.C``` files. 

In the first file (```launchBckgEstimateOnAll.py```): 
- You can set on which datasets you want to run, giving the path to the root file with all the needed histograms, in the ```confidatasetListg``` array. 
- The ```config``` array is used to indicate which kind of estimates are ran: nominal or the different systematics. 
- The ```nPe```variable set the number of pseudo-experiments done during the background estimate. 
- The ```odir``` array gives the directory where you can find fast produced plots. 

The code runs on 25 cores in parallel (in local) and it can be changed at line 798 of the file ```Regions.h```

In the second file (```step2_backgroundPrediction.C```):
- You can set which regions you want for estimation. Default: only SR1, SR2 and SR3.

After you ran the code, a new file is created with the histograms corresponding to the background estimate in the wanted regions, and labels are set to correspond to the different systematics cases. 

<!-- -->
**Important** Normalisation problem in CR
In the case you see normalisation problems (especially in control regions), be sure to apply the Ih > C cut at the preselection stage. 
<!-- -->

## Run the mass spectrum plotter code

This part concerns the run of mass spectrum plotter code, with nice style. 

The code runs with the command: 
```bash
python2.7 macroMass.py --ifile /opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/crab_Analysis_SingleMuon_Run2017_CodeVUnB_v1_v1_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_nPE200_test_v1.root --ofile test1 --region 999ias100 --odir test
```

Concerniing the options:
- option ```--ifile``` is for the input file; obtained at the previous step. 
- option ```--ofile``` is for the output file label. 
- option ```--region``` is for the region on which one wants to run. The possible values are: ```50ias60```, ```60ias70```, ```70ias80```, ```80ias90```, ```50ias90```, ```90ias100```, ```99ias100``` and ```999ias100```. 
- option ```--odir``` is for the output directory. 

The year on which runs is set directly in ```macroMass.py``` at line 204.
