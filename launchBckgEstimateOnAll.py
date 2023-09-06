import sys, os, time, re
import numpy as np
#from common_functions import *
from optparse import OptionParser
parser = OptionParser(usage="Usage: python %prog codeVersion")
(opt,args) = parser.parse_args()

datasetList=[
    #"\/opt\/sbg\/cms\/ui3_data1\/dapparu\/HSCP\/Production\/crab\_Analysis\_SingleMuon\_Run2017E\_CodeV76p0\_v1",
    #"\/opt\/sbg\/cms\/ui3_data1\/dapparu\/HSCP\/Production\/crab\_Analysis\_SingleMuon\_Run2017\_CodeV73p3\_v4",
    #"\/opt\/sbg\/cms\/ui3_data1\/dapparu\/HSCP\/Production\/crab\_Analysis\_SingleMuon\_Run2018\_CodeV73p3\_v4",
    #"\/opt\/sbg\/cms\/ui3_data1\/dapparu\/HSCP\/Production\/crab\_Analysis\_2018\_AllWJets\_CodeV73p0\_v1",
    #"\/opt\/sbg\/cms\/ui3_data1\/dapparu\/HSCP\/Production\/crab\_Analysis\_2018\_AllQCD\_CodeV73p0\_v1",
    #"\/opt\/sbg\/cms\/ui3_data1\/dapparu\/HSCP\/Production\/crab\_Analysis\_2018\_AllTTbar\_CodeV73p0\_v1",
    #"\/opt\/sbg\/cms\/ui3_data1\/dapparu\/HSCP\/Production\/crab\_Analysis\_2018\_AllZToMuMu\_CodeV73p0\_v1",
    #"\/opt\/sbg\/cms\/ui3_data1\/dapparu\/HSCP\/Production\/crab\_Analysis\_2018\_AllBackground\_CodeV73p0\_v1",
    #"\/opt\/sbg\/cms\/ui3_data1\/dapparu\/HSCP\/Production\/crab\_Analysis\_2018\_TTbarWjets\_CodeV73p0\_v1",
    #"\/opt\/sbg\/cms\/safe1\/cms\/dapparu\/HSCP\/CMSSW\_10\_6\_27\/src\/SUSYBSMAnalysis\/BackgroundPrediction\_backup\_20230122\/out\_muon\_pt200reg",
    "\/opt\/sbg\/cms\/ui3_data1\/dapparu\/HSCP\/Production\/crab\_Analysis\_SingleMuon\_Run2017\_CodeVUnB\_v1\_v1",
    "\/opt\/sbg\/cms\/ui3_data1\/dapparu\/HSCP\/Production\/crab\_Analysis\_SingleMuon\_Run2018\_CodeVUnB\_v1\_v1"
]
tagKC=[
    "data2017",
    "data2018",
    #"mc2018",
    #"mc2018",
    #"mc2018",
    #"mc2018",
    #"mc2018",
]
odir=[
    #"25april\_2017\_label\_testBinning\/",
    #"25april\_2018\_label\_testBinning\/",
    #"17april\_2018\_label\_wFit\/",
    #"17april\_2017\_label\_wFit\/",
    #"13april\_2017\_label\_ih0sigmapt2iso1\_wFit\/",
    #"13april\_2018\_label\_ih0sigmapt0iso1\_wFit\/",
    #"11april\_2017\_label\_VRnoBlinding\_wFit\_SR2\-3\/",
    #"11april\_2018\_label\_VRnoBlinding\_wFit\_SR2\-3\/",
    #"11april\_wjets\_label\_wFit\_\/",
    #"11april\_qcd\_label\_wFit\_\/",
    #"11april\_ttbar\_label\_wFit\_\/",
    #"12april\_ztomumu\_label\_wFit\_\/",
    #"11april\_allbckg\_label\_wFit\_\/",
    #"14april\_ttbarwjets\_label\_\/",
    #"testRaph\_pt200\/",
    #"21aug\_UnB\_v2\/"
    "31aug\_2017\_label\_UnB_v1_v2\/",
    "31aug\_2018\_label\_UnB_v1_v2\/",
]

nPE="200"

config=[
    #["nominal", "1", "1", "1", "0", "0", "1", "1"],
    #["eta60", "10", "4", "2", "0", "0", "1", "2"],
    #["eta50", "12", "4", "2", "0", "0", "1", "1"],
    #["eta40", "15", "4", "2", "0", "0", "1", "2"],
    #["eta30", "20", "4", "2", "0", "0", "1", "2"],
    #["eta24", "25", "4", "2", "0", "0", "1", "2"],
    #["eta20", "30", "4", "2", "0", "0", "1", "2"],
    #["eta12", "50", "4", "2", "0", "0", "1", "2"],
    #["eta10", "60", "4", "2", "0", "0", "1", "2"],
    
    #["nominal", "5", "4", "2", "0", "0", "1", "1"],
    #["etaup", "4", "4", "2", "0", "0", "1", "1"],
    #["etadown", "8", "4", "2", "0", "0", "1", "1"],
    #["ihup", "5", "2", "2", "0", "0", "1", "1"],
    #["ihdown", "5", "8", "2", "0", "0", "1", "1"],
    #["momup", "5", "4", "1", "0", "0", "1", "1"],
    #["momdown", "5", "4", "4", "0", "0", "1", "1"],
    #["corrIh", "5", "4", "2", "1", "0", "1", "1"],
    #["corrMom", "5", "4", "2", "0", "1", "1", "1"],
    #["FitIhUp", "5", "4", "2", "0", "0", "2", "1"],
    #["FitIhDown", "5", "4", "2", "0", "0", "0", "1"],
    #["FitMomUp", "5", "4", "2", "0", "0", "1", "2"],
    #["FitMomDown", "5", "4", "2", "0", "0", "1", "0"],

    ["nominal", "4", "4", "2", "0", "0", "1", "1"],
    ["etaup", "2", "4", "2", "0", "0", "1", "1"],
    ["etadown", "8", "4", "2", "0", "0", "1", "1"],
    ["ihup", "4", "2", "2", "0", "0", "1", "1"],
    ["ihdown", "4", "8", "2", "0", "0", "1", "1"],
    ["momup", "4", "4", "1", "0", "0", "1", "1"],
    ["momdown", "4", "4", "4", "0", "0", "1", "1"],
    ["corrIh", "4", "4", "2", "1", "0", "1", "1"],
    ["corrMom", "4", "4", "2", "0", "1", "1", "1"],
    ["FitIhUp", "4", "4", "2", "0", "0", "2", "1"],
    ["FitIhDown", "4", "4", "2", "0", "0", "0", "1"],
    ["FitMomUp", "4", "4", "2", "0", "0", "1", "2"],
    ["FitMomDown", "4", "4", "2", "0", "0", "1", "0"],

#    ["nominal", "5", "4", "2", "0", "0", "1", "1"],
#    ["etaup", "4", "4", "2", "0", "0", "1", "1"],
#    ["etadown", "8", "4", "2", "0", "0", "1", "1"],
#    ["ihup", "5", "2", "2", "0", "0", "1", "1"],
#    ["ihdown", "5", "8", "2", "0", "0", "1", "1"],
#    ["momup", "5", "4", "1", "0", "0", "1", "1"],
#    ["momdown", "5", "4", "4", "0", "0", "1", "1"],
#    ["corrIh", "5", "4", "2", "1", "0", "1", "1"],
#    ["corrMom", "5", "4", "2", "0", "1", "1", "1"],
#    ["FitIhUp", "5", "4", "2", "0", "0", "2", "1"],
#    ["FitIhDown", "5", "4", "2", "0", "0", "0", "1"],
#    ["FitMomUp", "5", "4", "2", "0", "0", "1", "2"],
#    ["FitMomDown", "5", "4", "2", "0", "0", "1", "0"],
]


i=0
for dataset in datasetList:
    print("Launch on dataset: "+dataset)
    os.system("cp configFile_readHist_template.txt configFile_readHisto_toLaunch_tmp.txt")
    os.system("sed -i 's/sample/"+dataset+"/g' configFile_readHisto_toLaunch_tmp.txt")
    os.system("sed -i 's/tag_KC/"+tagKC[i]+"/g' configFile_readHisto_toLaunch_tmp.txt")
    os.system("sed -i 's/dir/"+odir[i]+"/g' configFile_readHisto_toLaunch_tmp.txt")
    os.system("sed -i 's/nPE/"+nPE+"/g' configFile_readHisto_toLaunch_tmp.txt")
    for conf in config:
        os.system("cp configFile_readHisto_toLaunch_tmp.txt configFile_readHisto_toLaunch.txt")
        os.system("sed -i 's/label/"+conf[0]+"/g' configFile_readHisto_toLaunch.txt")
        os.system("sed -i 's/rebinEta/"+conf[1]+"/g' configFile_readHisto_toLaunch.txt")
        os.system("sed -i 's/rebinIh/"+conf[2]+"/g' configFile_readHisto_toLaunch.txt")
        os.system("sed -i 's/rebinMom/"+conf[3]+"/g' configFile_readHisto_toLaunch.txt")
        os.system("sed -i 's/corrTemplateIh/"+conf[4]+"/g' configFile_readHisto_toLaunch.txt")
        os.system("sed -i 's/corrTemplateMom/"+conf[5]+"/g' configFile_readHisto_toLaunch.txt")
        os.system("sed -i 's/fitIh/"+conf[6]+"/g' configFile_readHisto_toLaunch.txt")
        os.system("sed -i 's/fitMom/"+conf[7]+"/g' configFile_readHisto_toLaunch.txt")

        os.system("cat configFile_readHisto_toLaunch.txt")
        os.system("time root -l -q -b step2_backgroundPrediction.C")    
    i+=1

