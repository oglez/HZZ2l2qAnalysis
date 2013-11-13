#!/bin/bash -f
#
# This code performs some clean-up of directories to reduce the size
# of the sandbox that is sent to the Grid using crab
#
# USAGE HZZ2l2qAnalysis/Higgs2l2qCode/test/cleanforcrab.sh
#
#
# ./CMGTools/External/data
mv ./CMGTools/External/data/ ./CMGTools/External/data_NOTNEEDED
mkdir ./CMGTools/External/data/
mv ./CMGTools/External/data_NOTNEEDED/TMVAClassificationCategory_JetID_53X_Dec2012.weights.xml \
   ./CMGTools/External/data/
mv ./CMGTools/External/data_NOTNEEDED/TMVAClassificationCategory_JetID_MET_53X_Dec2012.weights.xml \
   ./CMGTools/External/data/
mv ./CMGTools/External/data_NOTNEEDED/dummy.txt ./CMGTools/External/data/

# ./EgammaAnalysis/ElectronTools/data
mv ./EgammaAnalysis/ElectronTools/data ./EgammaAnalysis/ElectronTools/data_NOTNEEDED
mkdir ./EgammaAnalysis/ElectronTools/data
mv ./EgammaAnalysis/ElectronTools/data_NOTNEEDED/eleEnergyRegWeights_WithSubClusters_VApr15.root \
   ./EgammaAnalysis/ElectronTools/data
mv ./EgammaAnalysis/ElectronTools/data_NOTNEEDED/linearityNewReg-May2013.csv \
   ./EgammaAnalysis/ElectronTools/data
mv ./EgammaAnalysis/ElectronTools/data_NOTNEEDED/scalesNewReg-May2013.csv \
   ./EgammaAnalysis/ElectronTools/data

# Other useless packages (I think)
mv ./CMGTools/Common/python/skims ./CMGTools/Common/python/skims_NOTNEEDED
mv ./CMGTools/Common/python/release_info ./CMGTools/Common/python/release_info_NOTNEEDED
mv ./CMGTools/H2TauTau/data ./CMGTools/H2TauTau/data_NOTNEEDED
mv ./CMGTools/H2TauTau/python ./CMGTools/H2TauTau/python_NOTNEEDED
mv ./CMGTools/Production/scripts ./CMGTools/Production/scripts_NOTNEEDED
mv ./CMGTools/Utilities/data ./CMGTools/Utilities/data_NOTNEEDED
mv ./DPGAnalysis/SiStripTools/test ./DPGAnalysis/SiStripTools/test_NOTNEEDED
mv ./DPGAnalysis/SiStripTools/python ./DPGAnalysis/SiStripTools/python_NOTNEEDED
mv ./PhysicsTools/PatUtils/data ./PhysicsTools/PatUtils/data_NOTNEEDED
mv ./ExoDiBosonResonances/GeneratorStudies/highmasspz_pku \
   ./ExoDiBosonResonances/GeneratorStudies/highmasspz_pku_NOTNEEDED
#
