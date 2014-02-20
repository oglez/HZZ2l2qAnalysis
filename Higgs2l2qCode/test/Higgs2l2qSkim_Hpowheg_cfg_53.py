#
# New way of defining the code: global code for the version, specific
# "header-like" files with parameters to separate the arguments
#

import HZZ2l2qAnalysis.Higgs2l2qCode.Hzz2l2qSetup_cfi as Hzz2l2qSetup

# Configuration for MC (Signal):

Hzz2l2qSetup.lhcEnergy = "8TeV"
Hzz2l2qSetup.higgsMass = 900

# Hzz2l2qSetup.usingevents = 5198

# Files to run over:

Hzz2l2qSetup.usingfiles = [
    'file:/data1/oglez/spool/Higgs900_gg_zzllqq_POWHEG_AOD532/A6B3E111-8602-E211-8B03-90E6BA442F12.root',
    'file:/data1/oglez/spool/Higgs900_gg_zzllqq_POWHEG_AOD532/26A46B19-8E02-E211-83F1-E0CB4E1A117D.root',
    'file:/data1/oglez/spool/Higgs900_gg_zzllqq_POWHEG_AOD532/48614CD2-7402-E211-88ED-485B39800BF3.root'
    ]

################
# Loading theglobal part for version
#
from HZZ2l2qAnalysis.Higgs2l2qCode.Higgs2l2qSkim_global53X_cfi import *

# For testing:

# To activate the debug of candidates:
# process.p.replace(process.allhcand,process.allhcand*process.debugCandidates)

#########################################################################

