#
# New way of defining the code: global code for the version, specific
# "header-like" files with parameters to separate the arguments
#

import HZZ2l2qAnalysis.Higgs2l2qCode.Hzz2l2qSetup_cfi as Hzz2l2qSetup

# Configuration for MC (Background):

# Default Hzz2l2qSetup.runOnMC = True

#Hzz2l2qSetup.usingevents = 4998
#Hzz2l2qSetup.usingevents = 98

# Files to run over:

Hzz2l2qSetup.usingfiles = [
    'file:/data1/oglez/spool/Zjets_MADGRAPH_v1_AOD532/FE414F4B-F6D2-E111-A4E9-003048674048.root'
    ,'file:/data1/oglez/spool/Zjets_MADGRAPH_v1_AOD532/EEBDBD8A-C4D1-E111-949E-003048673EBA.root'
    ,'file:/data1/oglez/spool/Zjets_MADGRAPH_v1_AOD532/DEEA6DD5-00D3-E111-B6F2-001E67397F3F.root'
    ,'file:/data1/oglez/spool/Zjets_MADGRAPH_v1_AOD532/AA4191D6-87D2-E111-A0B3-001E67397F71.root'
    ,'file:/data1/oglez/spool/Zjets_MADGRAPH_v1_AOD532/98C30DF6-5DD2-E111-853D-0025B3E05E10.root'
]

################
# Loading theglobal part for version
#
from HZZ2l2qAnalysis.Higgs2l2qCode.Higgs2l2qSkim_global53X_cfi import *

# For testing:

# To activate the debug of candidates:
# #process.p.replace(process.allhcand,process.allhcand*process.debugCandidates)

#########################################################################
