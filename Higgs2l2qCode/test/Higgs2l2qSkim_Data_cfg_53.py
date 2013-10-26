#
# New way of defining the code: global code for the version, specific
# "header-like" files with parameters to separate the arguments
#

import HZZ2l2qAnalysis.Higgs2l2qCode.Hzz2l2qSetup_cfi as Hzz2l2qSetup

# Configuration for DATA:

Hzz2l2qSetup.runOnMC = False

#Hzz2l2qSetup.usingevents = 4999
#Hzz2l2qSetup.usingevents = 99


# Files to run over:

Hzz2l2qSetup.usingfiles = [

    # data: ee
    'file:/data1/oglez/spool/data2012d_22jan2013_ee_AOD537/FEAA6EC6-DF8F-E211-B96F-00261894392D.root'
    ,'file:/data1/oglez/spool/data2012d_22jan2013_ee_AOD537/FC706CB1-1490-E211-9274-003048FFCB96.root'
    ,'file:/data1/oglez/spool/data2012d_22jan2013_ee_AOD537/F6ED4298-D18F-E211-BC8C-003048FFD770.root'
    ,'file:/data1/oglez/spool/data2012d_22jan2013_ee_AOD537/F0319D12-E78F-E211-9743-002618943856.root'
    ,'file:/data1/oglez/spool/data2012d_22jan2013_ee_AOD537/D6DE85AC-9A8F-E211-AAB7-00248C65A3EC.root'    

    # data: mumu
    ,'file:/data1/oglez/spool/data2012d_22jan2013_mumu_AOD537/FA55EDA0-3B84-E211-B731-00259073E406.root'
    ,'file:/data1/oglez/spool/data2012d_22jan2013_mumu_AOD537/C2E5C34A-5484-E211-8623-E0CB4E19F95A.root'
    ,'file:/data1/oglez/spool/data2012d_22jan2013_mumu_AOD537/A0280890-3184-E211-AC22-002590747E16.root'
    ,'file:/data1/oglez/spool/data2012d_22jan2013_mumu_AOD537/48170530-2B84-E211-AFAA-90E6BA442EEB.root'
]

################
# Loading theglobal part for version
#
from HZZ2l2qAnalysis.Higgs2l2qCode.Higgs2l2qSkim_global53X_cfi import *

# For testing:

# To activate the debug of candidates:
# # process.p.replace(process.allhcand,process.allhcand*process.debugCandidates)

#########################################################################
