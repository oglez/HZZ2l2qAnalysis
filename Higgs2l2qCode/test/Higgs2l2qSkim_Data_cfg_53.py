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

# # En CIEMAT:
#     'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleElectron/AOD/22Jan2013-v1/10002/B03E8149-2C90-E211-BB2E-003048678D86.root'
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleElectron/AOD/22Jan2013-v1/10002/B05B1B48-4990-E211-B54B-002618943935.root'
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleElectron/AOD/22Jan2013-v1/10002/B075F2E9-4290-E211-9A88-002618943811.root'
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleElectron/AOD/22Jan2013-v1/10002/B29CF47C-4D90-E211-85BD-003048679012.root'
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleElectron/AOD/22Jan2013-v1/10002/B2D02711-4390-E211-96A7-002618943918.root'
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleElectron/AOD/22Jan2013-v1/10002/B465052E-7590-E211-B351-003048678ED2.root'
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleElectron/AOD/22Jan2013-v1/10002/B473D278-2790-E211-BD66-0025905964C0.root'
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleElectron/AOD/22Jan2013-v1/10002/B4856543-2B90-E211-86E4-002618943894.root'
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleElectron/AOD/22Jan2013-v1/10002/B4E1280B-7690-E211-B3B8-003048678B3C.root'
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleElectron/AOD/22Jan2013-v1/10002/B648F0EB-9B90-E211-8D1F-003048D15D22.root'
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleElectron/AOD/22Jan2013-v1/10002/B6B8F071-9B90-E211-A8A9-002618943857.root'
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleElectron/AOD/22Jan2013-v1/10002/B6D6BF6C-2C90-E211-9254-00261894396D.root'
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleElectron/AOD/22Jan2013-v1/10002/B6DF33FB-5290-E211-9E88-00261894386E.root'
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleElectron/AOD/22Jan2013-v1/10002/B8B65BD9-3790-E211-B726-003048D15D22.root'
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleElectron/AOD/22Jan2013-v1/10002/B8BB7D40-3B90-E211-9B2C-002354EF3BE1.root'
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleElectron/AOD/22Jan2013-v1/10002/B8D29856-4890-E211-9853-00261894398D.root'
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleElectron/AOD/22Jan2013-v1/10002/BA57A596-3D90-E211-A087-002618FDA26D.root'
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleElectron/AOD/22Jan2013-v1/10002/BA69E1E8-5190-E211-ACFE-002590596486.root'
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleElectron/AOD/22Jan2013-v1/10002/BAB26259-3F90-E211-A826-003048D42D92.root'
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleElectron/AOD/22Jan2013-v1/10002/BC12219D-3090-E211-B2B4-00261894396A.root'
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleElectron/AOD/22Jan2013-v1/10002/BC3F3F07-5D90-E211-8D41-003048FFD728.root'
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleElectron/AOD/22Jan2013-v1/10002/BC6D03C3-5E90-E211-AE98-0025905938B4.root'
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleElectron/AOD/22Jan2013-v1/10002/BE31C1D2-3490-E211-ADCE-003048FFD71A.root'
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleElectron/AOD/22Jan2013-v1/10002/BE625CD9-5C90-E211-875C-002590596498.root'
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleElectron/AOD/22Jan2013-v1/10002/BE64C620-8190-E211-9743-002618943946.root'
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleElectron/AOD/22Jan2013-v1/10002/BEB82F9D-5190-E211-ACC8-002618943963.root'
    
   # data: mumu
    ,'file:/data1/oglez/spool/data2012d_22jan2013_mumu_AOD537/FA55EDA0-3B84-E211-B731-00259073E406.root'
    ,'file:/data1/oglez/spool/data2012d_22jan2013_mumu_AOD537/C2E5C34A-5484-E211-8623-E0CB4E19F95A.root'
    ,'file:/data1/oglez/spool/data2012d_22jan2013_mumu_AOD537/A0280890-3184-E211-AC22-002590747E16.root'
    ,'file:/data1/oglez/spool/data2012d_22jan2013_mumu_AOD537/48170530-2B84-E211-AFAA-90E6BA442EEB.root'

# En CIEMAT:
#      ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleMuParked/AOD/22Jan2013-v1/30002/A0280890-3184-E211-AC22-002590747E16.root'
#      ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleMuParked/AOD/22Jan2013-v1/30002/A04D4DA4-2684-E211-88A7-485B39800C13.root'
#      ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleMuParked/AOD/22Jan2013-v1/30002/A0523E1E-2584-E211-85FB-E0CB4E29C504.root'
#      ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleMuParked/AOD/22Jan2013-v1/30002/A05DDA7A-3584-E211-88EA-485B39897264.root'
#      ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleMuParked/AOD/22Jan2013-v1/30002/A095FEC5-3884-E211-A52E-20CF300E9EAD.root'
#      ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleMuParked/AOD/22Jan2013-v1/30002/A250486C-B384-E211-89FD-E0CB4E1A1177.root'
#      ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleMuParked/AOD/22Jan2013-v1/30002/A2A76DF3-4B84-E211-83E3-E0CB4E4408DF.root'
#      ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleMuParked/AOD/22Jan2013-v1/30002/A2B6A992-7084-E211-A41F-20CF300E9EB4.root'
#      ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleMuParked/AOD/22Jan2013-v1/30002/A4280D7C-2D84-E211-80FF-485B3989725B.root'
#      ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleMuParked/AOD/22Jan2013-v1/30002/A63815D8-3484-E211-BA59-00259073E3A6.root'
#      ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleMuParked/AOD/22Jan2013-v1/30002/A65B7F52-4684-E211-A1CC-00259073E3E6.root'
#      ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleMuParked/AOD/22Jan2013-v1/30002/A661115E-7084-E211-80E2-90E6BA19A25D.root'
#      ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleMuParked/AOD/22Jan2013-v1/30002/A696DAF0-4D84-E211-A7FD-90E6BA0D09E7.root'
#      ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleMuParked/AOD/22Jan2013-v1/30002/A86F964B-7084-E211-9BDA-E0CB4E29C4F1.root'
#      ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleMuParked/AOD/22Jan2013-v1/30002/A8E60EE9-5684-E211-95ED-E0CB4E1A1174.root'
#      ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleMuParked/AOD/22Jan2013-v1/30002/AA046A30-2B84-E211-BDBB-90E6BA442F3D.root'
#      ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleMuParked/AOD/22Jan2013-v1/30002/AA395B74-7084-E211-B0DC-20CF305B04D1.root'
#      ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleMuParked/AOD/22Jan2013-v1/30002/AAACCFA1-2984-E211-8265-002590747E1C.root'
#      ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleMuParked/AOD/22Jan2013-v1/30002/AAFE5369-4F84-E211-BA08-E0CB4E1A1180.root'
#      ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleMuParked/AOD/22Jan2013-v1/30002/AE056B33-2784-E211-AE5E-BCAEC54B302D.root'
#      ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleMuParked/AOD/22Jan2013-v1/30002/AE0D91A6-2684-E211-A6C2-00259074AE8C.root'
#      ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleMuParked/AOD/22Jan2013-v1/30002/AE40ACC5-3184-E211-AEF5-90E6BAE8CC1F.root'
#      ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleMuParked/AOD/22Jan2013-v1/30002/AE7AEB2A-4E84-E211-888B-20CF30561716.root'
#      ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleMuParked/AOD/22Jan2013-v1/30002/AECB2F36-2B84-E211-8FD0-00259074AE8C.root'
#      ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleMuParked/AOD/22Jan2013-v1/30002/AEF5B753-3384-E211-A650-E0CB4E19F9B8.root'
#      ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/data/Run2012D/DoubleMuParked/AOD/22Jan2013-v1/30002/AEFB6C1E-4984-E211-9EDB-00261834B526.root'
]

################
# Loading theglobal part for version
#
from HZZ2l2qAnalysis.Higgs2l2qCode.Higgs2l2qSkim_global53X_cfi import *

# For testing:

# To activate the debug of candidates:
# # process.p.replace(process.allhcand,process.allhcand*process.debugCandidates)

#########################################################################
