#
# New way of defining the code: global code for the version, specific
# "header-like" files with parameters to separate the arguments
#

import HZZ2l2qAnalysis.Higgs2l2qCode.Hzz2l2qSetup_cfi as Hzz2l2qSetup

# Configuration for MC (Background):

# Default Hzz2l2qSetup.runOnMC = True

#Hzz2l2qSetup.usingevents = 4998
#Hzz2l2qSetup.usingevents = 150

# Files to run over:

Hzz2l2qSetup.usingfiles = [

       'root://eoscms//eos/cms/store/caf/user/oglez/aod_532/ZZ_PYTHIA_v1/0000CAC5-D4DA-E111-8872-00A0D1EEF328.root'
    
#      'file:/data1/oglez/spool/Zjets_MADGRAPH_v1_AOD532/FE414F4B-F6D2-E111-A4E9-003048674048.root'
#      ,'file:/data1/oglez/spool/Zjets_MADGRAPH_v1_AOD532/EEBDBD8A-C4D1-E111-949E-003048673EBA.root'
#      ,'file:/data1/oglez/spool/Zjets_MADGRAPH_v1_AOD532/DEEA6DD5-00D3-E111-B6F2-001E67397F3F.root'
#      ,'file:/data1/oglez/spool/Zjets_MADGRAPH_v1_AOD532/AA4191D6-87D2-E111-A0B3-001E67397F71.root'
#      ,'file:/data1/oglez/spool/Zjets_MADGRAPH_v1_AOD532/98C30DF6-5DD2-E111-853D-0025B3E05E10.root'

# # En CIEMAT: Top sample
#    'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/mc/Summer12_DR53X/TTTo2L2Nu2B_8TeV-powheg-pythia6/AODSIM/PU_S10_START53_V7A-v1/0000/FEBDF4A1-DE03-E211-A20E-00266CFAE30C.root' 

# # En CIEMAT:
#  
#     'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0002/607F3D99-96D4-E111-8C09-001E67396D42.root'
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-rball/AODSIM/PU_S10_START53_V7A-v1/0002/6086366E-5DD4-E111-B42B-003048D45F54.root' 
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-rball/AODSIM/PU_S10_START53_V7A-v1/0002/60D1CF17-20D4-E111-8608-0025B3E05D62.root'
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-rball/AODSIM/PU_S10_START53_V7A-v1/0002/60DDE76F-5DD4-E111-9B19-002481E154CE.root'
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-rball/AODSIM/PU_S10_START53_V7A-v1/0002/6201B6B2-5DD4-E111-BE01-001E67398CA0.root'
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-rball/AODSIM/PU_S10_START53_V7A-v1/0002/6290BCAB-1AD4-E111-AED6-0025B3E05CAA.root'
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-rball/AODSIM/PU_S10_START53_V7A-v1/0002/62DD55E0-09D4-E111-B318-001E673984FD.root' 
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-rball/AODSIM/PU_S10_START53_V7A-v1/0002/62FCD042-35D4-E111-BDBB-001E673971CA.root'
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-rball/AODSIM/PU_S10_START53_V7A-v1/0002/640779EE-06D4-E111-9776-003048D4777A.root' 
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-rball/AODSIM/PU_S10_START53_V7A-v1/0002/644D01E9-06D4-E111-B7CB-001E673973D2.root'
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-rball/AODSIM/PU_S10_START53_V7A-v1/0002/64FB3C00-42D4-E111-BCDA-001E67396761.root'
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-rball/AODSIM/PU_S10_START53_V7A-v1/0002/66901438-50D4-E111-8A53-0025B3E05DDA.root' 
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-rball/AODSIM/PU_S10_START53_V7A-v1/0002/66BF98A6-1AD4-E111-9635-003048D460F4.root'
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-rball/AODSIM/PU_S10_START53_V7A-v1/0002/66C934D3-10D4-E111-89CF-001E67396FCC.root' 
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-rball/AODSIM/PU_S10_START53_V7A-v1/0002/66E4C923-38D4-E111-BF14-003048D47A42.root' 
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-rball/AODSIM/PU_S10_START53_V7A-v1/0002/683AFA8A-0CD4-E111-9F59-003048D45FE2.root' 
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-rball/AODSIM/PU_S10_START53_V7A-v1/0002/685E79E0-54D4-E111-99D2-0025B31E3C58.root' 
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-rball/AODSIM/PU_S10_START53_V7A-v1/0002/68630313-3ED4-E111-9AD4-003048D46046.root' 
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-rball/AODSIM/PU_S10_START53_V7A-v1/0002/6865CBD6-0AD4-E111-A331-003048D45FD2.root' 
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-rball/AODSIM/PU_S10_START53_V7A-v1/0002/6AA8E1DD-57D4-E111-8535-003048D45F72.root' 
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-rball/AODSIM/PU_S10_START53_V7A-v1/0002/6AD9A848-31D4-E111-B6E5-001E6739834A.root' 
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-rball/AODSIM/PU_S10_START53_V7A-v1/0002/6C398DC3-4CD4-E111-9545-003048673E9C.root' 
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-rball/AODSIM/PU_S10_START53_V7A-v1/0002/6C54ABFF-0ED4-E111-8594-003048D476AE.root' 
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-rball/AODSIM/PU_S10_START53_V7A-v1/0002/6CB0E9F8-02D4-E111-B542-001E67397CCE.root' 
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-rball/AODSIM/PU_S10_START53_V7A-v1/0002/6CC9FD37-50D4-E111-830D-001E67396874.root' 
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-rball/AODSIM/PU_S10_START53_V7A-v1/0002/6CEC2CFD-79D4-E111-9D51-001E67396D51.root' 
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-rball/AODSIM/PU_S10_START53_V7A-v1/0002/6E62BEDA-10D4-E111-B38C-002481E14F5C.root' 
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-rball/AODSIM/PU_S10_START53_V7A-v1/0002/6E96AB0E-3ED4-E111-A74E-003048D4609A.root' 
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-rball/AODSIM/PU_S10_START53_V7A-v1/0002/6EA93EC6-1CD4-E111-8231-001E67397CAB.root' 
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-rball/AODSIM/PU_S10_START53_V7A-v1/0002/6EB746C5-70D4-E111-96F1-003048D47934.root' 
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-rball/AODSIM/PU_S10_START53_V7A-v1/0002/6EBFB0D4-44D4-E111-9A8C-002481E14F74.root' 
#     ,'dcap://gaeds015.ciemat.es:22125/pnfs/ciemat.es/data/cms/prod/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-rball/AODSIM/PU_S10_START53_V7A-v1/0002/6EC3CD7A-23D4-E111-9D34-0025B31E3CBE.root'
]

################
# Loading theglobal part for version
#
from HZZ2l2qAnalysis.Higgs2l2qCode.Higgs2l2qSkim_global53X_cfi import *

# For testing:

# To activate the debug of candidates:
# #process.p.replace(process.allhcand,process.allhcand*process.debugCandidates)

# # Checking PAT objects using Oscar's validation tools:
# process.validatePatElectrons = cms.EDAnalyzer('OGPatElectronValidationModule'
#                                               #,src=cms.InputTag('calibratedPatElectrons')
#                                               ,src=cms.InputTag('selectedPatElectrons')
#                                               ,printInfo=cms.untracked.bool(True)
#                                               )
# 
# process.p.replace(process.selectedPatElectrons,process.selectedPatElectrons*process.validatePatElectrons)
# 
# # Global checks using Oscar's validation tools:
# process.globaldump = cms.EDAnalyzer('OGDumpEventsModule'
#                                     ,vertexCollection=cms.untracked.InputTag('offlinePrimaryVertices')
#                                     ,electronCollection=cms.untracked.InputTag('userDataSelectedElectrons')
#                                     ,muonCollection=cms.untracked.InputTag('userDataSelectedMuons')
#                                     ,jetCollections=cms.untracked.vstring('cleanPatJetsNoPUIsoLept'
#                                                                           ,'cleanCA8JetsNoPUIsoLept'
#                                                                           ,'patJetsCA8CHSprunedSubjetsOrig'
#                                                                           )
#                                     ,metProducer=cms.untracked.InputTag("metInfoProducer")
#                                     )
# 
# process.metInfoProducer = cms.EDProducer('MetVariablesProducer'
#                                           ,metTag = cms.InputTag("patMETsPFJetsAK5")
#                                           ,t1CorrMetTag = cms.InputTag("dummy")
#                                           )
# 
# process.p.replace(process.zee,process.metInfoProducer*process.globaldump*process.zee)
# 
# process.ciematMatching = cms.EDProducer('CiematJetVertexMatchingProducer'
#                                         ,jetCollection=cms.string('customPFJetsNoPUSub')
#                                         ,vertexCollection=cms.InputTag('goodvertices')
#                                     )
# process.ogJet = cms.EDProducer('OGPatJetStdVarsModule'
#                                ,src=cms.InputTag('customPFJetsNoPUSub')
#                                ,vertexJetMatching=cms.untracked.InputTag('ciematMatching')
#                                #,vertexCollection=cms.untracked.InputTag('goodvertices')
#                                )
# 
# process.cleanPatJetsNoPUIsoLept.src = cms.InputTag('ogJet')
# 
# process.goodvertices = cms.EDFilter('GoodVertexSelector',src=cms.InputTag('offlinePrimaryVertices'))
# 
# process.p.replace(process.customPFJetsNoPUSub,
#                   process.customPFJetsNoPUSub*process.goodvertices*process.ciematMatching*process.ogJet)

#########################################################################
