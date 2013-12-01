import HZZ2l2qAnalysis.Higgs2l2qCode.Hzz2l2qSetup_cfi as Hzz2l2qSetup

import os

## import skeleton process
from PhysicsTools.PatAlgos.patTemplate_cfg import *

#add the L2L3Residual corrections only for data
if Hzz2l2qSetup.runOnMC:#MC
    jetCorrections=['L1FastJet','L2Relative','L3Absolute']

    # Setting the global tag for MC
    process.GlobalTag.globaltag = 'START53_V27::All'

        
else:#Data
    jetCorrections=['L1FastJet','L2Relative','L3Absolute','L2L3Residual']

    # Setting the global tag for data
    #     + For 2012A/B/C/D Winter13 re-reco (= 22Jan2013 re-reco)
    process.GlobalTag.globaltag = 'FT_53_V21_AN5::All'


############ general options ####################
process.options.wantSummary = True
process.maxEvents.input = Hzz2l2qSetup.usingevents
process.MessageLogger.cerr.FwkReport.reportEvery = 100

#For 53x MC, the default Jet Probability Calibration from the
#GlobalTag is not optimal and needs to be replaced in the following way,
#when using CRAB: 
if Hzz2l2qSetup.runOnMC:#MC
    process.GlobalTag.toGet = cms.VPSet(
        cms.PSet(record = cms.string("BTagTrackProbability2DRcd"),
                 tag = cms.string("TrackProbabilityCalibration_2D_MC53X_v2"),
                 connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU")),
        cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
                 tag = cms.string("TrackProbabilityCalibration_3D_MC53X_v2"),
                 connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU"))
        )

# Needed for the smearing in the electron energy correction
# 
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                   #theSource = cms.PSet(
                                                   #    initialSeed = cms.untracked.uint32(123456789),
                                                   #   engineName = cms.untracked.string('HepJamesRandom')
                                                   #    )
                                                   calibratedPatElectrons = cms.PSet(
                                                       initialSeed = cms.untracked.uint32(975312468),
                                                       engineName = cms.untracked.string('TRandom3')
                                                       )
                                                   )
    
############ PRINTOUT ###################

fileWeight = "MMozer/powhegweight/data/mZZ_Higgs"+str(Hzz2l2qSetup.higgsMass)+"_"+Hzz2l2qSetup.lhcEnergy+"_Lineshape+Interference.txt"

sep_line = "-" * 50
print sep_line
print 'running the following PFBRECO+PAT sequences:'
print '\tAK5'
print 'run on MC           : ', Hzz2l2qSetup.runOnMC
print 'run on signal       : ', (Hzz2l2qSetup.higgsMass>0)
print sep_line
print 'Number of events    : ', Hzz2l2qSetup.usingevents
print sep_line
print 'Global tag          : ', process.GlobalTag.globaltag
print 'LHC Energy (signal) : ', Hzz2l2qSetup.lhcEnergy
print sep_line
if (Hzz2l2qSetup.higgsMass>0):
    print 'Processing signal for Higgs Mass: ', Hzz2l2qSetup.higgsMass
    print 'Lineshape Rew. using ', fileWeight
    print sep_line

######################################################
    
### INPUT COLLECTIONS ##########

process.source.fileNames = Hzz2l2qSetup.usingfiles

### DEFINITION OF THE PFBRECO+PAT SEQUENCES ##########
# load the PAT config
process.load("PhysicsTools.PatAlgos.patSequences_cff")
from PhysicsTools.PatAlgos.tools.coreTools import *

from PhysicsTools.PatAlgos.tools.pfTools import *

# ---------------- Sequence AK5 ----------------------
postfixAK5 ="PFJetsAK5"
jetAlgoAK5 ="AK5"
useType1 = True

# Configure PAT to use PFBRECO instead of AOD sources
# this function will modify the PAT sequences.
usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgoAK5, runOnMC=Hzz2l2qSetup.runOnMC, postfix=postfixAK5,
          jetCorrections=('AK5PF', jetCorrections), typeIMetCorrections=useType1)

### DO NOT APPLY PFnoPU ###
getattr(process,"pfNoPileUp"+postfixAK5).enable = False
### setup additional projections ###
getattr(process,"pfNoMuon"+postfixAK5).enable = False 
getattr(process,"pfNoElectron"+postfixAK5).enable = False 
getattr(process,"pfNoTau"+postfixAK5).enable = False 
getattr(process,"pfNoJet"+postfixAK5).enable = True

# removing default cuts on muons 	 
getattr(process,"pfMuonsFromVertex"+postfixAK5).dzCut = 99 	 
getattr(process,"pfMuonsFromVertex"+postfixAK5).d0Cut = 99 	 
getattr(process,"pfSelectedMuons"+postfixAK5).cut="pt()>3" 	 
getattr(process,"pfIsolatedMuons"+postfixAK5).isolationCut = 999999 	 

# removing default cuts on electrons 	 
getattr(process,"pfElectronsFromVertex"+postfixAK5).dzCut = 99 	 
getattr(process,"pfElectronsFromVertex"+postfixAK5).d0Cut = 99 	 
getattr(process,"pfSelectedElectrons"+postfixAK5).cut="pt()>5" 	 
getattr(process,"pfIsolatedElectrons"+postfixAK5).isolationCut = 999999 	 
# remove pfTau and pfPhoton from the sequence
getattr(process,"PFBRECO"+postfixAK5).remove( getattr(process,"pfTauSequence"+postfixAK5) )
getattr(process,"PFBRECO"+postfixAK5).remove( getattr(process,"pfNoTau"+postfixAK5) )
getattr(process,"PFBRECO"+postfixAK5).remove( getattr(process,"pfPhotonSequence"+postfixAK5) )

# make sure about patJets input
# There was a bug in the old version... incredible:
# # switchToPFJets(process, input=cms.InputTag('pfJetsAK5'), algo=jetAlgoAK5, postfix = postfixAK5, jetCorrections=('AK5PF', jetCorrections))
switchToPFJets(process, input=cms.InputTag('pfJetsPFJetsAK5'), algo=jetAlgoAK5, postfix = postfixAK5, jetCorrections=('AK5PF', jetCorrections), type1=useType1)
# It seems that it is completely irrelevant to set the type1 here or not... when it is set above (in usePF2PAT).

if useType1:
    process.patMETsPFJetsAK5.metSource = cms.InputTag("patType1CorrectedPFMetPFJetsAK5")
# In case we disable just the one in usePF2PAT   process.patMETsPFJetsAK5.metSource = cms.InputTag("pfMETPFJetsAK5")



### we use "classic" muons and electrons (see below)
removeSpecificPATObjects(process, ['Taus'], postfix = postfixAK5)
removeSpecificPATObjects(process, ['Electrons'], postfix = postfixAK5)
removeSpecificPATObjects(process, ['Muons'], postfix = postfixAK5)
removeSpecificPATObjects(process, ['Photons'], postfix = postfixAK5)

############### remove useless modules ####################
def removeUseless( modName ):
    getattr(process,"patDefaultSequence"+postfixAK5).remove(
        getattr(process, modName+postfixAK5)
        )

removeUseless( "produceCaloMETCorrections" )
#removeUseless( "pfCandsNotInJet" )
#removeUseless( "pfJetMETcorr" )
#removeUseless( "pfCandMETcorr" )
#removeUseless( "pfchsMETcorr" )
#removeUseless( "pfType1CorrectedMet" )
#removeUseless( "pfType1p2CorrectedMet" )
removeUseless( "electronMatch" )
removeUseless( "muonMatch" )
removeUseless( "patPFTauIsolation" )
removeUseless( "tauMatch" )
removeUseless( "tauGenJets" )
removeUseless( "tauGenJetsSelectorAllHadrons" )
removeUseless( "tauGenJetMatch" )
removeUseless( "patHPSPFTauDiscrimination" )


#########################################################

# curing a weird bug in PAT..
from CMGTools.Common.PAT.removePhotonMatching import removePhotonMatching
removePhotonMatching( process, postfixAK5 )

# #OLD ########## insert the PFMET significance calculation #############
# #OLD 
# #OLD process.load("CMGTools.Common.PAT.MetSignificance_cff")
# #OLD 
# #OLD setattr(process,"PFMETSignificance"+postfixAK5, process.pfMetSignificance.clone())
# #OLD getattr(process,"patDefaultSequence"+postfixAK5).insert(getattr(process,"patDefaultSequence"+postfixAK5).index(getattr(process,"patMETs"+postfixAK5)),getattr(process,"PFMETSignificance"+postfixAK5))
# #OLD 
# #OLD # Trying to add stuff on the PAT for Type1 corrections:
# #OLD # #process.load("CMGTools.Common.PAT.PATMet_cff")

####################################################################

########## add specific configuration for pat Jets ##############

getattr(process,"patJets"+postfixAK5).addTagInfos = True
getattr(process,"patJets"+postfixAK5).tagInfoSources  = cms.VInputTag(
    cms.InputTag("secondaryVertexTagInfosAOD"+postfixAK5),
    cms.InputTag("impactParameterTagInfosAOD"+postfixAK5)
    )
### non default embedding of AOD items for default patJets
getattr(process,"patJets"+postfixAK5).embedCaloTowers = False
getattr(process,"patJets"+postfixAK5).embedPFCandidates = True

### disable MC matching (will be done at analysis level)
getattr(process,"patJets"+postfixAK5).addGenPartonMatch = False

if not(Hzz2l2qSetup.runOnMC): # Data (in MC is needed for smearing)
    getattr(process,"patJets"+postfixAK5).addGenJetMatch = False

# Processing the jet configuration

process.load('HZZ2l2qAnalysis.Higgs2l2qCode.HZZ2l2qJetConfiguration_cfi')
#import HZZ2l2qAnalysis.Higgs2l2qCode.HZZ2l2qJetConfiguration_cfi as process

# Seleting jets: safe cut on 25 since it is before smearing and to be
# safe for energy scale studies. DO USE 30 IN ANALYSIS
getattr(process,"selectedPatJets"+postfixAK5).cut = cms.string('pt > 25.0')

############## "Classic" PAT Muons and Electrons ########################

#process.load('HZZ2l2qAnalysis.Higgs2l2qCode.HZZ2l2qMuonConfiguration_cfi')
from HZZ2l2qAnalysis.Higgs2l2qCode.HZZ2l2qMuonConfiguration_cfi import setupPatMuons
setupPatMuons(process)

from HZZ2l2qAnalysis.Higgs2l2qCode.HZZ2l2qElectronConfiguration_cfi import setupPatElectrons
setupPatElectrons(process)

#if not Hzz2l2qSetup.runOnMC:
process.stdMuonSeq.remove( process.muonMatch )
process.stdElectronSeq.remove( process.electronMatch )
process.patMuons.addGenMatch = False
process.patElectrons.addGenMatch = False
process.patMuons.embedGenMatch = False
process.patElectrons.embedGenMatch = False

# #OLDprocess.userDataStandardLeptonSequence = cms.Sequence(
# #OLD    process.stdMuonSeq *
# #OLD    process.stdElectronSeq *
# #OLD
# #OLD    process.userDataSelectedMuons *
# #OLD    process.userDataSelectedElectrons *
# #OLD    process.selectedIDMuons *
# #OLD    process.selectedIDElectrons *
# #OLD    process.selectedIsoMuons *
# #OLD    process.selectedIsoElectrons 
# #OLD    )

### PATH DEFINITION #############################################

# counters that can be used at analysis level to know the processed events
process.prePathCounter = cms.EDProducer("EventCountProducer")
process.postPathCounter = cms.EDProducer("EventCountProducer")

# trigger information (no selection)

process.p = cms.Path( process.prePathCounter )

# PFBRECO+PAT ---

# used by metsignificance module
process.p += process.pfNoJet

process.p += getattr(process,"patPF2PATSequence"+postfixAK5)

process.p += process.customPFJetsNoPUSub
# #OLD process.p += process.userDataStandardLeptonSequence
process.p += process.stdMuonSeq
process.p += process.stdElectronSeq
process.p += process.cleanPatJetsNoPUIsoLept

process.p += process.jetSubstructuresSequence

# #process.p.replace(process.selectedPatJetsAK5,process.selectedPatJetsAK5*process.puJetIdSqeuence)
process.p.replace(getattr(process,"selectedPatJets"+postfixAK5),getattr(process,"selectedPatJets"+postfixAK5)*process.puJetIdSqeuence)

#process.p.replace(getattr(process,"patMETs"+postfixAK5),
#                  getattr(process,"pfJetMETcorr"+postfixAK5)*getattr(process,"pfType1CorrectedMet"+postfixAK5)*getattr(process,"patMETs"+postfixAK5))

########## insert the PFMET significance calculation #############

process.load("CMGTools.Common.PAT.MetSignificance_cff")

setattr(process,"PFMETSignificance"+postfixAK5, process.pfMetSignificance.clone())
#getattr(process,"patDefaultSequence"+postfixAK5).insert(getattr(process,"patDefaultSequence"+postfixAK5).index(getattr(process,"patMETs"+postfixAK5))+1,getattr(process,"PFMETSignificance"+postfixAK5))
process.p.replace(getattr(process,"patMETs"+postfixAK5),
                  getattr(process,"patMETs"+postfixAK5)*getattr(process,"PFMETSignificance"+postfixAK5))



# Combinatorial process:

process.load('HZZ2l2qAnalysis.Higgs2l2qCode.HZZ2l2qCandidateCreation_cfi')

# Oscar: a module for checks
process.debugCandidates = cms.EDAnalyzer("Higgs2l2qDebugModule",
                                         collections = cms.vstring("hzzeejjBaseColl","hzzmmjjBaseColl","hzzemjjBaseColl",
                                                                   "hzzee1jBaseColl","hzzmm1jBaseColl","hzzem1jBaseColl")
                                         )

# Oscar: New combinatorial sequence that is much more efficient:

process.combinatorialSequence = cms.Sequence(
#        process.jetFilterNoPUSub *
#        
        process.zee *
        process.zmm *
        process.zem *
        process.zll *
        process.zllFilter *
        
        process.zjj *
        
        process.hzzeejjBaseColl *
        process.hzzmmjjBaseColl *
        process.hzzemjjBaseColl *

        process.hzzee1jBaseColl *
        process.hzzmm1jBaseColl *
        process.hzzem1jBaseColl *

        process.allhcand *
        # #OLD... not here process.debugCandidates *
        process.candFilter *

        process.hzzeejj *
        process.hzzmmjj *
        process.hzzemjj *

        process.hzzee1j *
        process.hzzmm1j *
        process.hzzem1j
        )

process.p += process.combinatorialSequence
process.p += getattr(process,"postPathCounter")

# We want to have some checking on the total-kinematics check
# performed on MC (for checking whther the sample needs the filter)

if Hzz2l2qSetup.runOnMC:
    process.ciematKinematicsTest = cms.EDProducer("CiematKinematicsTest",
                                                  src = cms.InputTag("genParticles")
                                                  )
    process.p += process.ciematKinematicsTest

# event cleaning (in tagging mode, no event rejected)
process.load('CMGTools.Common.PAT.addFilterPaths_cff')

# The clean-up are run as indepedent paths (and bits stored
# as trigger results in PAT (*_TriggerResults_*_PAT)
process.fullPath = cms.Schedule(
    process.p
    ,process.EcalDeadCellTriggerPrimitiveFilterPath
    ,process.hcalLaserEventFilterPath
    ,process.trackingFailureFilterPath
    ,process.CSCTightHaloFilterPath
    ,process.HBHENoiseFilterPath
    ,process.eeBadScFilterPath
    ,process.primaryVertexFilterPath
    ,process.noscrapingFilterPath
    ,process.ecalLaserFilterPath
    ,process.trkPOGFiltersPath
    )
# Following are not to be used:
#    ,process.metNoiseCleaningPath
#    ,process.trackIsolationMakerFilterPath
#   )
del process.metNoiseCleaningPath
del process.trackIsolationMakerFilterPath

if Hzz2l2qSetup.applyCleanUpFilters:
    # The clean-up filters are applied as part of the main path
# #OLD   del process.EcalDeadCellTriggerPrimitiveFilterPath
   process.p += process.EcalDeadCellTriggerPrimitiveFilter
# #OLD   del process.hcalLaserEventFilterPath
   process.p += process.hcalLaserEventFilter
# #OLD   del process.trackingFailureFilterPath
   process.p += process.trackingFailureSequence
# #OLD   del process.CSCTightHaloFilterPath
   process.p += process.CSCTightHaloFilter
# #OLD   del process.HBHENoiseFilterPath
   process.p += process.HBHENoiseFilter
# #OLD   del process.eeBadScFilterPath
   process.p += process.eeBadScFilter
# #OLD   del process.primaryVertexFilterPath
   process.p += process.primaryVertexFilter
# #OLD   del process.noscrapingFilterPath
   process.p += process.noscraping
# #OLD   del process.ecalLaserFilterPath
   process.p += process.ecalLaserCorrFilter
# #OLD   del process.trkPOGFiltersPath
   process.p += process.trkPOGFiltersSequence

# #OLD   del process.totalKinematicsFilterPath
# #OLD   del process.metNoiseCleaningPath
# #OLD   del process.trackIsolationMakerFilterPath

# This is needed only for Madgraph MC... but I think it is harmless
# for the rest.
# Please note that we run it always in path mode
if Hzz2l2qSetup.runOnMC:
    process.fullPath.append(process.totalKinematicsFilterPath)
    #process.p += process.totalKinematicsFilter
    # The default value does not work because it kills the massless problem
    # in Madgraph
    process.totalKinematicsFilter.tolerance = cms.double(5)

else:
    del process.totalKinematicsFilterPath

# For signal:

if (Hzz2l2qSetup.higgsMass>0):
    from MMozer.powhegweight.tongguang600 import *
    process.powWeightProducer = tongguangweights600.clone(
        filename = cms.FileInPath(fileWeight)
        )

    process.p += process.powWeightProducer
    
### OUTPUT DEFINITION #############################################

# PFBRECO+PAT ---

# Add PFBRECO output to the created file
#from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning, patTriggerEventContent, patTriggerStandAloneEventContent

process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string(os.getenv('OGCMS_ROOTFILE','h2l2qSkimData.root')),
                               SelectEvents = cms.untracked.PSet(
                                   SelectEvents = cms.vstring('p')
                                   ),
                               outputCommands =  cms.untracked.vstring(
                                   'drop *_*_*_*',
                                   )
                               )

process.out.dropMetaData = cms.untracked.string("DROPPED")

process.out.outputCommands.extend([
    'keep *_selectedIDMuons_*_PAT',
    'keep *_selectedPatElectrons_*_PAT',
#    'keep *_customPFJets_*_PAT',
    'keep *_selectedPatJets*_pfCandidates_PAT',
    'keep *_cleanPatJetsNoPUIsoLept_*_PAT',
    # rho variables
    'keep *_*_rho_PAT',
    'keep *_kt6PFJetsCentralNeutral_rho_*',
    'keep *_kt6PFJets*_rho_*',
    # PU jetID maps
    "keep *_puJetId*_*_*", # input variables
    "keep *_puJetMva*_*_*", # final MVAs and working point flags
    # ll, jj, lljj candidates
    'keep *_zee_*_PAT',
    'keep *_zmm_*_PAT',
    'keep *_zem_*_PAT',
    'keep *_zjj_*_PAT',
    'keep *_hzzeejj_*_PAT',
    'keep *_hzzmmjj_*_PAT',
    'keep *_hzzemjj_*_PAT',
    'keep *_hzzee1j_*_PAT',
    'keep *_hzzmm1j_*_PAT',
    'keep *_hzzem1j_*_PAT',
    ####
    'keep *_offlineBeamSpot_*_*',
    'keep *_offlinePrimaryVertices_*_*',
#    'keep *_secondaryVertexTagInfos*_*_*',
#    'keep *_impactParameterTagInfos*_*_*',
#    'keep *_*_*tagInfo*_*',
    # additional collections from AOD   
#    'keep *_generalTracks_*_*',
#    'keep *_electronGsfTracks_*_*',
#    'keep *_muons_*_*',
#    'keep *_globalMuons_*_*',
#    'keep *_standAloneMuons_*_*',
    'keep recoPFCandidates_particleFlow_*_*',
    # genParticles & genJets
    'keep *_genParticles_*_*',
    'keep recoGenJets_ak5GenJets_*_*',
#    'keep recoGenJets_kt6GenJets_*_*',
#    'keep *_powWeightProducer_*_*',
    # gen Info
    'keep PileupSummaryInfos_*_*_*',
    'keep GenEventInfoProduct_*_*_*',
    'keep GenRunInfoProduct_*_*_*',
    'keep LHEEventProduct_*_*_*',
    'keep *_genEventScale_*_*',
    # Signal weights for powheg
    'keep *_powWeightProducer_*_*',

    ###### MET products
    'keep *_patMETs*_*_*',
#    'keep *_patType1CorrectedPFMet_*_*', # NOT included for the moment
    'keep *_PFMETSignificance*_*_*',
    ### for HLT selection
    'keep edmTriggerResults_TriggerResults_*_HLT'])

# additional products for event cleaning
process.out.outputCommands.extend(['keep *_TriggerResults_*_PAT'])
process.out.outputCommands.extend(['keep edmMergeableCounter_*_*_*'])
process.out.outputCommands.extend(['keep *_ciematKinematicsTest_*_*'])


# Adding the CA8 information
process.out.outputCommands.extend(['keep *_ca8PFJetsCHSpruned_SubJets_*'])
process.out.outputCommands.extend(['keep *_ca8GenJetsNoNu_*_*'])
process.out.outputCommands.extend(['keep *_cleanCA8JetsNoPUIsoLept_*_*'])

process.out.outputCommands.extend(['keep *_btaggedPatJetsCA8CHSpruned_*_*'])
process.out.outputCommands.extend(['keep *_patJetsCA8CHSprunedSubjets_*_*'])

#
process.outPath = cms.EndPath(process.out)
#
