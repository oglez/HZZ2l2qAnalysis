
 ## import skeleton process
from PhysicsTools.PatAlgos.patTemplate_cfg import *

# turn on when running on MC
runOnMC = True

# turn on when running on powheg signal MC (-> to produce line-shape weights)
#isPowhegSignal = False

#add the L2L3Residual corrections only for data
if runOnMC:#MC
    jetCorrections=['L1FastJet','L2Relative','L3Absolute']
else:#Data
    jetCorrections=['L1FastJet','L2Relative','L3Absolute','L2L3Residual']

############ general options ####################
process.options.wantSummary = True
process.maxEvents.input = 10000
process.MessageLogger.cerr.FwkReport.reportEvery = 100
########### global tag ############################
#from CMGTools.Common.Tools.getGlobalTag import getGlobalTag
#process.GlobalTag.globaltag = cms.string(getGlobalTag(runOnMC))
process.GlobalTag.globaltag = 'START53_V27::All'
##################################################

#For 53x MC, the default Jet Probability Calibration from the
#GlobalTag is not optimal and needs to be replaced in the following way,
#when using CRAB: 

process.GlobalTag.toGet = cms.VPSet(
  cms.PSet(record = cms.string("BTagTrackProbability2DRcd"),
       tag = cms.string("TrackProbabilityCalibration_2D_MC53X_v2"),
       connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU")),
  cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
       tag = cms.string("TrackProbabilityCalibration_3D_MC53X_v2"),
       connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU"))
)

############ PRINTOUT ###################

sep_line = "-" * 50
print sep_line
print 'running the following PFBRECO+PAT sequences:'
print '\tAK5'
print 'run on MC        : ', runOnMC
print sep_line
print 'Global tag       : ', process.GlobalTag.globaltag
print sep_line

######################################################
    
### INPUT COLLECTIONS ##########

process.source.fileNames = [
#    'file:/data3/scratch/cms/mc/Summer12_DR53X/DYJetsToLL_M-50/FEF4E41A-46D4-E111-9594-0025B3E06424.root'
#    'file:/data3/scratch/cms/mc/Summer12/GluGluToHToZZTo2L2Q_M-700_8TeV/FE073BF1-CDB7-E111-8866-E0CB4E4408BE.root'
    'file:/data1/oglez/spool/Zjets_MADGRAPH_v1_AOD532/FE414F4B-F6D2-E111-A4E9-003048674048.root'
]

### DEFINITION OF THE PFBRECO+PAT SEQUENCES ##########
# load the PAT config
process.load("PhysicsTools.PatAlgos.patSequences_cff")
from PhysicsTools.PatAlgos.tools.coreTools import *

from PhysicsTools.PatAlgos.tools.pfTools import *

# ---------------- Sequence AK5 ----------------------
postfixAK5 ="AK5"
jetAlgoAK5 ="AK5"

# Configure PAT to use PFBRECO instead of AOD sources
# this function will modify the PAT sequences.
usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgoAK5, runOnMC=runOnMC, postfix=postfixAK5,
          jetCorrections=('AK5PF', jetCorrections))

### DO NOT APPLY PFnoPU ###
getattr(process,"pfNoPileUp"+postfixAK5).enable = False
### setup additional projections ###
getattr(process,"pfNoMuon"+postfixAK5).enable = False 
getattr(process,"pfNoElectron"+postfixAK5).enable = False 
getattr(process,"pfNoTau"+postfixAK5).enable = False 
getattr(process,"pfNoJet"+postfixAK5).enable = True

# removing default cuts on muons 	 
getattr(process,"pfMuonsFromVertexAK5").dzCut = 99 	 
getattr(process,"pfMuonsFromVertexAK5").d0Cut = 99 	 
getattr(process,"pfSelectedMuonsAK5").cut="pt()>3" 	 
getattr(process,"pfIsolatedMuons"+postfixAK5).isolationCut = 999999 	 

# removing default cuts on electrons 	 
getattr(process,"pfElectronsFromVertexAK5").dzCut = 99 	 
getattr(process,"pfElectronsFromVertexAK5").d0Cut = 99 	 
getattr(process,"pfSelectedElectronsAK5").cut="pt()>5" 	 
getattr(process,"pfIsolatedElectrons"+postfixAK5).isolationCut = 999999 	 
# remove pfTau and pfPhoton from the sequence
process.PFBRECOAK5.remove( process.pfTauSequenceAK5 )
process.PFBRECOAK5.remove( process.pfNoTauAK5 )
process.PFBRECOAK5.remove( process.pfPhotonSequenceAK5 )

# make sure about patJets input
switchToPFJets(process, input=cms.InputTag('pfJetsAK5'), algo=jetAlgoAK5, postfix = postfixAK5, jetCorrections=('AK5PF', jetCorrections))

### we use "classic" muons and electrons (see below)
removeSpecificPATObjects(process, ['Taus'], postfix = "AK5")
removeSpecificPATObjects(process, ['Electrons'], postfix = "AK5")
removeSpecificPATObjects(process, ['Muons'], postfix = "AK5")
removeSpecificPATObjects(process, ['Photons'], postfix = "AK5")

############### remove useless modules ####################
def removeUseless( modName ):
    getattr(process,"patDefaultSequence"+postfixAK5).remove(
        getattr(process, modName+postfixAK5)
        )

removeUseless( "produceCaloMETCorrections" )
removeUseless( "pfCandsNotInJet" )
removeUseless( "pfJetMETcorr" )
removeUseless( "pfCandMETcorr" )
removeUseless( "pfchsMETcorr" )
removeUseless( "pfType1CorrectedMet" )
removeUseless( "pfType1p2CorrectedMet" )
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

########## insert the PFMET significance calculation #############

process.load("CMGTools.Common.PAT.MetSignificance_cff")

setattr(process,"PFMETSignificance"+postfixAK5, process.pfMetSignificance.clone())
getattr(process,"patDefaultSequence"+postfixAK5).insert(getattr(process,"patDefaultSequence"+postfixAK5).index(getattr(process,"patMETs"+postfixAK5)),getattr(process,"PFMETSignificance"+postfixAK5))

####################################################################

########## add specific configuration for pat Jets ##############

getattr(process,"patJets"+postfixAK5).addTagInfos = True
getattr(process,"patJets"+postfixAK5).tagInfoSources  = cms.VInputTag(
    cms.InputTag("secondaryVertexTagInfosAODAK5"),
    cms.InputTag("impactParameterTagInfosAODAK5")
    )
### non default embedding of AOD items for default patJets
getattr(process,"patJets"+postfixAK5).embedCaloTowers = False
getattr(process,"patJets"+postfixAK5).embedPFCandidates = True

### disable MC matching (will be done at analysis level)
getattr(process,"patJets"+postfixAK5).addGenPartonMatch = False
getattr(process,"patJets"+postfixAK5).addGenJetMatch = False


##############################################################
#add user variables to PAT-jets 
#process.qglAK5PF   = cms.EDProducer(
#    "QuarkGluonTagger",
#    jets     = cms.InputTag("selectedPatJetsAK5"),
#    rho      = cms.InputTag("kt6PFJetsForIso:rho"),
#    jec      = cms.string('ak5PFL1FastL2L3'),
#    isPatJet = cms.bool(True),
#    )

process.customPFJetsNoPUSub = cms.EDProducer(
    'PFJetUserData',
    JetInputCollection=cms.untracked.InputTag("selectedPatJetsAK5"),
#    is2012Data=cms.untracked.bool(True),
#    qgMap=cms.untracked.InputTag("qglAK5PF"),
    Verbosity=cms.untracked.bool(False)
    )

#from  CMGTools.External.pujetidsequence_cff import puJetId, puJetMva
#process.puJetIdAK5 = puJetId.clone( jets = 'customPFJetsNoPUSub')
#process.puJetMvaAK5= puJetMva.clone(
#    jetids = cms.InputTag("puJetIdAK5"),
#    jets = 'customPFJetsNoPUSub',
#    )
#
#process.puJetIdSequenceAK5 = cms.Sequence(process.puJetIdAK5*process.puJetMvaAK5)

# central jets for filtering and Z->jj candidates

process.customPFJetsNoPUSubCentral = cms.EDFilter(
    "PATJetSelector",
    src = cms.InputTag("customPFJetsNoPUSub"),
    cut = cms.string("abs(eta) < 2.4")
    )

############## "Classic" PAT Muons and Electrons ########################
# (made from all reco muons, and all gsf electrons, respectively)
process.patMuons.embedTcMETMuonCorrs = False
process.patMuons.embedCaloMETMuonCorrs = False
process.patMuons.embedTrack = True

process.patElectrons.embedTrack = True
process.patElectrons.pfElectronSource = 'particleFlow'

# use PFIsolation
process.eleIsoSequence = setupPFElectronIso(process, 'gsfElectrons', 'PFIso')
process.muIsoSequence = setupPFMuonIso(process, 'muons', 'PFIso')
adaptPFIsoMuons( process, applyPostfix(process,"patMuons",""), 'PFIso')
adaptPFIsoElectrons( process, applyPostfix(process,"patElectrons",""), 'PFIso')


# setup recommended 0.3 cone for electron PF isolation
#process.pfIsolatedElectrons.isolationValueMapsCharged = cms.VInputTag(cms.InputTag("elPFIsoValueCharged03PFIdPFIso"))
#process.pfIsolatedElectrons.deltaBetaIsolationValueMap = cms.InputTag("elPFIsoValuePU03PFIdPFIso")
#process.pfIsolatedElectrons.isolationValueMapsNeutral = cms.VInputTag(cms.InputTag("elPFIsoValueNeutral03PFIdPFIso"), cms.InputTag("elPFIsoValueGamma03PFIdPFIso"))
process.patElectrons.isolationValues = cms.PSet(
    pfChargedHadrons = cms.InputTag("elPFIsoValueCharged03PFIdPFIso"),
    pfChargedAll = cms.InputTag("elPFIsoValueChargedAll03PFIdPFIso"),
    pfPUChargedHadrons = cms.InputTag("elPFIsoValuePU03PFIdPFIso"),
    pfNeutralHadrons = cms.InputTag("elPFIsoValueNeutral03PFIdPFIso"),
    pfPhotons = cms.InputTag("elPFIsoValueGamma03PFIdPFIso")
    )
process.patElectrons.isolationValuesNoPFId = cms.PSet(
    pfChargedHadrons = cms.InputTag("elPFIsoValueCharged03NoPFIdPFIso"),
    pfChargedAll = cms.InputTag("elPFIsoValueChargedAll03NoPFIdPFIso"),
    pfPUChargedHadrons = cms.InputTag("elPFIsoValuePU03NoPFIdPFIso"),
    pfNeutralHadrons = cms.InputTag("elPFIsoValueNeutral03NoPFIdPFIso"),
    pfPhotons = cms.InputTag("elPFIsoValueGamma03NoPFIdPFIso")
    )


process.stdMuonSeq = cms.Sequence(
    process.pfParticleSelectionSequence +
    process.muIsoSequence +
    process.makePatMuons +
    process.selectedPatMuons
    )
process.stdElectronSeq = cms.Sequence(
    process.pfParticleSelectionSequence +
    process.eleIsoSequence +
    process.makePatElectrons +
    process.selectedPatElectrons
    )

#if not runOnMC:
process.stdMuonSeq.remove( process.muonMatch )
process.stdElectronSeq.remove( process.electronMatch )
process.patMuons.addGenMatch = False
process.patElectrons.addGenMatch = False
process.patMuons.embedGenMatch = False
process.patElectrons.embedGenMatch = False

# Modules for Electron ID
# MVA Electron ID
process.load("EGamma.EGammaAnalysisTools.electronIdMVAProducer_cfi")
process.mvaeIdSequence = cms.Sequence(
    process.mvaTrigV0 +
    process.mvaNonTrigV0
)
# ElectronID in the VBTF prescription
#import ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi as vbtfid
# Switch to the official Electron VBTF Selection for 2011 Data (relax H/E cut in the Endcap):
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/VbtfEleID2011
import HZZ2l2qAnalysis.Higgs2l2qCode.simpleCutBasedElectronIDSummer11_cfi as vbtfid
process.eidVBTFRel95 = vbtfid.simpleCutBasedElectronID.clone( electronQuality = '95relIso' )
process.eidVBTFRel80 = vbtfid.simpleCutBasedElectronID.clone( electronQuality = '80relIso' )
process.eidVBTFCom95 = vbtfid.simpleCutBasedElectronID.clone( electronQuality = '95cIso'   )
process.eidVBTFCom80 = vbtfid.simpleCutBasedElectronID.clone( electronQuality = '80cIso'   )
        
process.vbtfeIdSequence = cms.Sequence(
        process.eidVBTFRel95 +
        process.eidVBTFRel80 +
        process.eidVBTFCom95 +
        process.eidVBTFCom80
)

process.eidSequence = cms.Sequence(
    process.mvaeIdSequence +
    process.vbtfeIdSequence
)

process.patElectrons.electronIDSources = cms.PSet(
        eidVBTFRel95 = cms.InputTag("eidVBTFRel95"),
        eidVBTFRel80 = cms.InputTag("eidVBTFRel80"),
        eidVBTFCom95 = cms.InputTag("eidVBTFCom95"),
        eidVBTFCom80 = cms.InputTag("eidVBTFCom80"),
        #MVA 
        mvaTrigV0 = cms.InputTag("mvaTrigV0"),
        mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0")
)

process.stdLeptonSequence = cms.Sequence(
    process.stdMuonSeq +
    process.eidSequence +
    process.stdElectronSeq 
    )

# Classic Electrons with UserData

process.userDataSelectedElectrons = cms.EDProducer(
    "Higgs2l2bElectronUserData",
    src = cms.InputTag("selectedPatElectrons"),
    rho = cms.InputTag("kt6PFJets:rho"),
    primaryVertices=cms.InputTag("offlinePrimaryVertices")
)

# ID+Isolated electrons: select electrons passing LOOSE
process.selectedIsoElectrons = cms.EDFilter(
    "PATElectronSelector",
    src = cms.InputTag("userDataSelectedElectrons"),
    cut = cms.string("(userFloat('cutIDCode') > 1) && (userFloat('passTriggerTight') > 0)")
#    src = cms.InputTag("selectedIDElectrons"),
#    cut = cms.string("electronID('eidVBTFCom95') == 7")
)

# Classic Muons with UserData
process.userDataSelectedMuons = cms.EDProducer(
    "Higgs2l2bMuonUserData",
    src = cms.InputTag("selectedPatMuons"),
    primaryVertices=cms.InputTag("offlinePrimaryVertices")
)

# ID-selected muons
process.selectedIDMuons = cms.EDFilter(
    "PATMuonSelector",
    src = cms.InputTag("userDataSelectedMuons"),
    cut = cms.string(
            "isGlobalMuon && isTrackerMuon && isPFMuon && globalTrack().normalizedChi2 < 10 && " +
            "globalTrack().hitPattern().numberOfValidMuonHits > 0 && "              +
            "numberOfMatchedStations > 1 && abs( userFloat('dzVtx') ) < 0.5 && "    +
            "innerTrack().hitPattern().numberOfValidPixelHits > 0 && "              +
            "track().hitPattern().trackerLayersWithMeasurement > 5 && abs(dB) < 0.2" )
)

# Isolated muons: standard isolation
process.selectedIsoMuons = cms.EDFilter(
    "PATMuonSelector",
    src = cms.InputTag("selectedIDMuons"),
#    cut = cms.string("trackIso + caloIso < 0.15 * pt")
# using DeltaBeta correction
    cut = cms.string("( max(0., (neutralHadronIso + photonIso - 0.5*puChargedHadronIso) ) + chargedHadronIso) < 0.12 * pt")
)

process.userDataStandardLeptonSequence = cms.Sequence(
    process.userDataSelectedMuons *
    process.userDataSelectedElectrons *
    process.selectedIDMuons *
#    process.selectedIDElectrons *
    process.selectedIsoMuons *
    process.selectedIsoElectrons 
    )


# Jet cleaning for patJets
process.cleanPatJetsNoPUIsoLept = cms.EDProducer(
    "PATJetCleaner",
    src = cms.InputTag("customPFJetsNoPUSubCentral"),
    preselection = cms.string(''),
    checkOverlaps = cms.PSet(
    muons = cms.PSet(
    src = cms.InputTag("selectedIsoMuons"),
    algorithm = cms.string("byDeltaR"),
    preselection = cms.string(""),
    deltaR = cms.double(0.5),
    checkRecoComponents = cms.bool(False),
    pairCut = cms.string(""),
    requireNoOverlaps = cms.bool(True),
    ),
    electrons = cms.PSet(
    src = cms.InputTag("selectedIsoElectrons"),
    algorithm = cms.string("byDeltaR"),
    preselection = cms.string(""),
    deltaR = cms.double(0.5),
    checkRecoComponents = cms.bool(False),
    pairCut = cms.string(""),
    requireNoOverlaps = cms.bool(True),
    )
    ),
    finalCut = cms.string('')
    )

#---------------- For Merged Jets, added by Oscar.

jetSubstructuresEventContent = cms.untracked.vstring()
jetSubstructuresSequence = cms.Sequence()

### Make correct pfNoPileUp ###########
from CommonTools.ParticleFlow.goodOfflinePrimaryVertices_cfi import goodOfflinePrimaryVertices
process.goodOfflinePrimaryVerticesForSubJets = goodOfflinePrimaryVertices.clone()

from CommonTools.ParticleFlow.pfPileUp_cfi  import *
from CommonTools.ParticleFlow.TopProjectors.pfNoPileUp_cfi import *

process.pfPileUpForSubJets = pfPileUp.clone(
        checkClosestZVertex = False,
            PFCandidates = 'particleFlow',
            Vertices = 'goodOfflinePrimaryVerticesForSubJets'
            )
process.pfNoPileUpForSubJets = pfNoPileUp.clone(
        topCollection = 'pfPileUpForSubJets',
            bottomCollection = 'particleFlow'
            )

jetSubstructuresSequence += process.goodOfflinePrimaryVerticesForSubJets
jetSubstructuresSequence += process.pfPileUpForSubJets
jetSubstructuresSequence += process.pfNoPileUpForSubJets

pfNoPileUpSrc = 'pfNoPileUpForSubJets'

#### Adding CA8 jets and CA8 pruned jets
process.load("ExoDiBosonResonances.PATtupleProduction.PAT_ca8jets_cff")
process.ca8PFJetsCHS.src = pfNoPileUpSrc
process.ca8PFJetsCHSpruned.src = pfNoPileUpSrc
jetSubstructuresSequence += process.PATCMGJetSequenceCA8CHS
jetSubstructuresSequence += process.PATCMGJetSequenceCA8CHSpruned
jetSubstructuresEventContent+=['keep *_ca8PFJetsCHSpruned_SubJets_*']
jetSubstructuresEventContent+=['keep *_ca8GenJetsNoNu_*_*']

## jetSubstructuresEventContent+=['keep patJets_selectedPatJetsCA8CHSwithNsub_*_*']

process.selectedPatJetsCA8CHS.cut = cms.string('pt > 25.0 && abs(eta) < 2.4')

jetSubstructuresSequence += process.PATCMGJetSequenceCA8CHS
jetSubstructuresSequence += process.PATCMGJetSequenceCA8CHSpruned
jetSubstructuresSequence += process.selectedPatJetsCA8CHSwithNsub
jetSubstructuresSequence += process.selectedPatJetsCA8CHSwithQjets
# OLD jetSubstructuresSequence += process.selectedPatJetsCA8CHSwithNsub

#### Jet cleaning for CA8 Jets
process.cleanCA8JetsNoPUIsoLept = cms.EDProducer(
        "PATJetCleaner",
        #src = cms.InputTag("selectedPatJetsCA8CHS"),
        src = cms.InputTag("selectedPatJetsCA8CHSwithQjets"),
        preselection = cms.string(''),
        checkOverlaps = cms.PSet(
                muons = cms.PSet(
                src = cms.InputTag("selectedIsoMuons"),
                algorithm = cms.string("byDeltaR"),
                preselection = cms.string(""),
                #deltaR = cms.double(0.5),
                deltaR = cms.double(0.7),
                checkRecoComponents = cms.bool(False),
                pairCut = cms.string(""),
                requireNoOverlaps = cms.bool(True),
                ),
                electrons = cms.PSet(
                src = cms.InputTag("selectedIsoElectrons"),
                algorithm = cms.string("byDeltaR"),
                preselection = cms.string(""),
                #deltaR = cms.double(0.5),
                deltaR = cms.double(0.7),
                checkRecoComponents = cms.bool(False),
                pairCut = cms.string(""),
                requireNoOverlaps = cms.bool(True),
                )
                ),
        finalCut = cms.string('')
        )

jetSubstructuresEventContent+=['keep *_cleanCA8JetsNoPUIsoLept_*_*']
jetSubstructuresSequence += process.cleanCA8JetsNoPUIsoLept

# ---------------- Common stuff ---------------

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

#process.p += process.qglAK5PF 
process.p += process.customPFJetsNoPUSub
#process.p += process.puJetIdSequenceAK5
process.p += process.customPFJetsNoPUSubCentral

process.p += process.stdLeptonSequence

process.p += process.userDataStandardLeptonSequence
process.p += process.cleanPatJetsNoPUIsoLept

# Select leptons
process.selectedPatMuons.cut = (
    "pt > 10 && abs(eta) < 2.4"
        )

process.selectedPatElectrons.cut = (
    "pt > 10.0 && abs(eta) < 2.5"
    )
# Select jets
process.selectedPatJetsAK5.cut = cms.string('pt > 25.0')

process.p += jetSubstructuresSequence


################# COMBINATORIAL ANALYSIS ###########################

process.zee = cms.EDProducer("CandViewShallowCloneCombiner",
                                 checkCharge = cms.bool(False),
                                 cut = cms.string('mass > 20 '),
                                 decay = cms.string("userDataSelectedElectrons@+ userDataSelectedElectrons@-")
                             )

process.zmm = cms.EDProducer("CandViewShallowCloneCombiner",
                                 checkCharge = cms.bool(False),
                                 cut = cms.string('mass > 20 '),
                                 decay = cms.string("userDataSelectedMuons@+ userDataSelectedMuons@-")
                             )

process.zem = cms.EDProducer("CandViewShallowCloneCombiner",
                                 checkCharge = cms.bool(False),
                                 cut = cms.string('mass > 20 '),
                                 decay = cms.string("userDataSelectedElectrons@+ userDataSelectedMuons@-")
                             )

process.zjj = cms.EDProducer("CandViewShallowCloneCombiner",
                             checkCharge = cms.bool(False),
                             checkOverlap = cms.bool(False),
                             cut = cms.string(''),
                             decay = cms.string("cleanPatJetsNoPUIsoLept cleanPatJetsNoPUIsoLept")
)

process.hzzeejjBaseColl = cms.EDProducer("CandViewCombiner",
                                             checkCharge = cms.bool(False),
                                             cut = cms.string(''),
                                             decay = cms.string("zee zjj")
                                         )

process.hzzmmjjBaseColl = cms.EDProducer("CandViewCombiner",
                                             checkCharge = cms.bool(False),
                                             cut = cms.string(''),
                                             decay = cms.string("zmm zjj")
                                         )

process.hzzemjjBaseColl = cms.EDProducer("CandViewCombiner",
                                             checkCharge = cms.bool(False),
                                             cut = cms.string(''),
                                             decay = cms.string("zem zjj")
                                         )

process.hzzeejj = cms.EDProducer("H2l2qCandidateData",
                                     higgs = cms.InputTag("hzzeejjBaseColl"),
                                     gensTag = cms.InputTag("genParticles"),
                                     PFCandidates = cms.InputTag("particleFlow"),
                                     primaryVertices = cms.InputTag("offlinePrimaryVertices"),
                                     dzCut = cms.double(0.1)
                                 )

process.hzzmmjj = cms.EDProducer("H2l2qCandidateData",
                                     higgs = cms.InputTag("hzzmmjjBaseColl"),
                                     gensTag = cms.InputTag("genParticles"),
                                     PFCandidates = cms.InputTag("particleFlow"),
                                     primaryVertices = cms.InputTag("offlinePrimaryVertices"),
                                     dzCut = cms.double(0.1)
                                     )

process.hzzemjj = cms.EDProducer("H2l2qCandidateData",
                                     higgs = cms.InputTag("hzzemjjBaseColl"),
                                     gensTag = cms.InputTag("genParticles"),
                                     PFCandidates = cms.InputTag("particleFlow"),
                                     primaryVertices = cms.InputTag("offlinePrimaryVertices"),
                                     dzCut = cms.double(0.1)
                                     )

# #OSCAR process.combinatorialSequence = cms.Sequence(
# #OSCAR     process.zee +
# #OSCAR     process.zmm +
# #OSCAR     process.zem +
# #OSCAR     process.zjj +
# #OSCAR     process.hzzeejjBaseColl +
# #OSCAR     process.hzzmmjjBaseColl +
# #OSCAR     process.hzzemjjBaseColl +
# #OSCAR     process.hzzeejj +
# #OSCAR     process.hzzmmjj +
# #OSCAR     process.hzzemjj
# #OSCAR )
# #OSCAR 
# #OSCAR process.p += process.combinatorialSequence

# #OSCAR #if runOnMC and isPowhegSignal:
# #OSCAR #    process.p += process.powWeightProducer
# #OSCAR 
# #OSCAR process.p += getattr(process,"postPathCounter") 

# We need to process the candidates for the merged jets

process.hzzee1jBaseColl = cms.EDProducer("CandViewCombiner",
                                         checkCharge = cms.bool(False),
                                         cut = cms.string(''),
                                         decay = cms.string("zee cleanCA8JetsNoPUIsoLept")
                                         )

process.hzzmm1jBaseColl = cms.EDProducer("CandViewCombiner",
                                         checkCharge = cms.bool(False),
                                         cut = cms.string(''),
                                         decay = cms.string("zmm cleanCA8JetsNoPUIsoLept")
                                         )
process.hzzem1jBaseColl = cms.EDProducer("CandViewCombiner",
                                         checkCharge = cms.bool(False),
                                         cut = cms.string(''),
                                         decay = cms.string("zem cleanCA8JetsNoPUIsoLept")
                                         )

process.hzzee1j = cms.EDProducer("H2l2qmergedCandidateData",
                                 higgs = cms.InputTag("hzzee1jBaseColl")
                                 #,gensTag = cms.InputTag("genParticles"),
                                 #,PFCandidates = cms.InputTag("particleFlow"),
                                 #,primaryVertices = cms.InputTag("offlinePrimaryVertices"),
                                 #,dzCut = cms.double(0.1)
                                 )

process.hzzmm1j = cms.EDProducer("H2l2qmergedCandidateData",
                                 higgs = cms.InputTag("hzzmm1jBaseColl")
                                 #,gensTag = cms.InputTag("genParticles"),
                                 #,PFCandidates = cms.InputTag("particleFlow"),
                                 #,primaryVertices = cms.InputTag("offlinePrimaryVertices"),
                                 #,dzCut = cms.double(0.1)
                                 )

process.hzzem1j = cms.EDProducer("H2l2qmergedCandidateData",
                                 higgs = cms.InputTag("hzzem1jBaseColl")
                                 #,gensTag = cms.InputTag("genParticles"),
                                 #,PFCandidates = cms.InputTag("particleFlow"),
                                 #,primaryVertices = cms.InputTag("offlinePrimaryVertices"),
                                 #,dzCut = cms.double(0.1)
                                 )

# Setup for a basic filtering
process.zll = cms.EDProducer("CandViewMerger",
                             src = cms.VInputTag("zee", "zmm", "zem")
)

process.zllFilter = cms.EDFilter("CandViewCountFilter",
                                 src = cms.InputTag("zll"),
                                 minNumber = cms.uint32(1)
                                 )

# Oscar modified the requirements to make them a bit cleaner: not a Z(ll)+2 jets but one
# Higgs candidate

process.allhcand = cms.EDProducer("CandViewMerger",
                                  src = cms.VInputTag("hzzeejjBaseColl","hzzmmjjBaseColl","hzzemjjBaseColl",
                                                      "hzzee1jBaseColl","hzzmm1jBaseColl","hzzem1jBaseColl")
                                  )

process.candFilter = cms.EDFilter("CandViewCountFilter",
                                  src = cms.InputTag("allhcand"),
                                  minNumber = cms.uint32(1)
                                  )

# process.jetFilterNoPUSub = cms.EDFilter("CandViewCountFilter",
#                                         src = cms.InputTag("customPFJetsNoPUSubCentral"),
#                                         minNumber = cms.uint32(2),
#                                         filter = cms.bool(True)
#                                         )

# #OSCAR process.filterPath1= cms.Path(
# #OSCAR     process.zll *
# #OSCAR     process.zllFilter *
# #OSCAR     #OSCAR process.jetFilterNoPUSub
# #OSCAR     process.allhcand *
# #OSCAR     process.candFilter
# #OSCAR     )

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

# event cleaning (in tagging mode, no event rejected)
process.load('CMGTools.Common.PAT.addFilterPaths_cff')

# process.fullPath = cms.Schedule(
#     process.p,
# # #OSCAR    process.filterPath1,
# #    process.filterPath2,
#     process.EcalDeadCellTriggerPrimitiveFilterPath,
#     process.hcalLaserEventFilterPath,
#     process.trackingFailureFilterPath,
#     process.CSCTightHaloFilterPath,
#     process.HBHENoiseFilterPath,
#     process.eeBadScFilterPath,
#     process.primaryVertexFilterPath,
#     process.noscrapingFilterPath,
#     process.ecalLaserFilterPath,
#     process.trkPOGFiltersPath
#     )

del process.EcalDeadCellTriggerPrimitiveFilterPath
process.p += process.EcalDeadCellTriggerPrimitiveFilter
del process.hcalLaserEventFilterPath
process.p += process.hcalLaserEventFilter
del process.trackingFailureFilterPath
process.p += process.trackingFailureSequence
del process.CSCTightHaloFilterPath
process.p += process.CSCTightHaloFilter
del process.HBHENoiseFilterPath
process.p += process.HBHENoiseFilter
del process.eeBadScFilterPath
process.p += process.eeBadScFilter
del process.primaryVertexFilterPath
process.p += process.primaryVertexFilter
del process.noscrapingFilterPath
process.p += process.noscraping
del process.ecalLaserFilterPath
process.p += process.ecalLaserCorrFilter
del process.trkPOGFiltersPath
process.p += process.trkPOGFiltersSequence

del process.totalKinematicsFilterPath
del process.metNoiseCleaningPath
del process.trackIsolationMakerFilterPath

#this is needed only for Madgraph MC:
if runOnMC:
    #process.fullPath.append(process.totalKinematicsFilterPath)
    process.p += process.totalKinematicsFilter
#else:
    #del process.totalKinematicsFilterPath

### OUTPUT DEFINITION #############################################

# PFBRECO+PAT ---

# Add PFBRECO output to the created file
#from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning, patTriggerEventContent, patTriggerStandAloneEventContent


process.out = cms.OutputModule("PoolOutputModule",
                 fileName = cms.untracked.string('h2l2qSkimData.root'),
                 SelectEvents = cms.untracked.PSet(
                     SelectEvents = cms.vstring('p')
# #OSCAR                         'filterPath1')#,
# #OSCAR #                        'filterPath2')
                     ),
                  outputCommands =  cms.untracked.vstring(
                  'drop *_*_*_*',
                  ),
)

process.out.dropMetaData = cms.untracked.string("DROPPED")

process.out.outputCommands.extend([
    'keep *_userDataSelectedElectrons_*_PAT',
    'keep *_userDataSelectedMuons_*_PAT',
#    'keep *_customPFJets_*_PAT',
    'keep *_selectedPatJetsAK5_pfCandidates_PAT',
    'keep *_customPFJetsNoPUSub_*_PAT',
    'keep *_cleanPatJetsNoPUIsoLept_*_PAT',
    # rho variables
    'keep *_*_rho_PAT',
    'keep *_kt6PFJetsCentralNeutral_rho_*',
    'keep *_kt6PFJets*_rho_*',
    # PU jetID maps
#    "keep *_puJetId*_*_*", # input variables
#    "keep *_puJetMva*_*_*", # final MVAs and working point flags
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
    ###### MET products
    'keep *_patMETs*_*_*',
#    'keep *_patType1CorrectedPFMet_*_*', # NOT included for the moment
    ### for HLT selection
    'keep edmTriggerResults_TriggerResults_*_HLT'])

# additional products for event cleaning
process.out.outputCommands.extend([
    'keep *_TriggerResults_*_PAT',
    ])

process.out.outputCommands.extend(['keep edmMergeableCounter_*_*_*'])

# Adding the CA8 information
process.out.outputCommands.extend(jetSubstructuresEventContent)

process.outPath = cms.EndPath(process.out)

#  LocalWords:  hzzeejjBaseColl
