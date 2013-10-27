import FWCore.ParameterSet.Config as cms

import HZZ2l2qAnalysis.Higgs2l2qCode.Hzz2l2qSetup_cfi as Hzz2l2qSetup

# This code configures the part related to the PAT Jets.
#
# In the main code one should add the configuration with:
#
#     process.load('HZZ2l2qAnalysis.Higgs2l2qCode.HZZ2l2qJetConfiguration_cfi')
#

# We put in the MVA PU since people for VBF want it...
# Using info under https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJetID
#process.load("CMGTools.External.pujetidsequence_cff")
from CMGTools.External.pujetidsequence_cff import *

# We do use selectedPatPFJetsAK5
puJetId.jets=cms.InputTag("selectedPatJetsPFJetsAK5")
puJetMva.jets = cms.InputTag("selectedPatJetsPFJetsAK5")
puJetIdChs.jets = cms.InputTag("selectedPatJetsPFJetsAK5")
puJetMvaChs.jets = cms.InputTag("selectedPatJetsPFJetsAK5")

# Adding user variables to PAT-jets:

customPFJetsNoPUSub = cms.EDProducer(
    'PFJetUserData',
    JetInputCollection=cms.untracked.InputTag("selectedPatJetsPFJetsAK5"),
    #    is2012Data=cms.untracked.bool(True),
    #    qgMap=cms.untracked.InputTag("qglAK5PF"),
    Verbosity=cms.untracked.bool(False)
    )

if Hzz2l2qSetup.runOnMC:
    customPFJetsNoPUSub.applySmearing=cms.untracked.bool(True)
    

# Jet cleaning for patJets
cleanPatJetsNoPUIsoLept = cms.EDProducer(
    "PATJetCleaner",
    #src = cms.InputTag("customPFJetsNoPUSubCentral"),
    src = cms.InputTag("customPFJetsNoPUSub"),
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

### For merged jets: Make correct pfNoPileUp ###########

# #OLD jetSubstructuresEventContent = cms.untracked.vstring()

jetSubstructuresSequence = cms.Sequence()

from CommonTools.ParticleFlow.goodOfflinePrimaryVertices_cfi import goodOfflinePrimaryVertices
goodOfflinePrimaryVerticesForSubJets = goodOfflinePrimaryVertices.clone()

from CommonTools.ParticleFlow.pfPileUp_cfi  import *
from CommonTools.ParticleFlow.TopProjectors.pfNoPileUp_cfi import *

pfPileUpForSubJets = pfPileUp.clone(
    checkClosestZVertex = False,
    PFCandidates = 'particleFlow',
    Vertices = 'goodOfflinePrimaryVerticesForSubJets'
    )
pfNoPileUpForSubJets = pfNoPileUp.clone(
    topCollection = 'pfPileUpForSubJets',
    bottomCollection = 'particleFlow'
    )

jetSubstructuresSequence += goodOfflinePrimaryVerticesForSubJets
jetSubstructuresSequence += pfPileUpForSubJets
jetSubstructuresSequence += pfNoPileUpForSubJets

pfNoPileUpSrc = 'pfNoPileUpForSubJets'
#### Adding CA8 jets and CA8 pruned jets
# load("ExoDiBosonResonances.PATtupleProduction.PAT_ca8jets_cff")
from ExoDiBosonResonances.PATtupleProduction.PAT_ca8jets_cff import *
ca8PFJetsCHS.src = pfNoPileUpSrc
ca8PFJetsCHSpruned.src = pfNoPileUpSrc
# #OLD jetSubstructuresEventContent+=['keep *_ca8PFJetsCHSpruned_SubJets_*']
# #OLD jetSubstructuresEventContent+=['keep *_ca8GenJetsNoNu_*_*']

selectedPatJetsCA8CHS.cut = cms.string('pt > 25.0 && abs(eta) < 2.4 && getPFConstituents().size > 1')

# According to Eiko and Eduardo, some information should come from the pruned
# jet but the original jet is the one that must be used to build the candidates.
# This is a clear issue, especially since the Qjet cannot be computed from the
# jets containing subjets...
#
# Anyway we need:

if (not Hzz2l2qSetup.runOnMC):
    PATCMGJetSequenceCA8CHS.remove(jetMCSequenceCA8CHS)
    PATCMGJetSequenceCA8CHSpruned.remove(jetMCSequenceCA8CHSpruned)

jetSubstructuresSequence += PATCMGJetSequenceCA8CHS
jetSubstructuresSequence += PATCMGJetSequenceCA8CHSpruned
# #OLD jetSubstructuresSequence += selectedPatJetsCA8CHSwithNsub
# #OLD jetSubstructuresSequence += selectedPatJetsCA8CHSwithQjets
#

# Adding subjet b-tagging  (by Matthias)

from PhysicsTools.PatAlgos.producersLayer1.jetProducer_cfi import patJets
from RecoJets.JetAssociationProducers.ak5JTA_cff import *
ca8CHSprunedSubjetsJetTracksAssociatorAtVertex=ak5JetTracksAssociatorAtVertex.clone()
ca8CHSprunedSubjetsJetTracksAssociatorAtVertex.jets=cms.InputTag('ca8PFJetsCHSpruned','SubJets')
from RecoBTag.Configuration.RecoBTag_cff import * # btagging sequence
impactParameterTagInfosCA8CHSprunedSubjets=impactParameterTagInfos.clone()
impactParameterTagInfosCA8CHSprunedSubjets.jetTracks='ca8CHSprunedSubjetsJetTracksAssociatorAtVertex'
secondaryVertexTagInfosCA8CHSprunedSubjets=secondaryVertexTagInfos.clone()
secondaryVertexTagInfosCA8CHSprunedSubjets.trackIPTagInfos='impactParameterTagInfosCA8CHSprunedSubjets'
combinedSecondaryVertexBJetTagsCA8CHSprunedSubjets=combinedSecondaryVertexBJetTags.clone()
combinedSecondaryVertexBJetTagsCA8CHSprunedSubjets.tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosCA8CHSprunedSubjets"),
                                                                                    cms.InputTag("secondaryVertexTagInfosCA8CHSprunedSubjets"))
btaggingCA8CHSprunedSubjets=cms.Sequence(ca8CHSprunedSubjetsJetTracksAssociatorAtVertex+impactParameterTagInfosCA8CHSprunedSubjets+secondaryVertexTagInfosCA8CHSprunedSubjets+combinedSecondaryVertexBJetTagsCA8CHSprunedSubjets)

##jet probability
jetProbabilityBJetTagsCA8CHSprunedSubjets=jetProbabilityBJetTags.clone()
jetProbabilityBJetTagsCA8CHSprunedSubjets.tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosCA8CHSprunedSubjets"))
jetBProbabilityBJetTagsCA8CHSprunedSubjets=jetBProbabilityBJetTags.clone()
jetBProbabilityBJetTagsCA8CHSprunedSubjets.tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosCA8CHSprunedSubjets"))

btaggingJPCA8CHSprunedSubjets=cms.Sequence(jetProbabilityBJetTagsCA8CHSprunedSubjets + jetBProbabilityBJetTagsCA8CHSprunedSubjets)

patJetsCA8CHSprunedSubjets = patJets.clone()
patJetsCA8CHSprunedSubjets.jetSource = cms.InputTag('ca8PFJetsCHSpruned','SubJets')
patJetsCA8CHSprunedSubjets.addGenJetMatch = False
patJetsCA8CHSprunedSubjets.addGenPartonMatch = False
patJetsCA8CHSprunedSubjets.addJetCharge = False
patJetsCA8CHSprunedSubjets.embedCaloTowers = False
patJetsCA8CHSprunedSubjets.embedPFCandidates = False
patJetsCA8CHSprunedSubjets.addAssociatedTracks = False
patJetsCA8CHSprunedSubjets.addBTagInfo = True
patJetsCA8CHSprunedSubjets.addDiscriminators = True
patJetsCA8CHSprunedSubjets.addJetID = False
patJetsCA8CHSprunedSubjets.tagInfoSources = cms.VInputTag(cms.InputTag("secondaryVertexTagInfosCA8CHSprunedSubjets"),cms.InputTag("impactParameterTagInfosCA8CHSprunedSubjets"))
patJetsCA8CHSprunedSubjets.trackAssociationSource = cms.InputTag("ca8CHSprunedSubjetsJetTracksAssociatorAtVertex")
patJetsCA8CHSprunedSubjets.discriminatorSources = cms.VInputTag(cms.InputTag("combinedSecondaryVertexBJetTagsCA8CHSprunedSubjets"),cms.InputTag("jetProbabilityBJetTagsCA8CHSprunedSubjets"),cms.InputTag("jetBProbabilityBJetTagsCA8CHSprunedSubjets"))
patJetsCA8CHSprunedSubjets.getJetMCFlavour = False
patJetsCA8CHSprunedSubjets.addJetCorrFactors = False

jetSubstructuresSequence += btaggingCA8CHSprunedSubjets
jetSubstructuresSequence += btaggingJPCA8CHSprunedSubjets
jetSubstructuresSequence += patJetsCA8CHSprunedSubjets
# #OLD jetSubstructuresEventContent+=['keep *_patJetsCA8CHSprunedSubjets_*_*']

btaggedPatJetsCA8CHSpruned = cms.EDProducer("BoostedJetMerger",
                                                    jetSrc=cms.InputTag("patJetsCA8CHSpruned"),
                                                    subjetSrc=cms.InputTag("patJetsCA8CHSprunedSubjets")
                                                    )
jetSubstructuresSequence += btaggedPatJetsCA8CHSpruned
# #OLD jetSubstructuresEventContent+=['keep *_btaggedPatJetsCA8CHSpruned_*_*']

#
# And the following changes... according to Eduardo and Eiko:
#
# Changed to include covnention selectedPatJetsCA8CHSwithNsub.src=cms.InputTag("selectedPatJetsCA8CHSpruned")
selectedPatJetsCA8CHSwithNsub.src=cms.InputTag("selectedPatJetsCA8CHS")
#cms.InputTag("btaggedPatJetsCA8CHSpruned")
selectedPatJetsCA8CHSwithQjets.src=cms.InputTag("selectedPatJetsCA8CHSwithNsub")
#
# Otras combinaciones no funcionan... de traca. 

jetSubstructuresSequence += selectedPatJetsCA8CHSwithNsub
jetSubstructuresSequence += selectedPatJetsCA8CHSwithQjets

#### Jet cleaning for CA8 Jets
cleanCA8JetsNoPUIsoLept = cms.EDProducer(
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

# #OLD jetSubstructuresEventContent+=['keep *_cleanCA8JetsNoPUIsoLept_*_*']
jetSubstructuresSequence += cleanCA8JetsNoPUIsoLept

##################################################
