import FWCore.ParameterSet.Config as cms
from PhysicsTools.PatAlgos.tools.pfTools import *

# This code configures the part related to the PAT Muons.
#
# In the main code one should add the configuration with:
#
#     from HZZ2l2qAnalysis.Higgs2l2qCode.HZZ2l2qMuonConfiguration_cfi import setupPatMuons
#     setupPatMuons(process)
#

def setupPatMuons (process):

    # made from all reco muons
    process.patMuons.embedTcMETMuonCorrs = False
    process.patMuons.embedCaloMETMuonCorrs = False
    process.patMuons.embedTrack = True

    # use PFIsolation
    process.muIsoSequence = setupPFMuonIso(process, 'muons', 'PFIso')
    adaptPFIsoMuons( process, applyPostfix(process,"patMuons",""), 'PFIso')

    #
    # MuscleFit for muons:
    #
    import HZZ2l2qAnalysis.Higgs2l2qCode.Hzz2l2qSetup_cfi as Hzz2l2qSetup

    # identifier of the MuScleFit is Data2012_53X_ReReco for data
    # and Summer12_DR53X_smearReReco for MC (to compare with ReReco data)
    muscleid = 'Data2012_53X_ReReco'
    if Hzz2l2qSetup.runOnMC: muscleid = 'Summer12_DR53X_smearReReco'
    
    process.MuScleFit = cms.EDProducer("MuScleFitPATMuonCorrector",
                                       src = cms.InputTag("patMuons"),
                                       debug = cms.bool(True),
                                       #identifier is Data2012_53X_ReReco for data
                                       # and Summer12_DR53X_smearReReco for MC (to compare with ReReco data)
                                       identifier = cms.string(muscleid),
                                       applySmearing = cms.bool(Hzz2l2qSetup.runOnMC),   # Must be false in data
                                       fakeSmearing = cms.bool(False)
                                    )

# #OLD     # Sequence for muons:
# #OLD 
# #OLD     process.stdMuonSeq = cms.Sequence (
# #OLD         process.pfParticleSelectionSequence +
# #OLD         process.muIsoSequence +
# #OLD         process.makePatMuons +
# #OLD         process.MuScleFit +
# #OLD         process.selectedPatMuons +
# #OLD         )

    # Kinematic cuts on electrons: tight to reduce ntuple size:
    process.selectedPatMuons.src = cms.InputTag("MuScleFit")

    process.selectedPatMuons.cut = (
        "pt > 18 && abs(eta) < 2.4"
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
            "isGlobalMuon && isPFMuon && globalTrack().normalizedChi2 < 10 && " +
            "globalTrack().hitPattern().numberOfValidMuonHits > 0 && "              +
            "numberOfMatchedStations > 1 &&" +
            #Not working        "abs( muonBestTrack()->dz(pv.position()) ) < 0.5 && "
            "abs( userFloat('dzVtx') ) < 0.5 && "    +
            #Not working        "abs( innerTrack()->dz(pv.position()))<0.5 &&" +
            "innerTrack().hitPattern().numberOfValidPixelHits > 0 && "              +
            "track().hitPattern().trackerLayersWithMeasurement > 5 && abs(dB) < 0.2" )
        )

    # Isolated muons: standard isolation
    process.selectedIsoMuons = cms.EDFilter(
        "PATMuonSelector",
        src = cms.InputTag("selectedIDMuons"),
        #    cut = cms.string("trackIso + caloIso < 0.15 * pt")
        # using DeltaBeta correction
        #OLD     cut = cms.string("( max(0., (neutralHadronIso + photonIso - 0.5*puChargedHadronIso) ) + chargedHadronIso) < 0.12 * pt")
        cut = cms.string("userInt('isIsolated')==1")
        )

    # Kinematic cuts on electrons: tight to reduce ntuple size:

    process.selectedPatMuons.cut = (
            "pt > 18 && abs(eta) < 2.4"
        )

    # Sequence for muons:

    process.stdMuonSeq = cms.Sequence (
        process.pfParticleSelectionSequence +
        process.muIsoSequence +
        process.makePatMuons +
        process.MuScleFit +
        process.selectedPatMuons +

        process.userDataSelectedMuons +
        process.selectedIDMuons +
        process.selectedIsoMuons
        )
        

##################################################

