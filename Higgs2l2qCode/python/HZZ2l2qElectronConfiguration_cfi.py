import FWCore.ParameterSet.Config as cms

# This code configures the part related to the PAT Electrons.
#
# In the main code one should add the configuration with:
#
#     from HZZ2l2qAnalysis.Higgs2l2qCode.HZZ2l2qElectronConfiguration_cfi import setupPatElectrons
#     setupPatElectrons(process)
#

def setupPatElectrons (process):

    from PhysicsTools.PatAlgos.tools.pfTools import *
    
    # We use electrons from GSF
    process.patElectrons.embedTrack = True
    process.patElectrons.pfElectronSource = 'particleFlow'

    # use PFIsolation
    process.eleIsoSequence = setupPFElectronIso(process, 'gsfElectrons', 'PFIso')
    adaptPFIsoElectrons( process, applyPostfix(process,"patElectrons",""), 'PFIso')

    #############
    # setup recommended 0.3 cone for electron PF isolation
    process.patElectrons.isolationValues = cms.PSet(
            pfChargedHadrons = cms.InputTag("elPFIsoValueCharged03PFIdPFIso"),
            pfChargedAll = cms.InputTag("elPFIsoValueChargedAll03PFIdPFIso"),
            pfPUChargedHadrons = cms.InputTag("elPFIsoValuePU03PFIdPFIso"),
            pfNeutralHadrons = cms.InputTag("elPFIsoValueNeutral03PFIdPFIso"),
            pfPhotons = cms.InputTag("elPFIsoValueGamma03PFIdPFIso")
            )
    # OLD: I think these are nor needed... do massive testing
    #process.patElectrons.isolationValuesNoPFId = cms.PSet(
    #    pfChargedHadrons = cms.InputTag("elPFIsoValueCharged03NoPFIdPFIso"),
    #    pfChargedAll = cms.InputTag("elPFIsoValueChargedAll03NoPFIdPFIso"),
    #    pfPUChargedHadrons = cms.InputTag("elPFIsoValuePU03NoPFIdPFIso"),
    #    pfNeutralHadrons = cms.InputTag("elPFIsoValueNeutral03NoPFIdPFIso"),
    #    pfPhotons = cms.InputTag("elPFIsoValueGamma03NoPFIdPFIso")
    #    )
    #############
    #
    # Using electron regression for momentum: TEST
    #
    process.load("EgammaAnalysis.ElectronTools.electronRegressionEnergyProducer_cfi")
    process.eleRegressionEnergy.inputElectronsTag = cms.InputTag('patElectrons')
    # Hay que cambiar para que selectedPatElectrons use estos.

    # Electron sequence:

    process.stdElectronSeq = cms.Sequence(
            process.pfParticleSelectionSequence +
            process.eleIsoSequence +
            process.makePatElectrons +
        #    process.eleRegressionEnergy +
            process.selectedPatElectrons
            )

    # Classic Electrons with UserData
    process.userDataSelectedElectrons = cms.EDProducer(
        "Higgs2l2bElectronUserData",
        src = cms.InputTag("selectedPatElectrons"),
        rho = cms.InputTag("kt6PFJets:rho"),
        primaryVertices=cms.InputTag("offlinePrimaryVertices")
        )

    # Basic electron selection user userdata:
    process.selectedIDElectrons = cms.EDFilter(
        "PATElectronSelector",
        src = cms.InputTag("userDataSelectedElectrons"),
        cut = cms.string("(userInt('passTriggerTight') > 0)")
        )

    # total ID+Isolated electrons: select electrons passing LOOSE
    process.selectedIsoElectrons = cms.EDFilter(
        "PATElectronSelector",
        src = cms.InputTag("selectedIDElectrons"),
        #    src = cms.InputTag("userDataSelectedElectrons"),
        # #OLD    cut = cms.string("(userFloat('cutIDCode') > 1) && (userFloat('passTriggerTight') > 0)")
        cut = cms.string("(userInt('isIsolated')==1)")
        # #OLD    cut = cms.string("electronID('eidVBTFCom95') == 7")
        )

    # Kinematic cuts on electrons: tight to reduce ntuple size:
    
    process.selectedPatElectrons.cut = (
        "pt > 18 && abs(eta) < 2.5"
        )

##################################################
