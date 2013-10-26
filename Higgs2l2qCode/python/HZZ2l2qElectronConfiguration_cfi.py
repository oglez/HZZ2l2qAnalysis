import FWCore.ParameterSet.Config as cms

# This code configures the part related to the PAT Electrons.
#
# In the main code one should add the configuration with:
#
#     from HZZ2l2qAnalysis.Higgs2l2qCode.HZZ2l2qElectronConfiguration_cfi import setupPatElectrons
#     setupPatElectrons(process)
#

from PhysicsTools.PatAlgos.tools.pfTools import *  

def setupPatElectrons (process):

    #from PhysicsTools.PatAlgos.tools.pfTools import *
    
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
    # process.eleRegressionEnergy.useRecHitCollections = cms.bool(True)
    #process.eleRegressionEnergy.correctionsType

    # It seems that the instructions were just to run the regression, but that is
    # not what we want (I think)... we want the following one:

    import HZZ2l2qAnalysis.Higgs2l2qCode.Hzz2l2qSetup_cfi as Hzz2l2qSetup
    
    process.load("EgammaAnalysis.ElectronTools.calibratedPatElectrons_cfi")
    if Hzz2l2qSetup.runOnMC:  # For MC need to change defauls
        process.calibratedPatElectrons.isMC = Hzz2l2qSetup.runOnMC
        process.calibratedPatElectrons.inputDataset = cms.string("Summer12_LegacyPaper")
        process.calibratedPatElectrons.lumiRatio = 0.607

# #OLD    # Electron sequence:
# #OLD
# #OLD    process.stdElectronSeq = cms.Sequence(
# #OLD            process.pfParticleSelectionSequence +
# #OLD            process.eleIsoSequence +
# #OLD            process.makePatElectrons +
# #OLD            process.eleRegressionEnergy +
# #OLD            process.selectedPatElectrons
# #OLD            )

    # Classic Electrons with UserData
    process.userDataSelectedElectrons = cms.EDProducer(
        "Higgs2l2bElectronUserData",
        #src = cms.InputTag("eleRegressionEnergy"),
        src = cms.InputTag("calibratedPatElectrons"),
        rho = cms.InputTag("kt6PFJets:rho"),
        primaryVertices=cms.InputTag("offlinePrimaryVertices"),
        #It is done before applyRegressionEnergy = cms.untracked.bool(True)
        )
    
    # Hay que cambiar para que selectedPatElectrons use estos.
    process.selectedPatElectrons.src = cms.InputTag('userDataSelectedElectrons')
    
    # Kinematic cuts on electrons: tight to reduce ntuple size:
    process.selectedPatElectrons.cut = (
        "pt > 18 && abs(eta) < 2.5 && (userInt('passTriggerTight') > 0)"
        )
        

    # Basic electron selection user userdata:
# #OLD    process.selectedIDElectrons = cms.EDFilter(
# #OLD        "PATElectronSelector",
# #OLD        src = cms.InputTag("userDataSelectedElectrons"),
# #OLD        cut = cms.string("(userInt('passTriggerTight') > 0)")
# #OLD        )

    # total ID+Isolated electrons: select electrons passing LOOSE
    process.selectedIsoElectrons = cms.EDFilter(
        "PATElectronSelector",
# #OLD        src = cms.InputTag("selectedIDElectrons"),
        src = cms.InputTag("selectedPatElectrons"),
        #src = cms.InputTag("userDataSelectedElectrons"),
        # #OLD    cut = cms.string("(userFloat('cutIDCode') > 1) && (userFloat('passTriggerTight') > 0)")
        cut = cms.string("(userInt('isIsolated')==1)")
        # #OLD    cut = cms.string("electronID('eidVBTFCom95') == 7")
        )

    # Electron sequence:

    process.stdElectronSeq = cms.Sequence(
        process.pfParticleSelectionSequence +
        process.eleIsoSequence +
        process.makePatElectrons +
        process.eleRegressionEnergy +
        process.calibratedPatElectrons +
        process.userDataSelectedElectrons +
        
## No longer needed        process.selectedIDElectrons +
        process.selectedPatElectrons +
        process.selectedIsoElectrons
        )



##################################################
