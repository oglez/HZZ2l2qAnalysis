import FWCore.ParameterSet.Config as cms

# This code creates the candidates modules and configures them
# to have some cleaning in the main python:
#
# In the main one one should add with
#
#    process.load('HZZ2l2qAnalysis.Higgs2l2qCode.HZZ2l2qCandidateCreation_cfi')
#

zee = cms.EDProducer("CandViewShallowCloneCombiner",
                     checkCharge = cms.bool(False),   # WE DO NOT CUT ON THE CHARGE
                     cut = cms.string('mass > 20 '),
                     #decay = cms.string("userDataSelectedElectrons@+ userDataSelectedElectrons@-")
                     decay = cms.string("selectedPatElectrons@+ selectedPatElectrons@-")
                     )

zmm = cms.EDProducer("CandViewShallowCloneCombiner",
                     checkCharge = cms.bool(False),   # WE DO NOT CUT ON THE CHARGE
                     cut = cms.string('mass > 20 '),
                     #decay = cms.string("userDataSelectedMuons@+ userDataSelectedMuons@-")
                     decay = cms.string("selectedIDMuons@+ selectedIDMuons@-")
                     )

zem = cms.EDProducer("CandViewShallowCloneCombiner",
                     checkCharge = cms.bool(False),   # WE DO NOT CUT ON THE CHARGE
                     cut = cms.string('mass > 20 '),
                     #                                 decay = cms.string("userDataSelectedElectrons@+ userDataSelectedMuons@-")
                     decay = cms.string("selectedIDMuons selectedPatElectrons")
                     #decay = cms.string("selectedPatElectrons@+ selectedIDMuons@-")
                     # It should be noted that we could safely cut on the charge of leptons since it
                     # takes the reverse combinations also... even when requesting same charge (++ gets also the --)
                     )

zjj = cms.EDProducer("CandViewShallowCloneCombiner",
                     checkCharge = cms.bool(False),
                     checkOverlap = cms.bool(False),
                     #OLD cut = cms.string(''),
                     cut = cms.string('abs( daughter(0).eta ) <2.4 && abs( daughter(1).eta ) <2.4'),
                     decay = cms.string("cleanPatJetsNoPUIsoLept cleanPatJetsNoPUIsoLept")
                     )

hzzeejjBaseColl = cms.EDProducer("CandViewCombiner",
                                 checkCharge = cms.bool(False),
                                 cut = cms.string(''),
                                 decay = cms.string("zee zjj")
                                )

hzzmmjjBaseColl = cms.EDProducer("CandViewCombiner",
                                 checkCharge = cms.bool(False),
                                 cut = cms.string(''),
                                 decay = cms.string("zmm zjj")
                                 )

hzzemjjBaseColl = cms.EDProducer("CandViewCombiner",
                                 checkCharge = cms.bool(False),
                                 cut = cms.string(''),
                                 decay = cms.string("zem zjj")
                                 )

import HZZ2l2qAnalysis.Higgs2l2qCode.Hzz2l2qSetup_cfi as Hzz2l2qSetup
    
usingGenTag = 'NotAvailable'
if Hzz2l2qSetup.runOnMC:
    usingGenTag = 'genParticles'
    
hzzeejj = cms.EDProducer("H2l2qCandidateData",
                         higgs = cms.InputTag("hzzeejjBaseColl"),
                         gensTag = cms.InputTag(usingGenTag),
                         PFCandidates = cms.InputTag("particleFlow"),
                         primaryVertices = cms.InputTag("offlinePrimaryVertices"),
                         dzCut = cms.double(0.1)
                         )

hzzmmjj = cms.EDProducer("H2l2qCandidateData",
                         higgs = cms.InputTag("hzzmmjjBaseColl"),
                         gensTag = cms.InputTag(usingGenTag),
                         PFCandidates = cms.InputTag("particleFlow"),
                         primaryVertices = cms.InputTag("offlinePrimaryVertices"),
                         dzCut = cms.double(0.1)
                         )

hzzemjj = cms.EDProducer("H2l2qCandidateData",
                         higgs = cms.InputTag("hzzemjjBaseColl"),
                         gensTag = cms.InputTag(usingGenTag),
                         PFCandidates = cms.InputTag("particleFlow"),
                         primaryVertices = cms.InputTag("offlinePrimaryVertices"),
                         dzCut = cms.double(0.1)
                         )

hzzee1jBaseColl = cms.EDProducer("CandViewCombiner",
                                 checkCharge = cms.bool(False),
                                 cut = cms.string(''),
                                 decay = cms.string("zee cleanCA8JetsNoPUIsoLept")
                                 )

hzzmm1jBaseColl = cms.EDProducer("CandViewCombiner",
                                 checkCharge = cms.bool(False),
                                 cut = cms.string(''),
                                 decay = cms.string("zmm cleanCA8JetsNoPUIsoLept")
                                 )

hzzem1jBaseColl = cms.EDProducer("CandViewCombiner",
                                 checkCharge = cms.bool(False),
                                 cut = cms.string(''),
                                 decay = cms.string("zem cleanCA8JetsNoPUIsoLept")
                                 )

hzzee1j = cms.EDProducer("H2l2qmergedCandidateData",
                         higgs = cms.InputTag("hzzee1jBaseColl"),
                         prunedjets = cms.InputTag("btaggedPatJetsCA8CHSpruned")                         
                         #,gensTag = cms.InputTag(usingGenTag),
                         #,PFCandidates = cms.InputTag("particleFlow"),
                         #,primaryVertices = cms.InputTag("offlinePrimaryVertices"),
                         #,dzCut = cms.double(0.1)
                         )

hzzmm1j = cms.EDProducer("H2l2qmergedCandidateData",
                         higgs = cms.InputTag("hzzmm1jBaseColl"),
                         prunedjets = cms.InputTag("btaggedPatJetsCA8CHSpruned")
                         #,gensTag = cms.InputTag(usingGenTag),
                         #,PFCandidates = cms.InputTag("particleFlow"),
                         #,primaryVertices = cms.InputTag("offlinePrimaryVertices"),
                         #,dzCut = cms.double(0.1)
                         )

hzzem1j = cms.EDProducer("H2l2qmergedCandidateData",
                         higgs = cms.InputTag("hzzem1jBaseColl"),
                         prunedjets = cms.InputTag("btaggedPatJetsCA8CHSpruned")                         
                         #,gensTag = cms.InputTag(usingGenTag),
                         #,PFCandidates = cms.InputTag("particleFlow"),
                         #,primaryVertices = cms.InputTag("offlinePrimaryVertices"),
                         #,dzCut = cms.double(0.1)
                         )

# Setup for a basic filtering
zll = cms.EDProducer("CandViewMerger",
                     src = cms.VInputTag("zee", "zmm", "zem")
                     )

zllFilter = cms.EDFilter("CandViewCountFilter",
                         src = cms.InputTag("zll"),
                         minNumber = cms.uint32(1)
                         )

# Oscar modified the requirements to make them a bit cleaner: not
# a Z(ll)+2 jets but one Higgs candidate

allhcand = cms.EDProducer("CandViewMerger",
                          src = cms.VInputTag("hzzeejjBaseColl","hzzmmjjBaseColl","hzzemjjBaseColl"
                                              ,"hzzee1jBaseColl","hzzmm1jBaseColl","hzzem1jBaseColl"
                                              )
                          )

candFilter = cms.EDFilter("CandViewCountFilter",
                          src = cms.InputTag("allhcand"),
                          minNumber = cms.uint32(1)
                          )

##################################################
