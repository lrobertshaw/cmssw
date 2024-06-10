import FWCore.ParameterSet.Config as cms

from L1Trigger.L1CaloTrigger.Phase1L1TJetSeedProducer_cfi import l1tPhase1JetSeedProducer
from L1Trigger.Phase2L1ParticleFlow.l1tDeregionizerProducer_cfi import l1tDeregionizerProducer as l1tLayer2Deregionizer
from L1Trigger.Phase2L1ParticleFlow.l1tSeedConePFJetProducer_cfi import l1tSeedConePFJetProducer, l1tSeedConePFJetEmulatorProducer

""" HSC4 JETS """
# SIM JETS
# Define seeds and jet producer for SW histo SC jets, so running on layer 1 puppi cands
l1tHSCPFL1PuppiSeedProducerTrimmed = l1tPhase1JetSeedProducer.clone(
    inputCollectionTag = cms.InputTag("l1tLayer1", "Puppi"),
    trimmedGrid = cms.bool(True),
    outputCollectionName = cms.string("HSCSWSEEDSTRIMMED"),
    )
# Define jet producer
l1tHSCPFL1PuppiTrimmed = l1tSeedConePFJetProducer.clone(
    useExternalSeeds = cms.bool(True),
    JetSeeds = cms.InputTag('l1tHSCPFL1PuppiSeedProducerTrimmed', 'HSCSWSEEDSTRIMMED')
)
# Define task
L1TPFHSCJetsSimTaskTrimmed = cms.Task(l1tLayer2Deregionizer, l1tHSCPFL1PuppiSeedProducerTrimmed, l1tHSCPFL1PuppiTrimmed)

# EMU JETS
# Define seeds and jet producer for HW histo SC jets, so running on deregionizer cands
l1tHSCPFL1PuppiEmulatorSeedProducerTrimmed = l1tPhase1JetSeedProducer.clone(
    inputCollectionTag = cms.InputTag('l1tLayer2Deregionizer:Puppi'),
    trimmedGrid = cms.bool(True),
    outputCollectionName = cms.string("HSCHWSEEDSTRIMMED"),
    )
# Defined jet producer
l1tHSCPFL1PuppiEmulatorTrimmed = l1tSeedConePFJetEmulatorProducer.clone(
    useExternalSeeds = cms.bool(True),
    JetSeeds = cms.InputTag('l1tHSCPFL1PuppiEmulatorSeedProducerTrimmed', 'HSCHWSEEDSTRIMMED'),
    L1PFObjects = cms.InputTag('l1tLayer2Deregionizer:Puppi')
)
# Define task
L1TPFHSCJetsEmuTaskTrimmed = cms.Task(
    l1tLayer2Deregionizer, l1tHSCPFL1PuppiEmulatorSeedProducerTrimmed, l1tHSCPFL1PuppiEmulatorTrimmed
)


""" NO TRIMMED HSC8 JETS YET """