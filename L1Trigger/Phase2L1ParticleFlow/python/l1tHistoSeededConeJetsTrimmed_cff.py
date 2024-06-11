import FWCore.ParameterSet.Config as cms

from L1Trigger.L1CaloTrigger.Phase1L1TJetSeedProducer_cfi import l1tPhase1JetSeedProducer
from L1Trigger.Phase2L1ParticleFlow.l1tDeregionizerProducer_cfi import l1tDeregionizerProducer as l1tLayer2Deregionizer
from L1Trigger.Phase2L1ParticleFlow.l1tSeedConePFJetProducer_cfi import l1tSeedConePFJetProducer, l1tSeedConePFJetEmulatorProducer

""" HSC4 JETS """
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
L1TPFHSCJetsEmuTaskTrimmed = cms.Task(l1tLayer2Deregionizer, l1tHSCPFL1PuppiEmulatorSeedProducerTrimmed, l1tHSCPFL1PuppiEmulatorTrimmed)


""" NO TRIMMED HSC8 JETS YET """
# EMU JETS
# Define seeds and jet producer for HW histo SC jets, so running on deregionizer cands
l1tHSCPFL1PuppiEmulatorWideSeedProducerTrimmed = l1tPhase1JetSeedProducer.clone(
    inputCollectionTag = cms.InputTag('l1tLayer2Deregionizer:Puppi'),
    jetIEtaSize = cms.uint32(17),
    jetIPhiSize = cms.uint32(17),
    trimmedGrid = cms.bool(True),
    outputCollectionName = cms.string("WIDEHSCHWSEEDSTRIMMED"),
    )
# Defined jet producer
l1tHSCPFL1PuppiWideEmulatorTrimmed = l1tSeedConePFJetEmulatorProducer.clone(
    coneSize = cms.double(0.8),
    useExternalSeeds = cms.bool(True),
    JetSeeds = cms.InputTag('l1tHSCPFL1PuppiEmulatorWideSeedProducerTrimmed', 'WIDEHSCHWSEEDSTRIMMED'),
    L1PFObjects = cms.InputTag('l1tLayer2Deregionizer:Puppi')
)
# Define task
L1TPFHSCWideJetsEmuTaskTrimmed = cms.Task(l1tLayer2Deregionizer, l1tHSCPFL1PuppiEmulatorWideSeedProducerTrimmed, l1tHSCPFL1PuppiWideEmulatorTrimmed)