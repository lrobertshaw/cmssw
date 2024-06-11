import FWCore.ParameterSet.Config as cms

from L1Trigger.L1CaloTrigger.Phase1L1TJetSeedProducer_cfi import l1tPhase1JetSeedProducer
from L1Trigger.Phase2L1ParticleFlow.l1tDeregionizerProducer_cfi import l1tDeregionizerProducer as l1tLayer2Deregionizer
from L1Trigger.Phase2L1ParticleFlow.l1tSeedConePFJetProducer_cfi import l1tSeedConePFJetEmulatorProducer


# 9x9
l1tHSCPFL1Puppi9x9SeedProducerTrimmed = l1tPhase1JetSeedProducer.clone(
    inputCollectionTag = cms.InputTag('l1tLayer2Deregionizer:Puppi'),
    outputCollectionName = cms.string("HSC9X9EmuSeedsTrimmed"),
    jetIEtaSize = cms.uint32(9),
    jetIPhiSize = cms.uint32(9),
    trimmedGrid = cms.bool(True)
    )
l1tHSCPFL1Puppi9x9EmuTrimmed = l1tSeedConePFJetEmulatorProducer.clone(
    useExternalSeeds = cms.bool(True),
    JetSeeds = cms.InputTag('l1tHSCPFL1Puppi9x9SeedProducerTrimmed', 'HSC9X9EmuSeedsTrimmed'),
    L1PFObjects = cms.InputTag('l1tLayer2Deregionizer:Puppi')

)
l1tHSCPFL1Puppi9x9EmuTaskTrimmed = cms.Task(l1tLayer2Deregionizer, l1tHSCPFL1Puppi9x9SeedProducerTrimmed, l1tHSCPFL1Puppi9x9EmuTrimmed)

# 7x7
l1tHSCPFL1Puppi7x7SeedProducerTrimmed = l1tPhase1JetSeedProducer.clone(
    inputCollectionTag = cms.InputTag('l1tLayer2Deregionizer:Puppi'),
    outputCollectionName = cms.string("HSC7X7EmuSeedsTrimmed"),
    jetIEtaSize = cms.uint32(7),
    jetIPhiSize = cms.uint32(7),
    trimmedGrid = cms.bool(True)
    )
l1tHSCPFL1Puppi7x7EmuTrimmed = l1tSeedConePFJetEmulatorProducer.clone(
    useExternalSeeds = cms.bool(True),
    JetSeeds = cms.InputTag('l1tHSCPFL1Puppi7x7SeedProducerTrimmed', 'HSC7X7EmuSeedsTrimmed'),
    L1PFObjects = cms.InputTag('l1tLayer2Deregionizer:Puppi')
)
l1tHSCPFL1Puppi7x7EmuTaskTrimmed = cms.Task(l1tLayer2Deregionizer, l1tHSCPFL1Puppi7x7SeedProducerTrimmed, l1tHSCPFL1Puppi7x7EmuTrimmed)

# 5x5
l1tHSCPFL1Puppi5x5SeedProducerTrimmed = l1tPhase1JetSeedProducer.clone(
    inputCollectionTag = cms.InputTag('l1tLayer2Deregionizer:Puppi'),
    outputCollectionName = cms.string("HSC5X5EmuSeedsTrimmed"),
    jetIEtaSize = cms.uint32(5),
    jetIPhiSize = cms.uint32(5),
    trimmedGrid = cms.bool(True)
    )
l1tHSCPFL1Puppi5x5EmuTrimmed = l1tSeedConePFJetEmulatorProducer.clone(
    useExternalSeeds = cms.bool(True),
    JetSeeds = cms.InputTag('l1tHSCPFL1Puppi5x5SeedProducerTrimmed', 'HSC5X5EmuSeedsTrimmed'),
    L1PFObjects = cms.InputTag('l1tLayer2Deregionizer:Puppi')
)
l1tHSCPFL1Puppi5x5EmuTaskTrimmed = cms.Task(l1tLayer2Deregionizer, l1tHSCPFL1Puppi5x5SeedProducerTrimmed, l1tHSCPFL1Puppi5x5EmuTrimmed)

# 3x3
l1tHSCPFL1Puppi3x3SeedProducerTrimmed = l1tPhase1JetSeedProducer.clone(
    inputCollectionTag = cms.InputTag('l1tLayer2Deregionizer:Puppi'),
    outputCollectionName = cms.string("HSC3X3EmuSeedsTrimmed"),
    jetIEtaSize = cms.uint32(3),
    jetIPhiSize = cms.uint32(3),
    trimmedGrid = cms.bool(True)
    )
l1tHSCPFL1Puppi3x3EmuTrimmed = l1tSeedConePFJetEmulatorProducer.clone(
    useExternalSeeds = cms.bool(True),
    JetSeeds = cms.InputTag('l1tHSCPFL1Puppi3x3SeedProducerTrimmed', 'HSC3X3EmuSeedsTrimmed'),
    L1PFObjects = cms.InputTag('l1tLayer2Deregionizer:Puppi')
)
l1tHSCPFL1Puppi3x3EmuTaskTrimmed = cms.Task(l1tLayer2Deregionizer, l1tHSCPFL1Puppi3x3SeedProducerTrimmed, l1tHSCPFL1Puppi3x3EmuTrimmed)