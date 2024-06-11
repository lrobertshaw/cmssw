import FWCore.ParameterSet.Config as cms

from L1Trigger.L1CaloTrigger.Phase1L1TJetSeedProducer_cfi import l1tPhase1JetSeedProducer
from L1Trigger.Phase2L1ParticleFlow.l1tDeregionizerProducer_cfi import l1tDeregionizerProducer as l1tLayer2Deregionizer
from L1Trigger.Phase2L1ParticleFlow.l1tSeedConePFJetProducer_cfi import l1tSeedConePFJetEmulatorProducer


# 9x9
l1tHSCPFL1Puppi9x9SeedProducer = l1tPhase1JetSeedProducer.clone(
    inputCollectionTag = cms.InputTag('l1tLayer2Deregionizer:Puppi'),
    outputCollectionName = cms.string("HSC9X9EmuSeeds"),
    jetIEtaSize = cms.uint32(9),
    jetIPhiSize = cms.uint32(9),
    )
l1tHSCPFL1Puppi9x9Emu = l1tSeedConePFJetEmulatorProducer.clone(
    useExternalSeeds = cms.bool(True),
    JetSeeds = cms.InputTag('l1tHSCPFL1Puppi9x9SeedProducer', 'HSC9X9EmuSeeds'),
    L1PFObjects = cms.InputTag('l1tLayer2Deregionizer:Puppi')
)
l1tHSCPFL1Puppi9x9EmuTask = cms.Task(l1tLayer2Deregionizer, l1tHSCPFL1Puppi9x9SeedProducer, l1tHSCPFL1Puppi9x9Emu)

# 7x7
l1tHSCPFL1Puppi7x7SeedProducer = l1tPhase1JetSeedProducer.clone(
    inputCollectionTag = cms.InputTag('l1tLayer2Deregionizer:Puppi'),
    outputCollectionName = cms.string("HSC7X7EmuSeeds"),
    jetIEtaSize = cms.uint32(7),
    jetIPhiSize = cms.uint32(7),
    )
l1tHSCPFL1Puppi7x7Emu = l1tSeedConePFJetEmulatorProducer.clone(
    useExternalSeeds = cms.bool(True),
    JetSeeds = cms.InputTag('l1tHSCPFL1Puppi7x7SeedProducer', 'HSC7X7EmuSeeds'),
    L1PFObjects = cms.InputTag('l1tLayer2Deregionizer:Puppi')
)
l1tHSCPFL1Puppi7x7EmuTask = cms.Task(l1tLayer2Deregionizer, l1tHSCPFL1Puppi7x7SeedProducer, l1tHSCPFL1Puppi7x7Emu)

# 5x5
l1tHSCPFL1Puppi5x5SeedProducer = l1tPhase1JetSeedProducer.clone(
    inputCollectionTag = cms.InputTag('l1tLayer2Deregionizer:Puppi'),
    outputCollectionName = cms.string("HSC5X5EmuSeeds"),
    jetIEtaSize = cms.uint32(5),
    jetIPhiSize = cms.uint32(5),
    )
l1tHSCPFL1Puppi5x5Emu = l1tSeedConePFJetEmulatorProducer.clone(
    useExternalSeeds = cms.bool(True),
    JetSeeds = cms.InputTag('l1tHSCPFL1Puppi5x5SeedProducer', 'HSC5X5EmuSeeds'),
    L1PFObjects = cms.InputTag('l1tLayer2Deregionizer:Puppi')
)
l1tHSCPFL1Puppi5x5EmuTask = cms.Task(l1tLayer2Deregionizer, l1tHSCPFL1Puppi5x5SeedProducer, l1tHSCPFL1Puppi5x5Emu)

# 3x3
l1tHSCPFL1Puppi3x3SeedProducer = l1tPhase1JetSeedProducer.clone(
    inputCollectionTag = cms.InputTag('l1tLayer2Deregionizer:Puppi'),
    outputCollectionName = cms.string("HSC3X3EmuSeeds"),
    jetIEtaSize = cms.uint32(3),
    jetIPhiSize = cms.uint32(3),
    )
l1tHSCPFL1Puppi3x3Emu = l1tSeedConePFJetEmulatorProducer.clone(
    useExternalSeeds = cms.bool(True),
    JetSeeds = cms.InputTag('l1tHSCPFL1Puppi3x3SeedProducer', 'HSC3X3EmuSeeds'),
    L1PFObjects = cms.InputTag('l1tLayer2Deregionizer:Puppi')
)
l1tHSCPFL1Puppi3x3EmuTask = cms.Task(l1tLayer2Deregionizer, l1tHSCPFL1Puppi3x3SeedProducer, l1tHSCPFL1Puppi3x3Emu)