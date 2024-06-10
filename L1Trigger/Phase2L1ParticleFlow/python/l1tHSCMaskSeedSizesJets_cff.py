import FWCore.ParameterSet.Config as cms

from L1Trigger.L1CaloTrigger.Phase1L1TJetSeedProducer_cfi import l1tPhase1JetSeedProducer
from L1Trigger.Phase2L1ParticleFlow.l1tDeregionizerProducer_cfi import l1tDeregionizerProducer as l1tLayer2Deregionizer
from L1Trigger.Phase2L1ParticleFlow.l1tSeedConePFJetProducer_cfi import l1tSeedConePFJetEmulatorProducer


# 9x9 MASK, 1X1 SEED, NOT TRIMMED
l1tHSCPFL1Puppi9x9Seed1x1Producer = l1tPhase1JetSeedProducer.clone(
    inputCollectionTag = cms.InputTag('l1tLayer2Deregionizer:Puppi'),
    outputCollectionName = cms.string("Mask9x9Seed1x1"),
    jetIEtaSize = cms.uint32(9),
    jetIPhiSize = cms.uint32(9),
    trimmedGrid = cms.bool(False),
    seedSize = cms.uint32(1),
    )
l1tHSCPFL1Puppi9x9Seed1x1Emu = l1tSeedConePFJetEmulatorProducer.clone(
    useExternalSeeds = cms.bool(True),
    JetSeeds = cms.InputTag('l1tHSCPFL1Puppi9x9Seed1x1Producer', 'Mask9x9Seed1x1'),
    L1PFObjects = cms.InputTag('l1tLayer2Deregionizer:Puppi')
)

# 9x9 MASK, 1X1 SEED, TRIMMED
l1tHSCPFL1Puppi9x9Seed1x1ProducerTrimmed = l1tPhase1JetSeedProducer.clone(
    inputCollectionTag = cms.InputTag('l1tLayer2Deregionizer:Puppi'),
    outputCollectionName = cms.string("Mask9x9Seed1x1Trimmed"),
    jetIEtaSize = cms.uint32(9),
    jetIPhiSize = cms.uint32(9),
    trimmedGrid = cms.bool(True)
    )
l1tHSCPFL1Puppi9x9Seed1x1EmuTrimmed = l1tSeedConePFJetEmulatorProducer.clone(
    useExternalSeeds = cms.bool(True),
    JetSeeds = cms.InputTag('l1tHSCPFL1Puppi9x9Seed1x1ProducerTrimmed', 'Mask9x9Seed1x1Trimmed'),
    L1PFObjects = cms.InputTag('l1tLayer2Deregionizer:Puppi')
)

# 9x9 MASK, 3X3 SEED, NOT TRIMMED
l1tHSCPFL1Puppi9x9Seed3x3Producer = l1tPhase1JetSeedProducer.clone(
    inputCollectionTag = cms.InputTag('l1tLayer2Deregionizer:Puppi'),
    outputCollectionName = cms.string("Mask9x9Seed3x3"),
    jetIEtaSize = cms.uint32(9),
    jetIPhiSize = cms.uint32(9),
    trimmedGrid = cms.bool(False),
    seedSize = cms.uint32(3),
    )
l1tHSCPFL1Puppi9x9Seed3x3Emu = l1tSeedConePFJetEmulatorProducer.clone(
    useExternalSeeds = cms.bool(True),
    JetSeeds = cms.InputTag('l1tHSCPFL1Puppi9x9Seed3x3Producer', 'Mask9x9Seed3x3'),
    L1PFObjects = cms.InputTag('l1tLayer2Deregionizer:Puppi')
)

# 9x9 MASK, 3X3 SEED, TRIMMED
l1tHSCPFL1Puppi9x9Seed3x3ProducerTrimmed = l1tPhase1JetSeedProducer.clone(
    inputCollectionTag = cms.InputTag('l1tLayer2Deregionizer:Puppi'),
    outputCollectionName = cms.string("Mask9x9Seed3x3Trimmed"),
    jetIEtaSize = cms.uint32(9),
    jetIPhiSize = cms.uint32(9),
    trimmedGrid = cms.bool(True),
    seedSize = cms.uint32(3),
    )
l1tHSCPFL1Puppi9x9Seed3x3EmuTrimmed = l1tSeedConePFJetEmulatorProducer.clone(
    useExternalSeeds = cms.bool(True),
    JetSeeds = cms.InputTag('l1tHSCPFL1Puppi9x9Seed3x3ProducerTrimmed', 'Mask9x9Seed3x3Trimmed'),
    L1PFObjects = cms.InputTag('l1tLayer2Deregionizer:Puppi')
)

""" 7X7 MASK """
# 7x7 MASK, 1X1 SEED, NOT TRIMMED
l1tHSCPFL1Puppi7x7Seed1x1Producer = l1tPhase1JetSeedProducer.clone(
    inputCollectionTag = cms.InputTag('l1tLayer2Deregionizer:Puppi'),
    outputCollectionName = cms.string("Mask7x7Seed1x1"),
    jetIEtaSize = cms.uint32(7),
    jetIPhiSize = cms.uint32(7),
    trimmedGrid = cms.bool(False),
    seedSize = cms.uint32(1),
    )
l1tHSCPFL1Puppi7x7Seed1x1Emu = l1tSeedConePFJetEmulatorProducer.clone(
    useExternalSeeds = cms.bool(True),
    JetSeeds = cms.InputTag('l1tHSCPFL1Puppi7x7Seed1x1Producer', 'Mask7x7Seed1x1'),
    L1PFObjects = cms.InputTag('l1tLayer2Deregionizer:Puppi')
)

# 7x7 MASK, 1X1 SEED, TRIMMED
l1tHSCPFL1Puppi7x7Seed1x1ProducerTrimmed = l1tPhase1JetSeedProducer.clone(
    inputCollectionTag = cms.InputTag('l1tLayer2Deregionizer:Puppi'),
    outputCollectionName = cms.string("Mask7x7Seed1x1Trimmed"),
    jetIEtaSize = cms.uint32(7),
    jetIPhiSize = cms.uint32(7),
    trimmedGrid = cms.bool(True),
    seedSize = cms.uint32(1),
    )
l1tHSCPFL1Puppi7x7Seed1x1EmuTrimmed = l1tSeedConePFJetEmulatorProducer.clone(
    useExternalSeeds = cms.bool(True),
    JetSeeds = cms.InputTag('l1tHSCPFL1Puppi7x7Seed1x1ProducerTrimmed', 'Mask7x7Seed1x1Trimmed'),
    L1PFObjects = cms.InputTag('l1tLayer2Deregionizer:Puppi')
)

# 7x7 MASK, 3X3 SEED, NOT TRIMMED
l1tHSCPFL1Puppi7x7Seed3x3Producer = l1tPhase1JetSeedProducer.clone(
    inputCollectionTag = cms.InputTag('l1tLayer2Deregionizer:Puppi'),
    outputCollectionName = cms.string("Mask7x7Seed3x3"),
    jetIEtaSize = cms.uint32(7),
    jetIPhiSize = cms.uint32(7),
    trimmedGrid = cms.bool(False),
    seedSize = cms.uint32(3),
    )
l1tHSCPFL1Puppi7x7Seed3x3Emu = l1tSeedConePFJetEmulatorProducer.clone(
    useExternalSeeds = cms.bool(True),
    JetSeeds = cms.InputTag('l1tHSCPFL1Puppi7x7Seed3x3Producer', 'Mask7x7Seed3x3'),
    L1PFObjects = cms.InputTag('l1tLayer2Deregionizer:Puppi')
)

# 7x7 MASK, 3X3 SEED, TRIMMED
l1tHSCPFL1Puppi7x7Seed3x3ProducerTrimmed = l1tPhase1JetSeedProducer.clone(
    inputCollectionTag = cms.InputTag('l1tLayer2Deregionizer:Puppi'),
    outputCollectionName = cms.string("Mask7x7Seed3x3Trimmed"),
    jetIEtaSize = cms.uint32(7),
    jetIPhiSize = cms.uint32(7),
    trimmedGrid = cms.bool(True),
    seedSize = cms.uint32(3),
    )
l1tHSCPFL1Puppi7x7Seed3x3EmuTrimmed = l1tSeedConePFJetEmulatorProducer.clone(
    useExternalSeeds = cms.bool(True),
    JetSeeds = cms.InputTag('l1tHSCPFL1Puppi7x7Seed3x3ProducerTrimmed', 'Mask7x7Seed3x3Trimmed'),
    L1PFObjects = cms.InputTag('l1tLayer2Deregionizer:Puppi')
)

""" TASKS """
l1tHSCMaskSeedSizesJetsTask = cms.Task(
    l1tLayer2Deregionizer,
    # 9x9 MASK
    l1tHSCPFL1Puppi9x9Seed1x1Producer, l1tHSCPFL1Puppi9x9Seed1x1Emu,
    l1tHSCPFL1Puppi9x9Seed3x3Producer, l1tHSCPFL1Puppi9x9Seed3x3Emu,
    # 7x7 MASK
    l1tHSCPFL1Puppi7x7Seed1x1Producer, l1tHSCPFL1Puppi7x7Seed1x1Emu,
    l1tHSCPFL1Puppi7x7Seed3x3Producer, l1tHSCPFL1Puppi7x7Seed3x3Emu
)

l1tHSCMaskSeedSizesJetsTaskTrimmed = cms.Task(
    l1tLayer2Deregionizer,
    # 9x9 MASK
    l1tHSCPFL1Puppi9x9Seed1x1ProducerTrimmed, l1tHSCPFL1Puppi9x9Seed1x1EmuTrimmed,
    l1tHSCPFL1Puppi9x9Seed3x3ProducerTrimmed, l1tHSCPFL1Puppi9x9Seed3x3EmuTrimmed,
    # 7x7 MASK
    l1tHSCPFL1Puppi7x7Seed1x1ProducerTrimmed, l1tHSCPFL1Puppi7x7Seed1x1EmuTrimmed,
    l1tHSCPFL1Puppi7x7Seed3x3ProducerTrimmed, l1tHSCPFL1Puppi7x7Seed3x3EmuTrimmed
)