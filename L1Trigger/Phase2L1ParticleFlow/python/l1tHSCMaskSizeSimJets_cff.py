import FWCore.ParameterSet.Config as cms

from L1Trigger.L1CaloTrigger.Phase1L1TJetSeedProducer_cfi import l1tPhase1JetSeedProducer
from L1Trigger.Phase2L1ParticleFlow.l1tDeregionizerProducer_cfi import l1tDeregionizerProducer as l1tLayer2Deregionizer
from L1Trigger.Phase2L1ParticleFlow.l1tSeedConePFJetProducer_cfi import l1tSeedConePFJetProducer

"""
This config file defines sim (not emu) HSC jets with different mask sizes.
The seeds and jets run on l1tLayer1 Puppi candidates, not on l1tLayer2 Deregionizer PUPPI candidates.
"""

# 9x9
l1tHSCPFL1Puppi9x9SeedProducer = l1tPhase1JetSeedProducer.clone(
    inputCollectionTag = cms.InputTag("l1tLayer1", "Puppi"),
    outputCollectionName = "HSC9X9SimSeeds",
    jetIEtaSize = cms.uint32(9),
    jetIPhiSize = cms.uint32(9),
    )
l1tHSCPFL1Puppi9x9Sim = l1tSeedConePFJetProducer.clone(
    useExternalSeeds = True,
    JetSeeds = ('l1tHSCPFL1Puppi9x9SeedProducer', 'HSC9X9SimSeeds')
)
l1tHSCPFL1Puppi9x9SimTask = cms.Task(l1tLayer2Deregionizer, l1tHSCPFL1Puppi9x9SeedProducer, l1tHSCPFL1Puppi9x9Sim)

# 7x7
l1tHSCPFL1Puppi7x7SeedProducer = l1tPhase1JetSeedProducer.clone(
    inputCollectionTag = cms.InputTag("l1tLayer1", "Puppi"),
    outputCollectionName = "HSC7X7SimSeeds",
    jetIEtaSize = cms.uint32(7),
    jetIPhiSize = cms.uint32(7),
    )
l1tHSCPFL1Puppi7x7Sim = l1tSeedConePFJetProducer.clone(
    useExternalSeeds = True,
    JetSeeds = ('l1tHSCPFL1Puppi7x7SeedProducer', 'HSC7X7SimSeeds')
)
l1tHSCPFL1Puppi7x7SimTask = cms.Task(l1tLayer2Deregionizer, l1tHSCPFL1Puppi7x7SeedProducer, l1tHSCPFL1Puppi7x7Sim)

# 5x5
l1tHSCPFL1Puppi5x5SeedProducer = l1tPhase1JetSeedProducer.clone(
    inputCollectionTag = cms.InputTag("l1tLayer1", "Puppi"),
    outputCollectionName = "HSC5X5SimSeeds",
    jetIEtaSize = cms.uint32(5),
    jetIPhiSize = cms.uint32(5),
    )
l1tHSCPFL1Puppi5x5Sim = l1tSeedConePFJetProducer.clone(
    useExternalSeeds = True,
    JetSeeds = ('l1tHSCPFL1Puppi5x5SeedProducer', 'HSC5X5SimSeeds')
)
l1tHSCPFL1Puppi5x5SimTask = cms.Task(l1tLayer2Deregionizer, l1tHSCPFL1Puppi5x5SeedProducer, l1tHSCPFL1Puppi5x5Sim)

# 3x3
l1tHSCPFL1Puppi3x3SeedProducer = l1tPhase1JetSeedProducer.clone(
    inputCollectionTag = cms.InputTag("l1tLayer1", "Puppi"),
    outputCollectionName = "HSC3X3SimSeeds",
    jetIEtaSize = cms.uint32(3),
    jetIPhiSize = cms.uint32(3),
    )
l1tHSCPFL1Puppi3x3Sim = l1tSeedConePFJetProducer.clone(
    useExternalSeeds = True,
    JetSeeds = ('l1tHSCPFL1Puppi3x3SeedProducer', 'HSC3X3SimSeeds')
)
l1tHSCPFL1Puppi3x3SimTask = cms.Task(l1tLayer2Deregionizer, l1tHSCPFL1Puppi3x3SeedProducer, l1tHSCPFL1Puppi3x3Sim)