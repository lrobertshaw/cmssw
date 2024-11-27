import FWCore.ParameterSet.Config as cms

from L1Trigger.L1CaloTrigger.Phase1L1TJetSeedProducer_cfi import l1tPhase1JetSeedProducer
from L1Trigger.Phase2L1ParticleFlow.l1tDeregionizerProducer_cfi import l1tDeregionizerProducer as l1tLayer2Deregionizer
from L1Trigger.Phase2L1ParticleFlow.l1SeedConePFJetEmulatorProducer_cfi import l1SeedConePFJetEmulatorProducer

# Define seeds and jet producer for HW histo SC jets, so running on deregionizer cands
l1tHSC8PFL1PuppiEmuSeedProducer = l1tPhase1JetSeedProducer.clone(
    jetIEtaSize = cms.uint32(17),
    jetIPhiSize = cms.uint32(17),
    trimmedGrid = cms.bool(False),
    seedSize = cms.uint32(1),
    inputCollectionTag = cms.InputTag('l1tLayer2Deregionizer:Puppi'),
    outputCollectionName = cms.string("WIDEHSCEMUSEEDS")
    )
# Defined jet producer
l1tHSC8PFL1PuppiEmu = l1SeedConePFJetEmulatorProducer.clone(
    coneSize = cms.double(0.8),
    useExternalSeeds = cms.bool(True),
    JetSeeds = cms.InputTag('l1tHSC8PFL1PuppiEmuSeedProducer', 'WIDEHSCEMUSEEDS'),
    L1PFObjects = cms.InputTag('l1tLayer2Deregionizer:Puppi')
)
# Define task
L1TPFHSC8JetsEmuTask = cms.Task(l1tLayer2Deregionizer, l1tHSC8PFL1PuppiEmuSeedProducer, l1tHSC8PFL1PuppiEmu)

""" Trimmed """
# Define seeds and jet producer for HW histo SC jets, so running on deregionizer cands
l1tHSC8PFL1PuppiEmuSeedProducerTrimmed = l1tPhase1JetSeedProducer.clone(
    inputCollectionTag = cms.InputTag('l1tLayer2Deregionizer:Puppi'),
    jetIEtaSize = cms.uint32(9),
    jetIPhiSize = cms.uint32(9),
    
    nBinsEta = cms.uint32(72),
    etaLow = cms.double(-3),
    etaUp = cms.double(3),
    
    nBinsPhi = cms.uint32(72),
    phiLow = cms.double(-3.1415926535897931),
    phiUp = cms.double(3.1415926535897931),

    fatJet = cms.bool(True),
    seedPtThreshold = cms.double(1),    # GeV

    trimmedGrid = cms.bool(True),
    seedSize = cms.uint32(1),
    outputCollectionName = cms.string("WIDEHSCEMUSEEDSTRIMMED")
    )
# Defined jet producer
l1tHSC8PFL1PuppiEmuTrimmed = l1SeedConePFJetEmulatorProducer.clone(
    coneSize = cms.double(0.8),
    useExternalSeeds = cms.bool(True),
    JetSeeds = cms.InputTag('l1tHSC8PFL1PuppiEmuSeedProducerTrimmed', 'WIDEHSCEMUSEEDSTRIMMED'),
    L1PFObjects = cms.InputTag('l1tLayer2Deregionizer:Puppi')
)
# Define task
L1TPFHSC8JetsEmuTaskTrimmed = cms.Task(l1tLayer2Deregionizer, l1tHSC8PFL1PuppiEmuSeedProducerTrimmed, l1tHSC8PFL1PuppiEmuTrimmed)