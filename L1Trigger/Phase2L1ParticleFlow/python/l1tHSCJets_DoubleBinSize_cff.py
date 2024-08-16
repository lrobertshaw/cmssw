import FWCore.ParameterSet.Config as cms

from L1Trigger.L1CaloTrigger.Phase1L1TJetSeedProducer_cfi import l1tPhase1JetSeedProducer
from L1Trigger.Phase2L1ParticleFlow.l1tDeregionizerProducer_cfi import l1tDeregionizerProducer as l1tLayer2Deregionizer
from L1Trigger.Phase2L1ParticleFlow.l1SeedConePFJetEmulatorProducer_cfi import l1SeedConePFJetEmulatorProducer

l1tHSC8PFL1PuppiEmuSeedProducerDoubleBinSize = l1tPhase1JetSeedProducer.clone(
    inputCollectionTag = cms.InputTag('l1tLayer2Deregionizer:Puppi'),
    jetIEtaSize = cms.uint32(17),
    jetIPhiSize = cms.uint32(17),
    
    nBinsEta = cms.uint32(72),
    etaLow = cms.double(-3),
    etaUp = cms.double(3),
    
    nBinsPhi = cms.uint32(72),
    phiLow = cms.double(-3.1415926535897931),
    phiUp = cms.double(3.1415926535897931),

    fatJet = cms.bool(False),
    trimmedGrid = cms.bool(True),    # Trim the grid
    seedSize = cms.uint32(1),    # Standard seed size
    seedPtThreshold = cms.double(1),    # GeV
    outputCollectionName = cms.string("HSCseedsDoubleBinSize")
    )

l1tHSC8PFL1PuppiEmuDoubleBinSize = l1SeedConePFJetEmulatorProducer.clone(
    coneSize = cms.double(0.8),
    useExternalSeeds = cms.bool(True),
    JetSeeds = cms.InputTag('l1tHSC8PFL1PuppiEmuSeedProducerDoubleBinSize', 'HSCseedsDoubleBinSize'),
    L1PFObjects = cms.InputTag('l1tLayer2Deregionizer:Puppi')
)

L1TPFHSC8JetsEmuTaskDoubleBinSize = cms.Task(l1tLayer2Deregionizer, l1tHSC8PFL1PuppiEmuSeedProducerDoubleBinSize, l1tHSC8PFL1PuppiEmuDoubleBinSize)