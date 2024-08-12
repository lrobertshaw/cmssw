import FWCore.ParameterSet.Config as cms

from L1Trigger.L1CaloTrigger.Phase1L1TJetSeedProducer_cfi import l1tPhase1JetSeedProducer, caloEtaSegmentation
from L1Trigger.Phase2L1ParticleFlow.l1tDeregionizerProducer_cfi import l1tDeregionizerProducer as l1tLayer2Deregionizer
from L1Trigger.Phase2L1ParticleFlow.l1SeedConePFJetEmulatorProducer_cfi import l1SeedConePFJetEmulatorProducer

newEtaBins = caloEtaSegmentation
l1tHSC8PFL1PuppiEmuSeedProducerDoubleBinSize = l1tPhase1JetSeedProducer.clone(
    inputCollectionTag = cms.InputTag('l1tLayer2Deregionizer:Puppi'),
    jetIEtaSize = cms.uint32(9),
    jetIPhiSize = cms.uint32(9),
    etaBinning = newEtaBins,    # New array with double the eta bin width
    nBinsPhi = cms.uint32(72),    # Half the num of phi bins
    trimmedGrid = cms.bool(False),    # Trim the grid
    seedSize = cms.uint32(1),    # Standard seed size
    seedPtThreshold = cms.double(1),    # GeV
    # etaRegions = cms.vdouble( -3, -2.5, -1.5, -1.0, -0.5, 0, 0.5, 1, 1.5, 2.5, 3 ),
    # phiRegions = cms.vdouble( -3.15, -2.45, -1.75, -1.05, -0.35, 0.35, 1.05, 1.75, 2.45, 3.15 ),#, 4.2, 4.9, 5.6, 6.3 ),
    outputCollectionName = cms.string("HSCseedsDoubleBinSize")
    )

l1tHSC8PFL1PuppiEmuDoubleBinSize = l1SeedConePFJetEmulatorProducer.clone(
    coneSize = cms.double(0.8),
    useExternalSeeds = cms.bool(True),
    JetSeeds = cms.InputTag('l1tHSC8PFL1PuppiEmuSeedProducerDoubleBinSize', 'HSCseedsDoubleBinSize'),
    L1PFObjects = cms.InputTag('l1tLayer2Deregionizer:Puppi')
)

L1TPFHSC8JetsEmuTaskDoubleBinSize = cms.Task(l1tLayer2Deregionizer, l1tHSC8PFL1PuppiEmuSeedProducerDoubleBinSize, l1tHSC8PFL1PuppiEmuDoubleBinSize)