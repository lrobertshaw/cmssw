import FWCore.ParameterSet.Config as cms

from L1Trigger.L1CaloTrigger.Phase1L1TJetSeedProducer_cfi import l1tPhase1JetSeedProducer
from L1Trigger.Phase2L1ParticleFlow.l1tDeregionizerProducer_cfi import l1tDeregionizerProducer as l1tLayer2Deregionizer
from L1Trigger.Phase2L1ParticleFlow.phase2L1TJetSeedReductionProducer_cfi import phase2L1TJetSeedReductionProducer

l1tPuppiSeedProducer = l1tPhase1JetSeedProducer.clone(
    jetIEtaSize = cms.uint32(9),
    jetIPhiSize = cms.uint32(9),
    trimmedGrid = cms.bool(True),
    seedSize = cms.uint32(1),
    inputCollectionTag = cms.InputTag('l1tLayer2Deregionizer:Puppi'),
    outputCollectionName = cms.string("seeds")
    )

l1tPuppiSeedPostprocessor = phase2L1TJetSeedReductionProducer.clone(
    seeds = cms.InputTag("l1tPuppiSeedProducer", "seeds"),
    outputCollectionName = cms.string("reducedSeeds")
    )

Phase2L1TSeedReductionTask = cms.Task(l1tLayer2Deregionizer, l1tPuppiSeedProducer, l1tPuppiSeedPostprocessor)