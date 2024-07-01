import FWCore.ParameterSet.Config as cms

from L1Trigger.L1CaloTrigger.l1tPhase1JetProducer_cfi import l1tPhase1JetProducer, caloEtaSegmentation
from L1Trigger.L1CaloTrigger.l1tPhase1JetCalibrator9_cfi import l1tPhase1JetCalibrator9x9
from L1Trigger.L1CaloTrigger.l1tPhase1JetSumsProducer_cfi import l1tPhase1JetSumsProducer
from L1Trigger.Phase2L1ParticleFlow.l1tDeregionizerProducer_cfi import l1tDeregionizerProducer as l1tLayer2Deregionizer

""" SINGLE BINNING 18 X 18"""
l1tPhase1WideJetProducer18x18 = l1tPhase1JetProducer.clone(
    inputCollectionTag = cms.InputTag("l1tLayer2Deregionizer:Puppi"),
    jetIEtaSize = 18,
    jetIPhiSize = 18,
    etaBinning = caloEtaSegmentation,
    nBinsPhi = cms.uint32(72),
    trimmedGrid = cms.bool(True),
    #   etaRegions = cms.vdouble( -5., -4.5, -4., -3.5, -3, -2.5, -1.5, -1.0, -0.5, 0, 0.5, 1, 1.5, 2.5, 3, 3.5, 4., 4.5, 5. ),
    #   phiRegions = cms.vdouble( -3.15, -2.45, -1.75, -1.05, -0.35, 0.35, 1.05, 1.75, 2.45, 3.15 ),#, 4.2, 4.9, 5.6, 6.3 ),
    outputCollectionName = "Uncalibratedl1tPhase1WideJet18x18FromPfCandidates"
    )

l1tPhase1WideJetCalibrator18x18 = l1tPhase1JetCalibrator9x9.clone(
    inputCollectionTag = cms.InputTag("l1tPhase1WideJetProducer18x18", "Uncalibratedl1tPhase1WideJet18x18FromPfCandidates", ""),
    outputCollectionName = cms.string("l1tPhase1WideJet18x18FromPfCandidates")
    )

l1tPhase1WideJetSumsProducer18x18 = l1tPhase1JetSumsProducer.clone(
    inputJetCollectionTag = cms.InputTag("l1tPhase1WideJetCalibrator18x18", "l1tPhase1WideJet18x18FromPfCandidates"),
    )

l1tPhase1WideJetCalibrator18x18Unsorted = l1tPhase1JetCalibrator9x9.clone(
    inputCollectionTag = cms.InputTag("l1tPhase1WideJetProducer18x18", "Uncalibratedl1tPhase1WideJet18x18FromPfCandidatesUnsorted", ""),
    outputCollectionName = cms.string("l1tPhase1WideJet18x18FromPfCandidates")
    )

l1tPhase1WideJetsSequence18x18 = cms.Sequence(
    l1tLayer2Deregionizer *
    l1tPhase1WideJetProducer18x18 *
    l1tPhase1WideJetCalibrator18x18 *
    l1tPhase1WideJetCalibrator18x18Unsorted *
    l1tPhase1WideJetSumsProducer18x18
    )

L1TPFWideJetsPhase1Task_18x18 = cms.Task(
    l1tLayer2Deregionizer,
    l1tPhase1WideJetProducer18x18,
    l1tPhase1WideJetCalibrator18x18,
    l1tPhase1WideJetCalibrator18x18Unsorted,
    l1tPhase1WideJetSumsProducer18x18
    )


""" DOUBLE BINNING 9 X 9 """
l1tPhase1WideJetProducer9x9 = l1tPhase1JetProducer.clone(
    inputCollectionTag = cms.InputTag("l1tLayer2Deregionizer:Puppi"),
    jetIEtaSize = 9,
    jetIPhiSize = 9,
    etaBinning = caloEtaSegmentation[::1],
    nBinsPhi = cms.uint32(36),
    trimmedGrid = cms.bool(True),
    #   etaRegions = cms.vdouble( -5., -4.5, -4., -3.5, -3, -2.5, -1.5, -1.0, -0.5, 0, 0.5, 1, 1.5, 2.5, 3, 3.5, 4., 4.5, 5. ),
    #   phiRegions = cms.vdouble( -3.15, -2.45, -1.75, -1.05, -0.35, 0.35, 1.05, 1.75, 2.45, 3.15 ),#, 4.2, 4.9, 5.6, 6.3 ),
    outputCollectionName = "Uncalibratedl1tPhase1WideJet9x9FromPfCandidates"
    )

l1tPhase1WideJetCalibrator9x9 = l1tPhase1JetCalibrator9x9.clone(
    inputCollectionTag = cms.InputTag("l1tPhase1WideJetProducer9x9", "Uncalibratedl1tPhase1WideJet9x9FromPfCandidates", ""),
    outputCollectionName = cms.string("l1tPhase1WideJet9x9FromPfCandidates")
    )

l1tPhase1WideJetSumsProducer9x9 = l1tPhase1JetSumsProducer.clone(
    inputJetCollectionTag = cms.InputTag("l1tPhase1WideJetCalibrator9x9", "l1tPhase1WideJet9x9FromPfCandidates"),
    )

l1tPhase1WideJetCalibrator9x9Unsorted = l1tPhase1JetCalibrator9x9.clone(
    inputCollectionTag = cms.InputTag("l1tPhase1WideJetProducer9x9", "Uncalibratedl1tPhase1WideJet9x9FromPfCandidatesUnsorted", ""),
    outputCollectionName = cms.string("l1tPhase1WideJet9x9FromPfCandidates")
    )

l1tPhase1WideJetsSequence9x9 = cms.Sequence(
    l1tLayer2Deregionizer *
    l1tPhase1WideJetProducer9x9 *
    l1tPhase1WideJetCalibrator9x9 *
    l1tPhase1WideJetCalibrator9x9Unsorted *
    l1tPhase1WideJetSumsProducer9x9
    )

L1TPFWideJetsPhase1Task_9x9 = cms.Task(
    l1tLayer2Deregionizer,
    l1tPhase1WideJetProducer9x9,
    l1tPhase1WideJetCalibrator9x9,
    l1tPhase1WideJetCalibrator9x9Unsorted,
    l1tPhase1WideJetSumsProducer9x9
    )