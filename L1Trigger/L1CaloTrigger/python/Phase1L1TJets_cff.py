import FWCore.ParameterSet.Config as cms

from L1Trigger.L1CaloTrigger.l1tPhase1JetProducer_cfi import l1tPhase1JetProducer
from L1Trigger.L1CaloTrigger.l1tPhase1JetCalibrator9_cfi import l1tPhase1JetCalibrator9x9
from L1Trigger.L1CaloTrigger.l1tPhase1JetCalibrator_cfi import l1tPhase1JetCalibrator
from L1Trigger.L1CaloTrigger.l1tPhase1JetSumsProducer_cfi import l1tPhase1JetSumsProducer


""" EMULATOR JETS - run over layer2 deregionizer puppi """

l1tPhase1PuppiJetProducer9x9 = l1tPhase1JetProducer.clone(
  inputCollectionTag = cms.InputTag("l1tLayer2Deregionizer:Puppi"),
  jetIEtaSize = 9,
  jetIPhiSize = 9,
  outputCollectionName = cms.string("Uncalibratedl1tPhase1JetFromPuppiCandidates")
)

l1tPhase1PuppiJetCalibrator9x9 = l1tPhase1JetCalibrator9x9.clone(
    inputCollectionTag = cms.InputTag("l1tPhase1PuppiJetProducer9x9", "Uncalibratedl1tPhase1JetFromPuppiCandidates", ""),
    outputCollectionName = cms.string("l1tPhase1JetFromPuppiCandidates")
)

l1tPhase1PuppiJetSumsProducer9x9 = l1tPhase1JetSumsProducer.clone(
  inputJetCollectionTag = cms.InputTag("l1tPhase1PuppiJetCalibrator9x9", "l1tPhase1JetFromPuppiCandidates"),
)

l1tPhase1PuppiJetCalibrator9x9Unsorted = l1tPhase1JetCalibrator9x9.clone(
    inputCollectionTag = cms.InputTag("l1tPhase1PuppiJetProducer9x9", "Uncalibratedl1tPhase1JetFromPuppiCandidatesUnsorted", ""),
    outputCollectionName = cms.string("l1tPhase1JetFromPuppiCandidates")
)

l1tPhase1PuppiJetsSequence9x9 = cms.Sequence(
  l1tPhase1PuppiJetProducer9x9 *
  l1tPhase1PuppiJetCalibrator9x9 * 
  l1tPhase1PuppiJetCalibrator9x9Unsorted *
  l1tPhase1PuppiJetSumsProducer9x9
)

L1TPuppiJetsPhase1Task_9x9 = cms.Task( 
    l1tPhase1PuppiJetProducer9x9,
    l1tPhase1PuppiJetCalibrator9x9,
    l1tPhase1PuppiJetCalibrator9x9Unsorted,
    l1tPhase1PuppiJetSumsProducer9x9
  )


""" 7 X 7 """
l1tPhase1PuppiJetProducer7x7 = l1tPhase1JetProducer.clone(
  inputCollectionTag = cms.InputTag("l1tLayer2Deregionizer:Puppi"),
  outputCollectionName = cms.string("Uncalibratedl1tPhase1Jet7x7FromPuppiCandidates")
)

l1tPhase1PuppiJetCalibrator7x7 = l1tPhase1JetCalibrator.clone(
    inputCollectionTag = cms.InputTag("l1tPhase1PuppiJetProducer7x7", "Uncalibratedl1tPhase1Jet7x7FromPuppiCandidates", ""),
    outputCollectionName = cms.string("l1tPhase1Jet7x7FromPuppiCandidates")
)

l1tPhase1PuppiJetSumsProducer7x7 = l1tPhase1JetSumsProducer.clone(
  inputJetCollectionTag = cms.InputTag("l1tPhase1PuppiJetCalibrator7x7", "l1tPhase1Jet7x7FromPuppiCandidates"),
)

l1tPhase1PuppiJetCalibrator7x7Unsorted = l1tPhase1JetCalibrator.clone(
    inputCollectionTag = cms.InputTag("l1tPhase1PuppiJetProducer7x7", "Uncalibratedl1tPhase1Jet7x7FromPuppiCandidates", ""),
    outputCollectionName = cms.string("l1tPhase1Jet7x7FromPuppiCandidates")
)

l1tPhase1PuppiJetsSequence7x7 = cms.Sequence(
  l1tPhase1PuppiJetProducer7x7 *
  l1tPhase1PuppiJetCalibrator7x7 * 
  l1tPhase1PuppiJetCalibrator7x7Unsorted *
  l1tPhase1PuppiJetSumsProducer7x7
)

L1TPuppiJetsPhase1Task_7x7 = cms.Task( 
    l1tPhase1PuppiJetProducer7x7,
    l1tPhase1PuppiJetCalibrator7x7,
    l1tPhase1PuppiJetCalibrator7x7Unsorted,
    l1tPhase1PuppiJetSumsProducer7x7
  )