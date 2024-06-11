import FWCore.ParameterSet.Config as cms

from L1Trigger.L1CaloTrigger.l1tPhase1JetProducer_cfi import l1tPhase1JetProducer
from L1Trigger.L1CaloTrigger.l1tPhase1JetCalibrator_9x9trimmed_cfi import l1tPhase1JetCalibrator9x9trimmed
from L1Trigger.L1CaloTrigger.l1tPhase1JetCalibrator_cfi import l1tPhase1JetCalibrator
from L1Trigger.L1CaloTrigger.l1tPhase1JetSumsProducer_cfi import l1tPhase1JetSumsProducer


""" Trimmed emu histojets """

l1tPhase1PuppiJetProducer9x9trimmed = l1tPhase1JetProducer.clone(
  inputCollectionTag = cms.InputTag("l1tLayer2Deregionizer:Puppi"),
  jetIEtaSize = 9,
  jetIPhiSize = 9,
  trimmedGrid = True,
  outputCollectionName = "Uncalibratedl1tPhase1TrimmedJetFromPuppiCandidates"  
)

l1tPhase1PuppiJetCalibrator9x9trimmed = l1tPhase1JetCalibrator9x9trimmed.clone(
    inputCollectionTag = cms.InputTag("l1tPhase1PuppiJetProducer9x9trimmed", "Uncalibratedl1tPhase1TrimmedJetFromPuppiCandidates", ""),
    outputCollectionName = cms.string("l1tPhase1TrimmedJetFromPuppiCandidates")
)

l1tPhase1PuppiJetSumsProducer9x9trimmed = l1tPhase1JetSumsProducer.clone(
  inputJetCollectionTag = cms.InputTag("l1tPhase1PuppiJetCalibrator9x9trimmed", "l1tPhase1TrimmedJetFromPuppiCandidates"),
)

l1tPhase1PuppiJetCalibrator9x9trimmedUnsorted = l1tPhase1JetCalibrator9x9trimmed.clone(
    inputCollectionTag = cms.InputTag("l1tPhase1PuppiJetProducer9x9trimmed", "Uncalibratedl1tPhase1TrimmedJetFromPuppiCandidatesUnsorted", ""),
    outputCollectionName = cms.string("l1tPhase1TrimmedJetFromPuppiCandidates")
)

l1tPhase1PuppiJetsSequence9x9trimmed = cms.Sequence(
  l1tPhase1PuppiJetProducer9x9trimmed *
  l1tPhase1PuppiJetCalibrator9x9trimmed * 
  l1tPhase1PuppiJetCalibrator9x9trimmedUnsorted *
  l1tPhase1PuppiJetSumsProducer9x9trimmed
)

L1TPuppiJetsPhase1Task_9x9trimmed = cms.Task( 
    l1tPhase1PuppiJetProducer9x9trimmed,
    l1tPhase1PuppiJetCalibrator9x9trimmed,
    l1tPhase1PuppiJetCalibrator9x9trimmedUnsorted,
    l1tPhase1PuppiJetSumsProducer9x9trimmed
  )



""" 7x7 """
l1tPhase1PuppiJetProducer7x7trimmed = l1tPhase1JetProducer.clone(
  inputCollectionTag = cms.InputTag("l1tLayer2Deregionizer:Puppi"),
  jetIEtaSize = 7,
  jetIPhiSize = 7,
  trimmedGrid = True,
  outputCollectionName = "Uncalibratedl1tPhase1Trimmed7x7JetFromPuppiCandidates"  
)

l1tPhase1PuppiJetCalibrator7x7trimmed = l1tPhase1JetCalibrator.clone(
    inputCollectionTag = cms.InputTag("l1tPhase1PuppiJetProducer7x7trimmed", "Uncalibratedl1tPhase1Trimmed7x7JetFromPuppiCandidates", ""),
    outputCollectionName = cms.string("l1tPhase1Trimmed7x7JetFromPuppiCandidates")
)

l1tPhase1PuppiJetSumsProducer7x7trimmed = l1tPhase1JetSumsProducer.clone(
  inputJetCollectionTag = cms.InputTag("l1tPhase1PuppiJetCalibrator7x7trimmed", "l1tPhase1Trimmed7x7JetFromPuppiCandidates"),
)

l1tPhase1PuppiJetCalibrator7x7trimmedUnsorted = l1tPhase1JetCalibrator.clone(
    inputCollectionTag = cms.InputTag("l1tPhase1PuppiJetProducer7x7trimmed", "Uncalibratedl1tPhase1Trimmed7x7JetFromPuppiCandidatesUnsorted", ""),
    outputCollectionName = cms.string("l1tPhase1Trimmed7x7JetFromPuppiCandidates")
)

l1tPhase1PuppiJetsSequence7x7trimmed = cms.Sequence(
  l1tPhase1PuppiJetProducer7x7trimmed *
  l1tPhase1PuppiJetCalibrator7x7trimmed * 
  l1tPhase1PuppiJetCalibrator7x7trimmedUnsorted *
  l1tPhase1PuppiJetSumsProducer7x7trimmed
)

L1TPuppiJetsPhase1Task_7x7trimmed = cms.Task( 
    l1tPhase1PuppiJetProducer7x7trimmed,
    l1tPhase1PuppiJetCalibrator7x7trimmed,
    l1tPhase1PuppiJetCalibrator7x7trimmedUnsorted,
    l1tPhase1PuppiJetSumsProducer7x7trimmed
  )