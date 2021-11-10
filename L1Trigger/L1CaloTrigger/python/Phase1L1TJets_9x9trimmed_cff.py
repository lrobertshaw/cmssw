import FWCore.ParameterSet.Config as cms

from L1Trigger.L1CaloTrigger.l1tPhase1JetProducer_cfi import l1tPhase1JetProducer
from L1Trigger.L1CaloTrigger.l1tPhase1JetCalibrator_9x9trimmed_cfi import l1tPhase1JetCalibrator_9x9trimmed
from L1Trigger.L1CaloTrigger.l1tPhase1JetSumsProducer_cfi import l1tPhase1JetSumsProducer

l1tPhase1JetProducer9x9trimmed = l1tPhase1JetProducer.clone(
	  jetIEtaSize = 9,
	  jetIPhiSize = 9,
	  trimmedGrid = True,
	  outputCollectionName = "UncalibratedPhase1L1TJetFromPfCandidates"
)

Phase1L1TJetCalibrator9x9trimmed.inputCollectionTag = cms.InputTag("Phase1L1TJetProducer9x9trimmed", "UncalibratedPhase1L1TJetFromPfCandidates", "")
Phase1L1TJetCalibrator9x9trimmed.outputCollectionName = cms.string("Phase1L1TJetFromPfCandidates")

Phase1L1TJetCalibrator9x9trimmedUnsorted = Phase1L1TJetCalibrator9x9trimmed.clone(
	  inputCollectionTag = cms.InputTag("Phase1L1TJetProducer9x9trimmed", "UncalibratedPhase1L1TJetFromPfCandidatesUnsorted", ""),
    outputCollectionName = cms.string("Phase1L1TJetFromPfCandidates")
)

Phase1L1TJetSumsProducer9x9trimmed = Phase1L1TJetSumsProducer.clone(
  inputJetCollectionTag = cms.InputTag("Phase1L1TJetCalibrator9x9trimmedUnsorted", "Phase1L1TJetFromPfCandidates"),
)

Phase1L1TJetsSequence9x9trimmed = cms.Sequence(
  Phase1L1TJetProducer9x9trimmed +
  Phase1L1TJetCalibrator9x9trimmed + 
  Phase1L1TJetCalibrator9x9trimmedUnsorted +
  Phase1L1TJetSumsProducer9x9trimmed
)
