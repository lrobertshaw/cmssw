import FWCore.ParameterSet.Config as cms

from L1Trigger.L1CaloTrigger.l1tPhase1JetProducer_cfi import l1tPhase1JetProducer
from L1Trigger.L1CaloTrigger.l1tPhase1JetCalibrator9_cfi import l1tPhase1JetCalibrator9
from L1Trigger.L1CaloTrigger.l1tPhase1JetSumsProducer_cfi import l1tPhase1JetSumsProducer

l1tPhase1JetProducer9x9 = l1tPhase1JetProducer.clone(
	  jetIEtaSize = 9,
	  jetIPhiSize = 9,
	  outputCollectionName = "UncalibratedPhase1L1TJetFromPfCandidates"
)

Phase1L1TJetCalibrator9x9.inputCollectionTag = cms.InputTag("Phase1L1TJetProducer9x9", "UncalibratedPhase1L1TJetFromPfCandidates", "")
Phase1L1TJetCalibrator9x9.outputCollectionName = cms.string("Phase1L1TJetFromPfCandidates")

Phase1L1TJetSumsProducer9x9 = Phase1L1TJetSumsProducer.clone(
  inputJetCollectionTag = cms.InputTag("Phase1L1TJetCalibrator9x9Unsorted", "Phase1L1TJetFromPfCandidates"),
)

Phase1L1TJetCalibrator9x9Unsorted = Phase1L1TJetCalibrator9x9.clone(
	  inputCollectionTag = cms.InputTag("Phase1L1TJetProducer9x9", "UncalibratedPhase1L1TJetFromPfCandidatesUnsorted", ""),
    outputCollectionName = cms.string("Phase1L1TJetFromPfCandidates")
)


Phase1L1TJetsSequence9x9 = cms.Sequence(
  Phase1L1TJetProducer9x9 +
  Phase1L1TJetCalibrator9x9 + 
  Phase1L1TJetCalibrator9x9Unsorted +
  Phase1L1TJetSumsProducer9x9
)
