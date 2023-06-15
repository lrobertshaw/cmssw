import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils # ADDED

process = cms.Process("TEST")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

fileList = FileUtils.loadListFromFile('ttbar.txt')
readFiles = cms.untracked.vstring(*fileList)

process.source = process.source = cms.Source("PoolSource",
  # fileNames = readFiles,
  fileNames = cms.untracked.vstring(
   "file:test.root",
  ),
# skipEvents = cms.untracked.uint32(3)
)

# Loads jet seeding sequence
from L1Trigger.L1CaloTrigger.Phase1L1TJetSeedProducer_cfi import l1tPhase1JetSeedProducer9x9trimmed
process.seeds9x9Trimmed = l1tPhase1JetSeedProducer9x9trimmed
# Loads seeded cone sequence
from L1Trigger.Phase2L1ParticleFlow.l1tDeregionizerProducer_cfi import l1tDeregionizerProducer
process.deregionizer = l1tDeregionizerProducer
from L1Trigger.Phase2L1ParticleFlow.l1tSeedConePFJetProducer_cfi import l1tSeedConePFJetEmulatorProducer
process.seededcone = l1tSeedConePFJetEmulatorProducer.clone(
    debug = False,
    useExternalSeeds = True,
    JetSeeds = ('seeds9x9Trimmed', 'histoJetSeeds9x9trimmed')
)
process.seededconeNoSeeds = l1tSeedConePFJetEmulatorProducer.clone(
    debug = False,
    useExternalSeeds = False,
    JetSeeds = ('seeds9x9Trimmed', 'histoJetSeeds9x9trimmed')
)   
process.out = cms.OutputModule("PoolOutputModule",
  fileName = cms.untracked.string('myOutputFile_ttbar.root'),
  outputCommands = cms.untracked.vstring(
    "drop *",
    "keep *_seed*_*_*",
    "keep *_ak4GenJetsNoNu_*_*",
    "keep *_Phase1L1TJetSeedProducer9x9trimmed_*_*",
    "keep *_l1tLayer1_Puppi_*"
  ),
)

process.p = cms.Path(process.seeds9x9Trimmed * process.deregionizer * process.seededcone * process.seededconeNoSeeds )

process.e = cms.EndPath(process.out)
