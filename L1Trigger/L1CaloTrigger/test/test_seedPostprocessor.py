import FWCore.ParameterSet.Config as cms
process = cms.Process("TEST")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.source = cms.Source( "PoolSource", fileNames = cms.untracked.vstring("file:ttbar.root") )


# from L1Trigger.Phase2L1ParticleFlow.l1tDeregionizerProducer_cfi import l1tDeregionizerProducer as l1tLayer2Deregionizer
# process.deregionizer = l1tLayer2Deregionizer

# Produce jet seeds
from L1Trigger.L1CaloTrigger.Phase1L1TJetSeedProducer_cfi import l1tPhase1JetSeedProducer
process.histoSeedFinder = l1tPhase1JetSeedProducer.clone(
  # inputCollectionTag = cms.InputTag('l1tLayer2Deregionizer:Puppi'),
  inputCollectionTag = cms.InputTag("l1tLayer1", "Puppi"),
  outputCollectionName = cms.string("seeds"),
  jetIEtaSize = cms.uint32(9),
  jetIPhiSize = cms.uint32(9),
  trimmedGrid = cms.bool(True),
  seedSize = cms.uint32(1),
  seedPtThreshold = cms.double(5),
  nBinsEta = cms.uint32(72),
  etaLow = cms.double(-3),
  etaUp = cms.double(3),
  nBinsPhi = cms.uint32(72),
  phiLow = cms.double(-3.1415926535897931),
  phiUp = cms.double(3.1415926535897931)
)

# Postprocess jet seeds
from L1Trigger.L1CaloTrigger.phase2L1TJetSeedReductionProducer_cfi import phase2L1TJetSeedReductionProducer
process.seedReducer = phase2L1TJetSeedReductionProducer.clone(
  seeds = cms.InputTag("histoSeedFinder", "seeds"),
  outputCollectionName = cms.string("reducedSeeds")
)

# Regular seeded cone
from L1Trigger.Phase2L1ParticleFlow.l1SeedConePFJetEmulatorProducer_cfi import l1SeedConePFJetEmulatorProducer
process.seededcone = l1SeedConePFJetEmulatorProducer.clone(
  L1PFObjects = cms.InputTag("l1tLayer1", "Puppi"),
  # L1PFObjects = cms.InputTag('l1tLayer2Deregionizer:Puppi'),
  nJets = cms.uint32(16),
  coneSize = cms.double(0.8),
  allowDoubleCounting = cms.bool(False),
  useExternalSeeds = cms.bool(False),
  debug = cms.bool(False),
  doCorrections = cms.bool(False),
  correctorFile = cms.string(''),
  correctorDir = cms.string(''),
)

# Histo-seeded cone
process.histoseededcone = process.seededcone.clone(
  allowDoubleCounting = cms.bool(False),
  useExternalSeeds = cms.bool(True),
  JetSeeds = cms.InputTag('seedReducer', 'reducedSeeds'),
)

# process.p = cms.Path( process.deregionizer * process.histoSeedFinder * process.seedReducer * process.seededcone * process.histoseededcone )
# process.p = cms.Path(process.deregionizer)
process.p = cms.Path()

process.task = cms.Task(process.histoSeedFinder, process.seedReducer, process.seededcone, process.histoseededcone)
process.p.associate(process.task)

process.TFileService = cms.Service( "TFileService", fileName = cms.string('ntuples_test.root') )
process.out = cms.OutputModule( "PoolOutputModule", fileName = cms.untracked.string('test.root') )
process.end = cms.EndPath(process.out)