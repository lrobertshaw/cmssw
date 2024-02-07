import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process("IN", eras.Phase2C9)
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D49_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '123X_mcRun4_realistic_v3', '')

process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
process.load('CalibCalorimetry.CaloTPG.CaloTPGTranscoder_cfi')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('L1Trigger.TrackTrigger.TrackTrigger_cff')
process.load("L1Trigger.TrackFindingTracklet.L1HybridEmulationTracks_cff") 
process.load("L1Trigger.TrackerDTC.ProducerES_cff") 
process.load("L1Trigger.TrackerDTC.ProducerED_cff") 
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/mc/Phase2HLTTDRWinter20DIGI/TT_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW/PU200_110X_mcRun4_realistic_v3-v2/110000/005E74D6-B50E-674E-89E6-EAA9A617B476.root',)
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5))
process.options = cms.untracked.PSet( 
        wantSummary = cms.untracked.bool(True),
        numberOfThreads = cms.untracked.uint32(4),
        numberOfStreams = cms.untracked.uint32(4),
)

process.PFInputsTask = cms.Task(
    process.TTClustersFromPhase2TrackerDigis,
    process.TTStubsFromPhase2TrackerDigis,
    process.TrackerDTCProducer,
    process.offlineBeamSpot,
    process.l1tTTTracksFromTrackletEmulation,
    process.l1tTTTracksFromExtendedTrackletEmulation,
    process.TTTrackAssociatorFromPixelDigis,
    process.TTTrackAssociatorFromPixelDigisExtended,
    process.SimL1EmulatorTask
)
process.p = cms.Path(
        process.l1tLayer1 +
        process.l1tLayer2Deregionizer +
        process.l1tLayer2EG
)
process.p.associate(process.PFInputsTask)

process.out = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string("inputs110X.root"),
        outputCommands = cms.untracked.vstring("drop *",
            # --- GEN
            "keep *_genParticles_*_*",
            "keep *_ak4GenJetsNoNu_*_*",
            "keep *_genMetTrue_*_*",
            # --- Track TPs
            "keep *_l1tTTTracksFromTrackletEmulation_*_*",
            "keep *_l1tTTTracksFromExtendedTrackletEmulation_*_*",
            "keep *_TTTrackAssociatorFromPixelDigis_*_*",
            "keep *_TTTrackAssociatorFromPixelDigisExtended_*_*",
            # --- Calo TPs
            "keep *_simHcalTriggerPrimitiveDigis_*_*",
            "keep *_simCaloStage2Layer1Digis_*_*",
            "keep *_simCaloStage2Digis_*_*",
            # --- Muon TPs
            "keep *_simMuonRPCDigis_*_*",
            "keep *_simMuonGEMPadDigis_*_*",
            "keep *_simMuonGEMPadDigiClusters_*_*",
            "keep *_simDtTriggerPrimitiveDigis_*_*",
            "keep *_simCscTriggerPrimitiveDigis_*_*",
            "keep *_simTwinMuxDigis_*_*",
            "keep *_simBmtfDigis_*_*",
            "keep *_simKBmtfStubs_*_*",
            "keep *_simKBmtfDigis_*_*",
            "keep *_simEmtfDigis_*_*",
            "keep *_simOmtfDigis_*_*",
            "keep *_simGmtCaloSumDigis_*_*",
            "keep *_simGmtStage2Digis_*_*",
            "keep *_simEmtfShowers_*_*",
            "keep *_simGmtShowerDigis_*_*",
            "keep *_simCscTriggerPrimitiveDigisRun3_*_*",
            "keep *_simMuonME0PadDigis_*_*",
            "keep *_me0TriggerDigis_*_*",
            "keep *_simMuonME0PseudoReDigisCoarse_*_*",
            "keep *_me0RecHitsCoarse_*_*",
            "keep *_me0TriggerPseudoDigis_*_*",
            "keep *_me0RecHits_*_*",
            "keep *_me0Segments_*_*",
            "keep *_me0TriggerConvertedPseudoDigis_*_*",
            "keep *_simCscTriggerPrimitiveDigisPhase2_*_*",
            "keep *_simGtExtFakeStage2Digis_*_*",
            "keep *_simGtStage2Digis_*_*",
            "keep *_CalibratedDigis_*_*",
            "keep *_dtTriggerPhase2PrimitiveDigis_*_*",
            # --- HGCal TPs
            "keep l1tHGCalTriggerCellBXVector_l1tHGCalVFEProducer_*_*",
            #"keep l1tHGCalTriggerCellBXVector_l1tHGCalConcentratorProducer_*_*",
            "keep l1tHGCalMulticlusterBXVector_l1tHGCalBackEndLayer2Producer_*_*",
            "keep l1tHGCalTowerBXVector_l1tHGCalTowerProducer_*_*",
            # --- GCT reconstruction
            "keep *_l1tEGammaClusterEmuProducer_*_*",
            "keep *_l1tTowerCalibration_*_*",
            "keep *_l1tCaloJet_*_*",
            "keep *_l1tCaloJetHTT_*_*",
            # --- GTT reconstruction
            "keep *_l1tVertexFinder_*_*",
            "keep *_l1tVertexFinderEmulator_*_*",
            "keep *_l1tTrackJets_*_*",
            "keep *_l1tTrackJetsExtended_*_*",
            "keep *_l1tTrackFastJets_*_*",
            "keep *_l1tTrackerEtMiss_*_*",
            "keep *_l1tTrackerHTMiss_*_*",
            "keep *_l1tTrackJetsEmulation_*_*",
            "keep *_l1tTrackJetsExtendedEmulation_*_*",
            "keep *_l1tTrackerEmuEtMiss_*_*",
            "keep *_l1tTrackerEmuHTMiss_*_*",
            "keep *_l1tTrackerEmuHTMissExtended_*_*",
            # --- GMT reconstruction
            "keep *_l1tTkStubsGmt_*_*",
            "keep *_l1tTkMuonsGmt_*_*",
            "keep *_l1tSAMuonsGmt_*_*",
        ),
        compressionAlgorithm = cms.untracked.string('LZMA'),
        compressionLevel = cms.untracked.int32(4),
        dropMetaData = cms.untracked.string('ALL'),
        fastCloning = cms.untracked.bool(False),
        overrideInputFileSplitLevels = cms.untracked.bool(True),
        eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
        SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("p")),
)
process.e = cms.EndPath(process.out)

process.schedule = cms.Schedule([process.p,process.e])

def goSlim():
    process.out.outputCommands += [ "drop *_l1tHGCalVFEProducer_*_*", ]
