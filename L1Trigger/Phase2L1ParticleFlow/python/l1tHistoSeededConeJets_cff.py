import FWCore.ParameterSet.Config as cms

from L1Trigger.L1CaloTrigger.Phase1L1TJetSeedProducer_cfi import l1tPhase1JetSeedProducer
from L1Trigger.Phase2L1ParticleFlow.l1tDeregionizerProducer_cfi import l1tDeregionizerProducer as l1tLayer2Deregionizer
from L1Trigger.Phase2L1ParticleFlow.l1tSeedConePFJetProducer_cfi import l1tSeedConePFJetEmulatorProducer


""" 17x17 MASK, 1X1 SEED """
l1t17x17Phase1JetSeedProducer = l1tPhase1JetSeedProducer.clone(
    jetIEtaSize = 17,    # for wide cone jets check a 17x17 square for highest pt particle rather than default 9x9
    jetIPhiSize = 17,
    trimmedGrid = False,
    outputCollectionName = "histoJetSeeds17x17",
    seedSize = 1
    )
""" 17x17 HSC jet producer for seeds with no summing """
l1t17x17HistoSeedsSCPFL1PuppiEmulator = l1tSeedConePFJetEmulatorProducer.clone(
    useExternalSeeds = True,
    JetSeeds = ('l1t17x17Phase1JetSeedProducer', 'histoJetSeeds17x17'),
    nJets = cms.uint32(16),
    coneSize = cms.double(0.8),    # cone size 0.8 for wide cone
    L1PFObjects = 'l1tLayer2Deregionizer:Puppi',
    doCorrections = cms.bool(False),
    correctorFile = cms.string("L1Trigger/Phase2L1ParticleFlow/data/jecs/jecs_20220308.root"),
    correctorDir = cms.string('L1PuppiSC4EmuJets')
    )
L1TPFHistoSeedJetsTask = cms.Task(
    l1t17x17Phase1JetSeedProducer, l1tLayer2Deregionizer, l1t17x17HistoSeedsSCPFL1PuppiEmulator
)


""" 17x17 MASK, 3x3 SEED """
l1t17x17Phase1Jet3x3SeedProducer = l1tPhase1JetSeedProducer.clone(
    jetIEtaSize = 17,    # for wide cone jets check a 17x17 square for highest pt particle rather than default 9x9
    jetIPhiSize = 17,
    trimmedGrid = False,
    outputCollectionName = "histoJet3x3Seeds17x17",
    seedSize = 3
    )
""" 17x17 HSC jet producer for seeds with no summing """
l1t17x17Histo3x3SeedsSCPFL1PuppiEmulator = l1tSeedConePFJetEmulatorProducer.clone(
    useExternalSeeds = True,
    JetSeeds = ('l1t17x17Phase1Jet3x3SeedProducer', 'histoJet3x3Seeds17x17'),
    nJets = cms.uint32(16),
    coneSize = cms.double(0.8),    # cone size 0.8 for wide cone
    L1PFObjects = 'l1tLayer2Deregionizer:Puppi',
    doCorrections = cms.bool(False),
    correctorFile = cms.string("L1Trigger/Phase2L1ParticleFlow/data/jecs/jecs_20220308.root"),
    correctorDir = cms.string('L1PuppiSC4EmuJets')
    )
L1TPFHisto3x3SeedJetsTask = cms.Task(
    l1t17x17Phase1Jet3x3SeedProducer, l1tLayer2Deregionizer, l1t17x17Histo3x3SeedsSCPFL1PuppiEmulator
)


""" 17x17 MASK, 1X1 SEED TRIMMED """
l1t17x17Phase1TrimmedJetSeedProducer = l1tPhase1JetSeedProducer.clone(
    jetIEtaSize = 17,    # for wide cone jets check a 17x17 square for highest pt particle rather than default 9x9
    jetIPhiSize = 17,
    trimmedGrid = True,
    outputCollectionName = "histoTrimmedJetSeeds17x17",
    seedSize = 1
    )
""" 17x17 HSC jet producer for seeds with no summing """
l1t17x17TrimmedHistoSeedsSCPFL1PuppiEmulator = l1tSeedConePFJetEmulatorProducer.clone(
    useExternalSeeds = True,
    JetSeeds = ('l1t17x17Phase1TrimmedJetSeedProducer', 'histoTrimmedJetSeeds17x17'),
    nJets = cms.uint32(16),
    coneSize = cms.double(0.8),    # cone size 0.8 for wide cone
    L1PFObjects = 'l1tLayer2Deregionizer:Puppi',
    doCorrections = cms.bool(False),
    correctorFile = cms.string("L1Trigger/Phase2L1ParticleFlow/data/jecs/jecs_20220308.root"),
    correctorDir = cms.string('L1PuppiSC4EmuJets')
    )
L1TPFTrimmedHistoSeedJetsTask = cms.Task(
    l1t17x17Phase1TrimmedJetSeedProducer, l1tLayer2Deregionizer, l1t17x17TrimmedHistoSeedsSCPFL1PuppiEmulator
)


""" 17x17 MASK, 3x3 SEED """
l1t15x15Phase1Jet3x3SeedProducer = l1tPhase1JetSeedProducer.clone(
    jetIEtaSize = 15,    # for wide cone jets check a 17x17 square for highest pt particle rather than default 9x9
    jetIPhiSize = 15,
    trimmedGrid = False,
    outputCollectionName = "histoJet3x3Seeds15x15",
    seedSize = 3
    )
""" 17x17 HSC jet producer for seeds with no summing """
l1t15x15Histo3x3SeedsSCPFL1PuppiEmulator = l1tSeedConePFJetEmulatorProducer.clone(
    useExternalSeeds = True,
    JetSeeds = ('l1t15x15Phase1Jet3x3SeedProducer', 'histoJet3x3Seeds15x15'),
    nJets = cms.uint32(16),
    coneSize = cms.double(0.8),    # cone size 0.8 for wide cone
    L1PFObjects = 'l1tLayer2Deregionizer:Puppi',
    doCorrections = cms.bool(False),
    correctorFile = cms.string("L1Trigger/Phase2L1ParticleFlow/data/jecs/jecs_20220308.root"),
    correctorDir = cms.string('L1PuppiSC4EmuJets')
    )
L1TPF15x15Histo3x3SeedJetsTask = cms.Task(
    l1t15x15Phase1Jet3x3SeedProducer, l1tLayer2Deregionizer, l1t15x15Histo3x3SeedsSCPFL1PuppiEmulator
)