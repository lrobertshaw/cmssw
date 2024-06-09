import FWCore.ParameterSet.Config as cms

from L1Trigger.L1CaloTrigger.Phase1L1TJetSeedProducer_cfi import l1tPhase1JetSeedProducer
from L1Trigger.Phase2L1ParticleFlow.l1tDeregionizerProducer_cfi import l1tDeregionizerProducer as l1tLayer2Deregionizer
from L1Trigger.Phase2L1ParticleFlow.l1tSeedConePFJetProducer_cfi import l1tSeedConePFJetProducer, l1tSeedConePFJetEmulatorProducer

""" HSC4 JETS """
# Define seeds and jet producer for SW histo SC jets, so running on layer 1 puppi cands
l1tHSCPFL1PuppiSeedProducer = l1tPhase1JetSeedProducer.clone(
    inputCollectionTag = cms.InputTag("l1tLayer1", "Puppi"),
    outputCollectionName = cms.string("HSCSWSEEDS"),
    )
# Define jet producer
l1tHSCPFL1Puppi = l1tSeedConePFJetProducer.clone(
    useExternalSeeds = cms.bool(True),
    JetSeeds = cms.InputTag('l1tHSCPFL1PuppiSeedProducer', 'HSCSWSEEDS')
)
# Define task
L1TPFHSCJetsSimTask = cms.Task(l1tLayer2Deregionizer, l1tHSCPFL1PuppiSeedProducer, l1tHSCPFL1Puppi)

# Define seeds and jet producer for HW histo SC jets, so running on deregionizer cands
l1tHSCPFL1PuppiEmulatorSeedProducer = l1tPhase1JetSeedProducer.clone(
    inputCollectionTag = cms.InputTag('l1tLayer2Deregionizer:Puppi'),
    outputCollectionName = cms.string("HSCHWSEEDS"),
    )
# Defined jet producer
l1tHSCPFL1PuppiEmulator = l1tSeedConePFJetEmulatorProducer.clone(
    useExternalSeeds = cms.bool(True),
    JetSeeds = cms.InputTag('l1tHSCPFL1PuppiEmulatorSeedProducer', 'HSCHWSEEDS'),
    L1PFObjects = cms.InputTag('l1tLayer2Deregionizer:Puppi')
)
# Define corrected jet producer
l1tHSCPFL1PuppiCorrectedEmulator = l1tHSCPFL1PuppiEmulator.clone(
    doCorrections = cms.bool(True),
    correctorFile = cms.string("L1Trigger/Phase2L1ParticleFlow/data/jecs/jecs_20220308.root"),
    correctorDir = cms.string('L1PuppiSC4EmuJets')
)
# Define task
L1TPFHSCJetsEmuTask = cms.Task(
    l1tLayer2Deregionizer, l1tHSCPFL1PuppiEmulatorSeedProducer, l1tHSCPFL1PuppiEmulator, l1tHSCPFL1PuppiCorrectedEmulator
)


""" HSC8 JETS """
# Define seeds and jet producer for SW histo SC jets, so running on layer 1 puppi cands
l1tHSCPFL1PuppiSimWideSeedProducer = l1tPhase1JetSeedProducer.clone(
    jetIEtaSize = cms.uint32(17),
    jetIPhiSize = cms.uint32(17),
    
    inputCollectionTag = cms.InputTag("l1tLayer1", "Puppi"),
    outputCollectionName = cms.string("WIDEHSCSIMSEEDS")
    )
# Define jet producer
l1tHSCPFL1PuppiSimWide = l1tSeedConePFJetProducer.clone(
    coneSize = cms.double(0.1),    # CHANGED TO TRY AND OBVS BREAK JETS
    useExternalSeeds = cms.bool(True),
    JetSeeds = cms.InputTag('l1tHSCPFL1PuppiSimWideSeedProducer', 'WIDEHSCSIMSEEDS')
)
# Define task
L1TPFHSCWideJetsSimTask = cms.Task(l1tLayer2Deregionizer, l1tHSCPFL1PuppiSimWideSeedProducer, l1tHSCPFL1PuppiSimWide)

""" DEFINE EMU JETS """
# Define seeds and jet producer for HW histo SC jets, so running on deregionizer cands
l1tHSCPFL1PuppiEmuWideSeedProducer = l1tPhase1JetSeedProducer.clone(
    jetIEtaSize = cms.uint32(17),
    jetIPhiSize = cms.uint32(17),

    inputCollectionTag = cms.InputTag('l1tLayer2Deregionizer:Puppi'),
    outputCollectionName = cms.string("WIDEHSCEMUSEEDS")
    )
# Defined jet producer
l1tHSCPFL1PuppiEmuWide = l1tSeedConePFJetEmulatorProducer.clone(
    coneSize = cms.double(0.8),
    useExternalSeeds = cms.bool(True),
    JetSeeds = cms.InputTag('l1tHSCPFL1PuppiEmuWideSeedProducer', 'WIDEHSCEMUSEEDS'),
    L1PFObjects = cms.InputTag('l1tLayer2Deregionizer:Puppi')
)
# Define task
L1TPFHSCWideJetsEmuTask = cms.Task(l1tLayer2Deregionizer, l1tHSCPFL1PuppiEmuWideSeedProducer, l1tHSCPFL1PuppiEmuWide)