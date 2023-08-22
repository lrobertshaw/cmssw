import FWCore.ParameterSet.Config as cms

from L1Trigger.Phase2L1ParticleFlow.l1tSeedConePFJetProducer_cfi import l1tSeedConePFJetEmulatorProducer
from L1Trigger.Phase2L1ParticleFlow.l1tDeregionizerProducer_cfi import l1tDeregionizerProducer as l1tLayer2Deregionizer

from L1Trigger.L1CaloTrigger.Phase1L1TJetSeedProducer_cfi import l1tPhase1JetSeedProducer9x9trimmed
l1tHistoSeedsSCPFL1PuppiCorrectedEmulator = l1tSeedConePFJetEmulatorProducer.clone( useExternalSeeds = True,
                                                                                  JetSeeds = ('l1tPhase1JetSeedProducer9x9trimmed', 'histoJetSeeds9x9trimmed'),
                                                                                  nJets = 12,
                                                                                  L1PFObjects = 'l1tLayer2Deregionizer:Puppi',
                                                                                  doCorrections = cms.bool(True),
                                                                                  correctorFile = cms.string("L1Trigger/Phase2L1ParticleFlow/data/jecs/jecs_20220308.root"),
                                                                                  correctorDir = cms.string('L1PuppiSC4EmuJets')
                                                                                  )

L1TPFHistoSeedJetsTask = cms.Task(
    l1tPhase1JetSeedProducer9x9trimmed, l1tLayer2Deregionizer, l1tHistoSeedsSCPFL1PuppiCorrectedEmulator
)


l1tPhase1JetSeedProducer7x7 = l1tPhase1JetSeedProducer9x9trimmed.clone( jetIEtaSize = 7,
                                                                        jetIPhiSize = 7,
                                                                        trimmedGrid = False,
                                                                        outputCollectionName = "histoJetSeeds7x7"
                                                                        )

l1t7x7HistoSeedsSCPFL1PuppiCorrectedEmulator = l1tHistoSeedsSCPFL1PuppiCorrectedEmulator.clone( JetSeeds = ('l1tPhase1JetSeedProducer7x7', 'histoJetSeeds7x7') )

L1TPF7x7HistoSeedJetsTask = cms.Task(
    l1tPhase1JetSeedProducer7x7, l1tLayer2Deregionizer, l1t7x7HistoSeedsSCPFL1PuppiCorrectedEmulator
)