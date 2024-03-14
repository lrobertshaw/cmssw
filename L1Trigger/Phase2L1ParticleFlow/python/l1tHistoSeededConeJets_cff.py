import FWCore.ParameterSet.Config as cms

from L1Trigger.L1CaloTrigger.Phase1L1TJetSeedProducer_cfi import l1tPhase1JetSeedProducer
from L1Trigger.Phase2L1ParticleFlow.l1tDeregionizerProducer_cfi import l1tDeregionizerProducer as l1tLayer2Deregionizer
from L1Trigger.Phase2L1ParticleFlow.l1tSeedConePFJetProducer_cfi import l1tSeedConePFJetEmulatorProducer


""" 9X9 MASK, 1X1 SEED """
l1t9x9Phase1JetSeedProducer = l1tPhase1JetSeedProducer.clone(
    jetIEtaSize = 9,
    jetIPhiSize = 9,
    trimmedGrid = False,
    outputCollectionName = "histoJetSeeds9x9",
    seedSize = 1
    )
""" 9x9 HSC jet producer for seeds with no summing """
l1t9x9HistoSeedsSCPFL1PuppiCorrectedEmulator = l1tSeedConePFJetEmulatorProducer.clone(
    useExternalSeeds = True,
    JetSeeds = ('l1t9x9Phase1JetSeedProducer', 'histoJetSeeds9x9'),
    nJets = 16,
    L1PFObjects = 'l1tLayer2Deregionizer:Puppi',
    doCorrections = cms.bool(True),
    correctorFile = cms.string("L1Trigger/Phase2L1ParticleFlow/data/jecs/jecs_20220308.root"),
    correctorDir = cms.string('L1PuppiSC4EmuJets')
    )
L1TPFHistoSeedJetsTask = cms.Task(
    l1t9x9Phase1JetSeedProducer, l1tLayer2Deregionizer, l1t9x9HistoSeedsSCPFL1PuppiCorrectedEmulator
)


"""  9X9 TRIMMED MASK, 1X1 SEED  """
""" 9x9 HSC jet seeds with no summing """
l1t9x9Phase1TrimmedJetSeedProducer = l1tPhase1JetSeedProducer.clone(
    jetIEtaSize = 9,
    jetIPhiSize = 9,
    trimmedGrid = True,
    outputCollectionName = "histoTrimmedJetSeeds9x9",
    seedSize = 1
    )
""" 9x9 HSC jet producer for seeds with no summing """
l1t9x9HistoTrimmedSeedsSCPFL1PuppiCorrectedEmulator = l1tSeedConePFJetEmulatorProducer.clone(
    useExternalSeeds = True,
    JetSeeds = ('l1t9x9Phase1TrimmedJetSeedProducer', 'histoTrimmedJetSeeds9x9'),
    nJets = 16,
    L1PFObjects = 'l1tLayer2Deregionizer:Puppi',
    doCorrections = cms.bool(True),
    correctorFile = cms.string("L1Trigger/Phase2L1ParticleFlow/data/jecs/jecs_20220308.root"),
    correctorDir = cms.string('L1PuppiSC4EmuJets')
    )
L1TPFHistoSeedTrimmedJetsTask = cms.Task(
    l1t9x9Phase1TrimmedJetSeedProducer, l1tLayer2Deregionizer, l1t9x9HistoTrimmedSeedsSCPFL1PuppiCorrectedEmulator
)


"""  9X9 MASK, 3X3 SEED  """
""" 9x9 HSC jet seeds with 3x3 summing """
l1t9x9Phase1JetSeedProducer3x3 = l1tPhase1JetSeedProducer.clone(
    jetIEtaSize = 9,
    jetIPhiSize = 9,
    trimmedGrid = False,
    outputCollectionName = "histoJetSeeds9x93x3",
    seedSize = 3
    )
""" 9x9 HSC jet producer for seeds with 3x3 summing """
l1t9x9HistoSeedsSCPFL1PuppiCorrectedEmulator3x3 = l1tSeedConePFJetEmulatorProducer.clone(
    useExternalSeeds = True,
    JetSeeds = ('l1t9x9Phase1JetSeedProducer3x3', 'histoJetSeeds9x93x3'),
    nJets = 16,
    L1PFObjects = 'l1tLayer2Deregionizer:Puppi',
    doCorrections = cms.bool(True),
    correctorFile = cms.string("L1Trigger/Phase2L1ParticleFlow/data/jecs/jecs_20220308.root"),
    correctorDir = cms.string('L1PuppiSC4EmuJets')
    )

L1TPF9x9HistoSeedJetsTask3x3 = cms.Task(
    l1t9x9Phase1JetSeedProducer3x3, l1tLayer2Deregionizer, l1t9x9HistoSeedsSCPFL1PuppiCorrectedEmulator3x3
)


"""  9X9 TRIMMED MASK, 3X3 SEED  """
""" 9x9 HSC jet seeds with 3x3 summing """
l1t9x9Phase1TrimmedJetSeedProducer3x3 = l1tPhase1JetSeedProducer.clone(
    jetIEtaSize = 9,
    jetIPhiSize = 9,
    trimmedGrid = True,
    outputCollectionName = "histoTrimmedJetSeeds9x93x3",
    seedSize = 3
    )
""" 9x9 HSC jet producer for seeds with 3x3 summing """
l1t9x9HistoTrimmedSeedsSCPFL1PuppiCorrectedEmulator3x3 = l1tSeedConePFJetEmulatorProducer.clone(
    useExternalSeeds = True,
    JetSeeds = ('l1t9x9Phase1TrimmedJetSeedProducer3x3', 'histoTrimmedJetSeeds9x93x3'),
    nJets = 16,
    L1PFObjects = 'l1tLayer2Deregionizer:Puppi',
    doCorrections = cms.bool(True),
    correctorFile = cms.string("L1Trigger/Phase2L1ParticleFlow/data/jecs/jecs_20220308.root"),
    correctorDir = cms.string('L1PuppiSC4EmuJets')
    )

L1TPF9x9HistoSeedTrimmedJetsTask3x3 = cms.Task(
    l1t9x9Phase1TrimmedJetSeedProducer3x3, l1tLayer2Deregionizer, l1t9x9HistoTrimmedSeedsSCPFL1PuppiCorrectedEmulator3x3
)






"""  7X7 MASK, 1X1 SEED  """
""" 7x7 HSC jet seeds with no summing """
l1t7x7Phase1JetSeedProducer = l1tPhase1JetSeedProducer.clone(
    jetIEtaSize = 7,
    jetIPhiSize = 7,
    trimmedGrid = False,
    outputCollectionName = "histoJetSeeds7x7",
    seedSize = 1
    )
""" 9x9 HSC jet producer for seeds with no summing """
l1t7x7HistoSeedsSCPFL1PuppiCorrectedEmulator = l1tSeedConePFJetEmulatorProducer.clone(
    useExternalSeeds = True,
    JetSeeds = ('l1t7x7Phase1JetSeedProducer', 'histoJetSeeds7x7'),
    nJets = 16,
    L1PFObjects = 'l1tLayer2Deregionizer:Puppi',
    doCorrections = cms.bool(True),
    correctorFile = cms.string("L1Trigger/Phase2L1ParticleFlow/data/jecs/jecs_20220308.root"),
    correctorDir = cms.string('L1PuppiSC4EmuJets')
    )
L1TPF7x7HistoSeedJetsTask = cms.Task(
    l1t7x7Phase1JetSeedProducer, l1tLayer2Deregionizer, l1t7x7HistoSeedsSCPFL1PuppiCorrectedEmulator
)


"""  7X7 TRIMMED MASK, 1X1 SEED  """
""" 7x7 trimmed HSC jet seeds with no summing """
l1t7x7Phase1TrimmedJetSeedProducer = l1tPhase1JetSeedProducer.clone(
    jetIEtaSize = 7,
    jetIPhiSize = 7,
    trimmedGrid = True,
    outputCollectionName = "histoTrimmedJetSeeds7x7",
    seedSize = 1
    )
""" 9x9 HSC jet producer for seeds with no summing """
l1t7x7HistoTrimmedSeedsSCPFL1PuppiCorrectedEmulator = l1tSeedConePFJetEmulatorProducer.clone(
    useExternalSeeds = True,
    JetSeeds = ('l1t7x7Phase1TrimmedJetSeedProducer', 'histoTrimmedJetSeeds7x7'),
    nJets = 16,
    L1PFObjects = 'l1tLayer2Deregionizer:Puppi',
    doCorrections = cms.bool(True),
    correctorFile = cms.string("L1Trigger/Phase2L1ParticleFlow/data/jecs/jecs_20220308.root"),
    correctorDir = cms.string('L1PuppiSC4EmuJets')
    )
L1TPF7x7HistoSeedTrimmedJetsTask = cms.Task(
    l1t7x7Phase1TrimmedJetSeedProducer, l1tLayer2Deregionizer, l1t7x7HistoTrimmedSeedsSCPFL1PuppiCorrectedEmulator
)


"""  7X7 MASK, 3X3 SEED  """
l1t7x7Phase1JetSeedProducer3x3 = l1tPhase1JetSeedProducer.clone(
    jetIEtaSize = 7,
    jetIPhiSize = 7,
    trimmedGrid = False,
    outputCollectionName = "histoJetSeeds7x73x3",
    seedSize = 3
    )
""" 7x7 HSC jet producer for seeds with 3x3 summing """
l1t7x7HistoSeedsSCPFL1PuppiCorrectedEmulator3x3 = l1tSeedConePFJetEmulatorProducer.clone(
    useExternalSeeds = True,
    JetSeeds = ('l1t7x7Phase1JetSeedProducer3x3', 'histoJetSeeds7x73x3'),
    nJets = 16,
    L1PFObjects = 'l1tLayer2Deregionizer:Puppi',
    doCorrections = cms.bool(True),
    correctorFile = cms.string("L1Trigger/Phase2L1ParticleFlow/data/jecs/jecs_20220308.root"),
    correctorDir = cms.string('L1PuppiSC4EmuJets')
    )

L1TPF7x7HistoSeedJetsTask3x3 = cms.Task(
    l1t7x7Phase1JetSeedProducer3x3, l1tLayer2Deregionizer, l1t7x7HistoSeedsSCPFL1PuppiCorrectedEmulator3x3
)


"""  7X7 TRIMMED MASK, 3X3 SEED  """
l1t7x7Phase1TrimmedJetSeedProducer3x3 = l1tPhase1JetSeedProducer.clone(
    jetIEtaSize = 7,
    jetIPhiSize = 7,
    trimmedGrid = True,
    outputCollectionName = "histoTrimmedJetSeeds7x73x3",
    seedSize = 3
    )
""" 7x7 HSC jet producer for seeds with 3x3 summing """
l1t7x7HistoTrimmedSeedsSCPFL1PuppiCorrectedEmulator3x3 = l1tSeedConePFJetEmulatorProducer.clone(
    useExternalSeeds = True,
    JetSeeds = ('l1t7x7Phase1TrimmedJetSeedProducer3x3', 'histoTrimmedJetSeeds7x73x3'),
    nJets = 16,
    L1PFObjects = 'l1tLayer2Deregionizer:Puppi',
    doCorrections = cms.bool(True),
    correctorFile = cms.string("L1Trigger/Phase2L1ParticleFlow/data/jecs/jecs_20220308.root"),
    correctorDir = cms.string('L1PuppiSC4EmuJets')
    )
L1TPF7x7HistoSeedTrimmedJetsTask3x3 = cms.Task(
    l1t7x7Phase1TrimmedJetSeedProducer3x3, l1tLayer2Deregionizer, l1t7x7HistoTrimmedSeedsSCPFL1PuppiCorrectedEmulator3x3
)





# """ 9x9 HSC jet seeds with 3x3 summing """
# l1t9x9Phase1JetSeedProducer3x3 = l1tPhase1JetSeedProducer.clone(
#     jetIEtaSize = 9,
#     jetIPhiSize = 9,
#     trimmedGrid = False,
#     outputCollectionName = "histoJetSeeds9x93x3",
#     seedSize = 3
#     )
# """ 9x9 HSC jet producer for seeds with 3x3 summing """
# l1t9x9HistoSeedsSCPFL1PuppiCorrectedEmulator3x3 = l1tSeedConePFJetEmulatorProducer.clone(
#     useExternalSeeds = True,
#     JetSeeds = ('l1t9x9Phase1JetSeedProducer3x3', 'histoJetSeeds9x93x3'),
#     nJets = 16,
#     L1PFObjects = 'l1tLayer2Deregionizer:Puppi',
#     doCorrections = cms.bool(False),
#     correctorFile = cms.string("L1Trigger/Phase2L1ParticleFlow/data/jecs/jecs_20220308.root"),
#     correctorDir = cms.string('L1PuppiSC4EmuJets')
#     )

# L1TPF9x9HistoSeedJetsTask3x3 = cms.Task(
#     l1t9x9Phase1JetSeedProducer3x3, l1tLayer2Deregionizer, l1t9x9HistoSeedsSCPFL1PuppiCorrectedEmulator3x3
# )


# """ 9x9 HSC jet seeds with 5x5 summing """
# l1t9x9Phase1JetSeedProducer5x5 = l1tPhase1JetSeedProducer.clone(
#     jetIEtaSize = 9,
#     jetIPhiSize = 9,
#     trimmedGrid = False,
#     outputCollectionName = "histoJetSeeds9x95x5",
#     seedSize = 5
#     )
# """ 9x9 HSC jet producer for seeds with 5x5 summing """
# l1t9x9HistoSeedsSCPFL1PuppiCorrectedEmulator5x5 = l1tSeedConePFJetEmulatorProducer.clone(
#     useExternalSeeds = True,
#     JetSeeds = ('l1t9x9Phase1JetSeedProducer5x5', 'histoJetSeeds9x95x5'),
#     nJets = 16,
#     L1PFObjects = 'l1tLayer2Deregionizer:Puppi',
#     doCorrections = cms.bool(False),
#     correctorFile = cms.string("L1Trigger/Phase2L1ParticleFlow/data/jecs/jecs_20220308.root"),
#     correctorDir = cms.string('L1PuppiSC4EmuJets')
#     )

# L1TPF9x9HistoSeedJetsTask5x5 = cms.Task(
#     l1t9x9Phase1JetSeedProducer5x5, l1tLayer2Deregionizer, l1t9x9HistoSeedsSCPFL1PuppiCorrectedEmulator5x5
# )


# """ 9x9 HSC jet seeds with 7x7 summing """
# l1t9x9Phase1JetSeedProducer7x7 = l1tPhase1JetSeedProducer.clone(
#     jetIEtaSize = 9,
#     jetIPhiSize = 9,
#     trimmedGrid = False,
#     outputCollectionName = "histoJetSeeds9x97x7",
#     seedSize = 7
#     )
# """ 9x9 HSC jet producer for seeds with 7x7 summing """
# l1t9x9HistoSeedsSCPFL1PuppiCorrectedEmulator7x7 = l1tSeedConePFJetEmulatorProducer.clone(
#     useExternalSeeds = True,
#     JetSeeds = ('l1t9x9Phase1JetSeedProducer7x7', 'histoJetSeeds9x97x7'),
#     nJets = 16,
#     L1PFObjects = 'l1tLayer2Deregionizer:Puppi',
#     doCorrections = cms.bool(False),
#     correctorFile = cms.string("L1Trigger/Phase2L1ParticleFlow/data/jecs/jecs_20220308.root"),
#     correctorDir = cms.string('L1PuppiSC4EmuJets')
#     )

# L1TPF9x9HistoSeedJetsTask7x7 = cms.Task(
#     l1t9x9Phase1JetSeedProducer7x7, l1tLayer2Deregionizer, l1t9x9HistoSeedsSCPFL1PuppiCorrectedEmulator7x7
# )


# """ 9x9 HSC jet seeds with 9x9 summing """
# l1t9x9Phase1JetSeedProducer9x9 = l1tPhase1JetSeedProducer.clone(
#     jetIEtaSize = 9,
#     jetIPhiSize = 9,
#     trimmedGrid = False,
#     outputCollectionName = "histoJetSeeds9x99x9",
#     seedSize = 9
#     )
# """ 9x9 HSC jet producer for seeds with 9x9 summing """
# l1t9x9HistoSeedsSCPFL1PuppiCorrectedEmulator9x9 = l1tSeedConePFJetEmulatorProducer.clone(
#     useExternalSeeds = True,
#     JetSeeds = ('l1t9x9Phase1JetSeedProducer9x9', 'histoJetSeeds9x99x9'),
#     nJets = 16,
#     L1PFObjects = 'l1tLayer2Deregionizer:Puppi',
#     doCorrections = cms.bool(False),
#     correctorFile = cms.string("L1Trigger/Phase2L1ParticleFlow/data/jecs/jecs_20220308.root"),
#     correctorDir = cms.string('L1PuppiSC4EmuJets')
#     )

# L1TPF9x9HistoSeedJetsTask9x9 = cms.Task(
#     l1t9x9Phase1JetSeedProducer9x9, l1tLayer2Deregionizer, l1t9x9HistoSeedsSCPFL1PuppiCorrectedEmulator9x9
# )














# """ 7x7 HSC jet seeds with no summing """
# l1t7x7Phase1JetSeedProducer = l1tPhase1JetSeedProducer.clone(
#     jetIEtaSize = 7,
#     jetIPhiSize = 7,
#     trimmedGrid = False,
#     outputCollectionName = "histoJetSeeds7x7",
#     seedSize = 1
#     )
# """ 9x9 HSC jet producer for seeds with no summing """
# l1t7x7HistoSeedsSCPFL1PuppiCorrectedEmulator = l1tSeedConePFJetEmulatorProducer.clone(
#     useExternalSeeds = True,
#     JetSeeds = ('l1t7x7Phase1JetSeedProducer', 'histoJetSeeds7x7'),
#     nJets = 16,
#     L1PFObjects = 'l1tLayer2Deregionizer:Puppi',
#     doCorrections = cms.bool(False),
#     correctorFile = cms.string("L1Trigger/Phase2L1ParticleFlow/data/jecs/jecs_20220308.root"),
#     correctorDir = cms.string('L1PuppiSC4EmuJets')
#     )

# L1TPF7x7HistoSeedJetsTask = cms.Task(
#     l1t7x7Phase1JetSeedProducer, l1tLayer2Deregionizer, l1t7x7HistoSeedsSCPFL1PuppiCorrectedEmulator
# )


# """ 9x9 HSC jet seeds with 3x3 summing """
# l1t7x7Phase1JetSeedProducer3x3 = l1tPhase1JetSeedProducer.clone(
#     jetIEtaSize = 7,
#     jetIPhiSize = 7,
#     trimmedGrid = False,
#     outputCollectionName = "histoJetSeeds7x73x3",
#     seedSize = 3
#     )
# """ 7x7 HSC jet producer for seeds with 3x3 summing """
# l1t7x7HistoSeedsSCPFL1PuppiCorrectedEmulator3x3 = l1tSeedConePFJetEmulatorProducer.clone(
#     useExternalSeeds = True,
#     JetSeeds = ('l1t7x7Phase1JetSeedProducer3x3', 'histoJetSeeds7x73x3'),
#     nJets = 16,
#     L1PFObjects = 'l1tLayer2Deregionizer:Puppi',
#     doCorrections = cms.bool(False),
#     correctorFile = cms.string("L1Trigger/Phase2L1ParticleFlow/data/jecs/jecs_20220308.root"),
#     correctorDir = cms.string('L1PuppiSC4EmuJets')
#     )

# L1TPF7x7HistoSeedJetsTask3x3 = cms.Task(
#     l1t7x7Phase1JetSeedProducer3x3, l1tLayer2Deregionizer, l1t7x7HistoSeedsSCPFL1PuppiCorrectedEmulator3x3
# )


# """ 7x7 HSC jet seeds with 5x5 summing """
# l1t7x7Phase1JetSeedProducer5x5 = l1tPhase1JetSeedProducer.clone(
#     jetIEtaSize = 7,
#     jetIPhiSize = 7,
#     trimmedGrid = False,
#     outputCollectionName = "histoJetSeeds7x75x5",
#     seedSize = 5
#     )
# """ 7x7 HSC jet producer for seeds with 5x5 summing """
# l1t7x7HistoSeedsSCPFL1PuppiCorrectedEmulator5x5 = l1tSeedConePFJetEmulatorProducer.clone(
#     useExternalSeeds = True,
#     JetSeeds = ('l1t7x7Phase1JetSeedProducer5x5', 'histoJetSeeds7x75x5'),
#     nJets = 16,
#     L1PFObjects = 'l1tLayer2Deregionizer:Puppi',
#     doCorrections = cms.bool(False),
#     correctorFile = cms.string("L1Trigger/Phase2L1ParticleFlow/data/jecs/jecs_20220308.root"),
#     correctorDir = cms.string('L1PuppiSC4EmuJets')
#     )

# L1TPF7x7HistoSeedJetsTask5x5 = cms.Task(
#     l1t7x7Phase1JetSeedProducer5x5, l1tLayer2Deregionizer, l1t7x7HistoSeedsSCPFL1PuppiCorrectedEmulator5x5
# )


# """ 7x7 HSC jet seeds with 7x7 summing """
# l1t7x7Phase1JetSeedProducer7x7 = l1tPhase1JetSeedProducer.clone(
#     jetIEtaSize = 7,
#     jetIPhiSize = 7,
#     trimmedGrid = False,
#     outputCollectionName = "histoJetSeeds7x77x7",
#     seedSize = 7
#     )
# """ 7x7 HSC jet producer for seeds with 7x7 summing """
# l1t7x7HistoSeedsSCPFL1PuppiCorrectedEmulator7x7 = l1tSeedConePFJetEmulatorProducer.clone(
#     useExternalSeeds = True,
#     JetSeeds = ('l1t7x7Phase1JetSeedProducer7x7', 'histoJetSeeds7x77x7'),
#     nJets = 16,
#     L1PFObjects = 'l1tLayer2Deregionizer:Puppi',
#     doCorrections = cms.bool(False),
#     correctorFile = cms.string("L1Trigger/Phase2L1ParticleFlow/data/jecs/jecs_20220308.root"),
#     correctorDir = cms.string('L1PuppiSC4EmuJets')
#     )

# L1TPF7x7HistoSeedJetsTask7x7 = cms.Task(
#     l1t7x7Phase1JetSeedProducer7x7, l1tLayer2Deregionizer, l1t7x7HistoSeedsSCPFL1PuppiCorrectedEmulator7x7
# )


# """ 7x7 HSC jet seeds with 9x9 summing """
# l1t7x7Phase1JetSeedProducer9x9 = l1tPhase1JetSeedProducer.clone(
#     jetIEtaSize = 7,
#     jetIPhiSize = 7,
#     trimmedGrid = False,
#     outputCollectionName = "histoJetSeeds7x79x9",
#     seedSize = 9
#     )
# """ 7x7 HSC jet producer for seeds with 9x9 summing """
# l1t7x7HistoSeedsSCPFL1PuppiCorrectedEmulator9x9 = l1tSeedConePFJetEmulatorProducer.clone(
#     useExternalSeeds = True,
#     JetSeeds = ('l1t7x7Phase1JetSeedProducer9x9', 'histoJetSeeds7x79x9'),
#     nJets = 16,
#     L1PFObjects = 'l1tLayer2Deregionizer:Puppi',
#     doCorrections = cms.bool(False),
#     correctorFile = cms.string("L1Trigger/Phase2L1ParticleFlow/data/jecs/jecs_20220308.root"),
#     correctorDir = cms.string('L1PuppiSC4EmuJets')
#     )

# L1TPF7x7HistoSeedJetsTask9x9 = cms.Task(
#     l1t7x7Phase1JetSeedProducer9x9, l1tLayer2Deregionizer, l1t7x7HistoSeedsSCPFL1PuppiCorrectedEmulator9x9
# )




# # def HSC_jet_maker(maskSize=9, trimmed=False, numJets=16, seedSumMaskSize=1):
# #     assert maskSize % 2 = 0
# #     assert seedSumMaskSize % 2 = 0

# #     outputCollName = "histoJetSeeds{}x{}".format(str(maskSize), str(maskSize))
# #     NxN = str(maskSize)+ "+" +str(maskSize)
# #     MxM = str(seedSumMaskSize) + "+" + str(seedSumMaskSize)

# #     seedVarName = "l1tPhase1JetSeedProducer{}_{}".format(NxN, MxM)
# #     producerVarName = "l1t{}HistoSeedsSCPFL1PuppiCorrectedEmulator_{}".format(NxN, MxM)
# #     taskVarName = "L1TPF{}HistoSeedJetsTask{}".format(NxN, MxM)

# #     seedProducer = {
# #         seedVarName: l1tPhase1JetSeedProducer.clone(
# #                                                     jetIEtaSize = maskSize,
# #                                                     jetIPhiSize = maskSize,
# #                                                     trimmedGrid = trimmed,
# #                                                     outputCollectionName = outputCollName,
# #                                                     seedSize = seedSumMaskSize,
# #                                                     seedPtThreshold = 1
# #                                                     )
# #     }
# #     locals().update(seedProducer)

# #     jetProducer = {
# #         producerVarName: l1tSeedConePFJetEmulatorProducer.clone(
# #                                                                 useExternalSeeds = True,
# #                                                                 JetSeeds = (seedVarName, outputCollName),
# #                                                                 nJets = numJets,
# #                                                                 L1PFObjects = 'l1tLayer2Deregionizer:Puppi',
# #                                                                 doCorrections = cms.bool(True),
# #                                                                 correctorFile = cms.string("L1Trigger/Phase2L1ParticleFlow/data/jecs/jecs_20220308.root"),
# #                                                                 correctorDir = cms.string('L1PuppiSC4EmuJets')
# #                                                                 )
# #     }
# #     locals().update(jetProducer)

# #     task = {
# #         taskVarName: cms.Task(
# #                             l1tPhase1JetSeedProducer{}.format(NxN), l1tLayer2Deregionizer, l1t{}HistoSeedsSCPFL1PuppiCorrectedEmulator.format(NxN)
# #                             )
# #     }
# #     globals().update(task)

# #     return 0