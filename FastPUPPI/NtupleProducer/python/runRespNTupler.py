import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process("RESP", eras.Phase2C9)

process.load('Configuration.StandardSequences.Services_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True), allowUnscheduled = cms.untracked.bool(False) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:inputs110X.root'),
    inputCommands = cms.untracked.vstring("keep *", 
            "drop l1tPFClusters_*_*_*",
            "drop l1tPFTracks_*_*_*",
            "drop l1tPFCandidates_*_*_*")
)

process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D49_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff') # needed to read HCal TPs
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '111X_mcRun4_realistic_Candidate_2020_12_09_15_46_46', '')

process.load("L1Trigger.Phase2L1ParticleFlow.l1ParticleFlow_cff")
process.load("L1Trigger.Phase2L1ParticleFlow.pfClustersFromHGC3DClustersEM_cfi")

process.pfClustersFromL1EGClustersRaw    = process.pfClustersFromL1EGClusters.clone(corrector = "")
process.pfClustersFromHGC3DClustersRaw   = process.pfClustersFromHGC3DClusters.clone(corrector = "")
process.pfClustersFromHGC3DClustersEMRaw = process.pfClustersFromHGC3DClustersEM.clone(corrector = "")
process.pfClustersFromCombinedCaloHFOnly = process.pfClustersFromCombinedCaloHF.clone(hcalCandidates =[])


process.pfClustersFromCombinedCaloHCalUnclust = process.pfClustersFromCombinedCaloHCal.clone(
    ecalCandidates = []
)
process.runPF = cms.Sequence( 
    process.l1ParticleFlow_proper + # excludes the prerequisites (3D clusters and L1EG clusters)
    process.pfClustersFromL1EGClustersRaw +
    process.pfClustersFromHGC3DClustersRaw +
    process.pfClustersFromHGC3DClustersEMRaw
    + process.pfClustersFromCombinedCaloHCalUnclust
    + process.pfClustersFromHGC3DClustersEM 
    + process.pfClustersFromCombinedCaloHFOnly
    + process.pfTracksFromL1Tracks
)


process.ntuple = cms.EDAnalyzer("ResponseNTuplizer",
    genJets = cms.InputTag("ak4GenJetsNoNu"),
    genParticles = cms.InputTag("genParticles"),
    isParticleGun = cms.bool(False),
    doRandom = cms.bool(False),
    objects = cms.PSet(
        # -- inputs and PF --
        RawTK  = cms.VInputTag('pfTracksFromL1Tracks',),
        # for calibrations
        L1RawBarrelEcal   = cms.VInputTag('pfClustersFromL1EGClustersRaw' ),
        L1RawBarrelCalo   = cms.VInputTag('pfClustersFromCombinedCaloHCal:uncalibrated'),
        L1RawBarrelCaloEM = cms.VInputTag('pfClustersFromCombinedCaloHCal:emUncalibrated'),
        L1BareHGCal   = cms.VInputTag(cms.InputTag("hgcalBackEndLayer2Producer","HGCalBackendLayer2Processor3DClustering")),
        L1RawHGCal   = cms.VInputTag('pfClustersFromHGC3DClustersRaw'),
        L1RawHGCalEM = cms.VInputTag('pfClustersFromHGC3DClustersEMRaw'),
        L1RawHFCalo  = cms.VInputTag('pfClustersFromCombinedCaloHF:uncalibrated'),
        L1RawHFOnlyCalo  = cms.VInputTag('pfClustersFromCombinedCaloHFOnly:uncalibrated'),
        L1RawHFOnlyCells = cms.VInputTag('pfClustersFromCombinedCaloHFOnly:hcalCells'),
        L1BarrelEcal = cms.VInputTag('pfClustersFromL1EGClusters' ),
        L1BarrelCalo = cms.VInputTag('pfClustersFromCombinedCaloHCal:calibrated'),
        L1HGCal   = cms.VInputTag('pfClustersFromHGC3DClusters'),
        L1HGCalEM = cms.VInputTag('pfClustersFromHGC3DClustersEM', ),
        L1HFCalo  = cms.VInputTag('pfClustersFromCombinedCaloHF:calibrated'),
        # alternate: don't use L1EG clusters when making hadrons
        L1RawBarrelCaloUnclust   = cms.VInputTag('pfClustersFromCombinedCaloHCalUnclust:uncalibrated'),
        L1RawBarrelCaloEMUnclust = cms.VInputTag('pfClustersFromCombinedCaloHCalUnclust:emUncalibrated'),
        L1BarrelCaloUnclust = cms.VInputTag('pfClustersFromCombinedCaloHCalUnclust:calibrated'),
        # outputs
        L1Calo = cms.VInputTag("l1pfCandidates:Calo",),
        L1TK = cms.VInputTag("l1pfCandidates:TK",),
        L1TKV = cms.VInputTag("l1pfCandidates:TKVtx",),
        L1TKV5 = cms.VInputTag("l1pfCandidates:TKVtx",),
        L1TKV5_sel = cms.string("pfTrack.nStubs >= 5"),
        L1PF = cms.VInputTag("l1pfCandidates:PF",),
        L1Puppi = cms.VInputTag("l1pfCandidates:Puppi",),
    ),
    copyUInts = cms.VInputTag(),
    copyFloats = cms.VInputTag(),
    copyVecUInts = cms.VInputTag(),
)

process.p = cms.Path(process.runPF + process.ntuple)
process.TFileService = cms.Service("TFileService", fileName = cms.string("respTupleNew.root"))

# Below for more debugging
if True:
    process.genInAcceptance = cms.EDFilter("GenParticleSelector",
        src = cms.InputTag("genParticles"),
        cut = cms.string("status == 1 && (abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16) && "+
                         "(abs(eta) < 2.5 && pt > 2 && charge != 0 || "+
                         "abs(pdgId) == 22 && pt > 1 || "+
                         "charge == 0 && pt > 1 || "+
                         "charge != 0 && abs(eta) > 2.5 && pt > 2) ") # tracks below pT 2 bend by more than 0.4,
    )
    process.ntuple.objects.GenAcc = cms.VInputTag(cms.InputTag("genInAcceptance"))
    process.chGenInAcceptance = cms.EDFilter("GenParticleSelector",
        src = cms.InputTag("genParticles"),
        cut = cms.string("status == 1 && (abs(eta) < 2.5 && pt > 2 && charge != 0)")
    )
    process.phGenInAcceptance = cms.EDFilter("GenParticleSelector",
        src = cms.InputTag("genParticles"),
        cut = cms.string("status == 1 && pt > 1 && pdgId == 22")
    )
    process.ntuple.objects.ChGenAcc = cms.VInputTag(cms.InputTag("chGenInAcceptance"))
    process.ntuple.objects.PhGenAcc = cms.VInputTag(cms.InputTag("phGenInAcceptance"))
    process.p = cms.Path(process.genInAcceptance + process.chGenInAcceptance + process.phGenInAcceptance + process.p._seq)

    process.ntuple.objects.L1PFCharged = cms.VInputTag("l1pfCandidates:PF",)
    process.ntuple.objects.L1PFCharged_sel = cms.string("charge != 0")
    process.ntuple.objects.L1PFPhoton = cms.VInputTag("l1pfCandidates:PF",)
    process.ntuple.objects.L1PFPhoton_sel = cms.string("pdgId == 22")
    process.ntuple.objects.L1PFNeutral = cms.VInputTag("l1pfCandidates:PF",)
    process.ntuple.objects.L1PFNeutral_sel = cms.string("charge == 0")
    process.ntuple.objects.L1PFNeutral = cms.VInputTag("l1pfCandidates:PF",)
    process.ntuple.objects.L1PFNeutral_sel = cms.string("charge == 0 && pdgId != 22")
    process.ntuple.objects.L1PFMuon = cms.VInputTag("l1pfCandidates:PF",)
    process.ntuple.objects.L1PFMuon_sel = cms.string("abs(pdgId) == 13")
    process.ntuple.objects.L1PFElectron = cms.VInputTag("l1pfCandidates:PF",)
    process.ntuple.objects.L1PFElectron_sel = cms.string("abs(pdgId) == 11")
    process.ntuple.objects.L1PuppiCharged = cms.VInputTag("l1pfCandidates:Puppi",)
    process.ntuple.objects.L1PuppiCharged_sel = cms.string("charge != 0")
    process.ntuple.objects.L1PuppiPhoton = cms.VInputTag("l1pfCandidates:Puppi",)
    process.ntuple.objects.L1PuppiPhoton_sel = cms.string("pdgId == 22")
    process.ntuple.objects.L1PuppiNeutral = cms.VInputTag("l1pfCandidates:Puppi",)
    process.ntuple.objects.L1PuppiNeutral_sel = cms.string("charge == 0")
    process.ntuple.objects.L1PuppiNeutralHad = cms.VInputTag("l1pfCandidates:Puppi",)
    process.ntuple.objects.L1PuppiNeutralHad_sel = cms.string("charge == 0 && pdgId != 22")


def goGun():
    process.ntuple.isParticleGun = True
def goRandom():
    process.ntuple.doRandom = True
def goMT(nthreads=2):
    process.options.numberOfThreads = cms.untracked.uint32(nthreads)
    process.options.numberOfStreams = cms.untracked.uint32(0)
def goVerbose(v=3, point=None, R=0.7):
    process.l1pfProducer.debug = v
    if point:
        process.l1pfProducer.debugEta = cms.untracked.double(point[0])
        process.l1pfProducer.debugPhi = cms.untracked.double(point[1])
        process.l1pfProducer.debugR   = cms.untracked.double(R)
def dumpGen():
    process.dumpGen = cms.EDAnalyzer("ParticleListDrawer",
        maxEventsToPrint = cms.untracked.int32(100),
        printVertex = cms.untracked.bool(False),
        printOnlyHardInteraction = cms.untracked.bool(False), 
        src = cms.InputTag("genParticles")
    )
    process.p.replace(process.ntuple, process.dumpGen+process.ntuple)
    process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10))
    process.MessageLogger.cerr.FwkReport.reportEvery = 1
def saveOut():
    process.out = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string("debugPF.root"),
        SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("p"))
    )
    process.end = cms.EndPath(process.out)
def hgcAcc(pdgId,pt=10,eta=(1.6,2.6),prompt=False):
    if type(pdgId) == int: pdgId = [ pdgId ]
    pdgIdCut = "||".join("abs(pdgId) == %d" % p for p in pdgId)
    process.acceptance = cms.EDFilter("GenParticleSelector",
        src = cms.InputTag("genParticles"),
        cut = cms.string("(%s) && pt > %g && %g < abs(eta) && abs(eta) < %g%s" % (pdgIdCut, pt, eta[0], eta[1], "&& statusFlags.isPrompt" if prompt else "")),
        filter = cms.bool(True),
    )
    process.p.insert(0, process.acceptance)

def noPU():
    for X in "", "EM", "Raw", "EMRaw":
        pfc = getattr(process, 'pfClustersFromHGC3DClusters'+X)
        pfc.emVsPUID.wp = "-1.0"
