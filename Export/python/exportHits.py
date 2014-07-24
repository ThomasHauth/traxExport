import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("ExportTracks")
#keep the logging output to a nice level
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('RecoLocalTracker.Configuration.RecoLocalTracker_cff')
process.load('RecoTracker.TrackProducer.TrackRefitters_cff')

# setup 'analysis'  options
options = VarParsing.VarParsing ('analysis')

# setup any defaults you want
options.outputFile = 'mcTrackAnalysis.root'
options.inputFiles = 'file:RECO_DIGI_L1_DIGI2RAW_RAW2DIGI_RECO.root'
options.maxEvents = -1 # -1 means all events

# get and parse the command line arguments
options.parseArguments()

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:mc', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

process.load('SimGeneral.TrackingAnalysis.trackingParticles_cfi')
process.load('SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi')


process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(options.inputFiles)
)

process.TFileService = cms.Service("TFileService",
  fileName = cms.string(options.outputFile)
)

process.analyzeMCTracks = cms.EDAnalyzer("HitExporter",
  recoTracks = cms.untracked.InputTag("TrackRefitter"), #generalTracks has no position information
  mcTruth = cms.untracked.InputTag("mergedtruth","MergedTrackTruth"),
  matchedRecHits = cms.untracked.InputTag("siStripMatchedRecHits", 'matchedRecHit'),
  minPt = cms.untracked.double(0.3),
  hitFile = cms.untracked.string("hits.pb")
  
)

process.p = cms.Path(process.mix + process.siPixelRecHits + process.siStripMatchedRecHits + process.mergedtruth + 
			process.TrackRefitter + process.analyzeMCTracks) # + process.InitialStep)

