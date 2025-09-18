import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
process = cms.Process("OWNPARTICLES")


process.load('Configuration.StandardSequences.Services_cff')
process.load('JetMETCorrections.Configuration.JetCorrectors_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('RecoTracker.Configuration.RecoTracker_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
from Configuration.AlCa.autoCond import autoCond
process.GlobalTag.globaltag = cms.string("92X_upgrade2017_realistic_v7")#("92X_dataRun2_Prompt_v8")#("80X_dataRun2_Prompt_ICHEP16JEC_v0")#("76X_dataRun2_16Dec2015_v0")#( autoCond[ 'run2_data' ] )

#Jet Tracks Association
process.load('RecoJets.JetAssociationProducers.ak4JTA_cff')


#process.load("CondCore.DBCommon.CondDBCommon_cfi")
#from CondCore.DBCommon.CondDBSetup_cfi import *
process.load("CondCore.CondDB.CondDB_cfi")
from CondCore.DBCommon.CondDBSetup_cfi import *
 
process.jec = cms.ESSource("PoolDBESSource",
      DBParameters = cms.PSet(
        messageLevel = cms.untracked.int32(0)
        ),
      timetype = cms.string('runnumber'),
      toGet = cms.VPSet(
      cms.PSet(
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string('JetCorrectorParametersCollection_Spring16_23Sep2016AllV2_DATA_AK4Calo'),
            # tag    = cms.string('JetCorrectorParametersCollection_Fall15_25nsV2_MC_AK4PFchs'),
            label  = cms.untracked.string('AK4Calo')
            ),
     # cms.PSet(
     #       record = cms.string('JetCorrectionsRecord'),
     #       tag    = cms.string('JetCorrectorParametersCollection_Fall15_V2_DATA_AK8PFPuppi'),
     #       # tag    = cms.string('JetCorrectorParametersCollection_Fall15_25nsV2_MC_AK8PFPuppi'),
     #       label  = cms.untracked.string('AK8PFPuppi')
     #       ),
      ## here you add as many jet types as you need
      ## note that the tag name is specific for the particular sqlite file 
      ),
      connect = cms.string('sqlite:python/Spring16_23Sep2016AllV2_DATA.db')#('sqlite:/uscms_data/d3/jluo/work/JEC/Spring16_25nsV6_MC.db')
     # uncomment above tag lines and this comment to use MC JEC
     # connect = cms.string('sqlite:Fall15_25nsV2_DATA.db')
)

## add an es_prefer statement to resolve a possible conflict from simultaneous connection to a global tag
process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')


process.load('JetMETCorrections.Configuration.JetCorrectors_cff')
process.load('JetMETCorrections.Configuration.CorrectedJetProducers_cff')
from JetMETCorrections.Configuration.CorrectedJetProducers_cff import *


#PAT Trigger Info
from PhysicsTools.PatAlgos.patTemplate_cfg import *

process.out = cms.OutputModule("PoolOutputModule", 
    fileName = cms.untracked.string("pat.root"),
    outputCommands = cms.untracked.vstring('drop *' )
)
# Switch on "unscheduled" mode
## Options and Output Report
process.options.allowUnscheduled = cms.untracked.bool( True )
#process.Tracer = cms.Service( "Tracer" )

# Load default PAT
process.load( "PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff" )
process.load( "PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff" )
process.load( "PhysicsTools.PatAlgos.patSequences_cff")
#process.p = cms.Path(
#    process.selectedPatCandidates
#    )

### Get PAT trigger tools
from PhysicsTools.PatAlgos.tools.trigTools import *

# ------------------------------------------------------------------------------
# Depending on the purpose, comment/uncomment the following sections
# ------------------------------------------------------------------------------

# Add full trigger information
switchOnTrigger( process ) 

from PhysicsTools.PatAlgos.tools.coreTools import removeMCMatching
#removeMCMatching( process, ['All'])
# Add stand-alone trigger information
#switchOnTriggerStandAlone( process )

process.IsoMufilter = cms.EDFilter("HLTHighLevel", 
    TriggerResultsTag = cms.InputTag("TriggerResults", "", "HLT"),
    HLTPaths = cms.vstring("HLT_IsoMu27_v*"),
    eventSetupPathsKey = cms.string(""),
    andOr = cms.bool(True),
    throw = cms.bool(False)
)

process.TriggerNtuple = cms.EDAnalyzer('TriggerAnl'
     , trigHTFilter = cms.InputTag('hltHT430', '', 'HLT')#('hltCaloJetCollection40Filter', '', 'HLT' )
     , trigLowHTFilter = cms.InputTag('hltHT400', '', 'HLT')
     , trigJetFilter = cms.InputTag('hltDoubleCentralCaloJetpt40', '', 'HLT')
     , trigPromptJetFilter = cms.InputTag('hltL4PromptDisplacedDijetFullTracksHLTCaloJetTagFilterLowPt', '', 'HLT')
     , trigDisplacedJetFilter = cms.InputTag('hltL4DisplacedDijetFullTracksHLTCaloJetTagFilterLowPt', '', 'HLT') 
     , trigSummary = cms.InputTag('hltTriggerSummaryAOD', '', 'HLT')
     , jets = cms.untracked.InputTag('ak4CaloJets')
     , primaryVertices = cms.untracked.InputTag('offlinePrimaryVerticesWithBS')
     , jetCorr = cms.untracked.InputTag('ak4CaloL1FastL2L3ResidualCorrector')
     , triggerEvent = cms.untracked.InputTag('patTriggerEvent')
     , associatorVTX = cms.untracked.InputTag("ak4JTAatVX")
     , JECTag = cms.untracked.string("AK4Calo")
)

process.ak4JTAatVX = cms.EDProducer("JetTracksAssociatorAtVertex",
        tracks = cms.InputTag("generalTracks"),
        coneSize = cms.double(0.4),
        useAssigned = cms.bool(False),
        pvSrc = cms.InputTag("offlinePrimaryVerticesWithBS"),
        jets = cms.InputTag("ak4CaloJets")
)
   
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.TFileService = cms.Service("TFileService", 
        fileName = cms.string("test_output.root"),
)


process.MessageLogger = cms.Service("MessageLogger",
     destinations   = cms.untracked.vstring('messages'),
     debugModules   = cms.untracked.vstring('*'),
     messages       = cms.untracked.PSet(
                                 threshold =  cms.untracked.string('ERROR')
                                                    )

)
#process.ak4JTAatVX = cms.EDProducer("JetTracksAssociatorAtVertex",
#        tracks = cms.InputTag("generalTracks"),
#        coneSize = cms.double(0.4),
#        useAssigned = cms.bool(False),
#        pvSrc = cms.InputTag("offlinePrimaryVerticesWithBS"),
#        jets = cms.InputTag("ak4CaloJets")
#)
#process.ak4JTAatCAL = cms.EDProducer("JetTracksAssociatorAtCaloFace",
#        tracks = cms.InputTag("generalTracks"),
#        extrapolations = cms.InputTag("trackExtrapolator"),
#        trackQuality = cms.string("goodIterative"),
#        coneSize = cms.double(0.4),
#        useAssigned = cms.bool(False),
#        jets = cms.InputTag("ak4CaloJets")
#)
process.p = cms.Path(
        process.IsoMufilter*
        ak4CaloL1FastL2L3ResidualCorrectorChain*
        ak4CaloL1FastL2L3ResidualCorrector*
        process.ak4JTAatVX*
        process.TriggerNtuple
)

process.source = cms.Source('PoolSource', fileNames = cms.untracked.vstring('root://cmsxrootd.fnal.gov//store/data/Run2016G/DisplacedJet/AOD/23Sep2016-v1/100000/00CD619A-9491-E611-8F22-047D7BD6DE1E.root'))#('root://cmsxrootd.fnal.gov//store/data/Run2017C/SingleMuon/AOD/PromptReco-v1/000/299/368/00000/20C5E3B7-826D-E711-BF78-02163E014372.root'))#('root://cmsxrootd.fnal.gov//store/data/Run2017C/DisplacedJet/AOD/PromptReco-v1/000/299/368/00000/581E4AD0-956D-E711-B522-02163E0137FC.root'))#'root://cmsxrootd.fnal.gov//store/mc/RunIISpring15DR74/QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v2/70000/AE33914D-1606-E511-A826-0025905A48D0.root'))
