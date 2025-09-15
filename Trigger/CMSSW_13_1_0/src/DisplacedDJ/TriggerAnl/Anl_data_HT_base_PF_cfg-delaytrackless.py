import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
from Configuration.Eras.Modifier_stage2L1Trigger_cff import stage2L1Trigger

process = cms.Process("OWNPARTICLES", stage2L1Trigger)


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
process.GlobalTag.globaltag = autoCond["run3_data"]#cms.string("92X_dataRun2_Prompt_v8")#("80X_dataRun2_Prompt_ICHEP16JEC_v0")#("76X_dataRun2_16Dec2015_v0")#( autoCond[ 'run2_data' ] )

#Jet Tracks Association
process.load('RecoJets.JetAssociationProducers.ak4JTA_cff')


#process.load("CondCore.DBCommon.CondDBCommon_cfi")
#from CondCore.DBCommon.CondDBSetup_cfi import *
#process.load("CondCore.CondDB.CondDB_cfi")
#from CondCore.DBCommon.CondDBSetup_cfi import *
# 
#process.jec = cms.ESSource("PoolDBESSource",
#      DBParameters = cms.PSet(
#        messageLevel = cms.untracked.int32(0)
#        ),
#      timetype = cms.string('runnumber'),
#      toGet = cms.VPSet(
#      cms.PSet(
#            record = cms.string('JetCorrectionsRecord'),
#            tag    = cms.string('JetCorrectorParametersCollection_Spring16_23Sep2016AllV2_DATA_AK4Calo'),
#            # tag    = cms.string('JetCorrectorParametersCollection_Fall15_25nsV2_MC_AK4PFchs'),
#            label  = cms.untracked.string('AK4Calo')
#            ),
#     # cms.PSet(
#     #       record = cms.string('JetCorrectionsRecord'),
#     #       tag    = cms.string('JetCorrectorParametersCollection_Fall15_V2_DATA_AK8PFPuppi'),
#     #       # tag    = cms.string('JetCorrectorParametersCollection_Fall15_25nsV2_MC_AK8PFPuppi'),
#     #       label  = cms.untracked.string('AK8PFPuppi')
#     #       ),
#      ## here you add as many jet types as you need
#      ## note that the tag name is specific for the particular sqlite file 
#      ),
#      connect = cms.string('sqlite:python/Spring16_23Sep2016AllV2_DATA.db')#('sqlite:/uscms_data/d3/jluo/work/JEC/Spring16_25nsV6_MC.db')
#     # uncomment above tag lines and this comment to use MC JEC
#     # connect = cms.string('sqlite:Fall15_25nsV2_DATA.db')
#)
#
### add an es_prefer statement to resolve a possible conflict from simultaneous connection to a global tag
#process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')


process.load('JetMETCorrections.Configuration.JetCorrectors_cff')
process.load('JetMETCorrections.Configuration.CorrectedJetProducers_cff')
from JetMETCorrections.Configuration.CorrectedJetProducers_cff import *



#PAT Trigger Info
#from PhysicsTools.PatAlgos.patTemplate_cfg import *

process.out = cms.OutputModule("PoolOutputModule", 
    fileName = cms.untracked.string("pat.root"),
    outputCommands = cms.untracked.vstring('drop *' )
)

#from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning
#process.out = cms.OutputModule("PoolOutputModule",
#                               fileName = cms.untracked.string('patTuple.root'),
#                               ## save only events passing the full path
#                               #SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
#                               ## save PAT output; you need a '*' to unpack the list of commands
#                               ## 'patEventContent'
#                               outputCommands = cms.untracked.vstring('drop *', *patEventContentNoCleaning )
#                               )

#from PhysicsTools.PatAlgos.tools.helpers import getPatAlgosToolsTask
#patAlgosToolsTask = getPatAlgosToolsTask(process)
#process.outpath = cms.EndPath(process.out, patAlgosToolsTask)
## Switch on "unscheduled" mode
### Options and Output Report
#process.options.allowUnscheduled = cms.untracked.bool( True )
##process.Tracer = cms.Service( "Tracer" )
#
## Load default PAT
#process.load( "PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff" )
#process.load( "PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff" )
#process.load( "PhysicsTools.PatAlgos.patSequences_cff")
##process.p = cms.Path(
##    process.selectedPatCandidates
##    )
#
#### Get PAT trigger tools
#from PhysicsTools.PatAlgos.tools.trigTools import *
#
## ------------------------------------------------------------------------------
## Depending on the purpose, comment/uncomment the following sections
## ------------------------------------------------------------------------------
#
## Add full trigger information
#switchOnTrigger( process, outputModule="") 
#
#from PhysicsTools.PatAlgos.tools.coreTools import removeMCMatching
#removeMCMatching( process, ['All'])
# Add stand-alone trigger information
#switchOnTriggerStandAlone( process )


process.IsoMufilter = cms.EDFilter("HLTHighLevel", 
    TriggerResultsTag = cms.InputTag("TriggerResults", "", "HLT"),
    HLTPaths = cms.vstring("HLT_HT200_L1SingleLLPJet_DelayedJet40_DoubleDelay0p5nsTrackless_v*"),
    eventSetupPathsKey = cms.string(""),
    andOr = cms.bool(True),
    throw = cms.bool(False)
)

process.TriggerNtuple = cms.EDAnalyzer('TriggerAnlPF_miniAOD'
     , tracks = cms.untracked.InputTag('packedPFCandidates')
     , losttracks = cms.untracked.InputTag('lostTracks', '', "RECO")
     , trigHTFilter = cms.InputTag('hltHT240', '', 'HLT')#('hltCaloJetCollection40Filter', '', 'HLT' )
     , trigJetFilter = cms.InputTag('hltDoubleCentralCaloJetpt40', '', 'HLT')
     , trigPromptJetFilter = cms.InputTag('hltL4PromptDisplacedDijetFullTracksHLTCaloJetTagFilterLowPt', '', 'HLT')
     , trigDisplacedJetFilter = cms.InputTag('hltL4DisplacedDijetFullTracksHLTCaloJetTagFilterLooseLowPt1PtrkPt0p5ShortSig5', '', 'HLT')#hltDisplacedHLTCaloJetCollectionProducerLowPt1PtrkShortSig5', '', 'HLT')
     , bits = cms.InputTag("TriggerResults", "", "HLT")
     , objects = cms.InputTag("slimmedPatTrigger")
     #, trigSummary = cms.InputTag('hltTriggerSummaryAOD', '', 'HLT')
     , primaryVertices = cms.untracked.InputTag('offlineSlimmedPrimaryVertices')
     , jets = cms.untracked.InputTag('slimmedJets')
     , jetCorr = cms.untracked.InputTag('ak4PFCHSL1FastL2L3ResidualCorrector')
     #, triggerEvent = cms.untracked.InputTag('patTriggerEvent')
     , associatorVTX = cms.untracked.InputTag("ak4JTAatVX")
     #, JECTag = cms.untracked.string("AK4Calo")
)



    
process.TFileService = cms.Service("TFileService", 
        fileName = cms.string("test_output.root"),
)


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.MessageLogger = cms.Service("MessageLogger",
#     destinations   = cms.untracked.vstring('messages'),
#     debugModules   = cms.untracked.vstring('*'),
#     messages       = cms.untracked.PSet(
#                                 threshold =  cms.untracked.string('ERROR')
#                                                    )
#
#)
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

#process.out = cms.OutputModule("")

process.p = cms.Path(
        process.IsoMufilter*
        ak4PFCHSL1FastL2L3ResidualCorrectorChain*
        ak4PFCHSL1FastL2L3ResidualCorrector*
        #process.ak4JTAatVX*
        #process.patDefaultSequence*
        process.TriggerNtuple
)


process.source = cms.Source('PoolSource', fileNames = cms.untracked.vstring('root://cmsxrootd.fnal.gov//store/data/Run2023C/DisplacedJet/MINIAOD/PromptReco-v1/000/367/095/00000/1a2378fe-1d64-4c64-b067-984c53c02d19.root'))
