import FWCore.ParameterSet.Config as cms

from Configuration.AlCa.GlobalTag import GlobalTag
from Configuration.StandardSequences.Eras import eras

process = cms.Process('PLOT', eras.Run2_2016)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag = GlobalTag(process.GlobalTag, '90X_dataRun2_HLT_Candidate_HCAL_payloads_tests_2017_03_24', '')

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring('/store/data/Commissioning2017/HcalNZS/RAW/v1/000/290/198/00000/0E82AAF0-8C13-E711-8328-02163E0138F1.root'))

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.RawToDigi_Data_cff")
process.load("SimCalorimetry.Configuration.hcalDigiSequence_cff")
process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')

process.TFileService = cms.Service("TFileService",
                                   closeFileFast=cms.untracked.bool(True),
                                   fileName=cms.string('analyze_290198.root'))

process.emulTP = process.simHcalTriggerPrimitiveDigis.clone()
process.emulTP.upgradeHF = cms.bool(True)
process.emulTP.upgradeHE = cms.bool(True)
process.emulTP.inputLabel = cms.VInputTag("hcalDigis", "hcalDigis")
process.emulTP.inputUpgradeLabel = cms.VInputTag("hcalDigis", "hcalDigis")
process.emulTP.numberOfSamples = cms.int32(3)
process.emulTP.numberOfPresamples = cms.int32(1)

# process.hcalDigis.InputLabel = cms.InputTag("source")
process.analyzeRAW = cms.EDAnalyzer("AnalyzeTP",
                                    triggerPrimitives=cms.InputTag("hcalDigis", "", ""))
process.analyzeSIM = cms.EDAnalyzer("AnalyzeTP",
                                    triggerPrimitives=cms.InputTag("emulTP", "", ""))
process.compare = cms.EDAnalyzer("CompareTP",
                                 triggerPrimitives=cms.InputTag("hcalDigis"),
                                 emulTriggerPrimitives=cms.InputTag("emulTP"),
                                 swapIphi=cms.bool(False))

process.emulTP2016 = process.simHcalTriggerPrimitiveDigis.clone()
process.emulTP2016.upgradeHF = cms.bool(False)
process.emulTP2016.upgradeHE = cms.bool(True)
process.emulTP2016.inputLabel = cms.VInputTag("hcalDigis", "hcalDigis")
process.emulTP2016.inputUpgradeLabel = cms.VInputTag("hcalDigis", "hcalDigis")
process.emulTP2016.numberOfSamples = cms.int32(3)
process.emulTP2016.numberOfPresamples = cms.int32(1)

process.compare2016 = cms.EDAnalyzer("CompareTP",
                                     triggerPrimitives=cms.InputTag("hcalDigis"),
                                     emulTriggerPrimitives=cms.InputTag("emulTP2016"),
                                     swapIphi=cms.bool(False))

process.dump = cms.EDAnalyzer("EventContentAnalyzer")

process.p = cms.Path(
    process.hcalDigis *
    # process.dump *
    process.emulTP *
    process.analyzeRAW *
    process.analyzeSIM *
    process.compare *
    process.emulTP2016 *
    process.compare2016)

# print process.dumpPython()
