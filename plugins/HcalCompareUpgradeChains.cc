//#if 0
// -*- C++ -*-
//
// Package:    HcalCompareUpgradeChains
// Class:      HcalCompareUpgradeChains
// 
/**\class HcalCompareUpgradeChains HcalCompareUpgradeChains.cc HcalDebug/CompareChans/src/HcalCompareUpgradeChains.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  matthias wolf
//         Created:  Fri Aug 26 11:37:21 CDT 2013
// $Id$
//
// Updates by: georgia karapostoli [2019]

// system include files
#include <memory>
#include <array>
#include <vector>
#include <algorithm>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CalibFormats/CaloTPG/interface/CaloTPGTranscoder.h"
#include "CalibFormats/CaloTPG/interface/CaloTPGRecord.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"

#include "CondFormats/DataRecord/interface/HcalChannelQualityRcd.h"
#include "CondFormats/DataRecord/interface/L1CaloGeometryRecord.h"
#include "CondFormats/HcalObjects/interface/HcalChannelQuality.h"
#include "CondFormats/L1TObjects/interface/L1CaloGeometry.h"

#include "CondFormats/L1TObjects/interface/L1RCTParameters.h"
#include "CondFormats/DataRecord/interface/L1RCTParametersRcd.h"
#include "CondFormats/L1TObjects/interface/L1CaloHcalScale.h"
#include "CondFormats/DataRecord/interface/L1CaloHcalScaleRcd.h"

#include "DataFormats/Common/interface/SortedCollection.h"
#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalTriggerPrimitiveDigi.h"
#include "DataFormats/HcalDigi/interface/HcalUpgradeTriggerPrimitiveDigi.h"  
#include "DataFormats/HcalDetId/interface/HcalTrigTowerDetId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "DataFormats/HcalRecHit/interface/HFRecHit.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/HcalTowerAlgo/interface/HcalGeometry.h"
#include "Geometry/HcalTowerAlgo/interface/HcalTrigTowerGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "RecoLocalCalo/HcalRecAlgos/interface/HcalSeverityLevelComputer.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalSeverityLevelComputerRcd.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TTree.h"

class HcalCompareUpgradeChains : public edm::EDAnalyzer {
   public:
      explicit HcalCompareUpgradeChains(const edm::ParameterSet&);
      ~HcalCompareUpgradeChains();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void beginLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&) override;

      double get_cosh(const HcalDetId&);

      // ----------member data ---------------------------
      bool first_;

      std::vector<edm::InputTag> frames_;
      edm::InputTag digis_;
      std::vector<edm::InputTag> rechits_;

      edm::ESHandle<CaloGeometry> gen_geo_; 
      edm::ESHandle<CaloTPGTranscoder> decoder_;
      edm::ESHandle<HcalTrigTowerGeometry> tpd_geo_h_;
      edm::ESHandle<HcalDbService> conditions_;

      //TH2D *df_multiplicity_;
      //TH2D *tp_multiplicity_;

      //TTree *tps_;
      //TTree *tpsplit_;
      TTree *events_;
      TTree *matches_;

      double tp_energy_;
      int tp_ieta_;
      int tp_iphi_;
      int tp_depth_max_;
      int tp_depth_start_;
      int tp_depth_end_;
      int tp_depth_;
      int tp_event_;
      std::vector<int> tp_ts_adc_;

      int tp_soi_;
      double tpsplit_energy_;
  
      double tpsplit_oot_;
      
      int tpsplit_ieta_;
      int tpsplit_iphi_;
      int tpsplit_depth_;
      double tpsplit_ettot_;
      int tpsplit_bx_;
      int tpsplit_event_;

      double tpsplit_rise_avg_;
      double tpsplit_rise_rms_;
      double tpsplit_fall_avg_;
      double tpsplit_fall_rms_;

      double ev_rh_energy0_;    
      double ev_rh_energy_;
      double ev_tp_energy_;
      int ev_rh_unmatched_;
      int ev_tp_unmatched_;
      int ev_tp_event_;
      std::vector<int> ev_tp_ieta_;
      std::vector<int> ev_tp_iphi_;
      std::vector<int> ev_tp_depth_;
      std::vector<double> ev_tp_et_;
      std::vector<int> ev_tp_soi_;
      std::vector<int> ev_tp_ts0_;
      std::vector<int> ev_tp_ts1_;
      std::vector<int> ev_tp_ts2_;
      std::vector<int> ev_tp_ts3_;
      std::vector<int> ev_tp_ts4_;
      std::vector<int> ev_tp_ts5_;
      std::vector<int> ev_tp_ts6_;
      std::vector<int> ev_tp_ts7_;

      double mt_rh_energy0_;
      double mt_rh_energy_;
      double mt_tp_energy_;
      std::vector<double> mt_rh_energy0_depth_;
      std::vector<double> mt_rh_energy_depth_;
      std::vector<double> mt_tp_energy_depth_;

      int mt_ieta_;
      int mt_iphi_;
      int mt_depth_;
      int mt_tp_soi_;
      int mt_bx_;
      int mt_event_;
      int mt_ts0_;
      int mt_ts1_;
      int mt_ts2_;
      int mt_ts3_;
      int mt_ts4_;
      int mt_ts5_;
      int mt_ts6_;
      int mt_ts7_;

      int mt_d0_;
      int mt_d1_;
      int mt_d2_;
      int mt_d3_;
      int mt_d4_;
      int mt_d5_;
      int mt_d6_;
      int mt_d7_;

      bool swap_iphi_;

      int max_severity_;
      const HcalChannelQuality* status_;
      const HcalSeverityLevelComputer* comp_;
};

HcalCompareUpgradeChains::HcalCompareUpgradeChains(const edm::ParameterSet& config) :
    edm::EDAnalyzer(),
    first_(true),
    frames_(config.getParameter<std::vector<edm::InputTag>>("dataFrames")),
    digis_(config.getParameter<edm::InputTag>("triggerPrimitives")),
    rechits_(config.getParameter<std::vector<edm::InputTag>>("recHits")),
    swap_iphi_(config.getParameter<bool>("swapIphi")),
    max_severity_(config.getParameter<int>("maxSeverity"))

{

    tp_ts_adc_.resize(8, 0);
    mt_rh_energy0_depth_.resize(8, 0.0);
    mt_rh_energy_depth_.resize(8, 0.0);
    mt_tp_energy_depth_.resize(8, 0.0);

    consumes<HcalUpgradeTrigPrimDigiCollection>(digis_);
    consumes<QIE11DigiCollection>(frames_[0]);
    consumes<QIE10DigiCollection>(frames_[1]);
    consumes<edm::SortedCollection<HBHERecHit>>(rechits_[0]);
    consumes<edm::SortedCollection<HFRecHit>>(rechits_[1]);

    edm::Service<TFileService> fs;

    //df_multiplicity_ = fs->make<TH2D>("df_multiplicity", "DataFrame multiplicity;ieta;iphi", 65, -32.5, 32.5, 72, 0.5, 72.5);
    //tp_multiplicity_ = fs->make<TH2D>("tp_multiplicity", "TrigPrim multiplicity;ieta;iphi", 65, -32.5, 32.5, 72, 0.5, 72.5);

    //tps_ = fs->make<TTree>("tps", "Trigger primitives");
    //tps_->Branch("et", &tp_energy_);
    //tps_->Branch("ieta", &tp_ieta_);
    //tps_->Branch("iphi", &tp_iphi_);
    //tps_->Branch("depth_max", &tp_depth_max_);
    //tps_->Branch("depth_start", &tp_depth_start_);
    //tps_->Branch("depth_end", &tp_depth_end_);
    //tps_->Branch("soi", &tp_soi_);
    //tps_->Branch("ts_adc", &tp_ts_adc_, 32000, 0);
    //tps_->Branch("depth", &tp_depth_);
    //tps_->Branch("event", &tp_event_);

    //tpsplit_ = fs->make<TTree>("tpsplit", "Trigger primitives");
    //tpsplit_->Branch("et", &tpsplit_energy_);
    //tpsplit_->Branch("oot", &tpsplit_oot_);
    //tpsplit_->Branch("ieta", &tpsplit_ieta_);
    //tpsplit_->Branch("iphi", &tpsplit_iphi_);
    //tpsplit_->Branch("depth", &tpsplit_depth_);
    //tpsplit_->Branch("etsum", &tpsplit_ettot_);
    //tpsplit_->Branch("soi", &tp_soi_);
    //tpsplit_->Branch("event", &tpsplit_event_);
    //tpsplit_->Branch("bx", &tpsplit_bx_);
    //tpsplit_->Branch("rise_avg", &tpsplit_rise_avg_);
    //tpsplit_->Branch("rise_rms", &tpsplit_rise_rms_);
    //tpsplit_->Branch("fall_avg", &tpsplit_fall_avg_);
    //tpsplit_->Branch("fall_rms", &tpsplit_fall_rms_);

    events_ = fs->make<TTree>("events", "Event quantities");
    events_->Branch("RH_energy", &ev_rh_energy_);
    events_->Branch("TP_energy", &ev_tp_energy_);
    events_->Branch("RH_unmatched", &ev_rh_unmatched_);
    events_->Branch("TP_unmatched", &ev_tp_unmatched_);
    events_->Branch("ieta", &ev_tp_ieta_);
    events_->Branch("iphi", &ev_tp_iphi_);
    events_->Branch("depth", &ev_tp_depth_);
    events_->Branch("soi", &ev_tp_soi_);
    events_->Branch("et", &ev_tp_et_);
    events_->Branch("ts0", &ev_tp_ts0_);
    events_->Branch("ts1", &ev_tp_ts1_);
    events_->Branch("ts2", &ev_tp_ts2_);
    events_->Branch("ts3", &ev_tp_ts3_);
    events_->Branch("ts4", &ev_tp_ts4_);
    events_->Branch("ts5", &ev_tp_ts5_);
    events_->Branch("ts6", &ev_tp_ts6_);
    events_->Branch("ts7", &ev_tp_ts7_);
    events_->Branch("event", &ev_tp_event_);

    matches_ = fs->make<TTree>("matches", "Matched RH and TP");
    matches_->Branch("RH_energyM0", &mt_rh_energy0_);
    matches_->Branch("RH_energy", &mt_rh_energy_);
    matches_->Branch("TP_energy", &mt_tp_energy_);
    matches_->Branch("RH_energyM0_depth", &mt_rh_energy0_depth_, 32000, 0);
    matches_->Branch("RH_energy_depth", &mt_rh_energy_depth_, 32000, 0);
    matches_->Branch("TP_energy_depth", &mt_tp_energy_depth_, 32000, 0);
    matches_->Branch("ieta", &mt_ieta_);
    matches_->Branch("iphi", &mt_iphi_);
    matches_->Branch("depth", &mt_depth_);
    matches_->Branch("tp_soi", &mt_tp_soi_);
    matches_->Branch("bx", &mt_bx_);
    matches_->Branch("event", &mt_event_);
    matches_->Branch("ts0", &mt_ts0_);
    matches_->Branch("ts1", &mt_ts1_);
    matches_->Branch("ts2", &mt_ts2_);
    matches_->Branch("ts3", &mt_ts3_);
    matches_->Branch("ts4", &mt_ts4_);
    matches_->Branch("ts5", &mt_ts5_);
    matches_->Branch("ts6", &mt_ts6_);
    matches_->Branch("ts7", &mt_ts7_);
    matches_->Branch("d0", &mt_d0_);
    matches_->Branch("d1", &mt_d1_);
    matches_->Branch("d2", &mt_d2_);
    matches_->Branch("d3", &mt_d3_);
    matches_->Branch("d4", &mt_d4_);
    matches_->Branch("d5", &mt_d5_);
    matches_->Branch("d6", &mt_d6_);
    matches_->Branch("d7", &mt_d7_);
}

HcalCompareUpgradeChains::~HcalCompareUpgradeChains() {}


double
HcalCompareUpgradeChains::get_cosh(const HcalDetId& id)
{

  const auto *sub_geo = dynamic_cast<const HcalGeometry*>(gen_geo_->getSubdetectorGeometry(id));
  auto eta = sub_geo->getPosition(id).eta();
  return cosh(eta);
}

void
HcalCompareUpgradeChains::beginLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup& setup)
{
  edm::ESHandle<HcalChannelQuality> status;
  setup.get<HcalChannelQualityRcd>().get("withTopo", status);
  status_ = status.product();
  edm::ESHandle<HcalSeverityLevelComputer> comp;
  setup.get<HcalSeverityLevelComputerRcd>().get(comp);
  comp_ = comp.product();
}


void
HcalCompareUpgradeChains::analyze(const edm::Event& event, const edm::EventSetup& setup)
{

    using namespace edm;

    setup.get<CaloGeometryRecord>().get(tpd_geo_h_);
    const HcalTrigTowerGeometry& tpd_geo = *tpd_geo_h_;

    // ==========
    // Dataframes
    // ==========

    mt_bx_ = event.bunchCrossing(); tpsplit_bx_ = mt_bx_;
    mt_event_ = event.id().event(); tpsplit_event_ = mt_event_; tp_event_ = mt_event_; ev_tp_event_ = mt_event_;
    
    Handle<QIE11DigiCollection> frames;
    Handle<QIE10DigiCollection> hfframes;
    if (first_ && event.getByLabel(frames_[0], frames) && event.getByLabel(frames_[1], hfframes)) {
       //std::set<HcalTrigTowerDetId> ids;

       //for (const auto& frame: *(frames.product())) {
       //   auto mapped = tpd_geo_h_->towerIds(frame.id());

       //   for (const auto& id: mapped) {
       //      df_multiplicity_->Fill(id.ieta(), id.iphi());
       //      ids.insert(id);
       //   }
       //}

       //for (const auto& frame: *(hfframes.product())) {
       //   auto mapped = tpd_geo_h_->towerIds(frame.id());

       //   for (const auto& id: mapped) {
       //      df_multiplicity_->Fill(id.ieta(), id.iphi());
       //      ids.insert(id);
       //   }
       //}

       //for (const auto& id: ids) {
       //   tp_multiplicity_->Fill(id.ieta(), id.iphi());
       //}

       first_ = false;
    }

    // ==============
    // Matching stuff
    // ==============

    ev_rh_energy0_ = 0.;
    ev_rh_energy_ = 0.;
    ev_rh_unmatched_ = 0.;
    ev_tp_energy_ = 0.;
    ev_tp_unmatched_ = 0.;

    std::map<HcalTrigTowerDetId, std::vector<HBHERecHit>> rhits;
    //std::map<HcalTrigTowerDetId, std::vector<HFRecHit>> fhits;
    std::map<HcalTrigTowerDetId, std::vector<HcalUpgradeTriggerPrimitiveDigi>> tpdigis;

    Handle<HcalUpgradeTrigPrimDigiCollection> digis;
    if (!event.getByLabel(digis_, digis)) {
       LogError("HcalTrigPrimDigiCleaner") <<
          "Can't find hcal trigger primitive digi collection with tag '" <<
          digis_ << "'" << std::endl;
       return;
    }

    edm::Handle< edm::SortedCollection<HBHERecHit> > hits;
    if (!event.getByLabel(rechits_[0], hits)) {
       edm::LogError("HcalCompareUpgradeChains") <<
          "Can't find rec hit collection with tag '" << rechits_[0] << "'" << std::endl;
       /* return; */
    }

    //edm::Handle< edm::SortedCollection<HFRecHit> > hfhits;
    //if (!event.getByLabel(rechits_[1], hfhits)) {
    //  edm::LogError("HcalCompareUpgradeChains") <<
    //    "Can't find rec hit collection with tag '" << rechits_[1] << "'" << std::endl;
    //  /* return; */
    //}

    setup.get<CaloGeometryRecord>().get(gen_geo_);

    auto isValid = [&](const auto& hit) {
        HcalDetId id(hit.id());
        auto s = status_->getValues(id);
        int level = comp_->getSeverityLevel(id, 0, s->getValue());
        return level <= max_severity_;
    };

    if (hits.isValid()) {
        for (auto& hit: *(hits.product())) {
            HcalDetId id(hit.id());
            if (not isValid(hit))  continue;
            ev_rh_energy0_ += hit.eraw() / get_cosh(id);
            ev_rh_energy_ += hit.energy() / get_cosh(id); //cosh(local_geo->getPosition().eta());

            auto tower_ids = tpd_geo.towerIds(id);
            //if (id.depth() == 2)
            //    std::cout << "HIT ID IS: " << hit.id() << " AND TOWER ID IS: " << tower_ids[0] << std::endl;

            for (auto& tower_id: tower_ids) {
                tower_id = HcalTrigTowerDetId(tower_id.ieta(), tower_id.iphi(), 0);
                rhits[tower_id].push_back(hit);
            }
        }
    }

    //if (hfhits.isValid()) {
    //    for (auto& hit: *(hfhits.product())) {
    //        HcalDetId id(hit.id());
    //        if (not isValid(hit)) continue;
    //        ev_rh_energy0_ += hit.energy() / get_cosh(id);
    //        ev_rh_energy_ += hit.energy() / get_cosh(id);

    //        auto tower_ids = tpd_geo.towerIds(id);
    //        for (auto& tower_id: tower_ids) {
    //            tower_id = HcalTrigTowerDetId(tower_id.ieta(), tower_id.iphi(), id.depth(), tower_id.version());
    //            fhits[tower_id].push_back(hit);
    //        }
    //    }
    //}

    setup.get<CaloTPGRecord>().get(decoder_);

    ev_tp_et_.clear();     ev_tp_et_.shrink_to_fit();  
    ev_tp_soi_.clear();    ev_tp_soi_.shrink_to_fit();
    ev_tp_ieta_.clear();   ev_tp_ieta_.shrink_to_fit();
    ev_tp_iphi_.clear();   ev_tp_iphi_.shrink_to_fit();
    ev_tp_depth_.clear();  ev_tp_depth_.shrink_to_fit();
    ev_tp_ts0_.clear();    ev_tp_ts0_.shrink_to_fit();
    ev_tp_ts1_.clear();    ev_tp_ts1_.shrink_to_fit();
    ev_tp_ts2_.clear();    ev_tp_ts2_.shrink_to_fit();
    ev_tp_ts3_.clear();    ev_tp_ts3_.shrink_to_fit();
    ev_tp_ts4_.clear();    ev_tp_ts4_.shrink_to_fit();
    ev_tp_ts5_.clear();    ev_tp_ts5_.shrink_to_fit();
    ev_tp_ts6_.clear();    ev_tp_ts6_.shrink_to_fit();
    ev_tp_ts7_.clear();    ev_tp_ts7_.shrink_to_fit();

    for (const auto& digi: *digis) {

        std::fill(tp_ts_adc_.begin(), tp_ts_adc_.end(), 0);

        HcalTrigTowerDetId id = HcalTrigTowerDetId(digi.id().ieta(), digi.id().iphi(), digi.id().depth(), digi.id().version());

        ev_tp_energy_ += decoder_->hcaletValue(id,digi.SOI_compressedEt());  

        tpdigis[id].push_back(digi);

        tp_energy_ = decoder_->hcaletValue(id, digi.SOI_compressedEt());
        tp_ieta_ = id.ieta();
        tp_iphi_ = id.iphi();

        tp_depth_start_ = -1;
        tp_depth_end_ = -1;
        tp_depth_max_ = -1;
        int et_max = 0;
        int et_sum = 0;

        std::vector<int> energy_depth = digi.getDepthData();
        for (int i = 0; i < static_cast<int>(energy_depth.size()); ++i) {
           int depth = energy_depth[i];
           if (depth > 0) {
              et_sum += depth;
              tp_depth_end_ = i;
              if (tp_depth_start_ < 0)
                 tp_depth_start_ = i;
              if (depth > et_max) {
                 tp_depth_max_ = i;
                 et_max = depth;
              }
           }
        }

        tp_soi_ = digi.SOI_compressedEt();
        tp_depth_ = 1;

        ev_tp_et_.push_back(decoder_->hcaletValue(id, digi.SOI_compressedEt()));
        ev_tp_soi_.push_back(digi.SOI_compressedEt());
        ev_tp_ieta_.push_back(id.ieta());
        ev_tp_iphi_.push_back(id.iphi());
        ev_tp_depth_.push_back(id.depth());

        std::vector<int> ts_adc = digi.getSampleData();
        for (int i = 0; i < static_cast<int>(ts_adc.size()); i++) { tp_ts_adc_[i] = ts_adc[i]; }
        ev_tp_ts0_.push_back(ts_adc[0]);
        ev_tp_ts1_.push_back(ts_adc[1]);
        ev_tp_ts2_.push_back(ts_adc[2]);
        ev_tp_ts3_.push_back(ts_adc[3]);
        ev_tp_ts4_.push_back(ts_adc[4]);
        ev_tp_ts5_.push_back(ts_adc[5]);
        ev_tp_ts6_.push_back(ts_adc[6]);
        ev_tp_ts7_.push_back(ts_adc[7]);

        //tps_->Fill();

        if (et_sum > 0) {
            for (int i = 0; i < static_cast<int>(energy_depth.size()); ++i) {

                int depth = energy_depth[i];
                tpsplit_energy_ = tp_energy_ * float(depth) / et_sum;
                tpsplit_oot_ = tp_energy_ * float(digi.SOI_oot_linear(i)) / et_sum;

                tpsplit_ieta_ = tp_ieta_;
                tpsplit_iphi_ = tp_iphi_;
                tpsplit_depth_ = i;
                tpsplit_ettot_ = tp_energy_;

                tpsplit_rise_avg_ = digi.SOI_rising_avg(i);
                tpsplit_rise_rms_ = digi.SOI_rising_rms(i);
                tpsplit_fall_avg_ = digi.SOI_falling_avg(i);
                tpsplit_fall_rms_ = digi.SOI_falling_rms(i);

                //tpsplit_->Fill();
            }
        }
    }

    for (const auto& pair: tpdigis) {

        std::fill(mt_rh_energy0_depth_.begin(), mt_rh_energy0_depth_.end(), 0.0);
        std::fill(mt_rh_energy_depth_.begin(),  mt_rh_energy_depth_.end(),  0.0);
        std::fill(mt_tp_energy_depth_.begin(),  mt_tp_energy_depth_.end(),  0.0);

        auto id = pair.first;

        auto new_id(id);
        if (swap_iphi_ and id.version() == 1 and id.ieta() > 28 and id.ieta() < 40) {
            if (id.iphi() % 4 == 1)
                new_id = HcalTrigTowerDetId(id.ieta(), (id.iphi() + 70) % 72, id.depth(), id.version());
            else
                new_id = HcalTrigTowerDetId(id.ieta(), (id.iphi() + 2) % 72 , id.depth(), id.version());
        }

        auto rh = rhits.find(new_id); //pair.first);
        if (rh != rhits.end()) {

            mt_ieta_ = new_id.ieta();//pair.first.ieta();
            mt_iphi_ = new_id.iphi(); //pair.first.iphi();
            mt_depth_ = new_id.depth(); //pair.first.depth();
            mt_tp_energy_ = 0;
            mt_tp_soi_ = 0;
            for (const auto& tp: pair.second) {
                 
                std::vector<int> sampleData = tp.getSampleData();
                mt_ts0_ = sampleData[0];
                mt_ts1_ = sampleData[1];
                mt_ts2_ = sampleData[2];
                mt_ts3_ = sampleData[3];
                mt_ts4_ = sampleData[4];
                mt_ts5_ = sampleData[5];
                mt_ts6_ = sampleData[6];
                mt_ts7_ = sampleData[7];

                std::vector<int> depthData = tp.getDepthData();
                mt_d0_ = depthData[0];
                mt_d1_ = depthData[1];
                mt_d2_ = depthData[2];
                mt_d3_ = depthData[3];
                mt_d4_ = depthData[4];
                mt_d5_ = depthData[5];
                mt_d6_ = depthData[6];
                mt_d7_ = depthData[7];

                mt_tp_energy_ += decoder_->hcaletValue( new_id, tp.SOI_compressedEt());
                for (int i = 0; i < static_cast<int>(depthData.size()); ++i) {  
                    mt_tp_energy_depth_[i] += depthData[i];
                }
                mt_tp_soi_ = tp.SOI_compressedEt();
            }

            double sum = 0;
            std::for_each(mt_tp_energy_depth_.begin(), mt_tp_energy_depth_.end(), [&](double &i){ sum += i; });
            std::for_each(mt_tp_energy_depth_.begin(), mt_tp_energy_depth_.end(), [=](double &i){ 
                if(sum != 0)
                    i = i / sum * mt_tp_energy_; 
            });

            mt_rh_energy0_ = 0.;
            mt_rh_energy_ = 0.;
            
            for (const auto& hit: rh->second) {
                HcalDetId id(hit.id());
                auto depth = id.depth();
                
                auto tower_ids = tpd_geo.towerIds(id);
                auto count = std::count_if(std::begin(tower_ids), std::end(tower_ids),
                  		     [&](const auto& t) { return t.version() == new_id.version(); });

                mt_rh_energy0_ += hit.eraw() / get_cosh(id) / count ; //cosh(local_geo->getPosition().eta());
                mt_rh_energy_ += hit.energy() / get_cosh(id) / count ; //cosh(local_geo->getPosition().eta());
                mt_rh_energy0_depth_[depth] += hit.eraw() / get_cosh(id) / count ; 
                mt_rh_energy_depth_[depth] += hit.energy() / get_cosh(id) / count ; 
            }
            
            matches_->Fill();
            rhits.erase(rh);
        } else {
            ++ev_tp_unmatched_;
        }
    }
    
    //for (const auto& pair: rhits) {
    //    ev_rh_unmatched_ += pair.second.size();
    //}
    
    events_->Fill();
}

void
HcalCompareUpgradeChains::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HcalCompareUpgradeChains);
