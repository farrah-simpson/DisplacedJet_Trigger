// -*- C++ -*-
//
// Package:    DisplacedDJ/TriggerAnl
// Class:      TriggerAnl
// 
/**\class TriggerAnl TriggerAnl.cc DisplacedDJ/TriggerAnl/plugins/TriggerAnl.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  jingyu luo
//         Created:  Tue, 05 Sep 2017 14:33:05 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DisplacedDJ/TriggerAnl/interface/TriggerAnl.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class TriggerAnl : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit TriggerAnl(const edm::ParameterSet&);
      ~TriggerAnl();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
  

      // ----------member data ---------------------------
      //
       
      const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> theTTBToken;
    //  const edm::ESGetToken<JetCorrector, JetCorrectionsRecord> theJECToken;
      edm::InputTag trigHTFilter_;
      edm::InputTag trigJetFilter_;
      edm::InputTag trigPromptJetFilter_;
      edm::InputTag trigDisplacedJetFilter_;
      edm::EDGetTokenT<trigger::TriggerEvent> trigsummaryToken_;
      edm::EDGetTokenT<edm::View<reco::CaloJet> > jet_collT_;
      edm::EDGetTokenT<reco::VertexCollection> PVCollT_;
      edm::EDGetTokenT<reco::JetCorrector> jetCorrT_;
      edm::EDGetTokenT<pat::TriggerEvent> triggerEvent_;
      edm::EDGetTokenT<reco::JetTracksAssociationCollection> assocVTXToken_;

      TTree* tree;
      float CaloHT;
      float unCorr_CaloHT;

      std::string JECTag_;
      bool HTfilter_fired;
  
      std::vector<float> calojet_pt;
      std::vector<float> calojet_eta;
      std::vector<float> calojet_phi;
      std::vector<float> calojet_energy;
      
      std::vector<bool> calojet_tagged;
 
      std::vector<float> offline_jetpt;
      std::vector<float> offline_jeteta;
      std::vector<float> offline_jetphi;
      std::vector<float> offline_jetenergy;

      std::vector<float> online_jetpt;
      std::vector<float> online_jeteta;
      std::vector<float> online_jetphi;
      std::vector<float> online_jetenergy;

      std::vector<int> calojet_nPromptTracks_1000;
      std::vector<int> calojet_nPromptTracks_500;
      std::vector<int> calojet_nDisplacedTracks_500;

      std::vector<int> nPromptTracks_1000;
      std::vector<int> nPromptTracks_500;
      std::vector<int> nDisplacedTracks_500;
    
      std::vector<bool> PromptTagged;
      std::vector<bool> DisplacedTagged;
 
       
      unsigned int run_;
      unsigned int lumi_;
      unsigned int evt_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TriggerAnl::TriggerAnl(const edm::ParameterSet& iConfig)
:
   theTTBToken(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
   //theJECToken(esConsumes(edm::ESInputTag("", "AK4Calo"))), 
   trigHTFilter_ (iConfig.getParameter<edm::InputTag>("trigHTFilter")),
   trigJetFilter_ (iConfig.getParameter<edm::InputTag>("trigJetFilter")),
   trigPromptJetFilter_ (iConfig.getParameter<edm::InputTag>("trigPromptJetFilter")),
   trigDisplacedJetFilter_ (iConfig.getParameter<edm::InputTag>("trigDisplacedJetFilter")),
   trigsummaryToken_ (consumes<trigger::TriggerEvent>(iConfig.getParameter<edm::InputTag>("trigSummary"))),
   jet_collT_ (consumes<edm::View<reco::CaloJet> >(iConfig.getUntrackedParameter<edm::InputTag>("jets"))),
   PVCollT_ (consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("primaryVertices"))),
   jetCorrT_ (consumes<reco::JetCorrector>(iConfig.getUntrackedParameter<edm::InputTag>("jetCorr"))),
   triggerEvent_ (consumes<pat::TriggerEvent>(iConfig.getUntrackedParameter<edm::InputTag >("triggerEvent"))),
   assocVTXToken_ (consumes<reco::JetTracksAssociationCollection>(iConfig.getUntrackedParameter<edm::InputTag>("associatorVTX")))
   //JECTag_ (iConfig.getUntrackedParameter<std::string>("JECTag"))

{
   //now do what ever initialization is needed
   usesResource("TFileService");
   edm::Service<TFileService> fs; 
   tree = fs->make<TTree>("tree", "tree");

}


TriggerAnl::~TriggerAnl()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TriggerAnl::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   run_ = iEvent.id().run();
   lumi_ = iEvent.luminosityBlock();
   evt_ = iEvent.id().event();

   CaloHT = 0;
   unCorr_CaloHT = 0;
   HTfilter_fired = false;

   calojet_pt.clear();
   calojet_eta.clear();
   calojet_phi.clear();
   calojet_energy.clear();
 
   calojet_tagged.clear();

   offline_jetpt.clear();
   offline_jeteta.clear();
   offline_jetphi.clear();
   offline_jetenergy.clear();

   online_jetpt.clear();
   online_jeteta.clear();
   online_jetphi.clear();
   online_jetenergy.clear();

   calojet_nPromptTracks_1000.clear();
   calojet_nPromptTracks_500.clear();
   calojet_nDisplacedTracks_500.clear();

   nPromptTracks_1000.clear();
   nPromptTracks_500.clear();
   nDisplacedTracks_500.clear();
   
   PromptTagged.clear();
   DisplacedTagged.clear();

   Handle<edm::View<reco::CaloJet> > jet_coll;
   Handle<reco::VertexCollection> pvHandle;
   Handle<reco::JetCorrector> jetCorr;
   Handle<trigger::TriggerEvent> trigsummary;
   Handle<pat::TriggerEvent> triggerEvent;
   //ESHandle<TransientTrackBuilder> theB; 
  // ESHandle<JetCorrectorParametersCollection> JetCorParColl; 
   Handle<reco::JetTracksAssociationCollection> JetTracksVTX;
   //ESHandle<TransientTrackBuilder> theB;
    

   iEvent.getByToken(jet_collT_, jet_coll);
   iEvent.getByToken(PVCollT_, pvHandle);
   iEvent.getByToken(jetCorrT_, jetCorr);
   //iSetup.get<JetCorrectionsRecord>().get(JECTag_, JetCorParColl);
   const auto& theB = &iSetup.getData(theTTBToken);
   //iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", theB);
   iEvent.getByToken(trigsummaryToken_, trigsummary);
   iEvent.getByToken(triggerEvent_, triggerEvent);
   iEvent.getByToken(assocVTXToken_, JetTracksVTX);

   reco::Vertex pv = (*pvHandle)[0];
   //Calculate the Offline CaloHT
   //if(jet_coll->size()!=0){
   for (auto const& ijet: *jet_coll){
       double corrfac = jetCorr->correction(ijet);
       calojet_pt.push_back(ijet.pt()*corrfac);
       calojet_eta.push_back(ijet.eta());
       calojet_phi.push_back(ijet.phi());
       calojet_energy.push_back(ijet.energy()*corrfac);

       if (ijet.pt()*corrfac>40 and fabs(ijet.eta())<2.5){
           CaloHT+=ijet.pt()*corrfac;
       }
       if (ijet.pt()>40 and fabs(ijet.eta())<2.5){
           unCorr_CaloHT+=ijet.pt();
       }

       reco::TrackRefVector jettrks = reco::JetTracksAssociation::getValue(*JetTracksVTX, (const reco::Jet&)ijet);
       int nPtrk_1000=0;
       int nPtrk_500=0;
       int nDtrk_500=0; 
       for (auto const& itrk: jettrks){
           if(!itrk->quality(reco::TrackBase::highPurity)) continue;
           if(itrk->pt()< 1.0) continue;
           reco::TransientTrack t_trk = (*theB).build(itrk);
           Measurement1D ip2d = IPTools::absoluteTransverseImpactParameter(t_trk, pv).second;
           if(fabs(ip2d.value())<0.1) nPtrk_1000+=1;
           if(fabs(ip2d.value())<0.05) nPtrk_500+=1;
           if(fabs(ip2d.value())>0.05 && ip2d.significance()>5.0) nDtrk_500+=1;
       }
       calojet_nPromptTracks_1000.push_back(nPtrk_1000);
       calojet_nPromptTracks_500.push_back(nPtrk_500);
       calojet_nDisplacedTracks_500.push_back(nDtrk_500);
   }

   //trigger::size_type trigHTFilterIndex = trigsummary->filterIndex(trigHTFilter_);
 
   const pat::TriggerFilter* trigHTFilter = triggerEvent->filter(trigHTFilter_.label());
   const pat::TriggerFilter* trigJetFilter = triggerEvent->filter(trigJetFilter_.label());
   //const pat::TriggerFilter* trigPromptJetFilter = triggerEvent->filter(trigPromptJetFilter_.label());
   //const pat::TriggerFilter* trigDisplacedJetFilter = triggerEvent->filter(trigDisplacedJetFilter_.label());

   pat::TriggerObjectRefVector JetFilterObjs = triggerEvent->filterObjects(trigJetFilter_.label());
   pat::TriggerObjectRefVector PromptJetFilterObjs = triggerEvent->filterObjects(trigPromptJetFilter_.label());
   pat::TriggerObjectRefVector DisplacedJetFilterObjs = triggerEvent->filterObjects(trigDisplacedJetFilter_.label()); 
   
   if (trigHTFilter!=0){
       if (trigHTFilter->status()==1) HTfilter_fired=true; 
   }
  
   if (trigJetFilter!=0){
       for(size_t ijet=0; ijet<calojet_eta.size(); ijet++){
           bool jet_tagged=false;
           for(size_t ojet=0; ojet<JetFilterObjs.size(); ojet++){
               float dR = deltaR(calojet_eta.at(ijet), calojet_phi.at(ijet), JetFilterObjs.at(ojet)->eta(), JetFilterObjs.at(ojet)->phi());
               if(dR<0.4) jet_tagged=true;
           }
           calojet_tagged.push_back(jet_tagged);
       }
    
   if(trigJetFilter->status()==1){
       //std::cout<<"Yest"<<std::endl;
       //std::cout<<JetFilterObjs.size()<<std::endl;
       //for(size_t i=0; i<JetFilterObjs.size(); i++){
       //    std::cout<<"Jet Eta:"<<JetFilterObjs.at(i)->eta()<<", Jet Phi:"<<JetFilterObjs.at(i)->phi()<<std::endl;
       //}
       //for(size_t j=0; j<PromptJetFilterObjs.size(); j++){
       //    std::cout<<"Prompt Jet Eta:"<<PromptJetFilterObjs.at(j)->eta()<<", Prompt Jet Phi:"<<PromptJetFilterObjs.at(j)->phi()<<std::endl;
       //}
       //   
       //for(size_t j=0; j<DisplacedJetFilterObjs.size(); j++){
       //    std::cout<<"Displaced Jet Eta:"<<DisplacedJetFilterObjs.at(j)->eta()<<", Displaced Jet Phi:"<<DisplacedJetFilterObjs.at(j)->phi()<<std::endl;
       //}

       for (size_t i = 0; i<JetFilterObjs.size(); i++){
           size_t jm = 0;
           float mindR = 1e6;
       //    size_t matched_ijet = 0;
           for(size_t j=0; j<calojet_eta.size(); j++){
               float dR = deltaR(calojet_eta.at(j), calojet_phi.at(j), JetFilterObjs.at(i)->eta(), JetFilterObjs.at(i)->phi());
               if(dR<0.4 && dR<mindR){
                   mindR = dR;
                   jm=j;
               }
           }
       //    std::cout<<"Min dR:"<<mindR<<std::endl;
           if(mindR<0.4){
               bool prompt=false;
               bool displaced=false;
               online_jetpt.push_back(JetFilterObjs.at(i)->pt());
               online_jeteta.push_back(JetFilterObjs.at(i)->eta());
               online_jetphi.push_back(JetFilterObjs.at(i)->phi());
               online_jetenergy.push_back(JetFilterObjs.at(i)->energy());

               offline_jetpt.push_back(calojet_pt.at(jm));
               offline_jeteta.push_back(calojet_eta.at(jm));
               offline_jetphi.push_back(calojet_phi.at(jm));
               offline_jetenergy.push_back(calojet_energy.at(jm));

               nPromptTracks_1000.push_back(calojet_nPromptTracks_1000.at(jm));
               nPromptTracks_500.push_back(calojet_nPromptTracks_500.at(jm));
               nDisplacedTracks_500.push_back(calojet_nDisplacedTracks_500.at(jm));

               for(size_t ip=0; ip<PromptJetFilterObjs.size(); ip++){
                   float dR = deltaR(JetFilterObjs.at(i)->eta(), JetFilterObjs.at(i)->phi(), PromptJetFilterObjs.at(ip)->eta(), PromptJetFilterObjs.at(ip)->phi());
                   if(dR<0.1) prompt=true; 
               }
               for(size_t id=0; id<DisplacedJetFilterObjs.size(); id++){
                   float dR = deltaR(JetFilterObjs.at(i)->eta(), JetFilterObjs.at(i)->phi(), DisplacedJetFilterObjs.at(id)->eta(), DisplacedJetFilterObjs.at(id)->phi());
                   if(dR<0.1) displaced=true;
               } 
               
               PromptTagged.push_back(prompt);
               DisplacedTagged.push_back(displaced);
           }
       }
       //    float mindR=1e6;
       //    size_t matched_jet = 1e4;
       //    for (auto const& ijet: *jet_coll){
       //        float dR = deltaR(ijet.eta(), ijet.phi(), JetFilterObjs.at(i)->eta(), JetFilterObjs.at(i)->pt());
       //        if(dR<0.4 && dR<mindR){
       //            matched_jet = 
       //        } 
       //    }
       //    
       ////for (auto const& iobj : JetFilterObjs){
       ////    std::cout<<iobj->origObjCand()->pt()<<std::endl;
       ////}

       //} 
   }
   }
   //trigger::size_type trigJetFilterIndex = trigsummary->filterIndex(trigJetFilter_);
   ////trigger::TriggerObjectCollection triggerObjects = trigsummary->getObjects();
   ////float HLT_HT = 0;
   //if(trigHTFilterIndex < trigsummary->sizeFilters() and trigJetFilterIndex< trigsummary->sizeFilters()) {
   //   
   //     HTfilter_fired=true;
   //     //const trigger::Keys& trigKeys = trigsummary->filterKeys(trigHTFilterIndex);
   //     //for(unsigned int ik = 0; ik < trigKeys.size(); ++ik){
   //     //  const trigger::TriggerObject& obj = triggerObjects[trigKeys[ik]];
   //     //  HLT_HT+=obj.pt();
   //     //  std::cout<<"objpt: "<<obj.pt()<<std::endl;
   //     //  std::cout<<"objeta:"<<obj.eta()<<std::endl;
   //     //  std::cout<<"objphi:"<<obj.phi()<<std::endl;
   //     //}
   //     //std::cout<<HLT_HT<<std::endl;

   //}
   
   //if(HLT_HT>430) HTfilter_fired=true;
   //std::cout<<"CaloHT: "<<CaloHT<<"; Fired: "<<HTfilter_fired<<std::endl; 

   //}
   tree->Fill(); 


#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
TriggerAnl::beginJob()
{

    tree->Branch("run", &run_, "run/i");
    tree->Branch("lumi", &lumi_, "lumi/i");
    tree->Branch("evt", &evt_, "evt/i");
    tree->Branch("CaloHT", &CaloHT);
    tree->Branch("unCorr_CaloHT", &unCorr_CaloHT);
    tree->Branch("HTfilter_fired", &HTfilter_fired);

    tree->Branch("calojet_pt", &calojet_pt);
    tree->Branch("calojet_eta", &calojet_eta);
    tree->Branch("calojet_phi", &calojet_phi);
    tree->Branch("calojet_energy", &calojet_energy);
    tree->Branch("calojet_tagged", &calojet_tagged);

    tree->Branch("offline_jetpt", &offline_jetpt);
    tree->Branch("offline_jeteta", &offline_jeteta);
    tree->Branch("offline_jetphi", &offline_jetphi); 
    tree->Branch("offline_jetenergy", &offline_jetenergy);

    tree->Branch("online_jetpt", &online_jetpt);
    tree->Branch("online_jeteta", &online_jeteta);
    tree->Branch("online_jetphi",  &online_jetphi);
    tree->Branch("online_jetenergy", &online_jetenergy);


    tree->Branch("nPromptTracks_1000", &nPromptTracks_1000);
    tree->Branch("nPromptTracks_500", &nPromptTracks_500);
    tree->Branch("nDisplacedTracks_500", &nDisplacedTracks_500);
    tree->Branch("PromptTagged", &PromptTagged);
    tree->Branch("DisplacedTagged", &DisplacedTagged);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TriggerAnl::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TriggerAnl::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TriggerAnl);
