// -*- C++ -*-
//
// Package:    DisplacedDJ/TriggerAnlPF_miniAOD
// Class:      TriggerAnlPF_miniAOD
// 
/**\class TriggerAnlPF_miniAOD TriggerAnlPF_miniAOD.cc DisplacedDJ/TriggerAnlPF_miniAOD/plugins/TriggerAnlPF_miniAOD.cc

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
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

typedef reco::JetTracksAssociation::Container JetTrks;

class TriggerAnlPF_miniAOD : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit TriggerAnlPF_miniAOD(const edm::ParameterSet&);
      ~TriggerAnlPF_miniAOD();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      void jettracks (JetTrks*, const edm::Handle<edm::View<reco::Jet> > &, const std::vector<reco::Track> &, float);
  

      // ----------member data ---------------------------
      //
       
      const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> theTTBToken;
    //  const edm::ESGetToken<JetCorrector, JetCorrectionsRecord> theJECToken;
      edm::EDGetTokenT<pat::PackedCandidateCollection> TrackCollT_;
      edm::EDGetTokenT<pat::PackedCandidateCollection> LostTrackCollT_;
      edm::InputTag trigHTFilter_;
      edm::InputTag trigJetFilter_;
      edm::InputTag trigPromptJetFilter_;
      edm::InputTag trigDisplacedJetFilter_;
      edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
      edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone> > triggerObjects_;
      //edm::EDGetTokenT<trigger::TriggerEvent> trigsummaryToken_;
      edm::EDGetTokenT<edm::View<reco::Jet> > jet_collT_;
      edm::EDGetTokenT<reco::VertexCollection> PVCollT_;
      edm::EDGetTokenT<reco::JetCorrector> jetCorrT_;
      //edm::EDGetTokenT<pat::TriggerEvent> triggerEvent_;
      edm::EDGetTokenT<reco::JetTracksAssociationCollection> assocVTXToken_;

      TTree* tree;
      float CaloHT;
      float unCorr_CaloHT;

      std::string JECTag_;
      bool HTfilter_fired;
      bool Promptfilter_fired;
  
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
TriggerAnlPF_miniAOD::TriggerAnlPF_miniAOD(const edm::ParameterSet& iConfig)
:
   theTTBToken(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
   TrackCollT_ (consumes<pat::PackedCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks"))),
   LostTrackCollT_ (consumes<pat::PackedCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("losttracks"))),
   //theJECToken(esConsumes(edm::ESInputTag("", "AK4Calo"))), 
   trigHTFilter_ (iConfig.getParameter<edm::InputTag>("trigHTFilter")),
   trigJetFilter_ (iConfig.getParameter<edm::InputTag>("trigJetFilter")),
   trigPromptJetFilter_ (iConfig.getParameter<edm::InputTag>("trigPromptJetFilter")),
   trigDisplacedJetFilter_ (iConfig.getParameter<edm::InputTag>("trigDisplacedJetFilter")),

   triggerBits_ (consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
   triggerObjects_ (consumes<std::vector<pat::TriggerObjectStandAlone> >(iConfig.getParameter<edm::InputTag>("objects"))), 
   //trigsummaryToken_ (consumes<trigger::TriggerEvent>(iConfig.getParameter<edm::InputTag>("trigSummary"))),
   jet_collT_ (consumes<edm::View<reco::Jet> >(iConfig.getUntrackedParameter<edm::InputTag>("jets"))),
   PVCollT_ (consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("primaryVertices"))),
   jetCorrT_ (consumes<reco::JetCorrector>(iConfig.getUntrackedParameter<edm::InputTag>("jetCorr"))),
   //triggerEvent_ (consumes<pat::TriggerEvent>(iConfig.getUntrackedParameter<edm::InputTag >("triggerEvent"))),
   assocVTXToken_ (consumes<reco::JetTracksAssociationCollection>(iConfig.getUntrackedParameter<edm::InputTag>("associatorVTX")))
   //JECTag_ (iConfig.getUntrackedParameter<std::string>("JECTag"))

{
   //now do what ever initialization is needed
   usesResource("TFileService");
   edm::Service<TFileService> fs; 
   tree = fs->make<TTree>("tree", "tree");

}


TriggerAnlPF_miniAOD::~TriggerAnlPF_miniAOD()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TriggerAnlPF_miniAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace pat;

   run_ = iEvent.id().run();
   lumi_ = iEvent.luminosityBlock();
   evt_ = iEvent.id().event();

   CaloHT = 0;
   unCorr_CaloHT = 0;
   HTfilter_fired = false;
   Promptfilter_fired = false;

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

   Handle<PackedCandidateCollection> patcan;
   Handle<PackedCandidateCollection> losttracks;
   Handle<edm::View<reco::Jet> > jet_coll;
   Handle<reco::VertexCollection> pvHandle;
   Handle<reco::JetCorrector> jetCorr;
   //Handle<trigger::TriggerEvent> trigsummary;
   //Handle<pat::TriggerEvent> triggerEvent;
   //ESHandle<TransientTrackBuilder> theB; 
  // ESHandle<JetCorrectorParametersCollection> JetCorParColl; 
   Handle<reco::JetTracksAssociationCollection> JetTracksVTX;
   Handle<edm::TriggerResults> triggerBits; 
   Handle<std::vector<pat::TriggerObjectStandAlone> > triggerObjects;
   //ESHandle<TransientTrackBuilder> theB;
    

  
   iEvent.getByToken(jet_collT_, jet_coll);
   iEvent.getByToken(PVCollT_, pvHandle);
   iEvent.getByToken(jetCorrT_, jetCorr);
   //iSetup.get<JetCorrectionsRecord>().get(JECTag_, JetCorParColl);
   const auto& theB = &iSetup.getData(theTTBToken);
   iEvent.getByToken(triggerBits_, triggerBits);
   iEvent.getByToken(triggerObjects_, triggerObjects);

   iEvent.getByToken(TrackCollT_, patcan);
   iEvent.getByToken(LostTrackCollT_, losttracks);
   
   const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
   //iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", theB);
   //iEvent.getByToken(trigsummaryToken_, trigsummary);
   //iEvent.getByToken(triggerEvent_, triggerEvent);
   //iEvent.getByToken(assocVTXToken_, JetTracksVTX);
   //
   //
   float onlineHT = 0;
   std::vector<float> jetobj_pt;
   std::vector<float> jetobj_eta;
   std::vector<float> jetobj_phi;
   std::vector<float> promptobj_pt;
   std::vector<float> promptobj_eta;
   std::vector<float> promptobj_phi;
   std::vector<float> dispobj_pt;
   std::vector<float> dispobj_eta;
   std::vector<float> dispobj_phi;
    
   
   //std::cout<< "\n TRIGGER OBJECTS "<<std::endl;
   
   for(pat::TriggerObjectStandAlone obj : *triggerObjects) {
       obj.unpackPathNames(names);
       //std::cout<< "\tTrigger object: pt "<<obj.pt() << ", eta "<<obj.eta() <<", phi "<<obj.phi() <<std::endl;
       //std::cout<< "\t Collection: "<<obj.collection()<<std::endl;
       //std::cout<< "\t Type IDs:  ";
       //  for (unsigned h = 0; h < obj.filterIds().size(); ++h) std::cout << " " << obj.filterIds()[h] ;
       // std::cout << std::endl;
            
        obj.unpackFilterLabels(iEvent, *triggerBits);
        //std::cout << "\t   Filters:    ";
        for (unsigned h = 0; h < obj.filterLabels().size(); ++h) {
        // std::cout << " " << obj.filterLabels()[h];
         if(obj.filterLabels()[h]==trigHTFilter_.label() && obj.phi()==0){
           onlineHT = obj.pt();
           break;
         } 
         if(obj.filterLabels()[h]==trigJetFilter_.label()){
           jetobj_pt.push_back(obj.pt());
           jetobj_eta.push_back(obj.eta());
           jetobj_phi.push_back(obj.phi());
          // std::cout<<"*******Jet PT "<<obj.pt()<<"*******"<<std::endl;
           break;
         }

         if(obj.filterLabels()[h]==trigPromptJetFilter_.label()){
           promptobj_pt.push_back(obj.pt());
           promptobj_eta.push_back(obj.eta());
           promptobj_phi.push_back(obj.phi());
          // std::cout<<"******Prompt Jet PT "<<obj.pt()<<"*******"<<std::endl;
           break;
        
         }
         
         if(obj.filterLabels()[h]==trigDisplacedJetFilter_.label()){
           dispobj_pt.push_back(obj.pt());
           dispobj_eta.push_back(obj.eta());
           dispobj_phi.push_back(obj.phi());
          // std::cout<<"******Displaced Jet PT "<<obj.pt()<<"*******"<<std::endl;
           break;
        
         }

        }
        //if(onlineHT>0) std::cout<<"****ONLINE HT "<<onlineHT<<"****"<<std::endl;
        //std::cout << std::endl;
        std::vector pathNamesAll  = obj.pathNames(false);
        std::vector pathNamesLast = obj.pathNames(true);

        //std::cout << "\t   Paths (" << pathNamesAll.size()<<"/"<<pathNamesLast.size()<<"):    ";
        //for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
        //    bool isBoth = obj.hasPathName( pathNamesAll[h], true, true );
        //    bool isL3   = obj.hasPathName( pathNamesAll[h], false, true );
        //    bool isLF   = obj.hasPathName( pathNamesAll[h], true, false );
        //    bool isNone = obj.hasPathName( pathNamesAll[h], false, false );
        //    std::cout << "   " << pathNamesAll[h];
        //    if (isBoth) std::cout << "(L,3)";
        //    if (isL3 && !isBoth) std::cout << "(*,3)";
        //    if (isLF && !isBoth) std::cout << "(L,*)";
        //    if (isNone && !isBoth && !isL3 && !isLF) std::cout << "(*,*)";
        //}
        //std::cout << std::endl;
        }
    //std::cout << std::endl;

   


   reco::Vertex pv = (*pvHandle)[0];

   std::vector<reco::Track> alltracks;

   for (auto const& itrack : *patcan){
       if (itrack.trackHighPurity() && itrack.hasTrackDetails()){
           reco::Track tmptrk = itrack.pseudoTrack();
           if (tmptrk.quality(reco::TrackBase::highPurity) && tmptrk.pt()> 1.0 && tmptrk.charge()!=0){
               alltracks.push_back(tmptrk);
           }
       }
   }

   for (auto const& itrack : *losttracks){
       if (itrack.trackHighPurity() && itrack.hasTrackDetails()){
           reco::Track tmptrk = itrack.pseudoTrack();
           if (tmptrk.quality(reco::TrackBase::highPurity) && tmptrk.pt()> 1.0 && tmptrk.charge()!=0){
               alltracks.push_back(itrack.pseudoTrack());
           }
       }
   }

   std::auto_ptr<JetTrks> jettracks_assoc (new JetTrks(reco::JetRefBaseProd(jet_coll)));
   jettracks(&*jettracks_assoc,  jet_coll, alltracks, 0.5);



   //Calculate the Offline CaloHT
   //if(jet_coll->size()!=0){
   for (auto const& ijet: *jet_coll){
       double corrfac = jetCorr->correction(ijet);
       calojet_pt.push_back(ijet.pt()*corrfac);
       calojet_pt.push_back(ijet.pt());
       calojet_eta.push_back(ijet.eta());
       calojet_phi.push_back(ijet.phi());
       calojet_energy.push_back(ijet.energy()*corrfac);
       calojet_energy.push_back(ijet.energy());

       if (ijet.pt()>40 and fabs(ijet.eta())<2.5){
           CaloHT+=ijet.pt()*corrfac;
       }
       if (ijet.pt()>40 and fabs(ijet.eta())<2.5){
           unCorr_CaloHT+=ijet.pt();
       }

       reco::TrackRefVector jettrks = reco::JetTracksAssociation::getValue(*jettracks_assoc, (const reco::Jet&)ijet);
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
 
   //const pat::TriggerFilter* trigHTFilter = triggerEvent->filter(trigHTFilter_.label());
   //const pat::TriggerFilter* trigJetFilter = triggerEvent->filter(trigJetFilter_.label());
   ////const pat::TriggerFilter* trigPromptJetFilter = triggerEvent->filter(trigPromptJetFilter_.label());
   ////const pat::TriggerFilter* trigDisplacedJetFilter = triggerEvent->filter(trigDisplacedJetFilter_.label());

   //pat::TriggerObjectRefVector JetFilterObjs = triggerEvent->filterObjects(trigJetFilter_.label());
   //pat::TriggerObjectRefVector PromptJetFilterObjs = triggerEvent->filterObjects(trigPromptJetFilter_.label());
   //pat::TriggerObjectRefVector DisplacedJetFilterObjs = triggerEvent->filterObjects(trigDisplacedJetFilter_.label()); 
   //
   //if (trigHTFilter!=0){
       if (onlineHT>200) HTfilter_fired=true; //fix
   //}
  
   //if (trigJetFilter!=0){
       for(size_t ijet=0; ijet<calojet_eta.size(); ijet++){
           bool jet_tagged=false;
           for(size_t ojet=0; ojet<jetobj_pt.size(); ojet++){
               float dR = deltaR(calojet_eta.at(ijet), calojet_phi.at(ijet), jetobj_eta.at(ojet), jetobj_phi.at(ojet));
               if(dR<0.4) jet_tagged=true;
           }
           calojet_tagged.push_back(jet_tagged);
       }
   // 
   if(jetobj_pt.size()>1){
//       std::cout<<"Yest"<<std::endl;
//       std::cout<<JetFilterObjs.size()<<std::endl;
//       for(size_t i=0; i<JetFilterObjs.size(); i++){
//           std::cout<<"Jet Eta:"<<JetFilterObjs.at(i)->eta()<<", Jet Phi:"<<JetFilterObjs.at(i)->phi()<<std::endl;
//       }
//       for(size_t j=0; j<PromptJetFilterObjs.size(); j++){
//           std::cout<<"Prompt Jet Eta:"<<PromptJetFilterObjs.at(j)->eta()<<", Prompt Jet Phi:"<<PromptJetFilterObjs.at(j)->phi()<<std::endl;
//       }
          
 //      for(size_t j=0; j<DisplacedJetFilterObjs.size(); j++){
 //          std::cout<<"Displaced Jet Eta:"<<DisplacedJetFilterObjs.at(j)->eta()<<", Displaced Jet Phi:"<<DisplacedJetFilterObjs.at(j)->phi()<<std::endl;
//       }
       if(promptobj_pt.size()>1) Promptfilter_fired=true;

       for (size_t i = 0; i<jetobj_pt.size(); i++){
           size_t jm = 0;
           float mindR = 1e6;
       //    size_t matched_ijet = 0;
           for(size_t j=0; j<calojet_eta.size(); j++){
               float dR = deltaR(calojet_eta.at(j), calojet_phi.at(j), jetobj_eta.at(i), jetobj_phi.at(i));
               if(dR<0.4 && dR<mindR){
                   mindR = dR;
                   jm=j;
               }
           }
       //    std::cout<<"Min dR:"<<mindR<<std::endl;
           if(mindR<0.4){
               bool prompt=false;
               bool displaced=false;
               online_jetpt.push_back(jetobj_pt.at(i));
               online_jeteta.push_back(jetobj_eta.at(i));
               online_jetphi.push_back(jetobj_phi.at(i));
               online_jetenergy.push_back(jetobj_eta.at(i));

               offline_jetpt.push_back(calojet_pt.at(jm));
               offline_jeteta.push_back(calojet_eta.at(jm));
               offline_jetphi.push_back(calojet_phi.at(jm));
               offline_jetenergy.push_back(calojet_energy.at(jm));

               nPromptTracks_1000.push_back(calojet_nPromptTracks_1000.at(jm));
               nPromptTracks_500.push_back(calojet_nPromptTracks_500.at(jm));
               nDisplacedTracks_500.push_back(calojet_nDisplacedTracks_500.at(jm));

               for(size_t ip=0; ip<promptobj_pt.size(); ip++){
                   float dR = deltaR(jetobj_eta.at(i), jetobj_phi.at(i), promptobj_eta.at(ip), promptobj_phi.at(ip));
		   std::cout<<"Distances prompt: "<<dR<<std::endl;
                   if(dR<0.01) prompt=true; 
               }
               for(size_t id=0; id<dispobj_pt.size(); id++){
                   float dR = deltaR(jetobj_eta.at(i), jetobj_phi.at(i), dispobj_eta.at(id), dispobj_phi.at(id));
		   std::cout<<"Distances displaced: "<<dR<<std::endl;
                   if(dR<0.01) displaced=true;
               } 
               
               PromptTagged.push_back(prompt);
               DisplacedTagged.push_back(displaced);
           }
       }
   //    //    float mindR=1e6;
   //    //    size_t matched_jet = 1e4;
   //    //    for (auto const& ijet: *jet_coll){
   //    //        float dR = deltaR(ijet.eta(), ijet.phi(), JetFilterObjs.at(i)->eta(), JetFilterObjs.at(i)->pt());
   //    //        if(dR<0.4 && dR<mindR){
   //    //            matched_jet = 
   //    //        } 
   //    //    }
   //    //    
   //    ////for (auto const& iobj : JetFilterObjs){
   //    ////    std::cout<<iobj->origObjCand()->pt()<<std::endl;
   //    ////}

   //    //} 
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
TriggerAnlPF_miniAOD::beginJob()
{

    tree->Branch("run", &run_, "run/i");
    tree->Branch("lumi", &lumi_, "lumi/i");
    tree->Branch("evt", &evt_, "evt/i");
    tree->Branch("CaloHT", &CaloHT);
    tree->Branch("unCorr_CaloHT", &unCorr_CaloHT);
    tree->Branch("HTfilter_fired", &HTfilter_fired);
    tree->Branch("Promptfilter_fired", &Promptfilter_fired);

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
TriggerAnlPF_miniAOD::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TriggerAnlPF_miniAOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void
TriggerAnlPF_miniAOD::jettracks (JetTrks* fAssociation, const edm::Handle<edm::View<reco::Jet> >&  jetcoll, const std::vector<reco::Track>& tracks, float cone){


    for (unsigned ijet=0; ijet<jetcoll->size(); ijet++){
        edm::RefToBase<reco::Jet> ijet_ref = jetcoll->refAt(ijet);
        reco::TrackRefVector  associated;
        for (unsigned itrk = 0; itrk < tracks.size(); itrk++){
            if (float(deltaR(ijet_ref->eta(), ijet_ref->phi(), tracks[itrk].eta(), tracks[itrk].phi()))<cone){
                associated.push_back(reco::TrackRef(&tracks, itrk));
            }
                
        }
                
        reco::JetTracksAssociation::setValue(fAssociation, ijet_ref, associated);
                
        }
                
}
                

//define this as a plug-in
DEFINE_FWK_MODULE(TriggerAnlPF_miniAOD);
