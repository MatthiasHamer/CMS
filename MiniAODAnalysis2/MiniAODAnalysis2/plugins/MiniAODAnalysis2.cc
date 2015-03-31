// -*- C++ -*-
//
// Package:    MiniAODAnalysis2/MiniAODAnalysis2
// Class:      MiniAODAnalysis2
// 
/**\class MiniAODAnalysis2 MiniAODAnalysis2.cc MiniAODAnalysis2/MiniAODAnalysis2/plugins/MiniAODAnalysis2.cc

 Description: 
 Simple analysis class to dump output in a plain rootfile

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matthias Hamer
//         Created:  Tue, 13 Jan 2015 17:58:12 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TFile.h"
#include "TTree.h"
//
//
// class declaration
//

class MiniAODAnalysis2 : public edm::EDAnalyzer {
   public:
      explicit MiniAODAnalysis2(const edm::ParameterSet&);
      ~MiniAODAnalysis2();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      edm::EDGetTokenT<pat::METCollection> metToken_;
      edm::EDGetTokenT<pat::JetCollection> jetToken_;
      edm::EDGetTokenT<reco::PFJetCollection> jetToken2_;
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
      edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> secVtxToken_;
      edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
      edm::EDGetTokenT<reco::GenJetCollection> genJetToken_;
      edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
      edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken_;
      edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
      edm::EDGetTokenT<pat::MuonCollection> muonToken_;
      edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
      
   
      // the output file and tree
      TFile *fOutput = new TFile("RecoOutput.root", "RECREATE");
      TTree *tOutput = new TTree("RecoData", "RecoData");
  

      // the variables used for output
      // missing transverse energy
      double met = 0.;
      
      // the jet variables
      std::vector<double> *jet_eta = new std::vector<double>;
      std::vector<double> *jet_phi = new std::vector<double>;
      std::vector<double> *jet_pt = new std::vector<double>;
     
      // the muon variables
      std::vector<double> *muon_px = new std::vector<double>;
      std::vector<double> *muon_py = new std::vector<double>;
      std::vector<double> *muon_pz = new std::vector<double>;
      std::vector<double> *muon_eta = new std::vector<double>;
      std::vector<double> *muon_phi = new std::vector<double>;
      
      // the electron variables
      std::vector<double> *electron_px = new std::vector<double>;
      std::vector<double> *electron_py = new std::vector<double>;
      std::vector<double> *electron_pz = new std::vector<double>;
      std::vector<double> *electron_eta = new std::vector<double>;
      std::vector<double> *electron_phi = new std::vector<double>;
     
      // the trigger bits and names
      std::vector<int> *triggerBits = new std::vector<int>;
      std::vector<std::string> *triggerNames = new std::vector<std::string>;

      // the jet constituents 
      std::vector<std::vector<double> >* jet_constVertex_x = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jet_constVertex_y = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jet_constVertex_z = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jet_const_pt      = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jet_const_eta     = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jet_const_phi     = new std::vector<std::vector<double> >;
      std::vector<std::vector<int> >*    jet_const_charge  = new std::vector<std::vector<int> >;
      std::vector<std::vector<double> >* jet_const_closestVertex_dxy = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jet_const_closestVertex_dz = new std::vector<std::vector<double> >;
      std::vector<std::vector<double> >* jet_const_closestVertex_d = new std::vector<std::vector<double> >;

      // the primary vertex information
      std::vector<double> *vertex_x = new std::vector<double>;
      std::vector<double> *vertex_y = new std::vector<double>;
      std::vector<double> *vertex_z = new std::vector<double>;
      std::vector<double> *vertex_nTracks = new std::vector<double>;
      std::vector<double> *vertex_pt = new std::vector<double>;

      // the secondary vertices
      std::vector<double> *secVertex_x = new std::vector<double>;
      std::vector<double> *secVertex_y = new std::vector<double>;
      std::vector<double> *secVertex_z = new std::vector<double>;


      // event metadata
      int RunNumber = 0;
      int EventNumber = 0;
      int LuminosityBlock = 0;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
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
MiniAODAnalysis2::MiniAODAnalysis2(const edm::ParameterSet& iConfig):
  metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
  jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
  jetToken2_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("jets2"))),
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  secVtxToken_(consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("secVertices"))),
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  genJetToken_(consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJets"))),
  prunedGenToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"))),
  packedGenToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packed"))),
  pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))),
  muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons")))
  {
   //now do what ever initialization is needed
   
   // define the triggers we (might) want to use
   // these seem to be interesting for us
   triggerNames->push_back( "HLT_PFJet260_v1" );
   triggerNames->push_back( "HLT_JetE30_NoBPTX_v1" );
   triggerNames->push_back( "HLT_JetE30_NoBPTX3BX_NoHalo_v1" );
   triggerNames->push_back( "HLT_JetE50_NoBPTX3BX_NoHalo_v1" );
   triggerNames->push_back( "HLT_JetE70_NoBPTX3BX_NoHalo_v1" );
   triggerNames->push_back( "HLT_PFMET170_NoiseCleaned_v1" );
   
   // set the output branches for the tree
   tOutput -> Branch("RecoJet_eta", &jet_eta );
   tOutput -> Branch("RecoJet_phi", &jet_phi );
   tOutput -> Branch("RecoJet_pt", &jet_pt );
   tOutput -> Branch("RecoJet_constVertex_x", &jet_constVertex_x );
   tOutput -> Branch("RecoJet_constVertex_y", &jet_constVertex_y );
   tOutput -> Branch("RecoJet_constVertex_z", &jet_constVertex_z );
   tOutput -> Branch("RecoJet_const_pt", &jet_const_pt );
   tOutput -> Branch("RecoJet_const_charge", &jet_const_charge );
   tOutput -> Branch("RecoJet_const_closestVertex_dxy", &jet_const_closestVertex_dxy );
   tOutput -> Branch("RecoJet_const_closestVertex_dz", &jet_const_closestVertex_dz );
   tOutput -> Branch("RecoJet_const_closestVertex_d", &jet_const_closestVertex_d );
   tOutput -> Branch("RecoJet_const_eta", &jet_const_eta );
   tOutput -> Branch("RecoJet_const_phi", &jet_const_phi );
   tOutput -> Branch("RecoMuon_px", &muon_px );
   tOutput -> Branch("RecoMuon_py", &muon_pz );
   tOutput -> Branch("RecoMuon_pz", &muon_py );
   tOutput -> Branch("RecoMuon_eta", &muon_eta );
   tOutput -> Branch("RecoMuon_phi", &muon_phi );
   tOutput -> Branch("RecoElectron_px", &electron_px );
   tOutput -> Branch("RecoElectron_py", &electron_pz );
   tOutput -> Branch("RecoElectron_pz", &electron_py );
   tOutput -> Branch("RecoElectron_eta", &electron_eta );
   tOutput -> Branch("RecoElectron_phi", &electron_phi );
   tOutput -> Branch("TriggerBits", &triggerBits );
   tOutput -> Branch("TriggerNames", &triggerNames );
   tOutput -> Branch("RecoVertex_x", &vertex_x );
   tOutput -> Branch("RecoVertex_y", &vertex_y );
   tOutput -> Branch("RecoVertex_z", &vertex_z );
   tOutput -> Branch("RecoSecVertex_x", &secVertex_x );
   tOutput -> Branch("RecoSecVertex_y", &secVertex_y );
   tOutput -> Branch("RecoSecVertex_z", &secVertex_z );
   tOutput -> Branch("RecoVertex_nTracks", &vertex_nTracks );
   tOutput -> Branch("RecoVertex_pt", &vertex_pt );
   tOutput -> Branch("MET", &met );
   tOutput -> Branch("EventNumber", &EventNumber );
   tOutput -> Branch("RunNumber", &RunNumber );
   tOutput -> Branch("LuminosityBlock", &LuminosityBlock );
}


MiniAODAnalysis2::~MiniAODAnalysis2()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MiniAODAnalysis2::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   // clear all variables
   met = 0.;
   muon_px->clear();
   muon_py->clear();
   muon_pz->clear();
   muon_eta->clear();
   muon_phi->clear();
   electron_px->clear();
   electron_py->clear();
   electron_pz->clear();
   electron_eta->clear();
   electron_phi->clear();
   triggerBits->clear();
   jet_eta->clear();
   jet_phi->clear();
   jet_pt->clear();
   jet_constVertex_x->clear();
   jet_constVertex_y->clear();
   jet_constVertex_z->clear();
   jet_const_pt->clear();
   jet_const_eta->clear();
   jet_const_phi->clear();
   jet_const_charge->clear();
   jet_const_closestVertex_dxy->clear();
   jet_const_closestVertex_dz->clear();
   jet_const_closestVertex_d->clear();
   vertex_x -> clear();
   vertex_y -> clear();
   vertex_z -> clear();
   secVertex_x -> clear();
   secVertex_y -> clear();
   secVertex_z -> clear();
   vertex_nTracks -> clear();
   vertex_pt -> clear();
  
   // get all the physics objects from the miniAOD
   edm::Handle<pat::ElectronCollection> electrons;
   iEvent.getByToken(electronToken_, electrons);
   
   edm::Handle<pat::MuonCollection> muons;
   iEvent.getByToken(muonToken_, muons);

   Handle<edm::View<reco::GenParticle> > pruned;
   iEvent.getByToken(prunedGenToken_, pruned);

   Handle<edm::View<pat::PackedGenParticle> > packed;
   iEvent.getByToken(packedGenToken_,packed);
   
   Handle<pat::JetCollection> jets;
   iEvent.getByToken( jetToken_, jets );
   
   Handle<reco::PFJetCollection> jets2;
   iEvent.getByToken( jetToken2_, jets2 );
  
   Handle<reco::VertexCollection> vertices;
   iEvent.getByToken( vtxToken_, vertices );
   
   Handle<reco::VertexCompositePtrCandidateCollection> secVertices;
   iEvent.getByToken( secVtxToken_, secVertices );

   Handle<pat::METCollection> mets;
   iEvent.getByToken( metToken_, mets );

   edm::Handle<edm::TriggerResults> evTriggerBits;
   iEvent.getByToken( triggerBits_, evTriggerBits );
  
   Handle<reco::GenJetCollection> genJets;
   iEvent.getByToken( genJetToken_, genJets );
  
   Handle<pat::PackedCandidateCollection> pfs;
   iEvent.getByToken(pfToken_, pfs);


   edm::EventAuxiliary aux = iEvent.eventAuxiliary();
   edm::EventID id = aux.id();
  
   EventNumber = id.event();
   RunNumber = id.run();
   LuminosityBlock = id.luminosityBlock();


    
   const edm::TriggerNames &names = iEvent.triggerNames(*evTriggerBits);
   //bool passTrigger = false;
   for(unsigned int i = 0; i < evTriggerBits->size(); ++i ) {
      for( unsigned int j = 0; j < triggerNames->size(); ++j ) {
        if( names.triggerName(i) == triggerNames->at(j) ) {
          triggerBits->push_back( evTriggerBits->accept(i) ? 1 : 0 );
          //if( evTriggerBits->accept(i) ) passTrigger = true;
        }
      }
   }

   
   // muons
   // currently using the muon id taken from here:
   // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Loose_Muon
   // using loose PF combined relative isolation cut (see https://twiki.cern.ch/twiki/bin/viewauth/CMS/SusyObjectExperts)
   // adding a pT cut of 15 GeV
   // not using a cut that constrains the muon to the pv
   for( const pat::Muon &m : *muons ) {
      if( !m.isPFMuon() ) continue;
      if( !( m.isGlobalMuon() || m.isTrackerMuon()) ) continue;
      if( !m.isLooseMuon() ) continue;
      if( m.pt() < 15. ) continue;
      if( fabs(m.eta()) > 2.5 ) continue;
      double pfRelIso = ( m.pfIsolationR04().sumChargedHadronPt 
                        + std::max(0., m.pfIsolationR04().sumNeutralHadronEt + m.pfIsolationR04().sumPhotonEt - 0.5*m.pfIsolationR04().sumPUPt ) ) 
                        / m.pt();
      if( pfRelIso > 0.2 ) continue;

      muon_px->push_back( m.px() );
      muon_py->push_back( m.pz() );
      muon_pz->push_back( m.py() );
      muon_phi->push_back( m.phi() );
      muon_eta->push_back( m.eta() );
   }
   
   // electrons
   // implemented the electron id taken from here:
   // https://twiki.cern.ch/twiki/bin/viewauth/CMS/Eg2012AnalysesSupportingMaterial#Electrons
   // and also
   // https://twiki.cern.ch/twiki/bin/viewauth/CMS/SusyObjectExperts
   // not using a cut that constraints the electron to the PV
   // calculate the isolation according to: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD
   for( const pat::Electron &e : *electrons ) {
      if(e.pt() < 15. ) continue;
      if( !(e.isEE() || e.isEB() ) ) continue;
      if( e.isEB() && e.hadronicOverEm() >= 0.15 ) continue;
      if( e.isEB() && (  fabs(e.deltaPhiSuperClusterTrackAtVtx()) >= 0.8 
                      || fabs(e.deltaEtaSuperClusterTrackAtVtx()) >= 0.007
                      || e.scSigmaIEtaIEta() >= 0.01 ) )
      continue;
      if( e.isEE() && ( fabs(e.deltaPhiSuperClusterTrackAtVtx()) >= 0.7 
                      || fabs(e.deltaEtaSuperClusterTrackAtVtx()) >= 0.01
                      || e.scSigmaIEtaIEta() >= 0.03 ) ) 
      continue;
       
      double charged = 0.;
      double neutral = 0.;
      double pileup = 0.;
      std::vector<reco::CandidatePtr> footprint;
      for( unsigned int i = 0; i < e.numberOfSourceCandidatePtrs(); ++i ) {
         footprint.push_back( e.sourceCandidatePtr(i) );
      }
      
      for( unsigned int i = 0; i < pfs->size(); ++i ) {
         const pat::PackedCandidate &pf = (*pfs)[i];
         
         if( deltaR( pf, e) < 0.3 ) {
            if( std::find(footprint.begin(), footprint.end(), reco::CandidatePtr(pfs,i)) != footprint.end() ) continue;
            if( pf.charge() == 0 ) {
               if( pf.pt() > 0.5 ) neutral += pf.pt();
            }
            else if ( pf.fromPV() >= 2 ) {
               charged += pf.pt();
            }
            else {
               if( pf.pt() > 0.5 ) pileup += pf.pt();
            }
         }  
      }
      double iso = charged + std::max(0.0, neutral-0.5*pileup) / e.pt();
      if( iso >= 0.15 ) continue;
    
      electron_px->push_back( e.px() );
      electron_py->push_back( e.pz() );
      electron_pz->push_back( e.py() );
      electron_phi->push_back( e.phi() );
      electron_eta->push_back( e.eta() );
   }

   // jets
   // using the tight selection from:
   // https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_13_TeV_data
   int ctrJet = -1;
   for( const pat::Jet &j : *jets ) {
     
     if( j.neutralHadronEnergyFraction() >= 0.90 ) continue;
     if( j.neutralEmEnergyFraction() >= 0.90 ) continue;
     if( j.numberOfDaughters() <= 1 ) continue;
     if( j.muonEnergyFraction() >= 0.8 ) continue;
     if( j.chargedEmEnergyFraction() >= 0.9 ) continue;
     
     if( fabs(j.eta()) > 2.4 ) {
        if( j.chargedHadronEnergyFraction() <= 0. ) continue;
        if( j.chargedMultiplicity() <= 0. ) continue;
        if( j.chargedEmEnergyFraction() >= 0.99 ) continue;
     }
     
     ctrJet += 1;
     
     std::vector<double> constVert_x;
     std::vector<double> constVert_y;
     std::vector<double> constVert_z;
     std::vector<double> const_pt;
     std::vector<double> const_eta;
     std::vector<double> const_phi;
     std::vector<int> const_charge;
     std::vector<double> constVert_closestVertex_dxy;
     std::vector<double> constVert_closestVertex_dz;
     std::vector<double> constVert_closestVertex_d;
     std::vector<double> average_distance( vertices->size(), 0. );
     
     // loop over the jet constituents to find the vertex closest to each:
     for( unsigned int iD = 0; iD < j.numberOfDaughters(); ++iD ) {
        
        const pat::PackedCandidate &dau1 = dynamic_cast<const pat::PackedCandidate &>(*j.daughter(iD));
        pat::PackedCandidate dau(dau1);
        
        // minimum distances (total, xy and z) for the jet constituent to any vertex
        double dMin = 100000.;
        double dxyMin = 10000.;
        double dzMin = 10000.;
        int ctr = -1;
        
        // first loop over the primary vertices
        for( const reco::Vertex &v : *vertices ) {
            ctr += 1;
            double d = sqrt( dau.dxy( v.position() ) * dau.dxy( v.position() ) + dau.dz( v.position() )*dau.dz( v.position() ) );
            
            // if the vertex is closer than the current reference vertex, set dMin, dxyMin, dzMin, and also change the vertex of reference
            if( d < dMin ) {
              dMin = d;
              dxyMin = dau.dxy( v.position() );
              dzMin = dau.dz( v.position() );
              const reco::VertexRef vref( vertices, ctr );
              dau.setVertexRef( vref );
              
            }
        }
  
        // now do the same for the secondary vertices
        // however, for some reason, I have to use additional variable here, jet constitutent won't accept a VertexCompositePtrCandidate as a new reference
        double jetVertex_x = -10000.;
        double jetVertex_y = -10000.;
        double jetVertex_z = -10000.;
        for( const reco::VertexCompositePtrCandidate &v : *secVertices ) {
            reco::Vertex::Point p;
            p.SetXYZ( v.vx(), v.vy(), v.vz() );
            double d = sqrt( dau.dxy( p ) * dau.dxy( p ) + dau.dz( p )*dau.dz( p ) );
            if( d < dMin ) {
              dMin = d;
              dxyMin = dau.dxy( p );
              dzMin = dau.dz( p );
              ctr = -1;  
              jetVertex_x = v.vx();
              jetVertex_y = v.vy();
              jetVertex_z = v.vz();
          }
        }
        


        // now fill all the variables for the jet constituent
        constVert_closestVertex_dxy.push_back( dxyMin );
        constVert_closestVertex_dz.push_back( dzMin );
        constVert_closestVertex_d.push_back( dMin );
        constVert_x.push_back( ctr != -1 ? dau.vertexRef().get()->position().x() : jetVertex_x );
        constVert_y.push_back( ctr != -1 ? dau.vertexRef().get()->position().y() : jetVertex_y );
        constVert_z.push_back( ctr != -1 ? dau.vertexRef().get()->position().z() : jetVertex_z );
        const_pt.push_back( dau.pt() );
        const_eta.push_back( dau.eta() );
        const_phi.push_back( dau.phi() );
        const_charge.push_back( dau.charge() );
     }


     // and fill allthe jet variables
     jet_pt->push_back( j.pt() );
     jet_eta->push_back( j.eta() );
     jet_phi->push_back( j.phi() );
     jet_constVertex_x->push_back( constVert_x ); 
     jet_constVertex_y->push_back( constVert_y ); 
     jet_constVertex_z->push_back( constVert_z ); 
     jet_const_closestVertex_dxy->push_back(constVert_closestVertex_dxy);
     jet_const_closestVertex_dz->push_back(constVert_closestVertex_dz);
     jet_const_closestVertex_d->push_back(constVert_closestVertex_d);
     jet_const_pt->push_back( const_pt );
     jet_const_eta->push_back( const_eta );
     jet_const_phi->push_back( const_phi );
     jet_const_charge->push_back( const_charge );
   }
  
  
   // now fill the primary vertex information
   for( const reco::Vertex &v : *vertices ) {
      vertex_x -> push_back( v.x() );
      vertex_y -> push_back( v.y() );
      vertex_z -> push_back( v.z() );
      vertex_nTracks -> push_back( v.nTracks() );
      vertex_pt -> push_back( v.p4().pt() );
   }
   
   // now fill the secondary vertex information
   for( const reco::VertexCompositePtrCandidate &v : *secVertices ) {
      secVertex_x -> push_back( v.vx() );
      secVertex_y -> push_back( v.vy() );
      secVertex_z -> push_back( v.vz() );
   }

   // and fill the met
   const pat::MET &themet = mets->front();
   met = themet.pt(); 

   // finally, write it to the tree
   tOutput->Fill(); 

}

// ------------ method called once each job just before starting event loop  ------------
void 
MiniAODAnalysis2::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MiniAODAnalysis2::endJob() 
{

  // save the output tree and file
  gDirectory = fOutput;
  tOutput->Write();
  fOutput->Close();

}

// ------------ method called when starting to processes a run  ------------
/*
void 
MiniAODAnalysis2::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
MiniAODAnalysis2::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
MiniAODAnalysis2::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
MiniAODAnalysis2::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MiniAODAnalysis2::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniAODAnalysis2);
