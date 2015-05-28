#include "LLGAnalysis.h"

void LLGAnalysis::SetupSignalRegion() {

    // setup the cutflow
    _cutFlow.insert(pair<string,int>("0_NoCut", 0) );
    _cutFlow.insert(pair<string,int>("1_Trigger", 0) );
    _cutFlow.insert(pair<string,int>("2_MuonVeto", 0) );
    _cutFlow.insert(pair<string,int>("3_ElectronVeto", 0) );
    _cutFlow.insert(pair<string,int>("4_HasPVJet", 0) );
    _cutFlow.insert(pair<string,int>("5_HasSV20", 0) );
    _cutFlow.insert(pair<string,int>("6_MET", 0) );
    _cutFlow.insert(pair<string,int>("7_BVeto", 0) );
    _cutFlow.insert(pair<string,int>("8_SVPVDistance", 0) );
    
    // and the histograms 
    makeHist( "selected_distances", 40, 0., 40., "Distance between PV and SV [mm]", "Number of PV-SV pairs" );
    makeHist( "selected_met", 50, 0., 500., "MET [GeV]", "Number of Events" );
    makeHist( "selected_nPVJet", 4, -0.5, 3.5, "# PV with at least 1 Jet > 75 GeV", "Number of Events" );
    makeHist( "selected_nJetsToSV", 7, -0.5, 6.5, "# Jets associated to SV", "Number of Vertices" ); 
    makeHist( "selected_nSV", 5, -0.5, 4.5, "# SV with at least 1 Jet", "Number of Events" ); 
    makeHist( "nBjetAtSV", 5, -0.5, 4.5, "Number of b-jets associated to SV", "Number of SV" );

    return;
}

void LLGAnalysis::SignalRegionSelection() {

    _cutFlow.at("0_NoCut") += 1;

    
    
    bool passTrigger = false;
    for( unsigned int iTrig = 0; iTrig < triggerNames->size(); ++iTrig ) {
        if( (triggerNames->at(iTrig) == "HLT_PFJet260_v1" || triggerNames->at(iTrig) == "HLT_PFMET170_NoiseCleaned_v1") && triggerBits->at(iTrig) == 1 ) passTrigger = true;
    }

    if( !passTrigger ) return; 
    _cutFlow.at("1_Trigger") += 1;

    // lepton veto:
    bool hasMuon = false;
    for( unsigned int im = 0; im < muon_px->size(); ++im ) {
        double pt = sqrt(muon_px->at(im)*muon_px->at(im) + muon_py->at(im)*muon_py->at(im));
        if( pt > MUON_PT_CUT ) hasMuon = true;
    }
    if( hasMuon ) return; 
    _cutFlow.at("2_MuonVeto") += 1;
        
    
    bool hasElectron = false;
    for( unsigned int im = 0; im < electron_px->size(); ++im ) {
        double pt = sqrt(electron_px->at(im)*electron_px->at(im) + electron_py->at(im)*electron_py->at(im));
        if( pt > ELECTRON_PT_CUT ) hasElectron = true;
    }
    if( hasElectron ) return;
    _cutFlow.at("3_ElectronVeto") += 1;




    // now assign jets to the vertices:
    vector<int> nJetsToPV( vertex_x->size(), 0 );
    vector<int> nJetsToSV( secVertex_x->size(), 0 );
    vector<vector<int> > idJetsToPV;
    vector<vector<int> > idJetsToSV;
    for( unsigned int iVtx = 0; iVtx < vertex_x->size(); ++iVtx ) {
        vector<int> idx;
        idJetsToPV.push_back( idx );
    }
    for( unsigned int iVtx = 0; iVtx < secVertex_x->size(); ++iVtx ) {
        vector<int> idx;
        idJetsToSV.push_back( idx );
    }
        
    for( unsigned int iJet = 0; iJet < recoJet_pt->size(); ++iJet ) {
        if( recoJet_pt->at(iJet) < JET_PT_CUT_SV ) continue;
        if( fabs(recoJet_eta->at(iJet)) > JET_ETA_CUT ) continue;
            
        //calculate jet vertex position:
        unsigned int nCons = 0;
        double weightednCons = 0.;
        vector<double> error(3,0.);
        vector<double> position = CalculateVertex( recoJet_constVertex_x->at(iJet), recoJet_constVertex_y->at(iJet), recoJet_constVertex_z->at(iJet), recoJet_const_pt->at(iJet), recoJet_const_charge->at(iJet), recoJet_const_closestVertex_d->at(iJet), nCons, weightednCons, error );
            
        int nMatch = 0;
        for( unsigned int iVtx = 0; iVtx < vertex_x->size(); ++iVtx ) {
            if( fabs(position.at(0) - vertex_x->at(iVtx) ) < 1.e-10 &&
                fabs(position.at(1) - vertex_y->at(iVtx) ) < 1.e-10 &&
                fabs(position.at(2) - vertex_z->at(iVtx) ) < 1.e-10 ) {
                nJetsToPV.at(iVtx) += 1;
                idJetsToPV.at(iVtx).push_back( iJet );
                nMatch += 1;
            }
        }
        for( unsigned int iVtx = 0; iVtx < secVertex_x->size(); ++iVtx ) {
            if( fabs(position.at(0) - secVertex_x->at(iVtx) ) < 1.e-10 &&
                fabs(position.at(1) - secVertex_y->at(iVtx) ) < 1.e-10 &&
                fabs(position.at(2) - secVertex_z->at(iVtx) ) < 1.e-10 ) {
                
                nJetsToSV.at(iVtx) += 1;
                idJetsToSV.at(iVtx).push_back( iJet );
                nMatch += 1;

            }
        }
        if( nMatch > 1 ) {
            cout << "WARNING! ASSOCIATED JET TO MORE THAN 1 VERTEX ?!" << endl;
        }
    }

    // now count the number of vertices with jets:
    vector<int> PVWithJet;
    vector<int> SVWithJets;
    vector<int> SVWith2Jets;
    for( unsigned int iPV = 0; iPV < vertex_x -> size(); ++iPV ) {
        bool hasJetPV = false;
        for( unsigned int iiJet = 0; iiJet < idJetsToPV.at(iPV).size(); ++iiJet ) {
            int iJet = idJetsToPV.at(iPV).at(iiJet);
            if( recoJet_pt->at(iJet) > JET_PT_CUT_PV ) hasJetPV = true;
        }
        if( hasJetPV ) PVWithJet.push_back( iPV );
    }
    for( unsigned int iSV = 0; iSV < secVertex_x->size(); ++iSV ) {
        if( idJetsToSV.at(iSV).size() > 0 ) SVWithJets.push_back( iSV );
        if( idJetsToSV.at(iSV).size() > 1 ) SVWith2Jets.push_back( iSV );
    }
        

    // and run the selection:
    if( PVWithJet.size() == 1 ) {
        _cutFlow.at("4_HasPVJet") += 1;

        if( SVWith2Jets.size() > 0 ) {
            _cutFlow.at("5_HasSV20") += 1;
            
            _histograms1D.at("selected_met").Fill( met, evtWeight );
            _histograms1D.at("selected_nPVJet").Fill( PVWithJet.size(), evtWeight  ); 
            _histograms1D.at("selected_nSV").Fill( SVWithJets.size(), evtWeight  );
                
            for( unsigned int iSV = 0; iSV < SVWithJets.size(); ++iSV ) {
                _histograms1D.at("selected_nJetsToSV").Fill( idJetsToSV.at(SVWithJets.at(iSV)).size(), evtWeight  );
            }
            vector<double> allDistances;

            for( unsigned int iPV = 0; iPV < PVWithJet.size(); ++iPV ) {
                double thispv_x = vertex_x->at(PVWithJet.at(iPV));
                double thispv_y = vertex_y->at(PVWithJet.at(iPV));
                double thispv_z = vertex_z->at(PVWithJet.at(iPV));
                for( unsigned int iSV = 0; iSV < SVWithJets.size(); ++iSV ) {
                    double thissv_x = secVertex_x->at(SVWithJets.at(iSV));
                    double thissv_y = secVertex_y->at(SVWithJets.at(iSV));
                    double thissv_z = secVertex_z->at(SVWithJets.at(iSV));
                    double dx = thissv_x - thispv_x;
                    double dy = thissv_y - thispv_y;
                    double dz = thissv_z - thispv_z;
                    double dist = 10.*sqrt( dx*dx + dy*dy + dz*dz );
                    _histograms1D.at("selected_distances").Fill( dist, evtWeight);
                    allDistances.push_back( dist );
                }
            }
            if( met > MET_CUT ) {
                _cutFlow.at("6_MET") += 1;
                   
                bool hasBjetFromSV = false;
                for( unsigned int iSV = 0; iSV < secVertex_x->size(); ++iSV ) {
                  if( idJetsToSV.at(iSV).size() <= 1 ) continue;
                  int nBjets = 0;
                  for( unsigned int iJToSV = 0; iJToSV < idJetsToSV.at(iSV).size(); ++iJToSV ) {
                    int jIdx = idJetsToSV.at(iSV).at(iJToSV);
                    if( recoJet_btag_combinedInclusiveSecondaryVertexV2BJetTags->at(jIdx) > 0.423 ) {
                      nBjets += 1;
                      hasBjetFromSV = true;
                    }
                  }
                  _histograms1D.at("nBjetAtSV").Fill(nBjets, evtWeight );
                }
                if( !hasBjetFromSV ) {
                  _cutFlow.at("7_BVeto") += 1;
                }
            }
        }
    }   
    return;
}
