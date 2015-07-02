#include "LLGAnalysis.h"
#include "TLorentzVector.h"

void LLGAnalysis::SetupSignalRegion() {

    // setup the cutflow
    _cutFlow.insert(pair<string,int>("0_NoCut", 0) );
    _cutFlow.insert(pair<string,int>("1_Trigger", 0) );
    _cutFlow.insert(pair<string,int>("2_MuonVeto", 0) );
    _cutFlow.insert(pair<string,int>("3_ElectronVeto", 0) );
    _cutFlow.insert(pair<string,int>("4_HasPVJet", 0) );
    _cutFlow.insert(pair<string,int>("5_HasSV20", 0) );
    _cutFlow.insert(pair<string,int>("6_MET", 0) );
    _cutFlow.insert(pair<string,int>("7_DiJetMass", 0 ) );
    _cutFlow.insert(pair<string,int>("8_BVeto", 0) );
    _cutFlow.insert(pair<string,int>("9_SVPVDistance", 0) );
    
    // and the histograms 
    makeHist( "selected_distances", 40, 0., 40., "Distance between PV and SV [mm]", "Number of PV-SV pairs" );
    makeHist( "selected_met", 50, 0., 500., "MET [GeV]", "Number of Events" );
    makeHist( "selected_nPVJet", 4, -0.5, 3.5, "# PV with at least 1 Jet > 75 GeV", "Number of Events" );
    makeHist( "selected_nJetsToSV", 7, -0.5, 6.5, "# Jets associated to SV", "Number of Vertices" ); 
    makeHist( "selected_nSV", 5, -0.5, 4.5, "# SV with at least 1 Jet", "Number of Events" ); 
    makeHist( "nBjetAtSV", 5, -0.5, 4.5, "Number of b-jets associated to SV", "Number of SV" );
    makeHist( "mJJSV", 100, 0., 500., "DiJet mass at SV", "Number of Jet Pairs" );
    makeHist( "nJetsSV", 7, -0.5, 6.5, "Number of Jets associated to SV", "Number of SV" );
    makeHist( "SVJet1Pt", 50, 0., 500., "SV Leading Jet pT [GeV]", "Number of SV" );
    makeHist( "SVJet2Pt", 50, 0., 500., "SV 2^{nd} Leading Jet pT [GeV]", "Number of SV" );
    makeHist( "SVJet3Pt", 50, 0., 500., "SV 3^{rd} Leading Jet pT [GeV]", "Number of SV" );
    makeHist( "SVJet4Pt", 50, 0., 500., "SV 4^{th} Leading Jet pT [GeV]", "Number of SV" );
    makeHist( "PVJet1Pt", 50, 0., 1000., "PV Leading Jet pT [GeV]", "Number of Events" );
    makeHist( "PVJet2Pt", 50, 0., 500., "PV 2^{nd} Leading Jet pT [GeV]", "Number of Events" );
    makeHist( "PVJet3Pt", 50, 0., 500., "PV 3^{rd} Leading Jet pT [GeV]", "Number of Events" );
    makeHist( "PVJet4Pt", 50, 0., 500., "PV 4^{th} Leading Jet pT [GeV]", "Number of Events" );
    makeHist( "nJetsTotal", 26, -0.5, 25.5, "Number of Jets with p_{T} > 30 GeV", "Number of Events" );
    makeHist( "BJet1Pt", 50, 0., 500., "Leading B-tagged Jet p_{T} [GeV]", "Number of Events" );
    makeHist( "BJet2Pt", 50, 0., 500., "Subleading B-tagged Jet p_{T} [GeV]", "Number of Events" );
    makeHist( "BJet3Pt", 50, 0., 500., "3^{rd} leading B-tagged Jet p_{T} [GeV]", "Number of Events" );
    makeHist( "BJet4Pt", 50, 0., 500., "4^{th} leading B-tagged Jet p_{T} [GeV]", "Number of Events" );
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
        if( muon_iso->at(im) / pt  > 0.2 ) continue;
        if( pt > MUON_PT_CUT ) hasMuon = true;
    }
    if( hasMuon ) return; 
    _cutFlow.at("2_MuonVeto") += 1;
        
    
    bool hasElectron = false;
    for( unsigned int im = 0; im < electron_px->size(); ++im ) {
        double pt = sqrt(electron_px->at(im)*electron_px->at(im) + electron_py->at(im)*electron_py->at(im));
        if( electron_iso->at(im) / pt >= 0.15 ) continue;
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

    int idxLeadingJetPV = -1;
    double ptLeadingJetPV = -1.;
    int idxSubLeadingJetPV = -1;
    double ptSubLeadingJetPV = -1;
    int idxThirdLeadingJetPV = -1;
    double ptThirdLeadingJetPV = -1;
    int idxFourthLeadingJetPV = -1;
    double ptFourthLeadingJetPV = -1;
  
    for( unsigned int iPV = 0; iPV < vertex_x -> size(); ++iPV ) {
        bool hasJetPV = false;
        for( unsigned int iiJet = 0; iiJet < idJetsToPV.at(iPV).size(); ++iiJet ) {
            int iJet = idJetsToPV.at(iPV).at(iiJet);
            if( recoJet_pt->at(iJet) > JET_PT_CUT_PV ) hasJetPV = true;
            
            if( recoJet_pt->at(iJet) > ptLeadingJetPV ) {
                idxFourthLeadingJetPV = idxThirdLeadingJetPV;
                ptFourthLeadingJetPV = ptThirdLeadingJetPV;
                idxThirdLeadingJetPV = idxSubLeadingJetPV;
                ptThirdLeadingJetPV = ptSubLeadingJetPV;
                idxSubLeadingJetPV = idxLeadingJetPV;
                ptSubLeadingJetPV = ptLeadingJetPV;
                idxLeadingJetPV = iJet;
                ptLeadingJetPV = recoJet_pt->at(iJet);
            }
            else if( recoJet_pt->at(iJet) > ptSubLeadingJetPV ) {
                idxFourthLeadingJetPV = idxThirdLeadingJetPV;
                ptFourthLeadingJetPV = ptThirdLeadingJetPV;
                idxThirdLeadingJetPV = idxSubLeadingJetPV;
                ptThirdLeadingJetPV = ptSubLeadingJetPV;
                idxSubLeadingJetPV = iJet;
                ptSubLeadingJetPV = recoJet_pt->at(iJet);
            }
            else if( recoJet_pt->at(iJet) > ptThirdLeadingJetPV ) {
                idxFourthLeadingJetPV = idxThirdLeadingJetPV;
                ptFourthLeadingJetPV = ptThirdLeadingJetPV;
                idxThirdLeadingJetPV = iJet;
                ptThirdLeadingJetPV = recoJet_pt->at(iJet);
            }
            else if( recoJet_pt->at(iJet) > ptThirdLeadingJetPV ) {
                idxFourthLeadingJetPV = iJet;
                ptFourthLeadingJetPV = recoJet_pt->at(iJet);
            }
        }
        if( hasJetPV ) PVWithJet.push_back( iPV );
    }
    _histograms1D.at("PVJet1Pt").Fill( ptLeadingJetPV, evtWeight );
    _histograms1D.at("PVJet2Pt").Fill( ptSubLeadingJetPV, evtWeight );
    _histograms1D.at("PVJet3Pt").Fill( ptThirdLeadingJetPV, evtWeight );
    _histograms1D.at("PVJet4Pt").Fill( ptFourthLeadingJetPV, evtWeight );

    for( unsigned int iSV = 0; iSV < secVertex_x->size(); ++iSV ) {
        if( idJetsToSV.at(iSV).size() > 0 ) SVWithJets.push_back( iSV );
        if( idJetsToSV.at(iSV).size() >= 2 ) SVWith2Jets.push_back( iSV );
    }
                  
    // and run the selection:
    if( PVWithJet.size() == 1 ) {
        _cutFlow.at("4_HasPVJet") += 1;
        for( unsigned int iSV = 0; iSV < secVertex_x->size(); ++iSV ) {
          if( met > MET_CUT ) {
            _histograms1D.at("nJetsSV").Fill( idJetsToSV.at(iSV).size(), evtWeight );
          }
        }
        
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
                 
                bool hasDiJetPair100 = false;
                
                int nJets30 = 0;
                double BJet1Pt = -1;
                int idxBJet1 = -1; 
                double BJet2Pt = -1;
                int idxBJet2 = -1; 
                double BJet3Pt = -1;
                int idxBJet3 = -1; 
                double BJet4Pt = -1;
                int idxBJet4 = -1; 
                for( unsigned int iJet = 0; iJet < recoJet_pt->size(); ++iJet ) {
                  if( recoJet_pt->at(iJet) < JET_PT_CUT_SV ) continue;
                  if( fabs(recoJet_eta->at(iJet)) > JET_ETA_CUT ) continue;
                  nJets30++;
                  if( recoJet_btag_combinedInclusiveSecondaryVertexV2BJetTags->at(iJet) > 0.814 ) {
                    if( recoJet_pt->at(iJet) > BJet1Pt ) {
                      BJet4Pt = BJet3Pt;
                      idxBJet4 = idxBJet3;
                      BJet3Pt = BJet2Pt;
                      idxBJet3 = idxBJet2;
                      BJet2Pt = BJet1Pt;
                      idxBJet2 = idxBJet1;
                      BJet1Pt = recoJet_pt->at(iJet);
                      idxBJet1 = iJet;
                    }
                    else if( recoJet_pt->at(iJet) > BJet2Pt ) {
                      BJet4Pt = BJet3Pt;
                      idxBJet4 = idxBJet3;
                      BJet3Pt = BJet2Pt;
                      idxBJet3 = idxBJet2;
                      BJet2Pt = recoJet_pt->at(iJet);
                      idxBJet2 = iJet; 
                    }
                    else if( recoJet_pt->at(iJet) > BJet3Pt ) {
                      BJet4Pt = BJet3Pt;
                      idxBJet4 = idxBJet3;
                      BJet3Pt = recoJet_pt->at(iJet);
                      idxBJet3 = iJet; 
                    }
                    else if( recoJet_pt->at(iJet) > BJet3Pt ) {
                      BJet4Pt = recoJet_pt->at(iJet);
                      idxBJet4 = iJet; 
                    }
                  }
                }
                _histograms1D.at("nJetsTotal").Fill( nJets30, evtWeight );
                _histograms1D.at("BJet1Pt").Fill( BJet1Pt, evtWeight );
                _histograms1D.at("BJet2Pt").Fill( BJet2Pt, evtWeight );
                _histograms1D.at("BJet3Pt").Fill( BJet3Pt, evtWeight );
                _histograms1D.at("BJet4Pt").Fill( BJet4Pt, evtWeight );

                for( unsigned int iSV = 0; iSV < secVertex_x->size(); ++iSV ) {
                  if( idJetsToSV.at(iSV).size() <= 1 ) continue;
                  _histograms1D.at("nJetsSV").Fill( idJetsToSV.at(iSV).size(), evtWeight );
                  
                  int idxLeadingJet = -1;
                  double ptLeadingJet = -1.;
                  int idxSubLeadingJet = -1;
                  double ptSubLeadingJet = -1;
                  int idxThirdLeadingJet = -1;
                  double ptThirdLeadingJet = -1;
                  int idxFourthLeadingJet = -1;
                  double ptFourthLeadingJet = -1;

                  for( unsigned int iJToSV = 0; iJToSV < idJetsToSV.at(iSV).size(); ++iJToSV ) {
                    int jIdx = idJetsToSV.at(iSV).at(iJToSV);
                    if( recoJet_pt->at(jIdx) > ptLeadingJet ) {
                      idxFourthLeadingJet = idxThirdLeadingJet;
                      ptFourthLeadingJet = ptThirdLeadingJet;
                      idxThirdLeadingJet = idxSubLeadingJet;
                      ptThirdLeadingJet = ptSubLeadingJet;
                      idxSubLeadingJet = idxLeadingJet;
                      ptSubLeadingJet = ptLeadingJet;
                      idxLeadingJet = jIdx;
                      ptLeadingJet = recoJet_pt->at(jIdx);
                    }
                    else if ( recoJet_pt->at(jIdx) > ptSubLeadingJet ) {
                      idxFourthLeadingJet = idxThirdLeadingJet;
                      ptFourthLeadingJet = ptThirdLeadingJet;
                      idxThirdLeadingJet = idxSubLeadingJet;
                      ptThirdLeadingJet = ptSubLeadingJet;
                      idxSubLeadingJet = jIdx;
                      ptSubLeadingJet = recoJet_pt->at(jIdx);
                    }
                    else if ( recoJet_pt->at(jIdx) > ptThirdLeadingJet ) {
                      idxFourthLeadingJet = idxThirdLeadingJet;
                      ptFourthLeadingJet = ptThirdLeadingJet;
                      idxThirdLeadingJet = jIdx;
                      ptThirdLeadingJet = recoJet_pt->at(jIdx);
                    }
                    else if ( recoJet_pt->at(jIdx) > ptFourthLeadingJet ) {
                      idxFourthLeadingJet = jIdx;
                      ptFourthLeadingJet = recoJet_pt->at(jIdx);
                    }
                  }
                 

                  _histograms1D.at("SVJet1Pt").Fill( ptLeadingJet, evtWeight );
                  _histograms1D.at("SVJet2Pt").Fill( ptSubLeadingJet, evtWeight );
                  _histograms1D.at("SVJet3Pt").Fill( ptThirdLeadingJet, evtWeight );
                  _histograms1D.at("SVJet4Pt").Fill( ptFourthLeadingJet, evtWeight );
                  TLorentzVector p4Jet1, p4Jet2;
                  p4Jet1.SetPtEtaPhiM( recoJet_pt->at(idxLeadingJet), recoJet_eta->at(idxLeadingJet), recoJet_phi->at(idxLeadingJet), 0. );
                  p4Jet2.SetPtEtaPhiM( recoJet_pt->at(idxSubLeadingJet), recoJet_eta->at(idxSubLeadingJet), recoJet_phi->at(idxSubLeadingJet), 0. );
                  TLorentzVector p4DiJet = p4Jet1 + p4Jet2;
                  _histograms1D.at("mJJSV").Fill( p4DiJet.M(), evtWeight );
                  
                  if( p4DiJet.M() > 100. ) hasDiJetPair100 = true;

                }
                
                if( hasDiJetPair100 ) return;
                _cutFlow.at("7_DiJetMass") += 1;

                bool hasBjetFromSV = false;
                for( unsigned int iSV = 0; iSV < secVertex_x->size(); ++iSV ) {
                  if( idJetsToSV.at(iSV).size() <= 1 ) continue;
                  int nBjets = 0;
                  for( unsigned int iJToSV = 0; iJToSV < idJetsToSV.at(iSV).size(); ++iJToSV ) {
                    int jIdx = idJetsToSV.at(iSV).at(iJToSV);
                    //if( recoJet_btag_combinedInclusiveSecondaryVertexV2BJetTags->at(jIdx) > 0.814 ) {
                    /*
                    if( recoJet_btag_jetProbabilityBJetTags -> at(jIdx) > 0.790 ) { 
                      nBjets += 1;
                      hasBjetFromSV = true;
                    */
                    //}
                    
                  }
                  _histograms1D.at("nBjetAtSV").Fill(nBjets, evtWeight );
                }
                if( !hasBjetFromSV ) {
                  _cutFlow.at("8_BVeto") += 1;
                }
            }
        }
    }   
    return;
}
