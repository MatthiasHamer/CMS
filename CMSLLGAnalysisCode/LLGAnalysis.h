#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <stdlib.h>
#include <TStyle.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TColor.h>
#include <math.h>
#include <THnSparse.h>
#include <TChain.h>
#include <map>
#include <string>
#include <vector>
#include <TRandom3.h>
#include "THnSparse.h"
#include "TF1.h"
#include "TSystem.h"

using namespace std;

class LLGAnalysis {
    public:
        static LLGAnalysis* GetInstance( char* configFileName );
        ~LLGAnalysis() {}
        
        vector<double> CalculateVertex( vector<double> x, vector<double> y, vector<double> z, vector<double> weight, vector<int> charge, vector<double> distance, unsigned int &nConsidered, double &weightednConsidered, vector<double> &error ); 
        
        bool Init();
        void RunEventLoop( int nEventsMax = -1);
        void FinishRun();
    
    private:
        static LLGAnalysis* _instance;

        void makeHist( string nametitle, int nbinsx, double xmin, double xmax, int nbinsy, double ymin, double ymax, string xtitle, string ytitle, string ztitle, string drawOption = "", double xAxisOffset = 1., double yAxisOffset = 1.2, double zAxisOffset = 1. ); 
        void makeHist( string nametitle, int nbins, double xmin, double xmax, string xtitle, string ytitle, string drawOption = "", double xAxisOffset = 1., double yAxisOffset = 1.2 );
        void setStyle(double ytoff = 1.0, bool marker = true, double left_margin = 0.15); 
        void MakeEfficiencyPlot( TH1D hpass, TH1D htotal, TCanvas *c, string triggerName = "");

    private:
        LLGAnalysis() {}
        LLGAnalysis( char* configFileName );
        
        map<string, int>        _cutFlow;
        map<string, TH1D>       _histograms1D;
        map<string, TH2D>       _histograms2D;
        map<string, string>     _histograms1DDrawOptions;
        map<string, string>     _histograms2DDrawOptions;
        vector<string>          _inputFileNames;
        string                  _inputTreeName;
        TChain                  *_inputTree;

        vector<double> *recoJet_pt; 
        vector<double> *recoJet_phi; 
        vector<double> *recoJet_eta; 
        vector<double> *recoJet_secVertex3D; 
        vector<double> *recoJet_secVertexSig; 
        vector<double> *recoJet_bJetTag; 
        vector<double> *recoJet_closestVertex_x; 
        vector<double> *recoJet_closestVertex_y; 
        vector<double> *recoJet_closestVertex_z; 
        vector<double> *recoJet_btag_combinedInclusiveSecondaryVertexV2BJetTags;
        vector<double> *muon_px; 
        vector<double> *muon_py; 
        vector<double> *muon_pz; 
        vector<double> *muon_phi; 
        vector<double> *muon_eta; 
        vector<double> *electron_px; 
        vector<double> *electron_py; 
        vector<double> *electron_pz; 
        vector<double> *electron_phi; 
        vector<double> *electron_eta; 
        vector<int> *triggerBits; 
        vector<string> *triggerNames; 
        vector<vector<double> > *recoJet_constVertex_x; 
        vector<vector<double> > *recoJet_constVertex_y; 
        vector<vector<double> > *recoJet_constVertex_z; 
        vector<vector<double> > *recoJet_constclosestVertex_x; 
        vector<vector<double> > *recoJet_constclosestVertex_y; 
        vector<vector<double> > *recoJet_constclosestVertex_z; 
        vector<vector<double> > *recoJet_const_pt; 
        vector<vector<double> > *recoJet_const_closestVertex_dxy; 
        vector<vector<double> > *recoJet_const_closestVertex_dz; 
        vector<vector<double> > *recoJet_const_closestVertex_d;
        vector<vector<int> > *recoJet_const_charge; 

        vector<double> *vertex_x;
        vector<double> *vertex_y;
        vector<double> *vertex_z; 
        vector<double> *secVertex_x;
        vector<double> *secVertex_y; 
        vector<double> *secVertex_z; 
        vector<double> *vertex_nTracks; 
        vector<double> *vertex_pt; 
    
        double met;

        double JET_PT_CUT;
        double JET_ETA_CUT;
        double MUON_PT_CUT;
        double ELECTRON_PT_CUT;
        double MET_CUT;

        vector<string> _plotFormats;

};

