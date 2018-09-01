#ifndef data_ana_hh
#define data_ana_hh

#include "includes.hh"

using namespace std;

class data_ana {
public:
    data_ana(TString _fileName);
    ~data_ana();
    void analysis(Double_t *N12, Double_t *N0, Double_t p0, Double_t p1, Bool_t _orient);
    void analysis2(Double_t &N12, Double_t &N0);
    void fillhisto(TH2D* histo, Double_t p0, Double_t p1);

    TChain * fChain;
    Long64_t nentries;
    TBranch * bEvent;
    TBranch * bTime;
    TBranch * bDate;
    TBranch * bGonioPos;
    TBranch * bMultiHit;
    TBranch * bMultiHits;
    TBranch * bSingleTrack;
    TBranch * bTracks;

    double  INx_multi;
    double  INy_multi;
    double  IMPx_multi;
    double  IMPy_multi;
    double  INx;
    double  INy;
    double  OUTx;
    double  OUTy;
    double  IMPx;
    double  IMPy;

    static const Int_t number_of_cutting_angel_center_position = 41;
    static constexpr Double_t theta_cut_center_min = -100.0e-6;
    static constexpr Double_t theta_cut_center_max =  100.0e-6;
    static constexpr Double_t theta_cut = 2.5e-6;
    static constexpr Double_t initial_gpos_x = 1570674.0;
    Double_t F12;
    Double_t errF12;
    Double_t max_d0x;
    Double_t min_d0x;
    Double_t max_d0y;
    Double_t min_d0y;

    Int_t cut_d0x;
    Int_t cut_d0y;

    Double_t theta_cut_center[number_of_cutting_angel_center_position];

    struct Event
    {
        int run;
        int evtnum;
        int nuclear;
    } evt;

    struct GonioPos
    {
        double x;
        double y;
        double z;
    } gPos;

    struct MultiHit
    {
        int p_nHits[5];
        double thetaIn_x;
        double thetaIn_y;
        double thetaInErr_x;
        double thetaInErr_y;
        double d0_x;
        double d0_y;
        double d0Err_x;
        double d0Err_y;
    } hits;

    struct Track
    {
        double thetaIn_x;
        double thetaIn_y;
        double thetaOut_x;
        double thetaOut_y;
        double thetaInErr_x;
        double thetaInErr_y;
        double thetaOutErr_x;
        double thetaOutErr_y;
        double d0_x;
        double d0_y;
        double d0Err_x;
        double d0Err_y;
        double chi2_x;
        double chi2_y;
    } tracks;

    int isHit;

    int isTrack;

    char time[80];
    char date[80];
};

#endif
