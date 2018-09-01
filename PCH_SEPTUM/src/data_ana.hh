#ifndef data_ana_hh
#define data_ana_hh

#include "includes.hh"

using namespace std;

class data_ana {
public:
    data_ana(TString data_file_name);
    ~data_ana();
    void plot_histogramms();
    void write_histo_to_file(TString outFileName);
    void find_inelastic_frequency1(Double_t *N12, Double_t *N0, Int_t signal1, Int_t signal2);
    void find_inelastic_frequency2(Double_t *N12, Double_t *N0, Int_t signal1, Int_t signal2);
    void find_inelastic_frequency3(Double_t *N12, Double_t *N0, Int_t signal1, Int_t signal2);
    void GetParameters(string filenamepath, double &_min_d0x, double &_max_d0x, double &_min_d0y, double &_max_d0y, double &_gpos, double &_k1, double &_k2, double &_a, double &_b);
    void get_counts(Double_t *N, Double_t &N_tot);
    void get_counts_pairs(Double_t *N, Double_t &N_tot);
    int bit2det(int bit);
    int det2bit(int det);

    Int_t nentries;
    TFile * file;
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

    vector<double> * Hits4x = 0;
    vector<double> * Hits4y = 0;

    TH1D *hist1;
    TH1D *hist2;
    TH1D *hist3;
    TH1D *hist4;
    TH1D *hist5;
    TH1D *hist6;
    TH2D *hist7;
    TH2D *hist8;
    TH2D *hist9;
    TH2D *hist10;
    TH2D *hist11;
    TH2D *hist12;
    TH2D *hist13;
    TH2D *hist14;
    TH1D *hist15;
    TH1D *hist16;
    TH1D *hist17;
    TH2D *hist18;
    TH2D *hist19;
    TH2D *hist20;
    TH2D *hist21;
    TH1D *hist22;
    TH1D *hist23;
    TH2D *hist24;
    TH2D *hist25;
    TH2D *hist26;
    TH1D *hist27;
    TH1D *hist28;
    TH1D *hist29;
    TH1D *hist30;
    TH2D *hist31;
    TProfile *profile11;
    TProfile *profile11_corr;

    static const Int_t number_of_cutting_angels = 32;        // number of the cut angles
    static constexpr Double_t theta_cut_min = 2.5e-6;       // rad
    static constexpr Double_t theta_cut_max = 80.0e-6;      // rad
    static const Int_t number_of_scintillators = 6;

    Double_t F12;
    Double_t errF12;
    Double_t max_d0x;
    Double_t min_d0x;
    Double_t max_d0x_sept;
    Double_t min_d0x_sept;
    Double_t max_d0y;
    Double_t min_d0y;
    Double_t initial_gpos_x;
    Double_t k1;
    Double_t k2;
    Double_t a;
    Double_t b;
    Double_t angular_cut;
    Double_t distL;


    Int_t cut_d0x;
    Int_t cut_d0y;

    Double_t theta_cut[number_of_cutting_angels];

    struct Event
    {
        int run;
        int evtnum;
        int nuclear;
        int nuclearRaw;
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

