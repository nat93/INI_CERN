#include <iostream>
#include <string>
#include "TNtuple.h"
#include "TFile.h"
#include "TTree.h"
#include <fstream>
#include <TROOT.h>
#include "TH1D.h"
#include <TStyle.h>
#include "TGaxis.h"
#include "TCanvas.h"
#include "vector"
#include "TProfile.h"
#include "time.h"
#include <TGraph.h>
#include <TF1.h>
#include "stdlib.h"
#include "TSystem.h"
#include <vector>
#include <fstream>
#include "TLine.h"
#include "TEllipse.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TPad.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TPolyLine.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLatex.h"
#include "TPaveStats.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TMath.h"
#include "THStack.h"
#include "TSystem.h"
#include "TBenchmark.h"
#include "TH3D.h"
#include "TMatrixD.h"

using namespace std;

void Function_1(TString fileName);
void Function_2(TString fileName);
double Function_3(Double_t *v, Double_t *par);

const Double_t thetaCut = 5.0e-6;

const Double_t Maxd0x =  0.0;
const Double_t Mind0x = -0.8;
const Double_t Maxd0y =  3.0;
const Double_t Mind0y = -1.0;
const Double_t MinFit_dThetaX = 50.0; // min value for ch peak searching
const Double_t MaxFit_dThetaX = 150.0; // max value for ch peak searching
const Double_t initial_gpos_x = 4822988.600; // initial position for the perfect CH orientation of the Gonio

int GetTorsionParameters()
{
    Function_1("fileName");
    Function_2("fileName");

    return 0;
}

void Function_1(TString fileName)
{
    TString outFileName = fileName;
    outFileName += "_torsion_histo.root";
    TString inFileName  = "/path/";
    fileName += ".root";
    inFileName += fileName;

    TGraph* gr = new TGraph();
    Long64_t gr_i = 0;

    cout<<"--> FileName: "<<inFileName<<endl;

    //struct declaration
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
        int plane[5];
        double thetaIn_x;
        double thetaIn_y;
        double thetaInErr_x;
        double thetaInErr_y;
        double d0_x;
        double d0_y;
        double d0Err_x;
        double d0Err_y;
   } mHits;

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


    int multiHit;
    int singleTrack;

    char time[80];
    char date[80];

    TFile * file = new TFile(inFileName.Data());
    TTree * tree = (TTree*)file->Get("simpleEvent");

    Long64_t nentries = (Int_t)tree->GetEntries();

    cout<<"--> nEntries = "<<nentries<<endl;

    TBranch * bEvent    = tree->GetBranch("Event");
    TBranch * bTime     = tree->GetBranch("Time");
    TBranch * bDate     = tree->GetBranch("Date");
    TBranch * bGonioPos = tree->GetBranch("GonioPos");
    TBranch * bnhits    = tree->GetBranch("MultiHit");
    TBranch * bHits     = tree->GetBranch("MultiHits");
    TBranch * bntracks  = tree->GetBranch("SingleTrack");
    TBranch * bTracks   = tree->GetBranch("Tracks");

    bEvent->SetAddress(&evt);
    bTime->SetAddress(&time);
    bDate->SetAddress(&date);
    bGonioPos->SetAddress(&gPos);
    bnhits->SetAddress(&multiHit);
    bHits->SetAddress(&mHits);
    bntracks->SetAddress(&singleTrack);
    bTracks->SetAddress(&tracks);

    TH1D * dTheta_x = new TH1D("dTheta_x","dTheta_x",40,-200.0,200.0);
    TH2D * dTheta_x_vs_ThetaIn_x = new TH2D("dTheta_x_vs_ThetaIn_x","dTheta_x_vs_ThetaIn_x",40,-200.0,200.0,40,-200.0,200.0);
    TH3D * dTheta_x_vs_ThetaIn_x_vs_Impact_y = new TH3D("dTheta_x_vs_ThetaIn_x_vs_Impact_y","dTheta_x_vs_ThetaIn_x_vs_Impact_y",40,-200.0,200.0,40,-200.0,200.0,100,-5.0,5.0);

    dTheta_x_vs_ThetaIn_x->GetXaxis()->SetTitle("Impact Angle #theta_{x} [#murad]");
    dTheta_x_vs_ThetaIn_x->GetYaxis()->SetTitle("#Delta#theta_{x} [#murad]");

    int single_track_num = 0;

    int cutd0y = 0, cutd0x = 0, angle_cut = 0;
    Int_t div = nentries*0.01;

    tree->LoadTree(0);
    tree->GetEntry(0);

    for (Long64_t i = 0; i < nentries; i++)
    {
        if(i%div == 0)
        {
            printf("\r--> Working: %3.1f %% | gPos.x: %10.3f", 100*(Double_t)i/nentries,gPos.x-initial_gpos_x);
            fflush(stdout);
        }

        gr->SetPoint(gr_i, i, gPos.x);
        gr_i++;

        tree->LoadTree(i);
        tree->GetEntry(i);

        cutd0x = 0;
        cutd0y = 0;
        angle_cut = 0;

        if (tracks.d0_y < Maxd0y && tracks.d0_y > Mind0y) {cutd0y = 1;}
        if (tracks.d0_x < Maxd0x && tracks.d0_x > Mind0x) {cutd0x = 1;}

        if(singleTrack == 1 && multiHit == 0 && cutd0x == 1 && cutd0y == 1)
        {
            single_track_num++;

            if (TMath::Abs(tracks.thetaIn_x-(gPos.x-initial_gpos_x)*1e-6) < thetaCut) {angle_cut = 1;}
            if (angle_cut) {dTheta_x->Fill((tracks.thetaOut_x - tracks.thetaIn_x)*1e6);}

            dTheta_x_vs_ThetaIn_x_vs_Impact_y->Fill((tracks.thetaOut_x - tracks.thetaIn_x)*1e6,(tracks.thetaIn_x*1e6)-(gPos.x-initial_gpos_x),tracks.d0_y);
            dTheta_x_vs_ThetaIn_x->Fill((tracks.thetaIn_x*1e6)-(gPos.x-initial_gpos_x),(tracks.thetaOut_x - tracks.thetaIn_x)*1e6);
        }
    }
    printf("\nDone!\n");

    cout<<"--> Number of the single track events: "<<single_track_num<<endl;

    TFile *outfile = new TFile(outFileName.Data(),"RECREATE");
    dTheta_x_vs_ThetaIn_x_vs_Impact_y->Write();
    dTheta_x->Write();
    dTheta_x_vs_ThetaIn_x->Write();
    gr->SetName("gr");
    gr->GetXaxis()->SetTitle("Event ID");
    gr->GetYaxis()->SetTitle("gPos.x [#murad]");
    gr->Write();
    outfile->Close();
    dTheta_x_vs_ThetaIn_x->Draw("colz");

    TCanvas* c_00 = new TCanvas();
    c_00->cd();
    gr->Draw("APL");
}

void Function_2(TString fileName)
{
    TString outFileName = fileName;
    outFileName += "_torsion_graph.root";
    fileName += "_torsion_histo.root";
    TFile * in_file = new TFile(fileName.Data());

    cout<<"--> FileName: "<<fileName<<endl;

    TH3D * h = (TH3D*) in_file->Get("dTheta_x_vs_ThetaIn_x_vs_Impact_y");
    TH1D * h1 = (TH1D*) in_file->Get("dTheta_x");

    Double_t par_ch[3] = {};

    TCanvas* c_21 = new TCanvas();
    c_21->cd();
    h1->SetLineColor(kBlack);
    TF1 *g1 = new TF1("g1","gaus",MinFit_dThetaX,MaxFit_dThetaX);
    g1->SetLineColor(kRed);
    h1->Fit(g1,"R");
    g1->GetParameters(&par_ch[0]);
    TF1 *g2 = new TF1("g2","gaus",par_ch[1],par_ch[1] + 5.0*par_ch[2]);
    g2->SetLineColor(kBlue);
    h1->Fit(g2,"R+");

    cout<<"Nch "<<"("<<h1->GetBinCenter(25)<<";"<<h1->GetBinCenter(40)<<") = "<<h1->Integral(25,40)<<
          " | Ntot "<<"("<<h1->GetBinCenter(1)<<";"<<h1->GetBinCenter(40)<<") = "<<h1->Integral(1,40)<<
          " | Efficiency = "<<100.0*h1->Integral(25,40)/h1->Integral(1,40)<<" %"<<endl;

    Int_t yn = h->GetNbinsY(); // ThetaIn_x
    Int_t zn = h->GetNbinsZ(); // d0_y

    TMatrixD int_am(yn+1,zn+1);
    TMatrixD int_ch(yn+1,zn+1);
    TMatrixD eff_int(yn+1,zn+1);

    for (int i = 0; i < yn; i++)
    {
        for (int j = 0; j < zn; j++)
        {
            TH1D *hx = h->ProjectionX("z slice",i,i,j,j);

            TAxis *axis = hx->GetXaxis();
            int bmin = 1;
            int bmax = axis->GetNbins();
            int_am[i][j] = hx->Integral(bmin,bmax);

            if (int_am[i][j] != 0)
            {
                int bmin_ch = axis->FindBin(par_ch[1]);
                int bmax_ch = axis->FindBin(par_ch[1] + 5.0*par_ch[2]);

                int_ch[i][j] = hx->Integral(bmin_ch,bmax_ch);
                eff_int[i][j] = 2.0*int_ch[i][j]/(int_am[i][j]);
            }
            else
            {
                eff_int[i][j] = 0.;
            }
        }
    }

    TH2D* efficiency_vs_ThetaIn_x_vs_Impact_y = new TH2D("efficiency_vs_ThetaIn_x_vs_Impact_y","efficiency_vs_ThetaIn_x_vs_Impact_y",100,-5.0,5.0,40,-200.0,200.0);
    efficiency_vs_ThetaIn_x_vs_Impact_y->SetTitle("Efficiency in Impact Angle x vs Impact Position y");
    efficiency_vs_ThetaIn_x_vs_Impact_y->GetXaxis()->SetTitle("Impact Position y [mm]");
    efficiency_vs_ThetaIn_x_vs_Impact_y->GetYaxis()->SetTitle("Impact Angle x [#murad]");

    for (int j = 0; j < zn; j++)
    {
        for (int i = 0; i < yn; i++)
        {
            if(eff_int[i][j] > 0.2 && i > 10 && i < 30)
                efficiency_vs_ThetaIn_x_vs_Impact_y->SetBinContent(j,i,eff_int[i][j]);
        }
    }

    TProfile *profile = efficiency_vs_ThetaIn_x_vs_Impact_y->ProfileX();
    profile->GetXaxis()->SetTitle("Impact Position y [mm]");
    profile->GetYaxis()->SetTitle("Impact Angle x [#murad]");

    TCanvas* c_22 = new TCanvas();
    c_22->cd();
    efficiency_vs_ThetaIn_x_vs_Impact_y->Draw("colz");

    TCanvas* c_23 = new TCanvas();
    c_23->cd();
    profile->Draw();


    Double_t par[2] = {};
    TF1 * func = new TF1("func","pol1",Mind0y,Maxd0y);
    func->SetLineColor(kRed);
    profile->Fit(func,"R");
    func->GetParameters(&par[0]);
    cout<<"--> Torsion parameters: p[0] = "<<par[0]*1e-6<<" [rad] | p[1] = "<<par[1]*1e-6<<" [rad/mm]"<<endl;

/*
    Double_t par[4] = {};
    TF1 *func = new TF1("fit",Function_3,Mind0y,Maxd0y,4);
    func->SetParameters(8.74423e-01,3.25249e-01,-3.13487e-01,3.99149e+00);
    func->SetParNames("k1","k2","a","b");
    profile->Fit("fit","R");
    func->GetParameters(&par[0]);

    printf("k1 \t\t k2 \t\t a \t\t b\n");
    printf("%3.3fe-6 \t %3.3fe-6 \t %3.3f \t %3.3fe-6\n",func->GetParameter(0),func->GetParameter(1),func->GetParameter(2),func->GetParameter(3));
*/
    TFile *outfile = new TFile(outFileName,"RECREATE");
    c_21->Write();
    c_22->Write();
    c_23->Write();

    outfile->Close();
}

double Function_3(Double_t *v, Double_t *par)
{
    Double_t k1 = par[0];
    Double_t k2 = par[1];
    Double_t a  = par[2];
    Double_t b  = par[3];

    return (k2*TMath::Power(v[0]-a,2) + k1*v[0] + b);
}
