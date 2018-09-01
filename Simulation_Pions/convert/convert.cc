#include <iostream>
#include <stdlib.h>
#include <assert.h>

#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TTree.h"
#include "TChain.h"
#include "TLegend.h"

using namespace std;

// int main(int argc, char *argv[])

const Int_t nPoint = 100;
const Int_t nWidth = 11;
const Int_t Width_init = 0;

void convert(TString fileName, Double_t *F12, Double_t *errF12);

int main()
{
    TMultiGraph *mg = new TMultiGraph();
    TGraphErrors* gr[nWidth];
    TLegend *leg = new TLegend(0.6,0.6,0.9,0.9,"","brNDC");
    Double_t F12[nWidth][nPoint], errF12[nWidth][nPoint];
    Double_t energyThreshold[nPoint];
    for(Int_t i = 0; i < nPoint; i++)
    {
        energyThreshold[i] = 0.1 + i*100.0;
    }

    for(Int_t i = 0; i < nWidth; i++)
    {
        TString fileName = "output_"; fileName += Width_init+i; fileName += ".root";
        cout<<"--> FileName: "<<fileName<<endl;
        convert(fileName,F12[i],errF12[i]);

        gr[i] = new TGraphErrors(nPoint,energyThreshold,F12[i],0,errF12[i]);
        TString title = "Gap: "; title += Width_init+i; title += " mm";
        gr[i]->SetTitle(title.Data());
        gr[i]->SetLineColor(i+1);
        gr[i]->SetLineWidth(2);
        gr[i]->SetMarkerColor(kBlack);
        gr[i]->SetMarkerStyle(20+i);
        mg->Add(gr[i]);
        leg->AddEntry(gr[i],title.Data(),"lp");
    }

    TFile* outputFile = new TFile("OUTPUT.root","RECREATE");
    TCanvas* canvas = new TCanvas();
    canvas->Clear();
    canvas->SetFillColor(kWhite);
    canvas->cd();
    mg->Draw("APL");
    mg->GetXaxis()->SetTitle("Energy threshold, keV");
    mg->GetYaxis()->SetTitle("F12");
    leg->Draw();
    canvas->Write();

    outputFile->Close();

    return 0;
}

void convert(TString fileName, Double_t *F12, Double_t *errF12)
{    
    Double_t energyThreshold[nPoint];
    for(Int_t i = 0; i < nPoint; i++)
    {
        energyThreshold[i] = 0.1 + i*100.0;
    }

    Int_t N12_sim[nPoint] = {}, N0_sim = 0, coincType = 0;

    TChain* fChain = new TChain("Tree");
    fChain->Add(fileName.Data());

    Int_t nentries  = (Int_t)fChain->GetEntries();

    Int_t           INI;
    Double_t        detEnergy[2];

    TBranch        *b_INI;
    TBranch        *b_detEnergy;

    if (!fChain) {cout<<"ERROR: A problem with input TREE."<<endl; assert(0);}
    fChain->SetBranchAddress("detEnergy", detEnergy,  &b_detEnergy);
    fChain->SetBranchAddress("INI",       &INI,       &b_INI);

    TH1D *h1 = new TH1D("h1_detEnergy_scintillator1","Deposition energy in the plastic scintillator 1",10000,-1.0,100000.0);
    TH1D *h2 = new TH1D("h2_detEnergy_scintillator2","Deposition energy in the plastic scintillator 2",10000,-1.0,100000.0);
    TH1D *h3 = new TH1D("h3_CrystalnNuclearPerEvent","Total number of nuclear interactions in the crystal per event",200,0.0,100.0);
    TH1D *h4 = new TH1D("h4_Coincidences","Type of the plastic scintillators coincidences",22,-1.0,10.0);

    N0_sim = 0;
    for(Int_t i = 0; i < nPoint; i++) N12_sim[i] = 0;

    Long64_t nINI = 0;
    for(Int_t i = 0; i < nentries; i++)
    {
        fChain->GetEntry(i);

        h1->Fill(detEnergy[0]);
        h2->Fill(detEnergy[1]);
        h3->Fill(INI);
        coincType = 0;
        if(detEnergy[0] > energyThreshold[0]) coincType = 1;
        if(detEnergy[1] > energyThreshold[0]) coincType = 2;
        if(detEnergy[0] > energyThreshold[0] && detEnergy[1] > energyThreshold[0]) coincType = 3;
        h4->Fill(coincType);

        if(INI > 0)
        {
            nINI++;
            N0_sim++;
            for(int ii = 0; ii < nPoint; ii++)
            {
                if (detEnergy[0] >= energyThreshold[ii] && detEnergy[1] >= energyThreshold[ii]) N12_sim[ii]++;
            }
        }

        if(i%10000 == 0)
        {
            printf("\r--> Working: %3.1f %%", 100*(Double_t)i/nentries);
            fflush(stdout);
        }
    }

    printf("\n");
    cout<<"--> nINI = "<<nINI<<endl;

    for(int jj = 0; jj < nPoint; jj++)
    {
        F12[jj] = (double)N12_sim[jj]/N0_sim;
        errF12[jj] = TMath::Sqrt((double)N12_sim[jj]/TMath::Power(N0_sim,2) + (double)TMath::Power(N12_sim[jj],2)*N0_sim/TMath::Power(N0_sim,4));
        //cout<<"--> F12 ("<<energyThreshold[jj]<<" keV) = "<<(double)N12_sim[jj]/N0_sim<<" +/- "<<TMath::Sqrt((double)N12_sim[jj]/TMath::Power(N0_sim,2) + (double)TMath::Power(N12_sim[jj],2)*N0_sim/TMath::Power(N0_sim,4))<<" | N12 = "<<N12_sim[jj]<<" | N0 = "<<N0_sim<<endl;
    }

    fileName += "_outGr.root";

    TFile *outfile = new TFile(fileName.Data(),"RECREATE");
    h1->Write();
    h2->Write();
    h3->Write();
    h4->Write();

    TCanvas* c1 = new TCanvas();
    c1->cd();
    TGraphErrors* gr = new TGraphErrors(nPoint,energyThreshold,F12,0,errF12);
    gr->GetXaxis()->SetTitle("Energy threshold, keV");
    gr->GetYaxis()->SetTitle("F12");
    gr->SetLineColor(kRed);
    gr->SetLineWidth(2);
    gr->SetMarkerColor(kBlack);
    gr->SetMarkerStyle(20);
    gr->Draw("AP");
    c1->Write();

    outfile->Close();

    h1->Delete();
    h2->Delete();
    h3->Delete();
    h4->Delete();
}

