#include "src/includes.hh"

using namespace std;

int plot()
{
    const Int_t number_of_files = 9;
    const Int_t number_of_scintillators = 6;
//    const Int_t number_of_scintillators = 6/2;

    //---------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------//

    TString filename[] =
    {
        "output/recoDataSimple_6363_4plane_offset0_0.root_COUNTS.txt",
        "output/recoDataSimple_6360_4plane_offset0_0.root_COUNTS.txt",
        "output/recoDataSimple_6355_4plane_offset0_0.root_COUNTS.txt",
        "output/recoDataSimple_6351_4plane_offset0_0.root_COUNTS.txt",
        "output/recoDataSimple_6354_4plane_offset0_0.root_COUNTS.txt",
        "output/recoDataSimple_6357_4plane_offset0_0.root_COUNTS.txt",
        "output/recoDataSimple_6361_4plane_offset0_0.root_COUNTS.txt",
        "output/recoDataSimple_6364_4plane_offset0_0.root_COUNTS.txt",
        "output/recoDataSimple_6365_4plane_offset0_0.root_COUNTS.txt"

        /*"output/recoDataSimple_6363_4plane_offset0_0.root_COUNTS_PAIRS.txt",
        "output/recoDataSimple_6360_4plane_offset0_0.root_COUNTS_PAIRS.txt",
        "output/recoDataSimple_6355_4plane_offset0_0.root_COUNTS_PAIRS.txt",
        "output/recoDataSimple_6351_4plane_offset0_0.root_COUNTS_PAIRS.txt",
        "output/recoDataSimple_6354_4plane_offset0_0.root_COUNTS_PAIRS.txt",
        "output/recoDataSimple_6357_4plane_offset0_0.root_COUNTS_PAIRS.txt",
        "output/recoDataSimple_6361_4plane_offset0_0.root_COUNTS_PAIRS.txt",
        "output/recoDataSimple_6364_4plane_offset0_0.root_COUNTS_PAIRS.txt",
        "output/recoDataSimple_6365_4plane_offset0_0.root_COUNTS_PAIRS.txt"*/
    };

    Double_t position [] =
    {
        66.970,
        67.020,
        67.045,
        67.070,
        67.095,
        67.120,
        67.170,
        67.270,
        67.570
    };

    Int_t det_id;
    TString word;
    Double_t counts[number_of_scintillators][number_of_files] = {};
    Double_t counts_err[number_of_scintillators][number_of_files] = {};

    Double_t alnm_counts[number_of_scintillators][number_of_files] = {};
    Double_t alnm_counts_err[number_of_scintillators][number_of_files] = {};

    for(Int_t i = 0; i < number_of_files; i++)
    {
        cout<<"--> Inputfile: "<<filename[i].Data()<<endl;
        ifstream inputfile(filename[i].Data());
        if(inputfile.is_open())
        {
            while(1)
            {
                inputfile>>det_id;
                inputfile>>word;
                inputfile>>counts[det_id][i];
                inputfile>>word;
                inputfile>>counts_err[det_id][i];

                if(inputfile.eof()) {break;}
            }
            inputfile.close();
        }
        else
        {
            cout<<"--> ERROR: Unable to open the input file!"<<endl;
            assert(0);
        }
    }

    //---------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------//
    // ALINMENT RUNS

    TString filename_alnm[] =
    {
        "output/recoDataSimple_6362_4plane_offset0_6363.root_COUNTS.txt",
        "output/recoDataSimple_6362_4plane_offset0_6360.root_COUNTS.txt",
        "output/recoDataSimple_6356_4plane_offset0_6355.root_COUNTS.txt",
        "output/recoDataSimple_6350_4plane_offset0_6351.root_COUNTS.txt",
        "output/recoDataSimple_6356_4plane_offset0_6354.root_COUNTS.txt",
        "output/recoDataSimple_6356_4plane_offset0_6357.root_COUNTS.txt",
        "output/recoDataSimple_6362_4plane_offset0_6361.root_COUNTS.txt",
        "output/recoDataSimple_6362_4plane_offset0_6364.root_COUNTS.txt",
        "output/recoDataSimple_6362_4plane_offset0_6365.root_COUNTS.txt"

        /*"output/recoDataSimple_6362_4plane_offset0_6363.root_COUNTS_PAIRS.txt",
        "output/recoDataSimple_6362_4plane_offset0_6360.root_COUNTS_PAIRS.txt",
        "output/recoDataSimple_6356_4plane_offset0_6355.root_COUNTS_PAIRS.txt",
        "output/recoDataSimple_6350_4plane_offset0_6351.root_COUNTS_PAIRS.txt",
        "output/recoDataSimple_6356_4plane_offset0_6354.root_COUNTS_PAIRS.txt",
        "output/recoDataSimple_6356_4plane_offset0_6357.root_COUNTS_PAIRS.txt",
        "output/recoDataSimple_6362_4plane_offset0_6361.root_COUNTS_PAIRS.txt",
        "output/recoDataSimple_6362_4plane_offset0_6364.root_COUNTS_PAIRS.txt",
        "output/recoDataSimple_6362_4plane_offset0_6365.root_COUNTS_PAIRS.txt"*/
    };

    for(Int_t i = 0; i < number_of_files; i++)
    {
        cout<<"--> Inputfile: "<<filename_alnm[i].Data()<<endl;
        ifstream inputfile(filename_alnm[i].Data());
        if(inputfile.is_open())
        {
            while(1)
            {
                inputfile>>det_id;
                inputfile>>word;
                inputfile>>alnm_counts[det_id][i];
                inputfile>>word;
                inputfile>>alnm_counts_err[det_id][i];

                if(inputfile.eof()) {break;}
            }
            inputfile.close();
        }
        else
        {
            cout<<"--> ERROR: Unable to open the input file!"<<endl;
            assert(0);
        }
    }
    //---------------------------------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------//

    TGraphErrors* gr[number_of_scintillators];
    TGraphErrors* gr2[number_of_scintillators];
    TGraphErrors* alnm_gr[number_of_scintillators];
    TMultiGraph* mg = new TMultiGraph();
    TMultiGraph* mg2 = new TMultiGraph();

    for(Int_t i = 0; i < number_of_scintillators; i++)
    {
        //if(i == 0 || i == 1){

        gr2[i] = new TGraphErrors(number_of_files,position,counts[i],0,counts_err[i]);

        for(Int_t yy = 0; yy < number_of_files; yy++)
        {
            counts[i][yy] = counts[i][yy] - alnm_counts[i][yy];
            counts_err[i][yy] = TMath::Sqrt(TMath::Power(counts_err[i][yy],2) + TMath::Power(alnm_counts_err[i][yy],2));
        }

        gr[i] = new TGraphErrors(number_of_files,position,counts[i],0,counts_err[i]);

        TString gr_name = "gr_";
        gr_name += i;
        gr[i]->SetName(gr_name.Data());
        gr2[i]->SetName(gr_name.Data());
        TString gr_title;
        if(number_of_scintillators == 6)
        {
            gr_title= "signal";
            gr_title += i;
        }
        else
        {
            gr_title = "INI station";
            gr_title += i+1;
        }
        gr[i]->SetTitle(gr_title.Data());
        gr2[i]->SetTitle(gr_title.Data());
        gr[i]->SetMarkerStyle(20+i);
        gr2[i]->SetMarkerStyle(20+i);
        gr[i]->SetLineColor(1+i);
        gr2[i]->SetLineColor(1+i);
        gr[i]->SetMarkerColor(1+i);
        gr2[i]->SetMarkerColor(1+i);
        gr[i]->SetLineWidth(2);
        gr2[i]->SetLineWidth(2);
        gr[i]->SetLineStyle(1);
        gr2[i]->SetLineStyle(1);
        gr[i]->SetFillColor(0);
        gr2[i]->SetFillColor(0);

        mg->Add(gr[i]);
        mg2->Add(gr2[i]);
        //}
    }

    for(Int_t i = 0; i < number_of_scintillators; i++)
    {
        alnm_gr[i] = new TGraphErrors(number_of_files,position,alnm_counts[i],0,alnm_counts_err[i]);

        TString alnm_gr_name = "alnm_gr_";
        alnm_gr_name += i;
        alnm_gr[i]->SetName(alnm_gr_name.Data());
        TString alnm_gr_title;
        if(number_of_scintillators == 6)
        {
            alnm_gr_title = "BKG INI station";
        }
        else
        {
            alnm_gr_title = "BKG signal";
        }
        alnm_gr_title += i+1;
        alnm_gr[i]->SetTitle(alnm_gr_title.Data());
        alnm_gr[i]->SetMarkerStyle(20+i);
        alnm_gr[i]->SetLineColor(1+i);
        alnm_gr[i]->SetMarkerColor(1+i);
        alnm_gr[i]->SetLineWidth(2);
        alnm_gr[i]->SetLineStyle(9);
        alnm_gr[i]->SetFillColor(0);

        mg2->Add(alnm_gr[i]);
    }

    TCanvas* c1 = new TCanvas();
    c1->cd();
    mg->SetTitle("Scintillator counts, normalized by the total number of incoming particles");
    mg->Draw("APL");
    c1->BuildLegend();

    TCanvas* c2 = new TCanvas();
    c2->cd();
    mg2->SetTitle("Scintillator counts, normalized by the total number of incoming particles");
    mg2->Draw("APL");
    c2->BuildLegend();

    return 0;
}

