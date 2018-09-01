#include "src/includes.hh"

using namespace std;

int plot_FINI()
{
    const Int_t number_of_files     = 3;
    const Int_t number_of_angles    = 32;

    TString filename[] =
    {
        "output/recoDataSimple_6381_4plane_offset0_0.root_FINI.txt",
        "output/recoDataSimple_6374_4plane_offset0_0.root_FINI.txt",
        "output/recoDataSimple_6382_4plane_offset0_0.root_FINI.txt"
    };

    TString word;
    Double_t angle[number_of_angles] = {};
    Double_t FINI[number_of_files][number_of_angles] = {};
    Double_t FINI_err[number_of_files][number_of_angles] = {};

    for(Int_t i = 0; i < number_of_files; i++)
    {
        cout<<"--> Inputfile: "<<filename[i].Data()<<endl;
        ifstream inputfile(filename[i].Data());
        Int_t angle_id = 0;
        if(inputfile.is_open())
        {
            while(1)
            {
                inputfile>>angle[angle_id];
                inputfile>>word;
                inputfile>>word;
                inputfile>>word;
                inputfile>>word;
                inputfile>>word;
                inputfile>>FINI[i][angle_id];
                inputfile>>word;
                inputfile>>FINI_err[i][angle_id];

                angle_id++;

                if(inputfile.eof()) {break;}
            }
            inputfile.close();
        }
        else
        {
            cout<<"--> ERROR: Unable to open the input file!"<<endl;
            assert(0);
        }
        if(angle_id-1 == number_of_angles) cout<<"--> Matched =)"<<endl;
        else
        {
            cout<<"--> ERROR: Something wrong with the number of angles: "<<angle_id<<" != "<<number_of_angles<<endl;
            assert(0);
        }
    }

    TGraphErrors* gr_fini[number_of_files];
    TMultiGraph* mg = new TMultiGraph();

    for(Int_t i = 0; i < number_of_files; i++)
    {
        gr_fini[i] = new TGraphErrors(10,angle,FINI[i],0,FINI_err[i]);

        TString gr_name = "gr_";
        gr_name += i;
        gr_fini[i]->SetName(gr_name.Data());
        if(i == 0) gr_fini[i]->SetTitle("F_{BKG}");
        if(i == 1) gr_fini[i]->SetTitle("F_{AM}");
        if(i == 2) gr_fini[i]->SetTitle("F_{CH}");
        gr_fini[i]->SetMarkerStyle(20+i);
        gr_fini[i]->SetLineColor(1+i);
        gr_fini[i]->SetMarkerColor(1+i);
        gr_fini[i]->SetLineWidth(2);
        gr_fini[i]->SetLineStyle(1);
        gr_fini[i]->SetFillColor(0);

        mg->Add(gr_fini[i]);
    }

    TCanvas* c1 = new TCanvas();
    c1->cd();
    mg->SetTitle("F_{INI}");
    mg->Draw("APL");
    c1->BuildLegend();

    Double_t PINI_CH[number_of_angles] = {};
    Double_t PINI_CH_err[number_of_angles] = {};
    Double_t PINI_AM[number_of_angles] = {};
    Double_t PINI_AM_err[number_of_angles] = {};
    Double_t F12 = 0.132; // from 2015 for 2 mm

    for(Int_t j = 0; j < number_of_angles; j++)
    {
        PINI_CH[j] = 100.0*(FINI[2][j] - FINI[0][j])/F12;
        PINI_AM[j] = 100.0*(FINI[1][j] - FINI[0][j])/F12;

        PINI_CH_err[j] = 100.0*TMath::Sqrt(TMath::Power(FINI_err[2][j]/F12,2) + TMath::Power(FINI_err[0][j]/F12,2));
        PINI_AM_err[j] = 100.0*TMath::Sqrt(TMath::Power(FINI_err[1][j]/F12,2) + TMath::Power(FINI_err[0][j]/F12,2));
    }
    TMultiGraph* mg2 = new TMultiGraph();
    TGraphErrors* gr_pini_ch = new TGraphErrors(number_of_angles,angle,PINI_CH,0,PINI_CH_err);
    TGraphErrors* gr_pini_am = new TGraphErrors(number_of_angles,angle,PINI_AM,0,PINI_AM_err);
    gr_pini_ch->SetTitle("PCH");
    gr_pini_am->SetTitle("AM");

    gr_pini_ch->SetMarkerStyle(22);
    gr_pini_am->SetMarkerStyle(21);
    gr_pini_ch->SetLineColor(2);
    gr_pini_am->SetLineColor(3);
    gr_pini_ch->SetMarkerColor(2);
    gr_pini_am->SetMarkerColor(3);
    gr_pini_ch->SetLineWidth(2);
    gr_pini_am->SetLineWidth(2);
    gr_pini_ch->SetLineStyle(1);
    gr_pini_am->SetLineStyle(1);
    gr_pini_ch->SetFillColor(0);
    gr_pini_am->SetFillColor(0);

    mg2->Add(gr_pini_ch);
    mg2->Add(gr_pini_am);

    TCanvas* c2 = new TCanvas();
    c2->cd();
    mg2->SetTitle("P_{INI}");
    mg2->SetMaximum(0.55);
    mg2->SetMinimum(0.0);
    mg2->Draw("APL");
    c2->BuildLegend();

    return 0;
}
