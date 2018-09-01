#include "src/data_ana.hh"
#include "src/includes.hh"

using namespace std;

int main()
{
    // MAY 2017 STF110

    TString name_ch = "./output/STF110/CH_FINI_"; name_ch += data_ana::theta_cut*1e6; name_ch += "_urad.txt";
    TString name_am = "./output/STF110/AM_FINI_"; name_am += data_ana::theta_cut*1e6; name_am += "_urad.txt";
    TString name_bg = "./output/STF110/BG_FINI_"; name_bg += data_ana::theta_cut*1e6; name_bg += "_urad.txt";
    TString name_ch2 = "./output/STF110/CH_PINI_"; name_ch2 += data_ana::theta_cut*1e6; name_ch2 += "_urad.txt";
    TString name_am2 = "./output/STF110/AM_PINI_"; name_am2 += data_ana::theta_cut*1e6; name_am2 += "_urad.txt";

    ofstream outfile_ch(name_ch.Data());
    ofstream outfile_am(name_am.Data());
    ofstream outfile_bg(name_bg.Data());
    ofstream outfile_ch2(name_ch2.Data());
    ofstream outfile_am2(name_am2.Data());

    time_t start_time, stop_time;
    start_time = time(NULL);

    const Int_t number_of_cutting_angel_center_position = data_ana::number_of_cutting_angel_center_position;

    //==================================================================================================================================================//
    cout<<endl<<"--> Channeling <--"<<endl;

    const Int_t nFiles_ch = 3;
    TString fileName_ch[nFiles_ch];
    fileName_ch[0] = "/run/media/anatochi/6094557b-7c6a-4bfd-9116-c0c405e03e59/anatochi/TB_22_07_17/recoDataSimple_4375_5plane_0_ini.root";
    fileName_ch[1] = "/run/media/anatochi/6094557b-7c6a-4bfd-9116-c0c405e03e59/anatochi/TB_22_07_17/recoDataSimple_4370_5plane_0_ini.root";
    fileName_ch[2] = "/run/media/anatochi/6094557b-7c6a-4bfd-9116-c0c405e03e59/anatochi/TB_22_07_17/recoDataSimple_4372_5plane_0_ini.root";

    TH2D* histo1 = new TH2D("histo1","#Delta#Theta_{X} vs #Theta^{in}_{X}",5000,-500.0,500.0,5000,-500.0,500.0);

    Double_t N12_ch_sum[number_of_cutting_angel_center_position]    = {};
    Double_t N0_ch_sum[number_of_cutting_angel_center_position]     = {};
    Double_t Fin_ch[number_of_cutting_angel_center_position]    = {};
    Double_t errFin_ch[number_of_cutting_angel_center_position] = {};

    data_ana *pointer_ch[nFiles_ch];

    for(Int_t i = 0; i < nFiles_ch; i++)
    {
        cout<<"--> File: "<<i+1<<"/"<<nFiles_ch<<endl;
        cout<<"--> File name: "<<fileName_ch[i]<<endl;
        pointer_ch[i] = new data_ana(fileName_ch[i]);

        Double_t* N12_ch    = new Double_t[number_of_cutting_angel_center_position];
        Double_t* N0_ch     = new Double_t[number_of_cutting_angel_center_position];

        pointer_ch[i]->analysis(N12_ch,N0_ch,-7.29589e-6,1.36759e-6,true);
        pointer_ch[i]->fillhisto(histo1,-7.29589e-6,1.36759e-6);

        for(Int_t j = 0; j < number_of_cutting_angel_center_position; j++)
        {
            N12_ch_sum[j] += N12_ch[j];
            N0_ch_sum[j] += N0_ch[j];
        }
        delete N12_ch;
        delete N0_ch;
    }
    //==================================================================================================================================================//
    cout<<endl<<"--> Amorphous <--"<<endl;

    const Int_t nFiles_am = 2;
    TString fileName_am[nFiles_am];
    fileName_am[0] = "/run/media/anatochi/6094557b-7c6a-4bfd-9116-c0c405e03e59/anatochi/TB_22_07_17/recoDataSimple_4373_5plane_0_ini.root";
    fileName_am[1] = "/run/media/anatochi/6094557b-7c6a-4bfd-9116-c0c405e03e59/anatochi/TB_22_07_17/recoDataSimple_4374_5plane_0_ini.root";

    Double_t N12_am_sum = 0;
    Double_t N0_am_sum = 0;
    Double_t Fin_am[number_of_cutting_angel_center_position]    = {};
    Double_t errFin_am[number_of_cutting_angel_center_position] = {};

    data_ana *pointer_am[nFiles_am];

    for(Int_t i = 0; i < nFiles_am; i++)
    {
        cout<<"--> File: "<<i+1<<"/"<<nFiles_am<<endl;
        cout<<"--> File name: "<<fileName_am[i]<<endl;
        pointer_am[i] = new data_ana(fileName_am[i]);

        Double_t N12_am = 0, N0_am = 0;

        pointer_am[i]->analysis2(N12_am,N0_am);

        for(Int_t j = 0; j < number_of_cutting_angel_center_position; j++)
        {
            N12_am_sum += N12_am;
            N0_am_sum += N0_am;
        }
    }
    //==================================================================================================================================================//
    cout<<endl<<"--> Background <--"<<endl;

    const Int_t nFiles_bg = 5;
    TString fileName_bg[nFiles_bg];
    fileName_bg[0] = "/run/media/anatochi/6094557b-7c6a-4bfd-9116-c0c405e03e59/anatochi/TB_22_07_17/recoDataSimple_4363_5plane_0_ini.root";
    fileName_bg[1] = "/run/media/anatochi/6094557b-7c6a-4bfd-9116-c0c405e03e59/anatochi/TB_22_07_17/recoDataSimple_4371_5plane_0_ini.root";
    fileName_bg[2] = "/run/media/anatochi/6094557b-7c6a-4bfd-9116-c0c405e03e59/anatochi/TB_22_07_17/recoDataSimple_4376_5plane_0_ini.root";
    fileName_bg[3] = "/run/media/anatochi/6094557b-7c6a-4bfd-9116-c0c405e03e59/anatochi/TB_22_07_17/recoDataSimple_4377_5plane_0_ini.root";
    fileName_bg[4] = "/run/media/anatochi/6094557b-7c6a-4bfd-9116-c0c405e03e59/anatochi/TB_22_07_17/recoDataSimple_4388_5plane_0_ini.root";

    Double_t N12_bg_sum = 0;
    Double_t N0_bg_sum = 0;
    Double_t Fin_bg[number_of_cutting_angel_center_position]    = {};
    Double_t errFin_bg[number_of_cutting_angel_center_position] = {};

    data_ana *pointer_bg[nFiles_bg];

    for(Int_t i = 0; i < nFiles_bg; i++)
    {
        cout<<"--> File: "<<i+1<<"/"<<nFiles_bg<<endl;
        cout<<"--> File name: "<<fileName_bg[i]<<endl;
        pointer_bg[i] = new data_ana(fileName_bg[i]);

        Double_t N12_bg = 0, N0_bg = 0;

        pointer_bg[i]->analysis2(N12_bg,N0_bg);

        for(Int_t j = 0; j < number_of_cutting_angel_center_position; j++)
        {
            N12_bg_sum  += N12_bg;
            N0_bg_sum   += N0_bg;
        }
    }
    //==================================================================================================================================================//
    Double_t Cutting_angle[number_of_cutting_angel_center_position] = {};
    Double_t errCutting_angle[number_of_cutting_angel_center_position] = {};

    cout<<"theta_cut_center | Fin_ch | errFin_ch | Fin_am | errFin_am | Fin_bg | errFin_bg"<<endl;
    for(Int_t j = 0; j < number_of_cutting_angel_center_position; j++)
    {
        Cutting_angle[j] = pointer_ch[0]->theta_cut_center[j]*1e6;
        errCutting_angle[j] = 1.0;

        Fin_ch[j] = (double)N12_ch_sum[j]/N0_ch_sum[j];
        Fin_am[j] = (double)N12_am_sum/N0_am_sum;
        Fin_bg[j] = (double)N12_bg_sum/N0_bg_sum;

        errFin_ch[j] = (TMath::Sqrt((double)N12_ch_sum[j]/TMath::Power(N0_ch_sum[j],2) + (double)TMath::Power(N12_ch_sum[j],2)*N0_ch_sum[j]/TMath::Power(N0_ch_sum[j],4)));
        errFin_am[j] = (TMath::Sqrt((double)N12_am_sum/TMath::Power(N0_am_sum,2) + (double)TMath::Power(N12_am_sum,2)*N0_am_sum/TMath::Power(N0_am_sum,4)));
        errFin_bg[j] = (TMath::Sqrt((double)N12_bg_sum/TMath::Power(N0_bg_sum,2) + (double)TMath::Power(N12_bg_sum,2)*N0_bg_sum/TMath::Power(N0_bg_sum,4)));

        outfile_ch<<pointer_ch[0]->theta_cut_center[j]*1e6<<" , "<<Fin_ch[j]<<" , "<<errFin_ch[j]<<"\n";
        outfile_am<<pointer_am[0]->theta_cut_center[j]*1e6<<" , "<<Fin_am[j]<<" , "<<errFin_am[j]<<"\n";
        outfile_bg<<pointer_bg[0]->theta_cut_center[j]*1e6<<" , "<<Fin_bg[j]<<" , "<<errFin_bg[j]<<"\n";

        cout<<pointer_ch[0]->theta_cut_center[j]*1e6<<"   "<<Fin_ch[j]<<"   "<<errFin_ch[j]<<"   "<<Fin_am[j]<<"   "<<errFin_am[j]<<"   "<<Fin_bg[j]<<"   "<<errFin_bg[j]<<endl;
    }

    //-------------------//
    TFile *resultsfile = new TFile("./output/STF110/results.root","recreate");

    // Plot of the interaction frequency for different orientation
    TCanvas *c1 = new TCanvas("c1","Canva",200,10,700,500);
    c1->cd();
    c1->SetGrid();
    TMultiGraph * mg = new TMultiGraph();

    TGraph *gr_CH = new TGraphErrors(number_of_cutting_angel_center_position,Cutting_angle,Fin_ch,errCutting_angle,errFin_ch);
    TGraph *gr_AM = new TGraphErrors(number_of_cutting_angel_center_position,Cutting_angle,Fin_am,errCutting_angle,errFin_am);
    TGraph *gr_BG = new TGraphErrors(number_of_cutting_angel_center_position,Cutting_angle,Fin_bg,errCutting_angle,errFin_bg);

    gr_CH->SetName("gr_CH");
    gr_CH->SetTitle("Channeling");
    gr_CH->SetMarkerStyle(20);
    gr_CH->SetMarkerColor(kBlack);
    gr_CH->SetLineColor(kBlue);
    gr_CH->SetLineWidth(2);
    gr_CH->SetFillStyle(0);

    gr_AM->SetName("gr_AM");
    gr_AM->SetTitle("Amorphous orientation");
    gr_AM->SetMarkerStyle(21);
    gr_AM->SetMarkerColor(kBlack);
    gr_AM->SetLineColor(kBlue);
    gr_AM->SetLineWidth(2);
    gr_AM->SetFillStyle(0);

    gr_BG->SetName("gr_BG");
    gr_BG->SetTitle("Hight statistic BG");
    gr_BG->SetMarkerStyle(23);
    gr_BG->SetMarkerColor(kBlack);
    gr_BG->SetLineColor(kBlue);
    gr_BG->SetLineWidth(2);
    gr_BG->SetFillStyle(0);

    mg->Add( gr_CH );
    mg->Add( gr_AM );
    mg->Add( gr_BG );

    mg->Draw("APC");
    mg->SetTitle("Interaction Frequency");
    mg->GetXaxis()->SetTitle("Cutting Angle (#murad)");
    mg->GetYaxis()->SetTitle("Interaction frequency");

    c1->BuildLegend();
    c1->Update();
    c1->Write();

    gr_CH->Write();
    gr_AM->Write();

    Double_t F12 = 1.0;
    Double_t errF12 = 0.0;

    // Calculation of the interaction probability with high statistical background
    Double_t IntProb_CH_BG[number_of_cutting_angel_center_position] = {};
    Double_t IntProb_AM_BG[number_of_cutting_angel_center_position] = {};
    Double_t errIntProb_CH_BG[number_of_cutting_angel_center_position] = {};
    Double_t errIntProb_AM_BG[number_of_cutting_angel_center_position] = {};

    for(Int_t i = 0; i < number_of_cutting_angel_center_position; i++)
    {
        IntProb_CH_BG[i] = 100.0*(Fin_ch[i] - Fin_bg[i])/F12;
        IntProb_AM_BG[i] = 100.0*(Fin_am[i] - Fin_bg[i])/F12;

        errIntProb_CH_BG[i] = 100.0*TMath::Sqrt(TMath::Power(errFin_ch[i]/F12,2) + TMath::Power(errFin_bg[i]/F12,2) +
                                       TMath::Power(errF12*(Fin_ch[i] - Fin_bg[i])/(F12*F12),2));
        errIntProb_AM_BG[i] = 100.0*TMath::Sqrt(TMath::Power(errFin_am[i]/F12,2) + TMath::Power(errFin_bg[i]/F12,2) +
                                       TMath::Power(errF12*(Fin_am[i] - Fin_bg[i])/(F12*F12),2));

        outfile_ch2<<pointer_ch[0]->theta_cut_center[i]*1e6<<" , "<<IntProb_CH_BG[i]<<" , "<<errIntProb_CH_BG[i]<<"\n";
        outfile_am2<<pointer_am[0]->theta_cut_center[i]*1e6<<" , "<<IntProb_AM_BG[i]<<" , "<<errIntProb_AM_BG[i]<<"\n";
    }
    TCanvas *c2 = new TCanvas("c2","Canva",200,10,700,500);
    c2->cd();
    c2->SetGrid();

    TMultiGraph * mg_prob_BG = new TMultiGraph();

    TGraph *gr_CH_prob_BG = new TGraphErrors(number_of_cutting_angel_center_position,Cutting_angle,IntProb_CH_BG,errCutting_angle,errIntProb_CH_BG);
    TGraph *gr_AM_prob_BG = new TGraphErrors(number_of_cutting_angel_center_position,Cutting_angle,IntProb_AM_BG,errCutting_angle,errIntProb_AM_BG);

    gr_CH_prob_BG->SetName("gr_CH_prob");
    gr_CH_prob_BG->SetTitle("Channeling");
    gr_CH_prob_BG->SetMarkerStyle(20);
    gr_CH_prob_BG->SetMarkerColor(kBlack);
    gr_CH_prob_BG->SetLineColor(kBlue);
    gr_CH_prob_BG->SetLineWidth(2);
    gr_CH_prob_BG->SetFillStyle(0);

    gr_AM_prob_BG->SetName("gr_AM_prob");
    gr_AM_prob_BG->SetTitle("Amorphous orientation");
    gr_AM_prob_BG->SetMarkerStyle(21);
    gr_AM_prob_BG->SetMarkerColor(kBlack);
    gr_AM_prob_BG->SetLineColor(kBlue);
    gr_AM_prob_BG->SetLineWidth(2);
    gr_AM_prob_BG->SetFillStyle(0);

    mg_prob_BG->Add( gr_CH_prob_BG );
    mg_prob_BG->Add( gr_AM_prob_BG );

    mg_prob_BG->Draw("APC");
    mg_prob_BG->SetTitle("Interaction Probability (high statistical background)");
    mg_prob_BG->GetXaxis()->SetTitle("Cutting Angle (#murad)");
    mg_prob_BG->GetYaxis()->SetTitle("Interaction probability (%)");

    c2->BuildLegend();
    //mg_prob_BG->SetMaximum(0.25);
    //mg_prob_BG->SetMinimum(0.0);
    c2->Update();
    c2->Write();

    gr_CH_prob_BG->Write();
    gr_AM_prob_BG->Write();

    // Calculation of the interaction probability ratio PCH/AM
    cout<<"--> Calculation of the interaction probability ratio PCH/AM"<<endl;
    Double_t Ratio_CH_AM[number_of_cutting_angel_center_position] = {};
    Double_t errRatio_CH_AM[number_of_cutting_angel_center_position] = {};

    for(Int_t i = 0; i < number_of_cutting_angel_center_position; i++)
    {
        Ratio_CH_AM[i] = IntProb_CH_BG[i]/IntProb_AM_BG[i];

        errRatio_CH_AM[i] = TMath::Sqrt(TMath::Power(errIntProb_CH_BG[i]/IntProb_AM_BG[i],2) + TMath::Power(IntProb_CH_BG[i]*errIntProb_AM_BG[i]/(IntProb_AM_BG[i]*IntProb_AM_BG[i]),2));
        cout<<Cutting_angle[i]<<"  "<<Ratio_CH_AM[i]<<"  "<<errRatio_CH_AM[i]<<endl;
    }
    TGraph *gr_ratio_CH_AM = new TGraphErrors(number_of_cutting_angel_center_position,Cutting_angle,Ratio_CH_AM,errCutting_angle,errRatio_CH_AM);
    gr_ratio_CH_AM->SetName("gr_ratio_CH_AM");
    gr_ratio_CH_AM->SetTitle("PCH/AM");
    gr_ratio_CH_AM->SetMarkerStyle(20);
    gr_ratio_CH_AM->SetMarkerColor(kBlack);
    gr_ratio_CH_AM->SetLineColor(kBlue);
    gr_ratio_CH_AM->SetLineWidth(2);
    gr_ratio_CH_AM->SetFillStyle(0);
    gr_ratio_CH_AM->Write();

    histo1->Write();

    resultsfile->Close();

    //-------------------//


    outfile_ch.close();
    outfile_am.close();
    outfile_bg.close();
    outfile_ch2.close();
    outfile_am2.close();

    stop_time = time(NULL);
    cout<<"--> Running time is : "<<stop_time - start_time<<" seconds"<<endl;
    return 0;
}
