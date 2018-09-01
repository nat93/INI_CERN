#include "src/data_ana.hh"
#include "src/includes.hh"

using namespace std;

int main()
{
    // PLOTST:    1 - histogramms  | 2 - graphics     | 3 - all together (1 & 2) |

    const Int_t plotST = 2;

    TString dir_name = "./output";
    cout<<endl<<"--> Output directory: "<<dir_name<<endl;
    struct stat st = {0};
    if (stat(dir_name.Data(), &st) == -1)
    {
        mkdir(dir_name.Data(), 0700);
    }

    TString name1 = dir_name; name1 += "/CH_FRQ.txt";
    TString name2 = dir_name; name2 += "/AM_FRQ.txt";
    TString name3 = dir_name; name3 += "/BG_FRQ.txt";
    TString name4 = dir_name; name4 += "/CH_PRB.txt";
    TString name5 = dir_name; name5 += "/AM_PRB.txt";
    TString name6 = dir_name; name6 += "/CH_AM_RT.txt";

    ofstream outfile_pl_frq(name1.Data());
    ofstream outfile_am_frq(name2.Data());
    ofstream outfile_bg_frq(name3.Data());
    ofstream outfile_pl_prb(name4.Data());
    ofstream outfile_am_prb(name5.Data());
    ofstream outfile_ch_am_rt(name6.Data());

    time_t start_time, stop_time;
    start_time = time(NULL);

    const Int_t number_of_cutting_angels = data_ana::number_of_cutting_angels;

    Double_t Cutting_angle[number_of_cutting_angels] = {};
    Double_t errCutting_angle[number_of_cutting_angels] = {};

    //==================================================================================================================================================//
    cout<<endl<<"--> The channeling orientation"<<endl;

    const Int_t nFiles_ch = 2;//7;
    string fileName_ch[nFiles_ch];

    // The channeling orientation
    fileName_ch[0] = "/media/andrii/NATOCHII_HDD/UA9_database/recoDataSimple_5733_5plane_offset0_0.root";
    fileName_ch[1] = "/media/andrii/NATOCHII_HDD/UA9_database/recoDataSimple_5738_5plane_offset0_0.root";
    /*fileName_ch[0] = "/media/andrii/NATOCHII_HDD/UA9_database/recoDataSimple_5743_5plane_offset0_0.root";
    fileName_ch[1] = "/media/andrii/NATOCHII_HDD/UA9_database/recoDataSimple_5745_5plane_offset0_0.root";
    fileName_ch[2] = "/media/andrii/NATOCHII_HDD/UA9_database/recoDataSimple_5750_5plane_offset0_0.root";
    fileName_ch[3] = "/media/andrii/NATOCHII_HDD/UA9_database/recoDataSimple_5752_5plane_offset0_0.root";
    fileName_ch[4] = "/media/andrii/NATOCHII_HDD/UA9_database/recoDataSimple_5752_5plane_offset0_1.root";
    fileName_ch[5] = "/media/andrii/NATOCHII_HDD/UA9_database/recoDataSimple_5752_5plane_offset0_2.root";
    fileName_ch[6] = "/media/andrii/NATOCHII_HDD/UA9_database/recoDataSimple_5757_5plane_offset0_0.root";*/

    Double_t N12_CH[number_of_cutting_angels] = {};
    Double_t N0_CH[number_of_cutting_angels] = {};

    Double_t Fin_CH[number_of_cutting_angels] = {};
    Double_t errFin_CH[number_of_cutting_angels] = {};

    data_ana *pointer_ch[nFiles_ch];

    if(plotST == 1 || plotST == 3)
    {
        for(Int_t i = 0; i < nFiles_ch; i++)
        {
            //-----------------------------------------------------------------------//
            size_t pos = fileName_ch[i].find_last_of('/');
            string filename = fileName_ch[i].substr(pos+1);
            TString name_file_tmp = dir_name;
            name_file_tmp += "/CH_TMP_";
            name_file_tmp += filename;
            name_file_tmp += "_.root";
            //-----------------------------------------------------------------------//
            cout<<"--> File: "<<i+1<<"/"<<nFiles_ch<<endl;
            cout<<"--> File name: "<<fileName_ch[i]<<endl;
            pointer_ch[i] = new data_ana(fileName_ch[i]);

            pointer_ch[i]->plot_histogramms();
            pointer_ch[i]->write_histo_to_file(name_file_tmp.Data());
        }
    }    

    if(plotST == 2 || plotST == 3)
    {
        for(Int_t i = 0; i < nFiles_ch; i++)
        {
            //-----------------------------------------------------------------------//
            size_t pos = fileName_ch[i].find_last_of('/');
            string filename = fileName_ch[i].substr(pos+1);
            TString name_file_tmp = dir_name;
            name_file_tmp += "/CH_FINI_TMP_";
            name_file_tmp += filename;
            name_file_tmp += "_.txt";
            ofstream outfile(name_file_tmp.Data());
            //-----------------------------------------------------------------------//

            cout<<"--> File: "<<i+1<<"/"<<nFiles_ch<<endl;
            cout<<"--> File name: "<<fileName_ch[i]<<endl;
            pointer_ch[i] = new data_ana(fileName_ch[i]);

            Double_t* N12_ch_tmp    = new Double_t[number_of_cutting_angels];
            Double_t* N0_ch_tmp     = new Double_t[number_of_cutting_angels];

            pointer_ch[i]->find_inelastic_frequency3(N12_ch_tmp,N0_ch_tmp);

            for(Int_t j = 0; j < number_of_cutting_angels; j++)
            {
                N12_CH[j] += N12_ch_tmp[j];
                N0_CH[j] += N0_ch_tmp[j];

                //-----------------------------------------------------------------------//
                outfile<<pointer_ch[i]->theta_cut[j]*1e6<<" , "<<(double)N12_ch_tmp[j]/N0_ch_tmp[j]<<" , "<<(TMath::Sqrt((double)N12_ch_tmp[j]/TMath::Power(N0_ch_tmp[j],2) + (double)TMath::Power(N12_ch_tmp[j],2)*N0_ch_tmp[j]/TMath::Power(N0_ch_tmp[j],4)))<<"\n";
                //-----------------------------------------------------------------------//
            }
            delete N12_ch_tmp;
            delete N0_ch_tmp;

            //-----------------------------------------------------------------------//
            outfile.close();
            //-----------------------------------------------------------------------//
        }
    }    

    //==================================================================================================================================================//
    cout<<endl<<"--> The amorphous orientation"<<endl;

    const Int_t nFiles_am = 5;
    string fileName_am[nFiles_am];

    // The amorphous orientation
    fileName_am[0] = "/media/andrii/NATOCHII_HDD/UA9_database/recoDataSimple_5736_5plane_offset0_0.root";
    fileName_am[1] = "/media/andrii/NATOCHII_HDD/UA9_database/recoDataSimple_5736_5plane_offset0_1.root";
    fileName_am[2] = "/media/andrii/NATOCHII_HDD/UA9_database/recoDataSimple_5736_5plane_offset0_2.root";
    fileName_am[3] = "/media/andrii/NATOCHII_HDD/UA9_database/recoDataSimple_5739_5plane_offset0_0.root";
    fileName_am[4] = "/media/andrii/NATOCHII_HDD/UA9_database/recoDataSimple_5739_5plane_offset0_1.root";
    //fileName_am[0] = "/media/andrii/NATOCHII_HDD/UA9_database/recoDataSimple_5748_5plane_offset0_0.root";
    //fileName_am[1] = "/media/andrii/NATOCHII_HDD/UA9_database/recoDataSimple_5755_5plane_offset0_0.root";
    //fileName_am[2] = "/media/andrii/NATOCHII_HDD/UA9_database/recoDataSimple_5755_5plane_offset0_1.root";


    Double_t N12_AM[number_of_cutting_angels] = {};
    Double_t N0_AM[number_of_cutting_angels] = {};

    Double_t Fin_AM[number_of_cutting_angels] = {};
    Double_t errFin_AM[number_of_cutting_angels] = {};

    data_ana *pointer_am[nFiles_am];

    if(plotST == 1 || plotST == 3)
    {
        for(Int_t i = 0; i < nFiles_am; i++)
        {
            //-----------------------------------------------------------------------//
            size_t pos = fileName_am[i].find_last_of('/');
            string filename = fileName_am[i].substr(pos+1);
            TString name_file_tmp = dir_name;
            name_file_tmp += "/AM_TMP_";
            name_file_tmp += filename;
            name_file_tmp += "_.root";
            //-----------------------------------------------------------------------//
            cout<<"--> File: "<<i+1<<"/"<<nFiles_am<<endl;
            cout<<"--> File name: "<<fileName_am[i]<<endl;
            pointer_am[i] = new data_ana(fileName_am[i]);

            pointer_am[i]->plot_histogramms();
            pointer_am[i]->write_histo_to_file(name_file_tmp.Data());
        }
    }

    if(plotST == 2 || plotST == 3)
    {
        for(Int_t i = 0; i < nFiles_am; i++)
        {
            //-----------------------------------------------------------------------//
            size_t pos = fileName_am[i].find_last_of('/');
            string filename = fileName_am[i].substr(pos+1);
            TString name_file_tmp = dir_name;
            name_file_tmp += "/AM_FINI_TMP_";
            name_file_tmp += filename;
            name_file_tmp += "_.txt";
            ofstream outfile(name_file_tmp.Data());
            //-----------------------------------------------------------------------//

            cout<<"--> File: "<<i+1<<"/"<<nFiles_am<<endl;
            cout<<"--> File name: "<<fileName_am[i]<<endl;
            pointer_am[i] = new data_ana(fileName_am[i]);

            Double_t* N12_am_tmp    = new Double_t[number_of_cutting_angels];
            Double_t* N0_am_tmp     = new Double_t[number_of_cutting_angels];

            pointer_am[i]->find_inelastic_frequency2(N12_am_tmp,N0_am_tmp);

            for(Int_t j = 0; j < number_of_cutting_angels; j++)
            {
                N12_AM[j] += N12_am_tmp[j];
                N0_AM[j] += N0_am_tmp[j];

                //-----------------------------------------------------------------------//
                outfile<<pointer_am[i]->theta_cut[j]*1e6<<" , "<<(double)N12_am_tmp[j]/N0_am_tmp[j]<<" , "<<(TMath::Sqrt((double)N12_am_tmp[j]/TMath::Power(N0_am_tmp[j],2) + (double)TMath::Power(N12_am_tmp[j],2)*N0_am_tmp[j]/TMath::Power(N0_am_tmp[j],4)))<<"\n";
                //-----------------------------------------------------------------------//
            }
            delete N12_am_tmp;
            delete N0_am_tmp;

            //-----------------------------------------------------------------------//
            outfile.close();
            //-----------------------------------------------------------------------//
        }
    }
    //==================================================================================================================================================//
    cout<<endl<<"--> Background"<<endl;

    const Int_t nFiles_bg = 5;
    string fileName_bg[nFiles_bg];

    // Background
    fileName_bg[0] = "/media/andrii/NATOCHII_HDD/UA9_database/recoDataSimple_5734_5plane_offset0_0.root";
    fileName_bg[1] = "/media/andrii/NATOCHII_HDD/UA9_database/recoDataSimple_5735_5plane_offset0_0.root";
    fileName_bg[2] = "/media/andrii/NATOCHII_HDD/UA9_database/recoDataSimple_5735_5plane_offset0_1.root";
    fileName_bg[3] = "/media/andrii/NATOCHII_HDD/UA9_database/recoDataSimple_5737_5plane_offset0_0.root";
    fileName_bg[4] = "/media/andrii/NATOCHII_HDD/UA9_database/recoDataSimple_5740_5plane_offset0_0.root";
    //fileName_bg[0] = "/media/andrii/NATOCHII_HDD/UA9_database/recoDataSimple_5742_5plane_offset0_0.root";
    //fileName_bg[1] = "/media/andrii/NATOCHII_HDD/UA9_database/recoDataSimple_5742_5plane_offset0_1.root";
    //fileName_bg[2] = "/media/andrii/NATOCHII_HDD/UA9_database/recoDataSimple_5756_5plane_offset0_0.root";
    //fileName_bg[3] = "/media/andrii/NATOCHII_HDD/UA9_database/recoDataSimple_5758_5plane_offset0_0.root";
    //fileName_bg[4] = "/media/andrii/NATOCHII_HDD/UA9_database/recoDataSimple_5758_5plane_offset0_1.root";


    Double_t N12_BG[number_of_cutting_angels] = {};
    Double_t N0_BG[number_of_cutting_angels] = {};

    Double_t Fin_BG[number_of_cutting_angels] = {};
    Double_t errFin_BG[number_of_cutting_angels] = {};

    data_ana *pointer_bg[nFiles_bg];

    if(plotST == 1 || plotST == 3)
    {
        for(Int_t i = 0; i < nFiles_bg; i++)
        {
            //-----------------------------------------------------------------------//
            size_t pos = fileName_bg[i].find_last_of('/');
            string filename = fileName_bg[i].substr(pos+1);
            TString name_file_tmp = dir_name;
            name_file_tmp += "/BG_TMP_";
            name_file_tmp += filename;
            name_file_tmp += "_.root";
            //-----------------------------------------------------------------------//
            cout<<"--> File: "<<i+1<<"/"<<nFiles_bg<<endl;
            cout<<"--> File name: "<<fileName_bg[i]<<endl;
            pointer_bg[i] = new data_ana(fileName_bg[i]);

            pointer_bg[i]->plot_histogramms();
            pointer_bg[i]->write_histo_to_file(name_file_tmp.Data());
        }
    }

    if(plotST == 2 || plotST == 3)
    {
        for(Int_t i = 0; i < nFiles_bg; i++)
        {
            //-----------------------------------------------------------------------//
            size_t pos = fileName_bg[i].find_last_of('/');
            string filename = fileName_bg[i].substr(pos+1);
            TString name_file_tmp = dir_name;
            name_file_tmp += "/BG_FINI_TMP_";
            name_file_tmp += filename;
            name_file_tmp += "_.txt";
            ofstream outfile(name_file_tmp.Data());
            //-----------------------------------------------------------------------//

            cout<<"--> File: "<<i+1<<"/"<<nFiles_bg<<endl;
            cout<<"--> File name: "<<fileName_bg[i]<<endl;
            pointer_bg[i] = new data_ana(fileName_bg[i]);

            Double_t* N12_bg_tmp    = new Double_t[number_of_cutting_angels];
            Double_t* N0_bg_tmp     = new Double_t[number_of_cutting_angels];

            pointer_bg[i]->find_inelastic_frequency2(N12_bg_tmp,N0_bg_tmp);

            for(Int_t j = 0; j < number_of_cutting_angels; j++)
            {
                N12_BG[j] += N12_bg_tmp[j];
                N0_BG[j] += N0_bg_tmp[j];

                //-----------------------------------------------------------------------//
                outfile<<pointer_bg[i]->theta_cut[j]*1e6<<" , "<<(double)N12_bg_tmp[j]/N0_bg_tmp[j]<<" , "<<(TMath::Sqrt((double)N12_bg_tmp[j]/TMath::Power(N0_bg_tmp[j],2) + (double)TMath::Power(N12_bg_tmp[j],2)*N0_bg_tmp[j]/TMath::Power(N0_bg_tmp[j],4)))<<"\n";
                //-----------------------------------------------------------------------//
            }
            delete N12_bg_tmp;
            delete N0_bg_tmp;

            //-----------------------------------------------------------------------//
            outfile.close();
            //-----------------------------------------------------------------------//
        }
    }

    //==================================================================================================================================================//
    if(plotST == 2 || plotST == 3)
    {
        Double_t F12 = pointer_ch[0]->F12;
        Double_t errF12 = pointer_ch[0]->errF12;
		
        TString output_file_name = dir_name;
        output_file_name += "/results.root";
        TFile *resultsfile = new TFile(output_file_name.Data(),"recreate");

        //==================================================================================================================================================//
        // Calculation of the interaction frequency for different orientation
        for(Int_t j = 0; j < number_of_cutting_angels; j++)
        {
            Cutting_angle[j] = 1e6*pointer_ch[0]->theta_cut[j];

        	Fin_CH[j] = (double)N12_CH[j]/N0_CH[j];
        	Fin_AM[j] = (double)N12_AM[j]/N0_AM[j];
        	Fin_BG[j] = (double)N12_BG[j]/N0_BG[j];

        	errFin_CH[j] = (TMath::Sqrt((double)N12_CH[j]/TMath::Power(N0_CH[j],2) + (double)TMath::Power(N12_CH[j],2)*N0_CH[j]/TMath::Power(N0_CH[j],4)));
        	errFin_AM[j] = (TMath::Sqrt((double)N12_AM[j]/TMath::Power(N0_AM[j],2) + (double)TMath::Power(N12_AM[j],2)*N0_AM[j]/TMath::Power(N0_AM[j],4)));
        	errFin_BG[j] = (TMath::Sqrt((double)N12_BG[j]/TMath::Power(N0_BG[j],2) + (double)TMath::Power(N12_BG[j],2)*N0_BG[j]/TMath::Power(N0_BG[j],4)));            
    	}
			
    
		
        TCanvas *c1 = new TCanvas("c1","Canva",200,10,700,500);
        c1->cd();
        c1->SetGrid();
        TMultiGraph * mg = new TMultiGraph();

        TGraph *gr_CH = new TGraphErrors(data_ana::number_of_cutting_angels,Cutting_angle,Fin_CH,errCutting_angle,errFin_CH);        
        TGraph *gr_AM = new TGraphErrors(data_ana::number_of_cutting_angels,Cutting_angle,Fin_AM,errCutting_angle,errFin_AM);
        TGraph *gr_BG = new TGraphErrors(data_ana::number_of_cutting_angels,Cutting_angle,Fin_BG,errCutting_angle,errFin_BG);

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

        //==================================================================================================================================================//
        // Calculation of the interaction probability with high statistical background
        Double_t IntProb_CH_BG[data_ana::number_of_cutting_angels] = {};
        Double_t IntProb_AM_BG[data_ana::number_of_cutting_angels] = {};
        Double_t errIntProb_CH_BG[data_ana::number_of_cutting_angels] = {};
        Double_t errIntProb_AM_BG[data_ana::number_of_cutting_angels] = {};

        for(Int_t i = 0; i < data_ana::number_of_cutting_angels; i++)
        {
            IntProb_CH_BG[i] = 100.0*(Fin_CH[i] - Fin_BG[i])/F12;
            IntProb_AM_BG[i] = 100.0*(Fin_AM[i] - Fin_BG[i])/F12;

            errIntProb_CH_BG[i] = 100.0*TMath::Sqrt(TMath::Power(errFin_CH[i]/F12,2) + TMath::Power(errFin_BG[i]/F12,2) +
                                           TMath::Power(errF12*(Fin_CH[i] - Fin_BG[i])/(F12*F12),2));
            errIntProb_AM_BG[i] = 100.0*TMath::Sqrt(TMath::Power(errFin_AM[i]/F12,2) + TMath::Power(errFin_BG[i]/F12,2) +
                                           TMath::Power(errF12*(Fin_AM[i] - Fin_BG[i])/(F12*F12),2));

            outfile_pl_frq<<Cutting_angle[i]<<" , "<<Fin_CH[i]<<" , "<<errFin_CH[i]<<"\n";
            outfile_am_frq<<Cutting_angle[i]<<" , "<<Fin_AM[i]<<" , "<<errFin_AM[i]<<"\n";
            outfile_bg_frq<<Cutting_angle[i]<<" , "<<Fin_BG[i]<<" , "<<errFin_BG[i]<<"\n";
            outfile_pl_prb<<Cutting_angle[i]<<" , "<<IntProb_CH_BG[i]<<" , "<<errIntProb_CH_BG[i]<<"\n";
            outfile_am_prb<<Cutting_angle[i]<<" , "<<IntProb_AM_BG[i]<<" , "<<errIntProb_AM_BG[i]<<"\n";

            //cout<<Cutting_angle[i]<<"       "<<IntProb_CH_BG[i]<<"        "<<errIntProb_CH_BG[i]<<"     "<<IntProb_AM_BG[i]<<"        "<<errIntProb_AM_BG[i]<<endl;
        }

        TCanvas *c2 = new TCanvas("c2","Canva",200,10,700,500);
        c2->cd();
        c2->SetGrid();

        TMultiGraph * mg_prob_BG = new TMultiGraph();

        TGraph *gr_CH_prob_BG = new TGraphErrors(data_ana::number_of_cutting_angels,Cutting_angle,IntProb_CH_BG,errCutting_angle,errIntProb_CH_BG);
        TGraph *gr_AM_prob_BG = new TGraphErrors(data_ana::number_of_cutting_angels,Cutting_angle,IntProb_AM_BG,errCutting_angle,errIntProb_AM_BG);

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
        mg_prob_BG->SetMinimum(0.0);
        c2->Update();
        c2->Write();

        // Calculation of the interaction probability ratio PCH/AM
        cout<<"--> Calculation of the interaction probability ratio PCH/AM"<<endl;
        Double_t Ratio_CH_AM[data_ana::number_of_cutting_angels] = {};
        Double_t errRatio_CH_AM[data_ana::number_of_cutting_angels] = {};

        for(Int_t i = 0; i < data_ana::number_of_cutting_angels; i++)
        {
            Ratio_CH_AM[i] = IntProb_CH_BG[i]/IntProb_AM_BG[i];

            errRatio_CH_AM[i] = TMath::Sqrt(TMath::Power(errIntProb_CH_BG[i]/IntProb_AM_BG[i],2) + TMath::Power(IntProb_CH_BG[i]*errIntProb_AM_BG[i]/(IntProb_AM_BG[i]*IntProb_AM_BG[i]),2));

            outfile_ch_am_rt<<Cutting_angle[i]<<" , "<<Ratio_CH_AM[i]<<" , "<<errRatio_CH_AM[i]<<"\n";

            cout<<Cutting_angle[i]<<"  "<<Ratio_CH_AM[i]<<"  "<<errRatio_CH_AM[i]<<endl;
        }
        TGraph *gr_ratio_CH_AM = new TGraphErrors(data_ana::number_of_cutting_angels,Cutting_angle,Ratio_CH_AM,0,errRatio_CH_AM);
        gr_ratio_CH_AM->SetName("gr_ratio_CH_AM");
        gr_ratio_CH_AM->SetTitle("PCH/AM");
        gr_ratio_CH_AM->SetMarkerStyle(20);
        gr_ratio_CH_AM->SetMarkerColor(kBlack);
        gr_ratio_CH_AM->SetLineColor(kBlue);
        gr_ratio_CH_AM->SetLineWidth(2);
        gr_ratio_CH_AM->SetFillStyle(0);
        gr_ratio_CH_AM->Write();

        resultsfile->Close();
    }

    outfile_pl_frq.close();
    outfile_am_frq.close();
    outfile_bg_frq.close();
    outfile_pl_prb.close();
    outfile_am_prb.close();
    outfile_ch_am_rt.close();

    stop_time = time(NULL);
    cout<<"--> Running time is : "<<stop_time - start_time<<" seconds"<<endl;
    return 0;
}
