#include "data_ana.hh"

using namespace std;

data_ana::data_ana(TString data_file_name)
{
    file = new TFile(data_file_name,"READ");
    cout<<"--> Input filename: "<<file->GetName()<<endl;
    TTree * tree = (TTree*)file->Get("simpleEvent");

    nentries = (Int_t)tree->GetEntries();

    cout<<"--> nEntries: "<<nentries<<endl;

    bEvent          = tree->GetBranch("Event");
    bTime           = tree->GetBranch("Time");
    bDate           = tree->GetBranch("Date");
    bGonioPos       = tree->GetBranch("GonioPos");
    bMultiHit       = tree->GetBranch("MultiHit");
    bMultiHits      = tree->GetBranch("MultiHits");
    bSingleTrack    = tree->GetBranch("SingleTrack");
    bTracks         = tree->GetBranch("Tracks");

    bEvent->SetAddress(&evt);
    bTime->SetAddress(&time);
    bDate->SetAddress(&date);
    bGonioPos->SetAddress(&gPos);
    bMultiHit->SetAddress(&isHit);
    bMultiHits->SetAddress(&hits);
    bSingleTrack->SetAddress(&isTrack);
    bTracks->SetAddress(&tracks);

    // Multi track
    // IN-coming
    TLeaf * In_x_m = bMultiHits->GetLeaf("thetaIn_x");
    In_x_m->SetAddress(&INx_multi);

    TLeaf * In_y_m = bMultiHits->GetLeaf("thetaIn_y");
    In_y_m->SetAddress(&INy_multi);

    // Impact parameter
    TLeaf * Imp_x_m = bMultiHits->GetLeaf("d0_x");
    Imp_x_m->SetAddress(&IMPx_multi);

    TLeaf * Imp_y_m = bMultiHits->GetLeaf("d0_y");
    Imp_y_m->SetAddress(&IMPy_multi);

    // Single track
    // IN-coming
    TLeaf * In_x = bTracks->GetLeaf("thetaIn_x");
    In_x->SetAddress(&INx);

    TLeaf * In_y = bTracks->GetLeaf("thetaIn_y");
    In_y->SetAddress(&INy);

    // OUT-coming
    TLeaf * Out_x = bTracks->GetLeaf("thetaOut_x");
    Out_x->SetAddress(&OUTx);

    TLeaf * Out_y = bTracks->GetLeaf("thetaOut_y");
    Out_y->SetAddress(&OUTy);

    // Impact parameter
    TLeaf * Imp_x = bTracks->GetLeaf("d0_x");
    Imp_x->SetAddress(&IMPx);

    TLeaf * Imp_y = bTracks->GetLeaf("d0_y");
    Imp_y->SetAddress(&IMPy);

    tree->SetBranchAddress("Hits4x",&Hits4x);
    tree->SetBranchAddress("Hits4y",&Hits4y);

    for(Int_t k = 0; k < number_of_cutting_angels; k++)
    {
        if(number_of_cutting_angels > 1)
            theta_cut[k] = theta_cut_min + k*(theta_cut_max - theta_cut_min)/(number_of_cutting_angels - 1);
        else
            theta_cut[k] = theta_cut_min;
    }

    //Default values
    F12     = 1.00;
    errF12  = 0.01;

    min_d0x_sept = -0.55;   // [mm]
    max_d0x_sept = -0.45;   // [mm]
    distL        =  4885.0; // [mm]

    GetParameters(data_file_name.Data(), min_d0x, max_d0x, min_d0y, max_d0y, initial_gpos_x, k1, k2, a, b);
    angular_cut = 5.0e-6;   //[rad]
}

data_ana::~data_ana()
{
    file->Delete();
}

void data_ana::plot_histogramms()
{
    cout<<endl;
    cout<<"--> Plotting histogramms ..."<<endl;
    cout<<"--> F12 = "<<F12<<endl;
    cout<<"--> errF12 = "<<errF12<<endl;
    cout<<"--> max_d0x = "<<max_d0x<<endl;
    cout<<"--> min_d0x = "<<min_d0x<<endl;
    cout<<"--> max_d0y = "<<max_d0y<<endl;
    cout<<"--> min_d0y = "<<min_d0y<<endl;
    cout<<"--> initial_gpos_x = "<<initial_gpos_x<<endl;
    cout<<"--> k1 = "<<k1<<endl;
    cout<<"--> k2 = "<<k2<<endl;
    cout<<"--> a = "<<a<<endl;
    cout<<"--> b = "<<b<<endl;
    cout<<"--> angular_cut = "<<angular_cut<<endl;
    cout<<endl;

    // ThetaX
    hist1 = new TH1D("hist1","#Theta^{in}_{X}",10000,-2000.0,2000.0);
    hist2 = new TH1D("hist2","#Theta^{out}_{X}",10000,-2000.0,2000.0);
    // ThetaY
    hist3 = new TH1D("hist3","#Theta^{in}_{Y}",10000,-2000.0,2000.0);
    hist4 = new TH1D("hist4","#Theta^{out}_{Y}",10000,-2000.0,2000.0);
    // dThetaX
    hist5 = new TH1D("hist5","#Delta#Theta_{X}",10000,-2000.0,2000.0);
    // dThetaX
    hist6 = new TH1D("hist6","#Delta#Theta_{Y}",10000,-2000.0,2000.0);
    // dThetaY vs dThetaX
    hist7 = new TH2D("hist7","#Delta#Theta_{Y} vs #Delta#Theta_{X}",2500,-500,500,2500,-500.0,500.0);
    // dThetaX vs ThetaX
    hist8 = new TH2D("hist8","#Delta#Theta_{X} vs #Theta^{in}_{X}",2500,-500.0,500.0,2500,-5000,5000);
    // dThetaY vs ThetaY
    hist9 = new TH2D("hist9","#Delta#Theta_{Y} vs #Theta^{in}_{Y}",2500,-500.0,500.0,2500,-500.0,500.0);
    // dThetaX vs d0X
    hist10 = new TH2D("hist10","#Delta#Theta_{X} vs d0_{X}",4000,-20.0,20.0,2500,-500,500);
    //hist10 = new TH2D("hist10","#Delta#Theta_{X} vs d0_{X}",200,-1.0,1.0,200,-100,100);
    // dThetaX vs d0Y
    hist11 = new TH2D("hist11","#Delta#Theta_{X} vs d0_{Y}",4000,-20.0,20.0,2500,-500,500);
    profile11 = new TProfile("profile11","#Delta#Theta_{X} vs d0_{Y}",2500,-20.0,20.0,-500,500);
    profile11_corr = new TProfile("profile11_corr","#Delta#Theta_{X} vs d0_{Y}",2500,-20.0,20.0,-500,500);
    // dThetaY vs d0X
    hist12 = new TH2D("hist12","#Delta#Theta_{Y} vs d0_{X}",2500,-20.0,20.0,2500,-500.0,500.0);
    // dThetaY vs d0Y
    hist13 = new TH2D("hist13","#Delta#Theta_{Y} vs d0_{Y}",2500,-20.0,20.0,2500,-500.0,500.0);
    // d0X vs d0Y
    hist14 = new TH2D("hist14","d0_{Y} vs d0_{X}",2500,-20.0,20.0,2500,-20.0,20.0);
    // d0X
    hist15 = new TH1D("hist15","d0_{X}",2500,-20.0,20.0);
    // d0Y
    hist16 = new TH1D("hist16","d0_{Y}",2500,-20.0,20.0);
    // dTime
    hist17 = new TH1D("hist17","#DeltaTime",10000,-0.0,0.01);
    // ThetaX vs d0Y
    hist18 = new TH2D("hist18","#Theta^{in}_{X} vs d0_{Y}",2500,-20.0,20.0,2500,-500.0,500.0);
    // ThetaY vs ThetaX
    hist19 = new TH2D("hist19","#Theta^{in}_{Y} vs #Theta^{in}_{X}",2500,-500.0,500.0,2500,-500.0,500.0);

    // dThetaX vs ThetaX with geometrical cut
    hist20 = new TH2D("hist20","#Delta#Theta_{X} vs #Theta^{in}_{X}",4000,-200.0,200.0,2500,-500,500);
    // dThetaX vs d0X with geometrical cut
    hist21 = new TH2D("hist21","#Delta#Theta_{X} vs d0_{X}",2500,-20.0,20.0,2500,-500,500);
    // dThetaX with geometrical cut
    hist22 = new TH1D("hist22","#Delta#Theta_{X}",2500,-500,500);
    // gPos.x
    hist23 = new TH1D("hist23","gPos.x",100000,0.0,5.0e6);
    // dThetaX vs d0Y with geometrical cut
    hist24 = new TH2D("hist24","#Delta#Theta_{X} vs d0_{Y}",2500,-20.0,20.0,2500,-500,500);

    // Hits4 yx (without cuts)
    hist25 = new TH2D("hist25","Hits4 YX",467,-14.01,14.01,467,-14.01,14.01);
    hist27 = new TH1D("hist27","Hits4 X",467,-14.01,14.01);
    hist28 = new TH1D("hist28","Hits4 Y",467,-14.01,14.01);

    // Hits4 yx (with cuts)
    hist26 = new TH2D("hist26","Hits4 YX (with cuts)",467,-14.01,14.01,467,-14.01,14.01);
    hist29 = new TH1D("hist29","Hits4 X (with cuts)",467,-14.01,14.01);
    hist30 = new TH1D("hist30","Hits4 Y (with cuts)",467,-14.01,14.01);

    // d0X vs d0Y cut
    hist31 = new TH2D("hist31","d0_{Y} vs d0_{X} (cut)",2500,-20.0,20.0,2500,-20.0,20.0);


    // Axis labels
    hist1->GetXaxis()->SetTitle("#murad");
    hist2->GetXaxis()->SetTitle("#murad");
    hist3->GetXaxis()->SetTitle("#murad");
    hist4->GetXaxis()->SetTitle("#murad");
    hist5->GetXaxis()->SetTitle("#murad");
    hist6->GetXaxis()->SetTitle("#murad");
    hist7->GetXaxis()->SetTitle("#Delta#Theta_{X} , #murad");
    hist7->GetYaxis()->SetTitle("#Delta#Theta_{Y} , #murad");
    hist8->GetXaxis()->SetTitle("#Theta^{in}_{X} , #murad");
    hist8->GetYaxis()->SetTitle("#Delta#Theta_{X} , #murad");
    hist9->GetXaxis()->SetTitle("#Theta^{in}_{Y} , #murad");
    hist9->GetYaxis()->SetTitle("#Delta#Theta_{Y} , #murad");
    hist10->GetXaxis()->SetTitle("d0_{X} , mm");
    hist10->GetYaxis()->SetTitle("#Delta#Theta_{X} , #murad");
    hist11->GetXaxis()->SetTitle("d0_{Y} , mm");
    hist11->GetYaxis()->SetTitle("#Delta#Theta_{X} , #murad");
    hist12->GetXaxis()->SetTitle("d0_{X} , mm");
    hist12->GetYaxis()->SetTitle("#Delta#Theta_{Y} , #murad");
    hist13->GetXaxis()->SetTitle("d0_{Y} , mm");
    hist13->GetYaxis()->SetTitle("#Delta#Theta_{Y} , #murad");
    hist14->GetXaxis()->SetTitle("d0_{X} , mm");
    hist14->GetYaxis()->SetTitle("d0_{Y} , mm");
    hist15->GetXaxis()->SetTitle("mm");
    hist16->GetXaxis()->SetTitle("mm");
    hist17->GetXaxis()->SetTitle("sec");
    hist18->GetXaxis()->SetTitle("d0_{Y} , mm");
    hist18->GetYaxis()->SetTitle("#Theta_{X} , #murad");
    hist19->GetXaxis()->SetTitle("#Theta_{X} , #murad");
    hist19->GetYaxis()->SetTitle("#Theta_{Y} , #murad");
    hist20->GetXaxis()->SetTitle("#Theta^{in}_{X} , #murad");
    hist20->GetYaxis()->SetTitle("#Delta#Theta_{X} , #murad");
    hist21->GetXaxis()->SetTitle("d0_{X} , mm");
    hist21->GetYaxis()->SetTitle("#Delta#Theta_{X} , #murad");
    hist22->GetXaxis()->SetTitle("#murad");
    hist23->GetXaxis()->SetTitle("#murad");
    hist31->GetXaxis()->SetTitle("d0_{X} , mm");
    hist31->GetYaxis()->SetTitle("d0_{Y} , mm");

    TTree * tree = (TTree*)file->Get("simpleEvent");

    Double_t newTime = 0, oldTime = 0;
    Int_t count = 0;

    tree->GetEvent(0);
    printf("--> gPos.x: %14.3f [urad]\n",gPos.x);
    for (Int_t i = 0; i < nentries; i++)
    {
        tree->GetEvent(i);

        hist23->Fill(gPos.x);

        if(isTrack == 1  && !isHit)
        {
            hist1->Fill(INx*1e6);
            hist2->Fill(OUTx*1e6);
            hist3->Fill(INy*1e6);
            hist4->Fill(OUTy*1e6);
            hist6->Fill((OUTy-INy)*1e6);
            hist9->Fill(INy*1e6,(OUTy-INy)*1e6);
            hist11->Fill(IMPy,(OUTx-INx)*1e6);
            profile11->Fill(IMPy,(OUTx-INx)*1e6,1);
            hist7->Fill((OUTx-INx)*1e6,(OUTy-INy)*1e6);
            hist12->Fill(IMPx,(OUTy-INy)*1e6);
            hist13->Fill(IMPy,(OUTy-INy)*1e6);
            hist14->Fill(IMPx,IMPy);
            hist15->Fill(IMPx);
            hist16->Fill(IMPy);
            hist18->Fill(IMPy,INx*1e6);
            hist19->Fill(INx*1e6,INy*1e6);
            hist5->Fill((OUTx-INx)*1e6);
            hist7->Fill((OUTx-INx)*1e6,(OUTy-INy)*1e6);
            hist8->Fill(INx*1e6-(gPos.x-initial_gpos_x),(OUTx-INx)*1e6);
            hist10->Fill(IMPx,(OUTx-INx)*1e6);

            if (IMPy < max_d0y && IMPy > min_d0y)
            {
                if (IMPx < max_d0x && IMPx > min_d0x)
                {
                    hist20->Fill(INx*1e6-(gPos.x-initial_gpos_x),(OUTx-INx)*1e6);
                    hist21->Fill(IMPx,(OUTx-INx)*1e6);
                    hist22->Fill((OUTx-INx)*1e6);
                }
            }
            if(IMPx < max_d0x && IMPx > min_d0x)
            {
                hist24->Fill(IMPy,(OUTx-INx)*1e6);
            }

            int Hits4xSize = Hits4x->size();
            int Hits4ySize = Hits4y->size();

            if(Hits4xSize > 0 && Hits4ySize > 0)
            {
                for(Int_t xi = 0; xi < Hits4xSize; xi++)
                {
                    for(Int_t yi = 0; yi < Hits4ySize; yi++)
                    {
                        hist25->Fill(Hits4x->at(xi),Hits4y->at(yi),1);
                        hist27->Fill(Hits4x->at(xi),1);
                        hist28->Fill(Hits4y->at(yi),1);
                        if (IMPy < max_d0y && IMPy > min_d0y)
                        {
                            if (IMPx < max_d0x && IMPx > min_d0x)
                            {
                                if (TMath::Abs(INx) < angular_cut)
                                {
                                    hist26->Fill(Hits4x->at(xi),Hits4y->at(yi),1);
                                    hist29->Fill(Hits4x->at(xi),1);
                                    hist30->Fill(Hits4y->at(yi),1);
                                }
                            }
                        }
                    }
                }
            }

            cut_d0x = 0;
            cut_d0y = 0;

            if (IMPy < max_d0y && IMPy > min_d0y && IMPx < max_d0x && IMPx > min_d0x)// Diffuser/Crystal
            {
                cut_d0x = 1;
                cut_d0y = 1;
            }
            else if (IMPy + distL*INy < max_d0y && IMPy + distL*INy > min_d0y &&
                     IMPx + distL*INx < max_d0x_sept && IMPx + distL*INx > min_d0x_sept)// Septum
            {
                cut_d0y = 1;
                cut_d0x = 1;
            }

            if(cut_d0x && cut_d0y) hist31->Fill(IMPx,IMPy);
        }
        else if(isHit == 1)
        {
            hist1->Fill(INx_multi*1e6);
            hist3->Fill(INy_multi*1e6);
            hist14->Fill(IMPx_multi,IMPy_multi);
            hist15->Fill(IMPx_multi);
            hist16->Fill(IMPy_multi);
            hist18->Fill(IMPy_multi,INx_multi*1e6);
            hist19->Fill(INx_multi*1e6,INy_multi*1e6);

            int Hits4xSize = Hits4x->size();
            int Hits4ySize = Hits4y->size();

            if(Hits4xSize > 0 && Hits4ySize > 0)
            {
                for(Int_t xi = 0; xi < Hits4xSize; xi++)
                {
                    for(Int_t yi = 0; yi < Hits4ySize; yi++)
                    {
                        hist25->Fill(Hits4x->at(xi),Hits4y->at(yi),1);
                        if (IMPy_multi < max_d0y && IMPy_multi > min_d0y)
                        {
                            if (IMPx_multi < max_d0x && IMPx_multi > min_d0x)
                            {
                                if (TMath::Abs(INx_multi) < angular_cut)
                                {
                                    hist26->Fill(Hits4x->at(xi),Hits4y->at(yi),1);
                                }
                            }
                        }
                    }
                }
            }

            cut_d0x = 0;
            cut_d0y = 0;

            if (IMPy_multi < max_d0y && IMPy_multi > min_d0y && IMPx_multi < max_d0x && IMPx_multi > min_d0x)// Diffuser/Crystal
            {
                cut_d0y = 1;
                cut_d0x = 1;
            }
            else if (IMPy_multi + distL*INy_multi < max_d0y && IMPy_multi + distL*INy_multi > min_d0y &&
                     IMPx_multi + distL*INx_multi < max_d0x_sept && IMPx_multi + distL*INx_multi > min_d0x_sept)// Septum
            {
                cut_d0y = 1;
                cut_d0x = 1;
            }

            if(cut_d0x && cut_d0y) hist31->Fill(IMPx_multi,IMPy_multi);
        }
        if(i%10000 == 0)
        {
            printf("\r--> Working: %3.1f %% | gPos.x: %14.3f [urad]", 100*(Double_t)i/nentries,gPos.x - initial_gpos_x);
            fflush(stdout);
        }

        newTime = atof(time);count++;
        if(count % 100 ==  0)
        {
            hist17->Fill((newTime-oldTime)/100);
            oldTime = newTime;
        }
    }
    printf("\nDone!\n");
}

void data_ana::write_histo_to_file(TString outFileName)
{
    cout<<"--> Writing output data to file <<"<<outFileName<<">>"<<endl;

    TFile *outfile = new TFile(outFileName,"recreate");

    hist1->Write();
    hist2->Write();
    hist3->Write();
    hist4->Write();
    hist5->Write();
    hist6->Write();
    hist7->Write();
    hist8->Write();
    hist9->Write();
    hist10->Write();
    hist11->Write();
    profile11->Write();
    profile11_corr->Write();
    hist12->Write();
    hist13->Write();
    hist14->Write();
    hist15->Write();
    hist16->Write();
    hist17->Write();
    hist18->Write();
    hist19->Write();
    hist20->Write();
    hist21->Write();
    hist22->Write();
    hist23->Write();
    hist24->Write();
    hist25->Write();
    hist26->Write();
    hist27->Write();
    hist28->Write();
    hist29->Write();
    hist30->Write();
    hist31->Write();

    outfile->Close();
}

void data_ana::find_inelastic_frequency1(Double_t *N12, Double_t *N0, Int_t signal1, Int_t signal2)
{
    cout<<endl;
    cout<<"--> Calculating of the inelastic frequency1 ..."<<endl;
    cout<<"--> F12 = "<<F12<<endl;
    cout<<"--> errF12 = "<<errF12<<endl;
    cout<<"--> max_d0x = "<<max_d0x<<endl;
    cout<<"--> min_d0x = "<<min_d0x<<endl;
    cout<<"--> max_d0y = "<<max_d0y<<endl;
    cout<<"--> min_d0y = "<<min_d0y<<endl;
    cout<<"--> initial_gpos_x = "<<initial_gpos_x<<endl;
    cout<<"--> k1 = "<<k1<<endl;
    cout<<"--> k2 = "<<k2<<endl;
    cout<<"--> a = "<<a<<endl;
    cout<<"--> b = "<<b<<endl;
    cout<<"--> angular_cut = "<<angular_cut<<endl;
    cout<<endl;

    TTree * tree = (TTree*)file->Get("simpleEvent");

    for(Int_t j = 0; j < number_of_cutting_angels; j++)
    {
        N0[j] = 0; N12[j] = 0;
    }

    tree->GetEvent(0);
    printf("--> gPos.x: %14.3f [urad]\n",gPos.x);

    for (Int_t i = 0; i < nentries; i++)
    {
        cut_d0x = 0;
        cut_d0y = 0;

        Double_t INx_corr = 0.0;

        tree->GetEvent(i);
        // Geometrical cuts
        if(isTrack == 1 && isHit == 0)
        {
            //INx_corr = INx - (gPos.x - initial_gpos_x)*1e-6 - (k2*IMPy + k1);
            INx_corr = INx - (k2*IMPy + k1);
            if (IMPy < max_d0y && IMPy > min_d0y) cut_d0y = 1;
            if (IMPx < max_d0x && IMPx > min_d0x) cut_d0x = 1;
        }
        else if(isHit == 1)
        {
            //INx_corr = INx_multi - (gPos.x - initial_gpos_x)*1e-6 - (k2*IMPy_multi + k1);
            INx_corr = INx_multi - (k2*IMPy_multi + k1);
            if (IMPy_multi < max_d0y && IMPy_multi > min_d0y) cut_d0y = 1;
            if (IMPx_multi < max_d0x && IMPx_multi > min_d0x) cut_d0x = 1;
        }
        if(isTrack == 1 || isHit == 1)
        {
            if(cut_d0x == 1 && cut_d0y == 1)
            {
                for(Int_t j = 0; j < number_of_cutting_angels; j++)
                {
                    if (TMath::Abs(INx_corr) < theta_cut[j])
                    {
                        N0[j]++;
                        Int_t incrm_0 = ((evt.nuclearRaw >> det2bit(signal1)) & 1);
                        Int_t incrm_1 = ((evt.nuclearRaw >> det2bit(signal2)) & 1);
                        if(incrm_0 == 1 && incrm_1 == 1)
                        {
                            N12[j]++;
                        }
                    }
                }
            }
        }
        if(i%10000 == 0)
        {
            printf("\r--> Working: %3.1f %% | gPos.x: %14.3f [urad]", 100*(Double_t)i/nentries,gPos.x - initial_gpos_x);
            fflush(stdout);
        }
    }
    printf("\nDone!\n");
}

void data_ana::find_inelastic_frequency2(Double_t *N12, Double_t *N0, Int_t signal1, Int_t signal2)
{
    cout<<endl;
    cout<<"--> Calculating of the inelastic frequency2 ..."<<endl;
    cout<<"--> F12 = "<<F12<<endl;
    cout<<"--> errF12 = "<<errF12<<endl;
    cout<<"--> max_d0x = "<<max_d0x<<endl;
    cout<<"--> min_d0x = "<<min_d0x<<endl;
    cout<<"--> max_d0y = "<<max_d0y<<endl;
    cout<<"--> min_d0y = "<<min_d0y<<endl;
    cout<<"--> initial_gpos_x = "<<initial_gpos_x<<endl;
    cout<<"--> k1 = "<<k1<<endl;
    cout<<"--> k2 = "<<k2<<endl;
    cout<<"--> a = "<<a<<endl;
    cout<<"--> b = "<<b<<endl;
    cout<<"--> angular_cut = "<<angular_cut<<endl;
    cout<<endl;

    TTree * tree = (TTree*)file->Get("simpleEvent");

    for(Int_t j = 0; j < number_of_cutting_angels; j++)
    {
        N0[j] = 0; N12[j] = 0;
    }

    tree->GetEvent(0);
    printf("--> gPos.x: %14.3f [urad]\n",gPos.x);

    for (Int_t i = 0; i < nentries; i++)
    {
        cut_d0x = 0;
        cut_d0y = 0;

        tree->GetEvent(i);
        // Geometrical cuts
        if(isTrack == 1 && isHit == 0)
        {
            if (IMPy < max_d0y && IMPy > min_d0y) cut_d0y = 1;
            if (IMPx < max_d0x && IMPx > min_d0x) cut_d0x = 1;
        }
        else if(isHit == 1)
        {
            if (IMPy_multi < max_d0y && IMPy_multi > min_d0y) cut_d0y = 1;
            if (IMPx_multi < max_d0x && IMPx_multi > min_d0x) cut_d0x = 1;
        }
        if(isTrack == 1 || isHit == 1)
        {
            if(cut_d0x == 1 && cut_d0y == 1)
            {
                for(Int_t j = 0; j < number_of_cutting_angels; j++)
                {
                    N0[j]++;
                    Int_t incrm_0 = ((evt.nuclearRaw >> det2bit(signal1)) & 1);
                    Int_t incrm_1 = ((evt.nuclearRaw >> det2bit(signal2)) & 1);
                    if(incrm_0 == 1 && incrm_1 == 1)
                    {
                        N12[j]++;
                    }
                }
            }
        }
        if(i%10000 == 0)
        {
            printf("\r--> Working: %3.1f %% | gPos.x: %14.3f [urad]", 100*(Double_t)i/nentries,gPos.x - initial_gpos_x);
            fflush(stdout);
        }
    }
    printf("\nDone!\n");
}

void data_ana::find_inelastic_frequency3(Double_t *N12, Double_t *N0, Int_t signal1, Int_t signal2)
{
    cout<<endl;
    cout<<"--> Calculating of the inelastic frequency3 ..."<<endl;
    cout<<"--> F12 = "<<F12<<endl;
    cout<<"--> errF12 = "<<errF12<<endl;
    cout<<"--> max_d0x = "<<max_d0x<<endl;
    cout<<"--> min_d0x = "<<min_d0x<<endl;
    cout<<"--> max_d0y = "<<max_d0y<<endl;
    cout<<"--> min_d0y = "<<min_d0y<<endl;
    cout<<"--> initial_gpos_x = "<<initial_gpos_x<<endl;
    cout<<"--> k1 = "<<k1<<endl;
    cout<<"--> k2 = "<<k2<<endl;
    cout<<"--> a = "<<a<<endl;
    cout<<"--> b = "<<b<<endl;
    cout<<"--> angular_cut = "<<angular_cut<<endl;
    cout<<endl;

    TTree * tree = (TTree*)file->Get("simpleEvent");

    for(Int_t j = 0; j < number_of_cutting_angels; j++)
    {
        N0[j] = 0; N12[j] = 0;
    }

    tree->GetEvent(0);
    printf("--> gPos.x: %14.3f [urad]\n",gPos.x);

    for (Int_t i = 0; i < nentries; i++)
    {
        cut_d0x = 0;
        cut_d0y = 0;

        Double_t INx_corr = 0.0;

        tree->GetEvent(i);
        // Geometrical cuts
        if(isTrack == 1 && isHit == 0)
        {
            INx_corr = INx - (gPos.x - initial_gpos_x)*1e-6 - (k2*TMath::Power(IMPy-a,2) + k1*IMPy + b);
            if (IMPy < max_d0y && IMPy > min_d0y) cut_d0y = 1;
            if (IMPx < max_d0x && IMPx > min_d0x) cut_d0x = 1;
        }
        else if(isHit == 1)
        {
            INx_corr = INx_multi - (gPos.x - initial_gpos_x)*1e-6 - (k2*TMath::Power(IMPy_multi-a,2) + k1*IMPy_multi + b);
            if (IMPy_multi < max_d0y && IMPy_multi > min_d0y) cut_d0y = 1;
            if (IMPx_multi < max_d0x && IMPx_multi > min_d0x) cut_d0x = 1;
        }
        if(isTrack == 1 || isHit == 1)
        {
            if(cut_d0x == 1 && cut_d0y == 1)
            {
                for(Int_t j = 0; j < number_of_cutting_angels; j++)
                {
                    if (TMath::Abs(INx_corr) < theta_cut[j])
                    {
                        N0[j]++;
                        Int_t incrm_0 = ((evt.nuclearRaw >> det2bit(signal1)) & 1);
                        Int_t incrm_1 = ((evt.nuclearRaw >> det2bit(signal2)) & 1);
                        if(incrm_0 == 1 && incrm_1 == 1)
                        {
                            N12[j]++;
                        }
                    }
                }
            }
        }
        if(i%10000 == 0)
        {
            printf("\r--> Working: %3.1f %% | gPos.x: %14.3f [urad]", 100*(Double_t)i/nentries,gPos.x - initial_gpos_x);
            fflush(stdout);
        }
    }
    printf("\nDone!\n");
}

void data_ana::GetParameters(string filenamepath, double &_min_d0x, double &_max_d0x, double &_min_d0y, double &_max_d0y, double &_gpos, double &_k1, double &_k2, double &_a, double &_b)
{
    size_t pos = filenamepath.find_last_of('/');
    string filename = filenamepath.substr(pos+1);

    string parameters_filename = "parameters.dat";

    ifstream inputfile(parameters_filename.data());
    TString word;
    bool status_read = false;

    if(inputfile.is_open())
    {
        while(1)
        {
            inputfile>>word;
            if(word == filename)
            {
                status_read = true;
                inputfile>>_min_d0x;
                inputfile>>_max_d0x;
                inputfile>>_min_d0y;
                inputfile>>_max_d0y;
                inputfile>>_gpos;
                inputfile>>_k1;
                inputfile>>_k2;
                inputfile>>_a;
                inputfile>>_b;
                break;
            }

            if(inputfile.eof()) {break;}
        }
        inputfile.close();
    }
    else
    {
        cout<<"--> ERROR: Unable to open the input file: "<<parameters_filename<<endl;
        assert(0);
    }
    if(!status_read) cout<<endl<<endl<<"--> Did not find file '"<<filename<<"' in the filelist '"<<parameters_filename<<"' !!!"<<endl<<endl;
}


void data_ana::get_counts(Double_t *N, Double_t &N_tot)
{
    cout<<endl;
    cout<<"--> Calculating of the scintillators counts ..."<<endl;
    cout<<"--> F12 = "<<F12<<endl;
    cout<<"--> errF12 = "<<errF12<<endl;
    cout<<"--> max_d0x = "<<max_d0x<<endl;
    cout<<"--> min_d0x = "<<min_d0x<<endl;
    cout<<"--> max_d0y = "<<max_d0y<<endl;
    cout<<"--> min_d0y = "<<min_d0y<<endl;
    cout<<"--> initial_gpos_x = "<<initial_gpos_x<<endl;
    cout<<"--> k1 = "<<k1<<endl;
    cout<<"--> k2 = "<<k2<<endl;
    cout<<"--> a = "<<a<<endl;
    cout<<"--> b = "<<b<<endl;
    cout<<"--> angular_cut = "<<angular_cut<<endl;
    cout<<endl;

    TTree * tree = (TTree*)file->Get("simpleEvent");
    for(Int_t j = 0; j < number_of_scintillators; j++)
    {
        N[j] = 0.0;
    }
    N_tot = 0.0;
    Double_t Nselected = 0;
    tree->GetEvent(0);
    printf("--> gPos.x: %14.3f [urad]\n",gPos.x);

    for (Int_t i = 0; i < nentries; i++)
    {
        cut_d0x = 0;
        cut_d0y = 0;

        tree->GetEvent(i);
        // Geometrical cuts
        if(isTrack == 1 && isHit == 0)
        {
            if (IMPy < max_d0y && IMPy > min_d0y && IMPx < max_d0x && IMPx > min_d0x)// Diffuser/Crystal
            {
                cut_d0x = 1;
                cut_d0y = 1;
            }
            else if (IMPy + distL*INy < max_d0y && IMPy + distL*INy > min_d0y &&
                     IMPx + distL*INx < max_d0x_sept && IMPx + distL*INx > min_d0x_sept)// Septum
            {
                cut_d0y = 1;
                cut_d0x = 1;
            }
        }
        else if(isHit == 1)
        {
            if (IMPy_multi < max_d0y && IMPy_multi > min_d0y && IMPx_multi < max_d0x && IMPx_multi > min_d0x)// Diffuser/Crystal
            {
                cut_d0y = 1;
                cut_d0x = 1;
            }
            else if (IMPy_multi + distL*INy_multi < max_d0y && IMPy_multi + distL*INy_multi > min_d0y &&
                     IMPx_multi + distL*INx_multi < max_d0x_sept && IMPx_multi + distL*INx_multi > min_d0x_sept)// Septum
            {
                cut_d0y = 1;
                cut_d0x = 1;
            }
        }
        if(isTrack == 1 || isHit == 1)
        {
            N_tot += 1;
            if(cut_d0x == 1 && cut_d0y == 1)
            {
                Nselected += 1;
                for (int bit_id = 0; bit_id < 8; bit_id++)
                {
                    Int_t incrm = ((evt.nuclearRaw >> bit_id) & 1);
                    if(incrm == 1)
                    {
                        N[bit2det(bit_id)] += 1;
                    }
                }
            }
        }

        if(i%10000 == 0)
        {
            printf("\r--> Working: %3.1f %% | gPos.x: %14.3f [urad]", 100*(Double_t)i/nentries,gPos.x - initial_gpos_x);
            fflush(stdout);
        }
    }
    printf("\nDone!\n");
    cout<<"--> N_SELECTED = "<<Nselected<<" +/- "<<TMath::Sqrt(Nselected)<<endl;
}

void data_ana::get_counts_pairs(Double_t *N, Double_t &N_tot)
{
    cout<<endl;
    cout<<"--> Calculating of the scintillators pairs counts ..."<<endl;
    cout<<"--> F12 = "<<F12<<endl;
    cout<<"--> errF12 = "<<errF12<<endl;
    cout<<"--> max_d0x = "<<max_d0x<<endl;
    cout<<"--> min_d0x = "<<min_d0x<<endl;
    cout<<"--> max_d0y = "<<max_d0y<<endl;
    cout<<"--> min_d0y = "<<min_d0y<<endl;
    cout<<"--> initial_gpos_x = "<<initial_gpos_x<<endl;
    cout<<"--> k1 = "<<k1<<endl;
    cout<<"--> k2 = "<<k2<<endl;
    cout<<"--> a = "<<a<<endl;
    cout<<"--> b = "<<b<<endl;
    cout<<"--> angular_cut = "<<angular_cut<<endl;
    cout<<endl;

    TTree * tree = (TTree*)file->Get("simpleEvent");
    for(Int_t j = 0; j < number_of_scintillators/2; j++)
    {
        N[j] = 0.0;
    }
    N_tot = 0.0;
    Double_t Nselected = 0;
    tree->GetEvent(0);
    printf("--> gPos.x: %14.3f [urad]\n",gPos.x);

    for (Int_t i = 0; i < nentries; i++)
    {
        cut_d0x = 0;
        cut_d0y = 0;

        tree->GetEvent(i);
        // Geometrical cuts
        if(isTrack == 1 && isHit == 0)
        {
            if (IMPy < max_d0y && IMPy > min_d0y && IMPx < max_d0x && IMPx > min_d0x)// Diffuser/Crystal
            {
                cut_d0x = 1;
                cut_d0y = 1;
            }
            else if (IMPy + distL*INy < max_d0y && IMPy + distL*INy > min_d0y &&
                     IMPx + distL*INx < max_d0x_sept && IMPx + distL*INx > min_d0x_sept)// Septum
            {
                cut_d0y = 1;
                cut_d0x = 1;
            }
        }
        else if(isHit == 1)
        {
            if (IMPy_multi < max_d0y && IMPy_multi > min_d0y && IMPx_multi < max_d0x && IMPx_multi > min_d0x)// Diffuser/Crystal
            {
                cut_d0y = 1;
                cut_d0x = 1;
            }
            else if (IMPy_multi + distL*INy_multi < max_d0y && IMPy_multi + distL*INy_multi > min_d0y &&
                     IMPx_multi + distL*INx_multi < max_d0x_sept && IMPx_multi + distL*INx_multi > min_d0x_sept)// Septum
            {
                cut_d0y = 1;
                cut_d0x = 1;
            }
        }
        if(isTrack == 1 || isHit == 1)
        {
            N_tot += 1;
            if(cut_d0x == 1 && cut_d0y == 1)
            {
                Nselected += 1;
                for(int det_id = 0; det_id < number_of_scintillators; det_id+=2)
                {
                    Int_t incrm_0 = ((evt.nuclearRaw >> det2bit(det_id+0)) & 1);
                    Int_t incrm_1 = ((evt.nuclearRaw >> det2bit(det_id+1)) & 1);
                    if(incrm_0 == 1 && incrm_1 == 1)
                    {
                        N[det_id/2] += 1;
                    }
                }
            }
        }
        if(i%10000 == 0)
        {
            printf("\r--> Working: %3.1f %% | gPos.x: %14.3f [urad]", 100*(Double_t)i/nentries,gPos.x - initial_gpos_x);
            fflush(stdout);
        }
    }
    printf("\nDone!\n");
    cout<<"--> N_SELECTED = "<<Nselected<<" +/- "<<TMath::Sqrt(Nselected)<<endl;
}

int data_ana::bit2det(int bit)
{
    int det = -1;
    switch (bit)
    {
    case 0:
        det = 5;
        break;
    case 1:
        det = 4;
        break;
    case 6:
        det = 3;
        break;
    case 7:
        det = 2;
        break;
    case 4:
        det = 1;
        break;
    case 5:
        det = 0;
        break;
    }
    return det;
}

int data_ana::det2bit(int det)
{
    int bit = -1;
    switch (det)
    {
    case 0:
        bit = 5;
        break;
    case 1:
        bit = 4;
        break;
    case 2:
        bit = 7;
        break;
    case 3:
        bit = 6;
        break;
    case 4:
        bit = 1;
        break;
    case 5:
        bit = 0;
        break;
    }
    return bit;
}

