#include "data_ana.hh"

using namespace std;

data_ana::data_ana(TString _fileName)
{    
    // STF110 MAY 2017
    max_d0x =  0.45;  // mm
    min_d0x =  0.05;  // mm
    max_d0y =   2.0;  // mm
    min_d0y =  -2.0;  // mm

    fChain = new TChain("simpleEvent");
    fChain->Add(_fileName.Data());

    nentries = fChain->GetEntries();

    bEvent          = fChain->GetBranch("Event");
    bTime           = fChain->GetBranch("Time");
    bDate           = fChain->GetBranch("Date");
    bGonioPos       = fChain->GetBranch("GonioPos");
    bMultiHit       = fChain->GetBranch("MultiHit");
    bMultiHits      = fChain->GetBranch("MultiHits");
    bSingleTrack    = fChain->GetBranch("SingleTrack");
    bTracks         = fChain->GetBranch("Tracks");

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

    for(Int_t k = 0; k < number_of_cutting_angel_center_position; k++)
    {
        if(number_of_cutting_angel_center_position > 1)
            theta_cut_center[k] = theta_cut_center_min + k*(theta_cut_center_max - theta_cut_center_min)/(number_of_cutting_angel_center_position - 1);
        else
            theta_cut_center[k] = theta_cut_center_min;

        //cout<<"--> Theta_cut_center["<<k<<"] = "<<theta_cut_center[k]*1e6<<" urad"<<endl;
    }
}

data_ana::~data_ana()
{}

void data_ana::analysis(Double_t *N12, Double_t *N0, Double_t p0, Double_t p1, Bool_t _orient)
{
    cout<<"--> nEntries: "<<nentries<<endl;

    for(Int_t j = 0; j < number_of_cutting_angel_center_position; j++)
    {
        N0[j] = 0; N12[j] = 0;
    }

    fChain->GetEvent(0);
    printf("--> gPos.x: %14.3f [urad]\n",gPos.x - initial_gpos_x);

    for (Long64_t i = 0; i < nentries; i++)
    {
        cut_d0x = 0;
        cut_d0y = 0;

        Double_t INx_corr = 0.0;

        fChain->GetEntry(i);
        // Geometrical cuts
        if(isTrack == 1 && isHit == 0)
        {
            if(_orient)
                INx_corr = INx - (gPos.x - initial_gpos_x)*1e-6 - (p1*IMPy + p0);
            else
                INx_corr = INx - (p1*IMPy + p0);

            if (IMPy < max_d0y && IMPy > min_d0y) cut_d0y = 1;
            if (IMPx < max_d0x && IMPx > min_d0x) cut_d0x = 1;
        }
        else if(isHit == 1)
        {
            if(_orient)
                INx_corr = INx_multi - (gPos.x - initial_gpos_x)*1e-6 - (p1*IMPy_multi + p0);
            else
                INx_corr = INx_multi - (p1*IMPy_multi + p0);

            if (IMPy_multi < max_d0y && IMPy_multi > min_d0y) cut_d0y = 1;
            if (IMPx_multi < max_d0x && IMPx_multi > min_d0x) cut_d0x = 1;
        }
        if(isTrack == 1 || isHit == 1)
        {
            if(cut_d0x == 1 && cut_d0y == 1)
            {
                for(Int_t j = 0; j < number_of_cutting_angel_center_position; j++)
                {
                    if (INx_corr <= (theta_cut_center[j] + theta_cut) && INx_corr >= (theta_cut_center[j] - theta_cut))
                    {
                        N0[j]++;
                        if (evt.nuclear == 3)
                        {
                            N12[j]++;
                        }
                    }
                }
            }
        }

        if(i%1000000 == 0)
        {
            printf("\r--> Working: %3.1f %%", 100*(Double_t)i/nentries);
            fflush(stdout);
        }
    }
    printf("\nDone!\n");
}

void data_ana::analysis2(Double_t &N12, Double_t &N0)
{
    cout<<"--> nEntries: "<<nentries<<endl;

    N0  = 0;
    N12 = 0;

    fChain->GetEvent(0);
    printf("--> gPos.x: %14.3f [urad]\n",gPos.x);

    for (Long64_t i = 0; i < nentries; i++)
    {
        cut_d0x = 0;
        cut_d0y = 0;

        fChain->GetEntry(i);
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
                N0++;
                if (evt.nuclear == 3)
                {
                    N12++;
                }
            }
        }

        if(i%1000000 == 0)
        {
            printf("\r--> Working: %3.1f %%", 100*(Double_t)i/nentries);
            fflush(stdout);
        }
    }
    printf("\nDone!\n");
}

void data_ana::fillhisto(TH2D* histo, Double_t p0, Double_t p1)
{
    cout<<"--> Filling the histogram ..."<<endl;
    for (Long64_t i = 0; i < nentries; i++)
    {
        cut_d0x = 0;
        cut_d0y = 0;

        Double_t INx_corr = 0.0;

        fChain->GetEntry(i);
        if(isTrack == 1 && isHit == 0)
        {
            INx_corr = INx - (gPos.x - initial_gpos_x)*1e-6 - (p1*IMPy + p0);

            if (IMPy < max_d0y && IMPy > min_d0y) cut_d0y = 1;
            if (IMPx < max_d0x && IMPx > min_d0x) cut_d0x = 1;

            if(cut_d0x == 1 && cut_d0y == 1)
            {
                histo->Fill(INx_corr*1e6,(OUTx-INx)*1e6);
            }
        }

        if(i%1000000 == 0)
        {
            printf("\r--> Working: %3.1f %%", 100*(Double_t)i/nentries);
            fflush(stdout);
        }
    }
    printf("\nDone!\n");
}
