#include "src/data_ana.hh"
#include "src/includes.hh"

using namespace std;

int main(int argc, char *argv[])
{
    cout<<endl<<endl;
    if(argc != 3)
    {
        cout<<endl;
        cout<<"--> ERROR:: Wrong number of input parameters!"<<endl;
        cout<<"--> [0]: script_name"<<endl;
        cout<<"--> [1]: intput file_name"<<endl;
        cout<<"--> [2]: run_mode (1 -- only historams, 2 -- only counts, 3 -- only frequency, 5 -- only 2 & 3, 6 -- only 1 & 2, 7 -- all together)"<<endl;
        cout<<endl;
        assert(0);
    }

    time_t start_time, stop_time;
    start_time = time(NULL);


    string filenamepath = argv[1];
    size_t pos = filenamepath.find_last_of('/');
    string filename = filenamepath.substr(pos+1);
    string output_dir = "./output/";
    TString output_filename;

    data_ana *pointer = new data_ana(filenamepath.data());
    Double_t N_COUNTS[pointer->number_of_scintillators] = {};
    Double_t N_COUNTS_PAIRS[pointer->number_of_scintillators/2] = {};
    Double_t N_COUNTS_TOTAL = 0.0;
    Double_t N_COUNTS_TOTAL_PAIRS = 0.0;

    Int_t run_mode = atoi(argv[2]);

    // HISTOGRAMS
    if(run_mode == 1 || run_mode == 6 || run_mode == 7)
    {
        pointer->plot_histogramms();
        output_filename = output_dir;
        output_filename += filename;
        output_filename += "_HISTOGRAMS.root";
        pointer->write_histo_to_file(output_filename.Data());
    }

    // COUNTS    
    if(run_mode == 2 || run_mode == 5 || run_mode == 6 || run_mode == 7)
    {
        pointer->get_counts(N_COUNTS,N_COUNTS_TOTAL);
        pointer->get_counts_pairs(N_COUNTS_PAIRS,N_COUNTS_TOTAL_PAIRS);
        cout<<"--> N_COUNTS_TOTAL = "<<(Double_t)N_COUNTS_TOTAL<<" +/- "<<TMath::Sqrt((Double_t)N_COUNTS_TOTAL)<<endl;
        cout<<"--> N_COUNTS_TOTAL_PAIRS = "<<(Double_t)N_COUNTS_TOTAL_PAIRS<<" +/- "<<TMath::Sqrt((Double_t)N_COUNTS_TOTAL_PAIRS)<<endl;

        output_filename = output_dir;
        output_filename += filename;
        output_filename += "_COUNTS.txt";
        ofstream outfile(output_filename.Data());
        cout<<endl<<"Detectors"<<endl;
        for(Int_t i = 0; i < pointer->number_of_scintillators; i++)
        {
            Double_t err_1 = TMath::Sqrt(N_COUNTS[i]);
            Double_t err_2 = TMath::Sqrt(N_COUNTS_TOTAL);

            outfile<<i<<" , "<<(Double_t)N_COUNTS[i]/N_COUNTS_TOTAL<<" , "<<TMath::Sqrt(TMath::Power((Double_t)err_1/N_COUNTS_TOTAL,2) +
                                                                                        TMath::Power((Double_t)N_COUNTS[i]*err_2/(N_COUNTS_TOTAL*N_COUNTS_TOTAL),2))<<"\n";
            cout<<i<<" , "<<N_COUNTS[i]<<" , "<<N_COUNTS_TOTAL<<endl;
        }
        outfile.close();

        output_filename = output_dir;
        output_filename += filename;
        output_filename += "_COUNTS_PAIRS.txt";
        ofstream outfile_pairs(output_filename.Data());
        cout<<endl<<"Pairs"<<endl;
        for(Int_t i = 0; i < pointer->number_of_scintillators/2; i++)
        {
            Double_t err_1 = TMath::Sqrt(N_COUNTS_PAIRS[i]);
            Double_t err_2 = TMath::Sqrt(N_COUNTS_TOTAL_PAIRS);

            outfile_pairs<<i<<" , "<<(Double_t)N_COUNTS_PAIRS[i]/N_COUNTS_TOTAL_PAIRS<<" , "<<TMath::Sqrt(TMath::Power((Double_t)err_1/N_COUNTS_TOTAL_PAIRS,2) +
                                                                                        TMath::Power((Double_t)N_COUNTS_PAIRS[i]*err_2/(N_COUNTS_TOTAL_PAIRS*N_COUNTS_TOTAL_PAIRS),2))<<"\n";
            cout<<i<<" , "<<N_COUNTS_PAIRS[i]<<" , "<<N_COUNTS_TOTAL_PAIRS<<endl;
        }
        outfile_pairs.close();
    }
    // FREQUENCY
    if(run_mode == 3 || run_mode == 5 || run_mode == 7)
    {
        Double_t N12[pointer->number_of_cutting_angels] = {};
        Double_t N0[pointer->number_of_cutting_angels]  = {};

        Double_t Fini[pointer->number_of_cutting_angels]    = {};
        Double_t errFini[pointer->number_of_cutting_angels] = {};

        Double_t Cutting_angle[pointer->number_of_cutting_angels] = {};

        if(filename == "recoDataSimple_6382_4plane_offset0_0.root")     // In case of crystal in CH
            pointer->find_inelastic_frequency1(N12,N0,0,1);
        else    // In case of diffuser of crystal in AM or BKG
            pointer->find_inelastic_frequency2(N12,N0,0,1);

        output_filename = output_dir;
        output_filename += filename;
        output_filename += "_FINI.txt";
        ofstream outfile_fini(output_filename.Data());

        for(Int_t j = 0; j < pointer->number_of_cutting_angels; j++)
        {
            Cutting_angle[j] = 1e6*pointer->theta_cut[j];
            Fini[j] = (double)N12[j]/N0[j];
            errFini[j] = (TMath::Sqrt((double)N12[j]/TMath::Power(N0[j],2) + (double)TMath::Power(N12[j],2)*N0[j]/TMath::Power(N0[j],4)));

            outfile_fini<<Cutting_angle[j]<<" , "<<N12[j]<<" , "<<N0[j]<<" , "<<Fini[j]<<" , "<<errFini[j]<<"\n";
        }
        outfile_fini.close();
    }

    stop_time = time(NULL);
    cout<<"--> Running time is : "<<stop_time - start_time<<" seconds"<<endl;
    return 0;
}
