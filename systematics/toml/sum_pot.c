void sum_source_pot() {
    const char* dir = "/pnfs/icarus/persistent/users/dcarber/spine/NuMI_CV_ext/v09_89_01_02p02_2";
    const char* pattern = "*.root";

    TString cmd = TString::Format("ls %s/%s 2>/dev/null", dir, pattern);
    TString file_list = gSystem->GetFromPipe(cmd.Data());
    TObjArray* files = file_list.Tokenize("\n");
    int n = files->GetEntries();
    std::cout << "Source files: " << n << std::endl;

    double total_pot = 0.0;
    double total_events_hist = 0.0;
    long long total_rec_entries = 0;
    int n_ok = 0;
    int n_missing_pot = 0;

    for (int i = 0; i < n; ++i) {
        TString fname = ((TObjString*)files->At(i))->GetString();
        TFile* f = TFile::Open(fname, "READ");
        if (!f || f->IsZombie()) { if (f) delete f; continue; }

        TH1D* h_pot = (TH1D*) f->Get("TotalPOT");
        TH1D* h_ev  = (TH1D*) f->Get("TotalEvents");
        TTree* t    = (TTree*) f->Get("recTree");

        double pot = h_pot ? h_pot->Integral() : 0.0;
        double evs = h_ev  ? h_ev->Integral()  : 0.0;
        long long rec_ent = t ? t->GetEntries() : 0;

        if (!h_pot || pot == 0) n_missing_pot++;

        total_pot += pot;
        total_events_hist += evs;
        total_rec_entries += rec_ent;
        n_ok++;

        if (i % 20 == 0) std::cout << "  " << i << "/" << n << " done" << std::endl;

        f->Close();
        delete f;
    }

    std::cout << std::endl;
    std::cout << "======================================================================" << std::endl;
    std::cout << "Opened OK                   : " << n_ok << " / " << n << std::endl;
    std::cout << "Files missing/zero POT      : " << n_missing_pot << std::endl;
    std::cout << Form("Sum POT (source CAFs)       : %.6e", total_pot) << std::endl;
    std::cout << Form("Sum TotalEvents histogram   : %.3e", total_events_hist) << std::endl;
    std::cout << Form("Sum recTree entries         : %lld", total_rec_entries) << std::endl;
    std::cout << std::endl;
    std::cout << "Compare to:" << std::endl;
    std::cout << "  Merged medulla file POT    : 7.78521e20" << std::endl;
    std::cout << Form("  Ratio source / merged      : %.4f", total_pot/7.78521e20) << std::endl;
    std::cout << "======================================================================" << std::endl;
}