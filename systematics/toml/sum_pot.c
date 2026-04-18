void sum_pot() {
    const char* dir = "/pnfs/icarus/persistent/users/dcarber/spine/NuMI_CV_ext/v09_89_01_02p02_2/combined_files";

    TString cmd = TString::Format("find %s -maxdepth 1 -name '*.flat.root' 2>/dev/null", dir);
    TString file_list = gSystem->GetFromPipe(cmd.Data());
    TObjArray* files = file_list.Tokenize("\n");
    int n = files->GetEntries();
    std::cout << "Combined files: " << n << std::endl << std::endl;

    double total_pot = 0.0;
    double total_events = 0.0;
    long long total_rec = 0;

    std::cout << Form("%-50s %-14s %-12s %-12s", "File", "POT", "Events", "recTree") << std::endl;
    std::cout << std::string(92, '-') << std::endl;

    for (int i = 0; i < n; ++i) {
        TString fname = ((TObjString*)files->At(i))->GetString();
        TFile* f = TFile::Open(fname, "READ");
        if (!f || f->IsZombie()) {
            std::cout << "FAILED: " << fname << std::endl;
            if (f) delete f;
            continue;
        }

        TH1D* h_pot = (TH1D*) f->Get("TotalPOT");
        TH1D* h_ev  = (TH1D*) f->Get("TotalEvents");
        TTree* t    = (TTree*) f->Get("recTree");

        double pot = h_pot ? h_pot->Integral() : 0.0;
        double evs = h_ev  ? h_ev->Integral()  : 0.0;
        long long rec_ent = t ? t->GetEntries() : 0;

        TString short_name = gSystem->BaseName(fname.Data());
        std::cout << Form("%-50s %.4e  %.3e  %lld",
                          short_name.Data(), pot, evs, rec_ent) << std::endl;

        total_pot += pot;
        total_events += evs;
        total_rec += rec_ent;

        f->Close();
        delete f;
    }

    std::cout << std::endl;
    std::cout << "======================================================================" << std::endl;
    std::cout << Form("Sum POT (combined)      : %.6e", total_pot) << std::endl;
    std::cout << Form("Sum TotalEvents         : %.3e", total_events) << std::endl;
    std::cout << Form("Sum recTree entries     : %lld", total_rec) << std::endl;
    std::cout << std::endl;
    std::cout << Form("Medulla file POT        : %.6e", 7.78521e20) << std::endl;
    std::cout << Form("Ratio combined / medulla: %.4f", total_pot/7.78521e20) << std::endl;
    std::cout << "======================================================================" << std::endl;
}