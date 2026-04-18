void sum_pot() {
    const char* dir = "/pnfs/icarus/scratch/users/rvizarr/1muNp0pi_uncontained_project/output";
    const char* pattern = "output_systematics_jobid*.root";

    TString cmd = TString::Format("ls %s/%s 2>/dev/null", dir, pattern);
    TString file_list = gSystem->GetFromPipe(cmd.Data());
    TObjArray* files = file_list.Tokenize("\n");
    int n = files->GetEntries();
    std::cout << "Files: " << n << std::endl;

    double total_pot = 0.0;
    long long total_entries = 0;
    int n_ok = 0;

    for (int i = 0; i < n; ++i) {
        TString fname = ((TObjString*)files->At(i))->GetString();
        TFile* f = TFile::Open(fname, "READ");
        if (!f || f->IsZombie()) { if (f) delete f; continue; }

        TH1D* h = (TH1D*) f->Get("events/nominal/POT");
        TTree* t = (TTree*) f->Get("events/nominal/signal");

        double pot = h ? h->Integral() : 0.0;
        long long ent = t ? t->GetEntries() : 0;

        total_pot += pot;
        total_entries += ent;
        n_ok++;

        if (i % 20 == 0) std::cout << "  " << i << "/" << n << " done" << std::endl;

        f->Close();
        delete f;
    }

    std::cout << std::endl;
    std::cout << "Opened OK   : " << n_ok << " / " << n << std::endl;
    std::cout << Form("Sum POT     : %.6e", total_pot) << std::endl;
    std::cout << Form("Sum entries : %lld", total_entries) << std::endl;
    std::cout << Form("Merged POT  : 7.78521e20") << std::endl;
    std::cout << Form("Ratio POT   : %.4f", total_pot/7.78521e20) << std::endl;
    std::cout << Form("Ratio ent.  : %.4f", (double)total_entries/113501.0) << std::endl;
}