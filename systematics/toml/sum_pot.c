#include <iostream>
#include <fstream>
#include <string>
#include "TString.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"

void sum_pot(TString inputFiles = "*.root") {
    // ==========================================
    // USER INPUT
    // ==========================================
    TString histName = "TotalPOT";
    // ==========================================

    double total_pot = 0.0;
    int valid_files = 0;

    TChain dummyChain("dummy");
    if (inputFiles.EndsWith(".txt") || inputFiles.EndsWith(".list")) {
        std::ifstream infile(inputFiles.Data());
        std::string line;
        while (std::getline(infile, line)) {
            if (!line.empty()) dummyChain.Add(line.c_str());
        }
    } else {
        dummyChain.Add(inputFiles);
    }

    TObjArray *fileElements = dummyChain.GetListOfFiles();
    if (!fileElements || fileElements->GetEntries() == 0) {
        std::cerr << "Error: No valid ROOT files found matching your input!" << std::endl;
        return;
    }

    int total_files = fileElements->GetEntries();
    std::cout << "Summing POT from " << total_files << " files..." << std::endl;

    int barWidth = 50;

    for (int i = 0; i < total_files; ++i) {
        TString fileName = fileElements->At(i)->GetTitle();

        TFile *f = TFile::Open(fileName, "READ");
        if (f && !f->IsZombie()) {
            TH1 *potHist = (TH1*)f->Get(histName);
            if (potHist) {
                total_pot += potHist->Integral();
                valid_files++;
            }
            f->Close();
            delete f;
        }

        // ==========================================
        // PROGRESS BAR LOGIC
        // ==========================================
        float progress = (float)(i + 1) / total_files;
        std::cout << "[";
        int pos = barWidth * progress;
        for (int p = 0; p < barWidth; ++p) {
            if (p < pos) std::cout << "=";
            else if (p == pos) std::cout << ">";
            else std::cout << " ";
        }
        // The \r at the end forces the terminal to overwrite the same line
        std::cout << "] " << int(progress * 100.0) << " %\r";
        std::cout.flush();
        // ==========================================
    }

    // Print a final newline so the results don't overwrite the progress bar
    std::cout << std::endl;

    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Files successfully read: " << valid_files << " / " << total_files << std::endl;
    printf("TOTAL POT:               %.4e\n", total_pot);
    std::cout << "----------------------------------------" << std::endl;
}