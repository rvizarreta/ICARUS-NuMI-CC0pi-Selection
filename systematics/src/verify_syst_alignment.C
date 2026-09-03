// Verify that each output table's systematic spline branches were built from
// that table's OWN systematic tree.
//
// Usage, on a machine with ROOT:
//   root -l -b -q 'verify_syst_alignment.C("/exp/icarus/data/users/rvizarr/medulla/icarus_numi_numu_mc_onbeam_offbeam_syst_gundam.root","/exp/icarus/data/users/rvizarr/medulla/icarus_numi_numu_mc_onbeam_offbeam_syst_ppfx.root")'
//
// Check 1 (entry counts) is the decisive one and needs only the output file.
// Before the fix, `sideband` reported 34420 entries while its spline branches
// held 98983 -- the mismatch that let the bug go unnoticed. After the fix all
// branch counts must equal their tree's count.
//
// Check 2 (values) proves WHICH source tree the splines came from, by printing
// the stored knots next to the same row of both candidate input trees.

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TGraph.h>
#include <TClonesArray.h>
#include <TObjArray.h>
#include <TString.h>
#include <iostream>
#include <vector>

// First branch on `t` whose name contains `pattern` (and is not a _sigma one).
TString first_branch(TTree * t, const char * pattern)
{
  TObjArray * br = t->GetListOfBranches();
  for(int i = 0; i < br->GetEntries(); ++i)
    {
      TString n = ((TBranch*)br->At(i))->GetName();
      if(n.Contains(pattern) && !n.EndsWith("_sigma")) return n;
    }
  return TString("");
}

void check_counts(TFile * fo)
{
  const char * tables[3] = {"events/full/selected", "events/full/sideband", "events/full/signal"};
  const char * pats[3]   = {"multisigma", "hysyst", "variation"};

  std::cout << "\n===== Check 1: tree entries vs branch entries =====\n" << std::endl;
  bool all_ok = true;

  for(int it = 0; it < 3; ++it)
    {
      TTree * t = (TTree*)fo->Get(tables[it]);
      if(t == nullptr) { std::cout << tables[it] << " : MISSING\n"; continue; }

      Long64_t nt = t->GetEntries();
      std::cout << tables[it] << "  tree entries = " << nt << std::endl;

      for(int ip = 0; ip < 3; ++ip)
        {
          TString bn = first_branch(t, pats[ip]);
          if(bn == "") { std::cout << "    (no '" << pats[ip] << "' branch)\n"; continue; }

          Long64_t nb = t->GetBranch(bn)->GetEntries();
          bool ok = (nb == nt);
          if(!ok) all_ok = false;
          std::cout << "    " << (ok ? "OK   " : "FAIL ") << bn << "  entries = " << nb;
          if(!ok) std::cout << "   <-- should be " << nt;
          std::cout << std::endl;
        }
    }
  std::cout << "\n  ==> " << (all_ok ? "PASS: every spline branch matches its tree."
                                     : "FAIL: a branch count differs -- the old binary/file is still in play.")
            << std::endl;
}

void check_values(TFile * fo, TFile * fi, int nev)
{
  TTree * out = (TTree*)fo->Get("events/full/sideband");
  TTree * sb  = (TTree*)fi->Get("events/full/sideband_multisigmaTree");
  TTree * sel = (TTree*)fi->Get("events/full/selected_multisigmaTree");
  if(out == nullptr || sb == nullptr || sel == nullptr)
    { std::cout << "\n(skipping check 2: a tree is missing)\n"; return; }

  TString dial = first_branch(out, "multisigma");
  if(dial == "") { std::cout << "\n(skipping check 2: no multisigma branch)\n"; return; }

  std::cout << "\n===== Check 2: which source tree did sideband splines come from? =====\n"
            << "dial: " << dial << "\n"
            << "The 'stored' row must match 'sideband_multisigmaTree', not 'selected_multisigmaTree'.\n"
            << std::endl;

  TClonesArray * arr = nullptr;
  std::vector<float> * w_sb  = nullptr;
  std::vector<float> * w_sel = nullptr;
  out->SetBranchAddress(dial, &arr);
  sb ->SetBranchAddress(dial, &w_sb);
  sel->SetBranchAddress(dial, &w_sel);

  for(int i = 0; i < nev; ++i)
    {
      out->GetEntry(i); sb->GetEntry(i); sel->GetEntry(i);
      TGraph * g = (arr && arr->GetEntries() > 0) ? (TGraph*)arr->At(0) : nullptr;

      std::cout << "row " << i << std::endl;
      std::cout << "   stored (sorted by nsigma) :";
      if(g) for(int k = 0; k < g->GetN(); ++k) printf(" %8.5f", g->GetY()[k]);
      std::cout << std::endl;

      std::cout << "   sideband_multisigmaTree   :";
      if(w_sb)  for(size_t k = 0; k < w_sb->size();  ++k) printf(" %8.5f", (*w_sb)[k]);
      std::cout << std::endl;

      std::cout << "   selected_multisigmaTree   :";
      if(w_sel) for(size_t k = 0; k < w_sel->size(); ++k) printf(" %8.5f", (*w_sel)[k]);
      std::cout << "\n" << std::endl;
    }
  std::cout << "(stored is sorted ascending in nsigma and carries an extra 1.0 at nsigma=0;\n"
            << " the input vectors are in raw [-1,1,-2,2,-3,3] order, so compare the value SET.)\n"
            << std::endl;
}

void verify_syst_alignment(const char * outfile, const char * infile = "", int nev = 3)
{
  TFile * fo = TFile::Open(outfile, "READ");
  if(fo == nullptr || fo->IsZombie()) { std::cout << "cannot open " << outfile << std::endl; return; }

  check_counts(fo);

  if(TString(infile) != "")
    {
      TFile * fi = TFile::Open(infile, "READ");
      if(fi && !fi->IsZombie()) check_values(fo, fi, nev);
      else std::cout << "\n(cannot open input file; skipping check 2)\n";
    }
}
