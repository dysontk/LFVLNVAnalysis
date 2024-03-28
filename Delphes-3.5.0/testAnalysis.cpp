{
  // Create chain of root trees g++ `root-config --cflags --libs` -o testAnalysis testAnalysis.cpp
  TChain chain("Delphes");
  chain.Add("delphes_output.root");

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  // Get pointers to branches used in this analysis
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");

  // Book histograms
  TH1 *histElectronPT = new TH1F("electron pt", "electron P_{T}", 50, 0.0, 100.0);

  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {

    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    // If event contains at least 1 electron
    if(branchElectron->GetEntries() > 0)
    {
      // Take first electron
      Electron *electron = (Electron*) branchElectron->At(0);

      // Plot electron transverse momentum
      histElectronPT->Fill(electron->PT);

      // Print electron transverse momentum
      cout << electron->PT << endl;
    }

  }

  // Show resulting histograms
  histElectronPT->Draw();
}