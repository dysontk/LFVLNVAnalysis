#include <iostream>
#include <string>
#include <complex>

#include "TCanvas.h"
#include "TChain.h"


#include "TROOT.h"
#include "TClonesArray.h"
#include "TRint.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TChain.h"

#include "TMath.h"
#include "TFile.h"
#include "TSystem.h"
#include "TGStatusBar.h"
#include "TSystem.h"
// #include "TXMLEngine.h"
#include "TTree.h"

#include "ExRootClasses.h"
#include "ExRootClassifier.h"

#include "ExRootFactory.h"
#include "ExRootFilter.h"
#include "ExRootLHEFReader.h"


#include "ExRootProgressBar.h"
#include "ExRootResult.h"
//#include "ExRootSTDHEPReader.h"
//#include "ExRootStream.h"

#include "ExRootTreeBranch.h"
#include "ExRootTreeReader.h"
#include "ExRootTreeWriter.h"
#include "ExRootUtilities.h"

// #include "fastjet/PseudoJet.hh"
// #include "fastjet/JetDefinition.hh"
// #include "fastjet/ClusterSequence.hh"


#include "classes/DelphesClasses.h"

//#include <random>
#include <cstdlib>
#include <time.h>
#include <unistd.h>
#include <iterator>
#include<stdio.h>
#include "TRandom3.h"

using namespace std;
using namespace TMath;
// using namespace fastjet;

#include "observables.h"
#include "ran.h"
// #include "mt2_bisect.h"
// #include "mt2w_bisect.h"

bool VERBOSE = false; // Allows for some troubleshooting or extra detail if true.
bool lowlepcut = false; //turns on and off Some sort of cut on the low energy lepton?????????????????????????

// For the simplification of simulated data in early stages,
// we do not have branches for MET nor Muons. The following Bools are to account for this. 
// When we use data that has MET and Muons, then we set these to true and it'll function properly 
// Used in lines 147-149, 234-237
bool hasMET = false;
bool hasMu = false;
// This is used if we want to post-simulate the LFV with ratios from CMS paper (true) or just leave as single flavor (false)
bool simLFV = false;

void Export_TH1F (TH1F * target_H1F, ofstream &out_file){
    double num_entries = target_H1F->GetEntries();
    for( int i = 1; i <= target_H1F->GetNbinsX(); i++){
        out_file<<target_H1F->GetXaxis()->GetBinCenter(i)<<" "<<(target_H1F->GetBinContent(i)/num_entries)<<endl;
    }//return normalized histogram
}


Bool Below_DeltaR_Same(vector<PseudoJet> par, double delta_R_min);
Bool Below_DeltaR_Same(vector<PseudoJet> par, double delta_R_min){
    if(par.size() < 2){
        return false;
    }
    else {
        for(int i=0; i<par.size()-1; i++){
            for (int j=i+1; j<par.size(); j++) {
                if( par[i].delta_R(par[j]) < delta_R_min ) return true;
            }
        }
    }
    return false;
}

Bool Below_DeltaR_Diff(vector<PseudoJet> par1, vector<PseudoJet> par2, double delta_R_min);
Bool Below_DeltaR_Diff(vector<PseudoJet> par1, vector<PseudoJet> par2, double delta_R_min){
    if(par1.size()==0 || par2.size()==0) {
        return false;
    }
    else {
        for(int i=0; i<par1.size(); i++){
            for (int j=0; j<par2.size(); j++) {
                if( par1[i].delta_R(par2[j]) < delta_R_min ) return true;
            }
        }
    }
    return false;
}


/*
Hello, 
This code is meant to isolate a dilepton dijet Lepton Number Violating process from the inputted
data. The purpose is to get rid of mostly jet fakes and diboson backgrounds. In commenting
the cuts made, I will reference https://arxiv.org/abs/1806.10905 
Search for heavy Majorana neutrinos in same-sign dilepton channels in proton-proton collisions at √s= 13 TeV.
Which was done by CMS collaboration. The idea of this code is to mostly replicate the analysis that they did on
experimental data. So any reference to CMS Sections or tables etc. are references to that paper.
*/



int main(int argc, const char * argv[])
{
    
    //Pulls and arranges data as needed.
    TChain chain("Delphes");

    for(int i=1; i<argc; i++)
    {
        cout << chain.Add(argv[i])<< endl;
    }

    ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
    Long64_t NumEntries = treeReader->GetEntries();
    cout << "There are "<< NumEntries <<" Entries." <<endl;
    TClonesArray *branchJet = treeReader->UseBranch("Jet"); 
    TClonesArray *branchElectron = treeReader->UseBranch("Electron");
    TClonesArray *branchMuon;
    if (hasMu) branchMuon = treeReader->UseBranch("Muon");
    TClonesArray *branchMET;
    if (hasMET) branchMET = treeReader->UseBranch("MET");


    TRootLHEFEvent *event;
    TRootLHEFParticle *particle;

    vector <PseudoJet> all_jets, w_jet_pairs; // all_jets = f_all_jets and w_jet_pairs = f_w_jet from Gang's main.cpp
    vector <PseudoJet> v_eM, v_eP, v_muM, v_muP, v_lep, v_e, v_mu, v_lepP, v_lepM;
    vector <PseudoJet> v_MET, b_jets;

    // all_jets = a vector with the every jet object for each entry
    // w_jet_pairs = a vector with pairs of jets that have inv. mass ~ M_W for each entry
    // v_... = {e=electron, mu=muon, P=Plus/+, M=Minus/-, lep=any lepton}
    //v_MET = vector for the Missing transverse energy
    // b_jets = vector of b-tagged jets. 

    MissingET *met;
    Jet *jet;

    // HISTOGRAMS---------------------------------------
    //                                                                  , bins, xlow, xhigh)
    // Obj. name; Vert. axis vs. Horiz. axis; images produced by this
    TH1F *MW2j = new TH1F("Inv_Mass_2Jets_close_to_W", "Inv. Mass 2 Jets", 100, 0, 140);
    // 
    TH1F *MW2j2l = new TH1F("Inv_Mass_2Jets_close_to_W_2l", "Inv. Mass 2 Jets 2 Lep",  70, 0, 700);
    TH1F *M2l;
    if (lowlepcut)  M2l = new TH1F("Inv_Mass_2l", "Inv. Mass 2 Lep above 1GeV",  70, 0, 300);
    else M2l = new TH1F("Inv_Mass_2l", "Inv. Mass 2 Lep",  70, 0, 300);
    TH1F *MW2j1l_0 = new TH1F("Inv_Mass_2Jets_close_to_W_1l_0", "Inv. Mass 2 Jets + Lep_0",  20, 0, 400);
    TH1F *MW2j1l_1  = new TH1F("Inv_Mass_2Jets_close_to_W_1l_1", "Inv. Mass 2 Jets + Lep_1",  20, 0, 250);
// Figure out how to best use this later based on how I do the cuts
    //    // event counts after each cut
    // vector<int> cntCuts;
    // cntCuts.clear();
    /*
    The four cuts are as follows:
    cut 0: signal definition -- dilepton + dijet
    cut 1: preselection -- ?
    cut 2: HMSR1 -- ?
    cut 3: Misc. -- ?
    numCutCats below counts how many events pass each cut. 
    */
    int numCutCats[4] = {0, 0, 0, 0}; //events that pass 0: signal def, 1: preselection, 2: HMSR1, 3: Misc.
    // cout << sizeof(numCutCats)<< endl;
    for (int i=0; i<sizeof(numCutCats)/sizeof(int); i++) cout << numCutCats[i] << ", ";
    cout << endl;
    // Okay, so all I want to do to start is to write the cut that 
    // Compares pairs of Jets to W boson mass

    if (NumEntries ==0) {

      cout<<"  ------------------  "<<endl;
      cout<<"  no event analyzed  "<<endl;
      cout<<"  ------------------  "<<endl;
    }

    else
    {
        cout << "Loading Events ..." << endl;
        if(VERBOSE) cout <<"There are "<<  NumEntries << " Entries." << endl;

        // Is there a reason to use Int_t as opposed to int
        for(int entry=0; entry < NumEntries; entry++)
        {
            if (VERBOSE) cout<< "I am on Entry " << entry<< endl;
            all_jets.clear();
            w_jet_pairs.clear();
            v_eM.clear();
            v_eP.clear();
            v_e.clear();
            v_mu.clear();
            v_muP.clear();
            v_muM.clear();
            v_lep.clear();
            v_lepP.clear();
            v_lepM.clear();
            v_MET.clear();
            b_jets.clear();

            // cout << "cleared vectors"<< endl;
            treeReader->ReadEntry(entry);
            int numJet = branchJet->GetEntries();
            int numEl = branchElectron->GetEntries();
            int numMu;
            int numMET; 
            if(hasMu) numMu = branchMuon->GetEntries(); 
            else numMu = 0;
            if(hasMET) numMET = branchMET->GetEntries();
            else numMET = 0;
            //I should diagram the rest of this out
            // I want to determine the pairs of jets that are closest to W mass

            // The following loops collects the data as PsuedoJet objects and organizes them in the
            // pre-defined vectors
            PseudoJet tempEvent;
            for(int jt=0; jt < numJet; jt++)
            {
                Jet *jettmp = (Jet *) branchJet->At(jt);
                //for some reason this .reset() causes it not to compile
                tempEvent.reset((jettmp->P4()).Px(), (jettmp->P4()).Py(), (jettmp->P4()).Pz(), (jettmp->P4()).E());
                // cout << jettmp->BTag << endl;
                if(jettmp->BTag == 1) b_jets.push_back(tempEvent);
                all_jets.push_back(tempEvent);
                // cout << "added a jet"<< endl;
            }
            
            //Time to collect some leptons
            //Electrons
            // cout << "Bout to make electrons"<< endl;
            // if (VERBOSE) cout<< "numEl = " << numEl << endl;
            for(int e=0; e < numEl; e++)
            {
                Electron *etmp = (Electron *) branchElectron->At(e);
                tempEvent.reset((etmp->P4()).Px(), (etmp->P4()).Py(), (etmp->P4()).Pz(), (etmp->P4()).E());

                v_lep.push_back(tempEvent);
                v_e.push_back(tempEvent);
                if (etmp->Charge == 1) 
                {
                    v_eP.push_back(tempEvent); // Positive 
                    v_lepP.push_back(tempEvent);
                }
                else 
                {
                    v_eM.push_back(tempEvent); // Negative
                    v_lepM.push_back(tempEvent);
                }

            }
            //Muons
            // cout << "bout to make muons"<< endl;
            // if (VERBOSE) cout<< "numMu = " << numMu << endl;

            for(int m=0; m < numMu; m++)
            {
                Muon *mutmp = (Muon *) branchMuon->At(m);
                tempEvent.reset((mutmp->P4()).Px(), (mutmp->P4()).Py(), (mutmp->P4()).Pz(), (mutmp->P4()).E());

                v_lep.push_back(tempEvent);
                v_mu.push_back(tempEvent);
                if (mutmp->Charge == 1) 
                {
                    v_muP.push_back(tempEvent); // Positive 
                    v_lepP.push_back(tempEvent);
                }
                else
                {
                    v_muM.push_back(tempEvent); // Negative
                    v_lepM.push_back(tempEvent);
                }

            }
            // cout << "i'm about to sort vectors" << endl;
            sort(v_lep.begin(), v_lep.end(), sort_by_pt());
            sort(v_eM.begin(), v_eM.end(), sort_by_pt());
            sort(v_eP.begin(), v_eP.end(), sort_by_pt());
            sort(v_muM.begin(), v_muM.end(), sort_by_pt());
            sort(v_muP.begin(), v_muP.end(), sort_by_pt());
            sort(v_e.begin(), v_e.end(), sort_by_pt());
            sort(v_mu.begin(), v_mu.end(), sort_by_pt());
            // cout << "I made it past the making of vectors"<< endl;
            // Signal definition

            // First make sure that there are s.s. dilep pairs
            if (v_lepP.size() < 2 && v_lepM.size() < 2) 
            {
                if (VERBOSE) cout << "No s.s. dilepton pair"<< endl;
                continue;
            } 
            // I don't know if I should make a vector for all + leptons and all - leptons
            // Which would allow for eμ pairs 

            if (VERBOSE) cout << "There is an s.s. dilepton pair"<< endl;
            //Now it is time to find the pair of jets that are close to W
            if (all_jets.size() < 1) 
            {   
                
                if (VERBOSE)cout << "Not enough jets"<< endl;
                continue;
            }
            numCutCats[0]++;

            // Preselection Criteria: CMS Sec 5.1
            
            //Third lepton too much energy
            if (v_lep.size()>2)
            {
                if (v_lep[2].pt() > 10) continue;//GeV
            }

            
            double m_ll = 0;
            // double m_mm = 0;
            // if (v_e.size()>2)
            // {
            //     m_ll = (v_e[0]+v_e[1]).m();
            // }
            // if (v_mu.size()>2)
            // {
            //     m_ll = (v_mu[0]+v_mu[1]).m();
            // }
            if (v_lep.size()>2)
            {
                m_ll = (v_lep[0]+v_lep[1]).m();
            }
            // m_ll must be at least 10 GeV
            if (m_ll < 10) continue;
            //Lepton pair from Z 
            if (abs(m_ll - 91.2) > 20) continue;

            numCutCats[1]++;

            // High Mass SR 1
            // The preceeding three cuts are the definition of High Mass Signal Region 1 (HM SR1)
            // CMS Table 1

            // Each jet must be above 25 GeV pt
            int numJet2 = all_jets.size();
            if (VERBOSE) cout << "Is the length of jet vector the same as branch size? " << (numJet2 == numJet);
            int htsum = 0;
            if (all_jets[0].pt()<=25) continue;
            for (int i=0; i < numJet2; i++)
            {
                double this_pt = all_jets[i].pt();
                htsum += this_pt;
                // if( this_pt < 25) continue;
            }

            for (int j=0; j < v_lep.size(); j++)
            {
                htsum += v_lep[j].pt();
            }
            // Ratio of missing Trans. momentum^2 and total PT less than 15 GeV
            if (pow(v_MET[0].pt(),2)/htsum > 15) continue;
            
            // For the next part, I need to find Δm_Wj which is the inv mass of
            // the pair of jets with inv mass closes to M_W

            int JLowPairInd[2] = {0,0};
            double LowDiff = 1000;
            for(int j1=0; j1<all_jets.size(); j1++)
            {
                for(int j2=j1; j2<all_jets.size(); j2++)
                {
                    // cout << "Current Pair (" << j1 << "," << j2<< ")"<< endl;
                    double thisDiff = abs((all_jets[j1]+all_jets[j2]).m()-80.4);
                    // cout << "Inv M is off from W by " << thisDiff << endl;
                    // cout << "Low diff: " << LowDiff << endl;
                    if (thisDiff < LowDiff)
                    {
                        LowDiff = thisDiff;
                        JLowPairInd[0] = j1;
                        JLowPairInd[1] = j2;
                    }
                    else continue;
                }
            }
            // if (!JLowPairInd[0] && !JLowPairInd[1])
            // {
            //     cout << "For some reason, I couldn't pick a smallest"<< endl;
            // }

            double M_Wj = (w_jet_pairs[0]+w_jet_pairs[1]).m();

            if (M_Wj < 30 || M_Wj > 150) continue;

            numCutCats[2]++;
            // We have concluded the HMSR1 cuts


            // Here come the Miscellaneous Cuts
            // First, we want angularly well separated events
            // ΔR information from CMS Sec. 4.1 & 4.2

            if (Below_DeltaR_Diff(v_lep, all_jets, 0.4)) continue;
            if (Below_DeltaR_Same(v_lep, 0.4)) continue;
            if (Below_DeltaR_Same(all_jets, 0.4)) continue;

            //Now we remove all the B_tagged events (This doesn't seem to come from CMS paper)
            //Not sure why we are doing it.
            if(!b_jets.size()) continue;
            /*
                Now We do a bit of a funky thing
                For now, we are not simulating any μ events. However according to CMS Sec 5 (paragraph 1)
                we want events with leptons at least loosely isolated. We use the "offline requirements".
                These requirements are different for ee, μμ, and eμ. 
                since we have no μ events, then we have to use ratios to split up ee events.
                I don't yet know where these ratios came from.
            */
           int lPairType = 0; // ee->0, μμ -> 1, eμ -> 2
           double r1 = 0.18;
           double r2 = 0.31;
           if (hasMu)
           {
                if (abs(v_lep[0].m() - 0.511*pow(10, 3)) < 0.01) // is 0.01 a good cut off for being an electron?
                {
                    if (abs(v_lep[1].m() - 0.511*pow(10,3)) < 0.01) lPairType = 0;
                    else lPairType = 2;
                }
                else
                {
                    if (abs(v_lep[1].m() - 0.511*pow(10,3)) < 0.01) lPairType = 2;
                    else lPairType = 1; 
                }
           }
            else if (simLFV)
           {
                double rs = gRandom->Uniform(); // random number 0-1
                if (rs < r1) lPairType = 0; // 0-0.18 => ee
                else lPairType = (r2 < rs) ? 2 : 1; // 0.18-0.31 =>μμ, 0.31-1 => eμ
           }
        //    CMS Section 5 paragraph 1. This seems to be simulating something about the detector
           int leadingThresh[3] = {25, 20, 25}; // GeV; corresponds to lPairType 0,1,2 indices
           int trailingThresh[3] = {15, 10, 10};
           if ((v_lep[0].pt() < leadingThresh[lPairType] || v_lep[1].pt() < trailingThresh[lPairType])) continue;
           numCutCats[3]++;

            
            w_jet_pairs.push_back(all_jets[JLowPairInd[0]]);
            w_jet_pairs.push_back(all_jets[JLowPairInd[1]]);
            if (VERBOSE) cout << "I'm filling the histograms now"<< endl;
            MW2j->Fill((w_jet_pairs[0]+w_jet_pairs[1]).m());
            MW2j2l->Fill((w_jet_pairs[0]+w_jet_pairs[1]+v_lep[0]+v_lep[1]).m());
            MW2j1l_0->Fill((w_jet_pairs[0]+w_jet_pairs[1]+v_lep[0]).m());
            MW2j1l_1->Fill((w_jet_pairs[0]+w_jet_pairs[1]+v_lep[1]).m());
            if (lowlepcut)
            {
                if ((v_lep[0]+v_lep[1]).m()>1) M2l->Fill((v_lep[0]+v_lep[1]).m());
            }
            else
            {
                M2l->Fill((v_lep[0]+v_lep[1]).m());
            }
            
            // TIME TO PLOT

            

            // return 0;
        }

        TCanvas *c1 = new TCanvas("c1", "ROOT Canvas", 900, 20, 540, 550);
        
        MW2j->GetXaxis()->SetTitle("GeV");
        MW2j->Draw();
        c1->SaveAs("Mass_2jW.png");

        MW2j2l->GetXaxis()->SetTitle("GeV");
        MW2j2l->Draw();
        c1->SaveAs("Mass_2jW2l.png");

        MW2j1l_0->GetXaxis()->SetTitle("GeV");
        MW2j1l_0->Draw();
        c1->SaveAs("Mass_2jW1l0.png");
x
        MW2j1l_1->GetXaxis()->SetTitle("GeV");
        MW2j1l_1->Draw();
        c1->SaveAs("Mass_2jW1l1.png");

        M2l->GetXaxis()->SetTitle("GeV");
        M2l->Draw();
        if (lowlepcut) c1->SaveAs("Mass_l2_above1GeV.png");
        else c1->SaveAs("Mass_l2.png");

        cout << "out of " << NumEntries << endl;
        for (int i=0; i < sizeof(numCutCats)/sizeof(int); i++) cout << numCutCats[i] << endl;
    }

}