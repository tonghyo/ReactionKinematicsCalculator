#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include <iostream>

void compare_proton() {
    cout << "Loading proton histograms for comparison..." << endl;
    
    // Open ROOT files
    TFile* file1 = TFile::Open("fusion_results_41Ti.root", "READ");
    TFile* file2 = TFile::Open("fusion_results_42V.root", "READ"); // You'll replace this with your second file
    
    if (!file1) {
        cout << "ERROR: Could not open first ROOT file!" << endl;
        return;
    }
    
    // List all histograms to find the correct names
    cout << "Available histograms in file1 (direct proton):" << endl;
    file1->ls();
    
    cout << "\nAvailable histograms in file2 (decay proton):" << endl;
    file2->ls();
    
    // Get proton angle histograms
    TH1D* his_proton1 = nullptr;
    TH1D* his_proton2 = nullptr;
    
    // Get the 2D histograms directly
    cout << "\nLoading 2D histograms..." << endl;
    TH2D* his_2d_proton1 = (TH2D*)file1->Get("his_product_4_Evsang");
    TH2D* his_2d_proton2 = (TH2D*)file2->Get("his_decay_1_Evsang");
    
    if (his_2d_proton1) {
        cout << "Found direct proton 2D histogram!" << endl;
    }
    if (his_2d_proton2) {
        cout << "Found decay proton 2D histogram!" << endl;
    }
    
    if (!his_2d_proton1 || !his_2d_proton2) {
        cout << "ERROR: Could not find proton histograms!" << endl;
        cout << "Please check the histogram names and update the code." << endl;
        file1->Close();
        file2->Close();
        return;
    }
    
    // Create comparison canvas
    TCanvas* c = new TCanvas("c", "Proton Theta vs Energy Comparison", 1200, 800);
    c->Divide(2, 2);
    
    // Plot 1: Direct fusion proton
    c->cd(1);
    his_2d_proton1->SetTitle("Direct Fusion Proton: Theta vs Energy");
    his_2d_proton1->GetXaxis()->SetTitle("Theta (degrees)");
    his_2d_proton1->GetYaxis()->SetTitle("Energy (MeV)");
    his_2d_proton1->Draw("COLZ");
    
    // Plot 2: Decay proton
    c->cd(2);
    his_2d_proton2->SetTitle("Decay Proton: Theta vs Energy");
    his_2d_proton2->GetXaxis()->SetTitle("Theta (degrees)");
    his_2d_proton2->GetYaxis()->SetTitle("Energy (MeV)");
    his_2d_proton2->Draw("COLZ");
    
    // Plot 3: Theta projection comparison
    c->cd(3);
    TH1D* his_theta1 = his_2d_proton1->ProjectionX("his_theta1");
    TH1D* his_theta2 = his_2d_proton2->ProjectionX("his_theta2");
    
    his_theta1->SetLineColor(kBlue);
    his_theta1->SetLineWidth(2);
    his_theta1->SetTitle("Theta Distribution Comparison");
    his_theta1->GetXaxis()->SetTitle("Theta (degrees)");
    his_theta1->GetYaxis()->SetTitle("Counts");
    
    // Set y-axis range for theta comparison
    double max_theta = TMath::Max(his_theta1->GetMaximum(), his_theta2->GetMaximum());
    his_theta1->GetYaxis()->SetRangeUser(0, max_theta * 1.1);
    his_theta1->Draw();
    
    his_theta2->SetLineColor(kRed);
    his_theta2->SetLineWidth(2);
    his_theta2->Draw("SAME");
    
    // Add legend
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(his_theta1, "Direct Fusion Proton", "l");
    legend->AddEntry(his_theta2, "Decay Proton", "l");
    legend->Draw();
    
    // Plot 4: Energy projection comparison
    c->cd(4);
    TH1D* his_energy1 = his_2d_proton1->ProjectionY("his_energy1");
    TH1D* his_energy2 = his_2d_proton2->ProjectionY("his_energy2");
    
    his_energy1->SetLineColor(kBlue);
    his_energy1->SetLineWidth(2);
    his_energy1->SetTitle("Energy Distribution Comparison");
    his_energy1->GetXaxis()->SetTitle("Energy (MeV)");
    his_energy1->GetYaxis()->SetTitle("Counts");
    
    // Set energy axis range: 0-50 MeV for decay proton
    his_energy1->GetXaxis()->SetRangeUser(0, 50);
    his_energy2->GetXaxis()->SetRangeUser(0, 50);
    
    // Set y-axis range for energy comparison
    double max_energy = TMath::Max(his_energy1->GetMaximum(), his_energy2->GetMaximum());
    his_energy1->GetYaxis()->SetRangeUser(0, max_energy * 1.1);
    his_energy1->Draw();
    
    his_energy2->SetLineColor(kRed);
    his_energy2->SetLineWidth(2);
    his_energy2->Draw("SAME");
    
    // Add legend
    TLegend* legend2 = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend2->AddEntry(his_energy1, "Direct Fusion Proton", "l");
    legend2->AddEntry(his_energy2, "Decay Proton", "l");
    legend2->Draw();
    
    // Print statistics
    cout << "\nDirect Proton Statistics:" << endl;
    cout << "  Mean theta: " << his_theta1->GetMean() << " degrees" << endl;
    cout << "  RMS theta: " << his_theta1->GetRMS() << " degrees" << endl;
    cout << "  Mean energy: " << his_energy1->GetMean() << " MeV" << endl;
    cout << "  RMS energy: " << his_energy1->GetRMS() << " MeV" << endl;
    cout << "  Total counts: " << his_2d_proton1->GetEntries() << endl;
    
    cout << "\nDecay Proton Statistics:" << endl;
    cout << "  Mean theta: " << his_theta2->GetMean() << " degrees" << endl;
    cout << "  RMS theta: " << his_theta2->GetRMS() << " degrees" << endl;
    cout << "  Mean energy: " << his_energy2->GetMean() << " MeV" << endl;
    cout << "  RMS energy: " << his_energy2->GetRMS() << " MeV" << endl;
    cout << "  Total counts: " << his_2d_proton2->GetEntries() << endl;
    
    // Save comparison plot
    c->SaveAs("proton_comparison.png");
    cout << "\nComparison plot saved as 'proton_comparison.png'" << endl;
    
}
