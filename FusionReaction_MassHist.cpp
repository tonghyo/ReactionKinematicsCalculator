#include "FusionReaction.h"
#include "TLegend.h"

// Read mass file and set particle masses
void FusionReaction::ReadMassFile(const char* filename) {
    ifstream mass_file(filename);
    int A_temp, Z_temp;
    double M_temp;
    
    cout << "Reading mass file: " << filename << endl;
    cout << "Looking for masses:" << endl;
    cout << "Beam: " << A_beam << " (Z=" << Z_beam << ")" << endl;
    cout << "Target: " << A_target << " (Z=" << Z_target << ")" << endl;
    for (int i = 0; i < product_A.size(); i++) {
        cout << "Product " << i+1 << ": " << product_names[i] << " (" << product_A[i] << ", Z=" << product_Z[i] << ")" << endl;
    }
    
    // Print decay products if enabled
    if (decay_enabled) {
        cout << "Decay products:" << endl;
        for (int i = 0; i < decay_A.size(); i++) {
            cout << "Decay " << i+1 << ": " << decay_names[i] << " (" << decay_A[i] << ", Z=" << decay_Z[i] << ")" << endl;
        }
    }
    
    bool beam_found = false, target_found = false, parent_found = false;
    vector<bool> products_found(product_A.size(), false);
    vector<bool> decay_found(decay_A.size(), false);
    
    while (mass_file >> A_temp >> Z_temp >> M_temp) {
        
        if (A_temp == A_beam && Z_temp == Z_beam && !beam_found) {
            M_beam = M_temp;
            beam_found = true;
            cout << "Found beam mass: " << A_beam << " (Z=" << Z_beam << ") = " << M_temp << " MeV" << endl;
        }
        if (A_temp == A_target && Z_temp == Z_target && !target_found) {
            M_target = M_temp;
            target_found = true;
            cout << "Found target mass: " << A_target << " (Z=" << Z_target << ") = " << M_temp << " MeV" << endl;
        }
        
        // Find parent particle mass for reconstruction
        if (enable_product_reconstruction && A_temp == parent_A && Z_temp == parent_Z && !parent_found) {
            parent_mass = M_temp;
            parent_found = true;
            cout << "Found parent mass: " << parent_name << " (" << parent_A << ", Z=" << parent_Z << ") = " << M_temp << " MeV" << endl;
        }
        
        // Find product masses
        for (int i = 0; i < product_A.size(); i++) {
            if (A_temp == product_A[i] && Z_temp == product_Z[i] && !products_found[i]) {
                product_masses[i] = M_temp;
                products[i].mass = M_temp;
                products_found[i] = true;
                cout << "Found product mass: " << product_names[i] << " (" << product_A[i] << ", Z=" << product_Z[i] << ") = " << M_temp << " MeV" << endl;
            }
        }
        
        // Find decay product masses
        for (int i = 0; i < decay_A.size(); i++) {
            if (A_temp == decay_A[i] && Z_temp == decay_Z[i] && !decay_found[i]) {
                decay_masses[i] = M_temp;
                decay_found[i] = true;
                cout << "Found decay product mass: " << decay_names[i] << " (" << decay_A[i] << ", Z=" << decay_Z[i] << ") = " << M_temp << " MeV" << endl;
            }
        }
        
        
        // Check if all masses found
        bool all_found = beam_found && target_found;
        for (int i = 0; i < products_found.size(); i++) {
            all_found = all_found && products_found[i];
        }
        for (int i = 0; i < decay_found.size(); i++) {
            all_found = all_found && decay_found[i];
        }
        // Also check parent mass if product reconstruction is enabled
        if (enable_product_reconstruction) {
            all_found = all_found && parent_found;
        }
        if (all_found) {
            cout << "All masses found successfully!" << endl;
            break;
        }
        
        if (A_temp == 0) break;
    }
    mass_file.close();
    
    // Check for missing masses
    if (!beam_found) {
        cout << "ERROR: Beam mass not found in mass.dat!" << endl;
    }
    if (!target_found) {
        cout << "ERROR: Target mass not found in mass.dat!" << endl;
    }
    for (int i = 0; i < product_A.size(); i++) {
        if (!products_found[i]) {
            cout << "ERROR: Product mass not found in mass.dat: " << product_names[i] << " (" << product_A[i] << ", Z=" << product_Z[i] << ")" << endl;
        }
    }
    for (int i = 0; i < decay_A.size(); i++) {
        if (!decay_found[i]) {
            cout << "ERROR: Decay product mass not found in mass.dat: " << decay_names[i] << " (" << decay_A[i] << ", Z=" << decay_Z[i] << ")" << endl;
        }
    }
    if (enable_product_reconstruction && !parent_found) {
        cout << "ERROR: Parent particle mass not found in mass.dat: " << parent_name << " (" << parent_A << ", Z=" << parent_Z << ")" << endl;
    }
    
    // Check conservation laws after all masses are loaded
    CheckConservation();
}

// Initialize all histograms
void FusionReaction::InitializeHistograms() {
    his_beam_E = new TH1D("his_beam_E", "Beam Energy", 4000, 0, 400);
    his_beam_pos = new TH2F("his_beam_pos", "Beam Position", 100, -10, 10, 100, -10, 10);
    his_beam_pos->SetOption("COL");
    
    his_multi_momentum = new TH2F("his_multi_momentum", "Multi-particle Momentum", 
                                 100, -1000, 1000, 100, -1000, 1000);
    his_multi_momentum->SetOption("COL");
    
    // Initialize energy reconstruction histograms (if enabled)
    if (enable_total_energy_reconstruction) {
        his_total_energy_initial = new TH1D("his_total_energy_initial", "Total Initial Energy", 1000, 0, 100);
        his_total_energy_final = new TH1D("his_total_energy_final", "Total Final Energy", 1000, 0, 100);
        his_energy_difference = new TH1D("his_energy_difference", "Energy Difference (Initial - Final)", 1000, -10, 10);
        his_total_momentum_mag = new TH1D("his_total_momentum_mag", "Total Momentum Magnitude in CM frame", 100, -10, 10);
    } else {
        his_total_energy_initial = nullptr;
        his_total_energy_final = nullptr;
        his_energy_difference = nullptr;
        his_total_momentum_mag = nullptr;
    }
    
    // Create histograms for each product
    his_product_angle.resize(products.size());
    his_product_energy.resize(products.size());
    his_product_Evsang.resize(products.size());
    his_product_theta_E_lab.resize(products.size());
    
    for (int i = 0; i < products.size(); i++) {
        char name[100], title[100];
        sprintf(name, "his_product_%d_angle", i);
        sprintf(title, "%s Angle", product_names[i].c_str());
        his_product_angle[i] = new TH1D(name, title, 1800, 0, 180);
        
        sprintf(name, "his_product_%d_energy", i);
        sprintf(title, "%s Energy", product_names[i].c_str());
        // Use wide energy range for auto-adjustment
        his_product_energy[i] = new TH1D(name, title, 2000, 0, 500); // Wide range
        
        sprintf(name, "his_product_%d_Evsang", i);
        sprintf(title, "%s E vs Angle (CM)", product_names[i].c_str());
        his_product_Evsang[i] = new TH2F(name, title, 180, 0, 180, 1000, 0, 500); // Wide range
        his_product_Evsang[i]->SetOption("COL");
        
        sprintf(name, "his_product_%d_theta_E_lab", i);
        sprintf(title, "%s Theta vs Energy (Lab)", product_names[i].c_str());
        his_product_theta_E_lab[i] = new TH2F(name, title, 180, 0, 180, 1000, 0, 500); // Wide range
        his_product_theta_E_lab[i]->SetOption("COL");
    }
    
    // Initialize product reconstruction histograms (if enabled)
    if (enable_product_reconstruction) {
        // Calculate appropriate mass range based on parent mass
        double mass_center = parent_mass;
        double mass_range = parent_mass * 0.1;  // ±10% range
        double mass_min = mass_center - mass_range;
        double mass_max = mass_center + mass_range;
        
        // Product 1 histograms
        his_product1_mass_reconstructed = new TH1D("his_product1_mass_reconstructed", 
            "Product 1 Reconstructed Mass", 200, mass_min, mass_max);
        his_product1_mass_actual = new TH1D("his_product1_mass_actual", 
            "Product 1 Actual Mass", 200, mass_min, mass_max);
        his_product1_mass_difference = new TH1D("his_product1_mass_difference", 
            "Product 1 Mass Reconstruction Error", 500, -50, 50);
        his_product1_energy_reconstructed = new TH1D("his_product1_energy_reconstructed", 
            "Product 1 Reconstructed Energy", 100, 0, 100);
        his_product1_energy_actual = new TH1D("his_product1_energy_actual", 
            "Product 1 Actual Energy", 100, 0, 100);
        his_product1_energy_difference = new TH1D("his_product1_energy_difference", 
            "Product 1 Energy Reconstruction Error", 100, -10, 10);
        
        // Product 2 histograms (same range as product 1)
        his_product2_mass_reconstructed = new TH1D("his_product2_mass_reconstructed", 
            "Product 2 Reconstructed Mass", 200, mass_min, mass_max);
        his_product2_mass_actual = new TH1D("his_product2_mass_actual", 
            "Product 2 Actual Mass", 200, mass_min, mass_max);
        his_product2_mass_difference = new TH1D("his_product2_mass_difference", 
            "Product 2 Mass Reconstruction Error", 500, -50, 50);
        his_product2_energy_reconstructed = new TH1D("his_product2_energy_reconstructed", 
            "Product 2 Reconstructed Energy", 100, 0, 100);
        his_product2_energy_actual = new TH1D("his_product2_energy_actual", 
            "Product 2 Actual Energy", 100, 0, 100);
        his_product2_energy_difference = new TH1D("his_product2_energy_difference", 
            "Product 2 Energy Reconstruction Error", 100, -10, 10);
    } else {
        // Set all product reconstruction histograms to nullptr
        his_product1_mass_reconstructed = nullptr;
        his_product1_mass_actual = nullptr;
        his_product1_mass_difference = nullptr;
        his_product1_energy_reconstructed = nullptr;
        his_product1_energy_actual = nullptr;
        his_product1_energy_difference = nullptr;
        
        his_product2_mass_reconstructed = nullptr;
        his_product2_mass_actual = nullptr;
        his_product2_mass_difference = nullptr;
        his_product2_energy_reconstructed = nullptr;
        his_product2_energy_actual = nullptr;
        his_product2_energy_difference = nullptr;
    }
    
    // Automatically initialize decay histograms if decay is enabled
    if (decay_enabled) {
        InitializeDecayHistograms();
    }
}

// Initialize decay histograms
void FusionReaction::InitializeDecayHistograms() {
    if (!decay_enabled) return;
    
    // Create histograms for each decay product with wide initial range
    his_decay_angle.resize(decay_A.size());
    his_decay_energy.resize(decay_A.size());
    his_decay_Evsang.resize(decay_A.size());
    his_decay_theta_E_lab.resize(decay_A.size());
    
    for (int i = 0; i < decay_A.size(); i++) {
        char name[100], title[100];
        sprintf(name, "his_decay_%d_angle", i);
        sprintf(title, "%s Decay Angle", decay_names[i].c_str());
        his_decay_angle[i] = new TH1D(name, title, 1800, 0, 180);
        
        sprintf(name, "his_decay_%d_energy", i);
        sprintf(title, "%s Decay Energy", decay_names[i].c_str());
        his_decay_energy[i] = new TH1D(name, title, 2000, 0, 500); // Very wide initial range
        
        sprintf(name, "his_decay_%d_Evsang", i);
        sprintf(title, "%s Decay E vs Angle (CM)", decay_names[i].c_str());
        his_decay_Evsang[i] = new TH2F(name, title, 180, 0, 180, 1000, 0, 500); // Very wide initial range
        his_decay_Evsang[i]->SetOption("COL");
        
        sprintf(name, "his_decay_%d_theta_E_lab", i);
        sprintf(title, "%s Decay Theta vs Energy (Lab)", decay_names[i].c_str());
        his_decay_theta_E_lab[i] = new TH2F(name, title, 180, 0, 180, 1000, 0, 500); // Very wide initial range
        his_decay_theta_E_lab[i]->SetOption("COL");
    }
    
    // Initialize parent particle reconstruction histograms (if enabled)
    if (enable_energy_reconstruction) {
        his_parent_energy_reconstructed = new TH1D("his_parent_energy_reconstructed", 
            "Reconstructed Parent Particle Energy", 2000, 0, 500); // Wide range
        his_parent_energy_actual = new TH1D("his_parent_energy_actual", 
            "Actual Parent Particle Energy", 2000, 0, 500); // Wide range
        his_parent_energy_difference = new TH1D("his_parent_energy_difference", 
            "Parent Energy Reconstruction Error", 100, -10, 10);
    } else {
        his_parent_energy_reconstructed = nullptr;
        his_parent_energy_actual = nullptr;
        his_parent_energy_difference = nullptr;
    }
    
    // Initialize parent particle mass reconstruction histograms (if enabled)
    if (enable_mass_reconstruction) {
        // Get parent particle mass for auto-range
        double parent_particle_mass = 0.0;
        if (decay_product_index >= 0 && decay_product_index < products.size()) {
            parent_particle_mass = products[decay_product_index].mass;
        }
        
        // Calculate appropriate mass range (±10% of parent mass)
        double mass_range = parent_particle_mass * 0.1;
        double mass_min = parent_particle_mass - mass_range;
        double mass_max = parent_particle_mass + mass_range;
        
        his_parent_mass_reconstructed = new TH1D("his_parent_mass_reconstructed", 
            "Reconstructed Parent Particle Mass", 200, mass_min, mass_max);
        his_parent_mass_actual = new TH1D("his_parent_mass_actual", 
            "Actual Parent Particle Mass", 200, mass_min, mass_max);
        his_parent_mass_difference = new TH1D("his_parent_mass_difference", 
            "Parent Mass Reconstruction Error", 500, -50, 50);
    } else {
        his_parent_mass_reconstructed = nullptr;
        his_parent_mass_actual = nullptr;
        his_parent_mass_difference = nullptr;
    }
    
    // Initialize product reconstruction histograms (if enabled)
    if (enable_product_reconstruction) {
        // Product 1 histograms
        his_product1_mass_reconstructed = new TH1D("his_product1_mass_reconstructed", 
            "Product 1 Reconstructed Mass", 100, 0, 1000);
        his_product1_mass_actual = new TH1D("his_product1_mass_actual", 
            "Product 1 Actual Mass", 100, 0, 1000);
        his_product1_mass_difference = new TH1D("his_product1_mass_difference", 
            "Product 1 Mass Reconstruction Error", 100, -50, 50);
        his_product1_energy_reconstructed = new TH1D("his_product1_energy_reconstructed", 
            "Product 1 Reconstructed Energy", 100, 0, 100);
        his_product1_energy_actual = new TH1D("his_product1_energy_actual", 
            "Product 1 Actual Energy", 100, 0, 100);
        his_product1_energy_difference = new TH1D("his_product1_energy_difference", 
            "Product 1 Energy Reconstruction Error", 100, -10, 10);
        
        // Product 2 histograms
        his_product2_mass_reconstructed = new TH1D("his_product2_mass_reconstructed", 
            "Product 2 Reconstructed Mass", 100, 0, 1000);
        his_product2_mass_actual = new TH1D("his_product2_mass_actual", 
            "Product 2 Actual Mass", 100, 0, 1000);
        his_product2_mass_difference = new TH1D("his_product2_mass_difference", 
            "Product 2 Mass Reconstruction Error", 100, -50, 50);
        his_product2_energy_reconstructed = new TH1D("his_product2_energy_reconstructed", 
            "Product 2 Reconstructed Energy", 100, 0, 100);
        his_product2_energy_actual = new TH1D("his_product2_energy_actual", 
            "Product 2 Actual Energy", 100, 0, 100);
        his_product2_energy_difference = new TH1D("his_product2_energy_difference", 
            "Product 2 Energy Reconstruction Error", 100, -10, 10);
    } else {
        // Set all product reconstruction histograms to nullptr
        his_product1_mass_reconstructed = nullptr;
        his_product1_mass_actual = nullptr;
        his_product1_mass_difference = nullptr;
        his_product1_energy_reconstructed = nullptr;
        his_product1_energy_actual = nullptr;
        his_product1_energy_difference = nullptr;
        
        his_product2_mass_reconstructed = nullptr;
        his_product2_mass_actual = nullptr;
        his_product2_mass_difference = nullptr;
        his_product2_energy_reconstructed = nullptr;
        his_product2_energy_actual = nullptr;
        his_product2_energy_difference = nullptr;
    }
}

// Draw all results on canvases
void FusionReaction::DrawResults() {
    // Auto-adjust histogram ranges based on actual data
    AutoAdjustHistogramRanges();
    
    // Canvas 1: Beam properties and momentum
    TCanvas* c1 = new TCanvas("c1", "Beam and Momentum", 1200, 400);
    c1->Divide(3, 1);
    
    c1->cd(1);
    his_beam_E->Draw();
    
    c1->cd(2);
    his_beam_pos->Draw("COL");
    
    c1->cd(3);
    his_multi_momentum->Draw("COL");
    
    // Canvas 2: Energy reconstruction (if enabled)
    if (enable_total_energy_reconstruction) {
        TCanvas* c2 = new TCanvas("c2", "Energy Reconstruction", 1200, 800);
        c2->Divide(2, 2);
        
        c2->cd(1);
        his_total_energy_initial->Draw();
        
        c2->cd(2);
        his_total_energy_final->Draw();
        
        c2->cd(3);
        his_energy_difference->Draw();
        
        c2->cd(4);
        his_total_momentum_mag->Draw();
    }
    
    // Canvas 3+: Lab Frame Distributions for all particles
    int n_products = his_product_theta_E_lab.size();
    if (n_products > 0) {
        // Calculate number of canvases needed (4 plots per canvas)
        int n_canvases = (n_products + 3) / 4;  // Round up
        
        for (int canvas_idx = 0; canvas_idx < n_canvases; canvas_idx++) {
            char canvas_name[50], canvas_title[100];
            sprintf(canvas_name, "c%d", canvas_idx + 3);
            sprintf(canvas_title, "Lab Frame Theta vs Energy (Products %d-%d)", 
                    canvas_idx * 4 + 1, min((canvas_idx + 1) * 4, n_products));
            
            TCanvas* c = new TCanvas(canvas_name, canvas_title, 1200, 800);
            
            int start_idx = canvas_idx * 4;
            int end_idx = min(start_idx + 4, n_products);
            int n_plots = end_idx - start_idx;
            
            // Arrange plots in 2x2 grid
            c->Divide(2, 2);
            
            for (int i = 0; i < n_plots; i++) {
                c->cd(i + 1);
                gPad->SetRightMargin(0.15);
                his_product_theta_E_lab[start_idx + i]->Draw("COLZ");
            }
        }
    }
    
    // Draw decay histograms if decay is enabled
    if (decay_enabled && his_decay_theta_E_lab.size() > 0) {
        int n_decay_products = his_decay_theta_E_lab.size();
        int n_decay_canvases = (n_decay_products + 3) / 4;  // Round up
        
        for (int canvas_idx = 0; canvas_idx < n_decay_canvases; canvas_idx++) {
            char canvas_name[50], canvas_title[100];
            sprintf(canvas_name, "c_decay_%d", canvas_idx + 1);
            sprintf(canvas_title, "Decay Products Theta vs Energy (Products %d-%d)", 
                    canvas_idx * 4 + 1, min((canvas_idx + 1) * 4, n_decay_products));
            
            TCanvas* c = new TCanvas(canvas_name, canvas_title, 1200, 800);
            
            int start_idx = canvas_idx * 4;
            int end_idx = min(start_idx + 4, n_decay_products);
            int n_plots = end_idx - start_idx;
            
            // Arrange plots in 2x2 grid
            c->Divide(2, 2);
            
            for (int i = 0; i < n_plots; i++) {
                c->cd(i + 1);
                gPad->SetRightMargin(0.15);
                his_decay_theta_E_lab[start_idx + i]->Draw("COLZ");
            }
        }
    }
    
    // Draw parent reconstruction results if decay and reconstruction are enabled
    if (decay_enabled && enable_energy_reconstruction && his_parent_energy_reconstructed) {
        TCanvas* c_parent = new TCanvas("c_parent", "Parent Energy Reconstruction", 1200, 400);
        c_parent->Divide(3, 1);
        
        c_parent->cd(1);
        his_parent_energy_reconstructed->SetLineColor(kBlack);
        his_parent_energy_reconstructed->SetLineWidth(2);
        his_parent_energy_reconstructed->GetXaxis()->SetTitle("Energy (MeV)");
        his_parent_energy_reconstructed->GetYaxis()->SetTitle("Counts");
        his_parent_energy_reconstructed->Draw();
        
        c_parent->cd(2);
        his_parent_energy_actual->SetLineColor(kBlack);
        his_parent_energy_actual->SetLineWidth(2);
        his_parent_energy_actual->GetXaxis()->SetTitle("Energy (MeV)");
        his_parent_energy_actual->GetYaxis()->SetTitle("Counts");
        his_parent_energy_actual->Draw();
        
        c_parent->cd(3);
        his_parent_energy_difference->SetLineColor(kBlack);
        his_parent_energy_difference->SetLineWidth(2);
        his_parent_energy_difference->GetXaxis()->SetTitle("Energy Difference (MeV)");
        his_parent_energy_difference->GetYaxis()->SetTitle("Counts");
        his_parent_energy_difference->Draw();
        
        // Add legend to first plot
        c_parent->cd(1);
        TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
        legend->AddEntry(his_parent_energy_reconstructed, "Reconstructed", "l");
        legend->AddEntry(his_parent_energy_actual, "Actual", "l");
        legend->Draw();
    }
    
    // Draw parent mass reconstruction results if decay and mass reconstruction are enabled
    if (decay_enabled && enable_mass_reconstruction && his_parent_mass_reconstructed) {
        TCanvas* c_mass = new TCanvas("c_mass", "Parent Mass Reconstruction", 1200, 400);
        c_mass->Divide(3, 1);
        
        c_mass->cd(1);
        his_parent_mass_reconstructed->SetLineColor(kBlack);
        his_parent_mass_reconstructed->SetLineWidth(2);
        his_parent_mass_reconstructed->GetXaxis()->SetTitle("Mass (MeV)");
        his_parent_mass_reconstructed->GetYaxis()->SetTitle("Counts");
        his_parent_mass_reconstructed->Draw();
        
        c_mass->cd(2);
        his_parent_mass_actual->SetLineColor(kBlack);
        his_parent_mass_actual->SetLineWidth(2);
        his_parent_mass_actual->GetXaxis()->SetTitle("Mass (MeV)");
        his_parent_mass_actual->GetYaxis()->SetTitle("Counts");
        his_parent_mass_actual->Draw();
        
        c_mass->cd(3);
        his_parent_mass_difference->SetLineColor(kBlack);
        his_parent_mass_difference->SetLineWidth(2);
        his_parent_mass_difference->GetXaxis()->SetTitle("Mass Difference (MeV)");
        his_parent_mass_difference->GetYaxis()->SetTitle("Counts");
        his_parent_mass_difference->Draw();
        
        // Add legend to first plot
        c_mass->cd(1);
        TLegend* legend_mass = new TLegend(0.7, 0.7, 0.9, 0.9);
        legend_mass->AddEntry(his_parent_mass_reconstructed, "Reconstructed", "l");
        legend_mass->AddEntry(his_parent_mass_actual, "Actual", "l");
        legend_mass->Draw();
    }
    
    // Draw product reconstruction results if enabled
    if (enable_product_reconstruction && his_product1_mass_reconstructed) {
        TCanvas* c_products = new TCanvas("c_products", "Product Reconstruction", 1200, 400);
        c_products->Divide(3, 1);
        
        // Product 1 plots
        c_products->cd(1);
        his_product1_mass_reconstructed->SetLineColor(kBlack);
        his_product1_mass_reconstructed->SetLineWidth(2);
        his_product1_mass_reconstructed->GetXaxis()->SetTitle("Mass (MeV)");
        his_product1_mass_reconstructed->GetYaxis()->SetTitle("Counts");
        his_product1_mass_reconstructed->SetTitle(("Parent Mass from " + selected_product1_name + " + " + selected_product2_name).c_str());
        his_product1_mass_reconstructed->Draw();
        
        c_products->cd(2);
        his_product1_mass_actual->SetLineColor(kBlack);
        his_product1_mass_actual->SetLineWidth(2);
        his_product1_mass_actual->GetXaxis()->SetTitle("Mass (MeV)");
        his_product1_mass_actual->GetYaxis()->SetTitle("Counts");
        his_product1_mass_actual->SetTitle(("Actual Sum Mass (" + selected_product1_name + " + " + selected_product2_name + ")").c_str());
        his_product1_mass_actual->Draw();
        
        c_products->cd(3);
        his_product1_mass_difference->SetLineColor(kBlack);
        his_product1_mass_difference->SetLineWidth(2);
        his_product1_mass_difference->GetXaxis()->SetTitle("Mass Difference (MeV)");
        his_product1_mass_difference->GetYaxis()->SetTitle("Counts");
        his_product1_mass_difference->SetTitle(("Mass Reconstruction Error (" + selected_product1_name + " + " + selected_product2_name + ")").c_str());
        his_product1_mass_difference->Draw();
    }
}

// Auto-adjust histogram ranges based on actual data
void FusionReaction::AutoAdjustHistogramRanges() {
    // Adjust fusion product histograms
    for (int i = 0; i < his_product_energy.size(); i++) {
        if (his_product_energy[i] && his_product_energy[i]->GetEntries() > 0) {
            // Find the actual maximum value in the data by scanning all bins
            double actual_max = 0.0;
            int n_bins = his_product_energy[i]->GetNbinsX();
            
            for (int bin = 1; bin <= n_bins; bin++) {
                if (his_product_energy[i]->GetBinContent(bin) > 0) {
                    double bin_center = his_product_energy[i]->GetXaxis()->GetBinCenter(bin);
                    actual_max = TMath::Max(actual_max, bin_center);
                }
            }
            
            // Also check 2D histograms for maximum energy
            if (his_product_theta_E_lab[i] && his_product_theta_E_lab[i]->GetEntries() > 0) {
                int n_y_bins = his_product_theta_E_lab[i]->GetNbinsY();
                for (int y_bin = 1; y_bin <= n_y_bins; y_bin++) {
                    for (int x_bin = 1; x_bin <= his_product_theta_E_lab[i]->GetNbinsX(); x_bin++) {
                        if (his_product_theta_E_lab[i]->GetBinContent(x_bin, y_bin) > 0) {
                            double y_center = his_product_theta_E_lab[i]->GetYaxis()->GetBinCenter(y_bin);
                            actual_max = TMath::Max(actual_max, y_center);
                        }
                    }
                }
            }
            
            if (actual_max > 0) {
                // Set new range with generous margin (50% more)
                double new_max = actual_max * 1.5;
                double min_val = his_product_energy[i]->GetXaxis()->GetXmin();
                
                cout << "Adjusting fusion product " << i << " (" << product_names[i] << ") range: " 
                     << min_val << " - " << new_max << " MeV (actual max: " << actual_max << " MeV)" << endl;
                
                his_product_energy[i]->GetXaxis()->SetRangeUser(min_val, new_max);
                
                // Also adjust 2D histograms
                if (his_product_Evsang[i]) {
                    his_product_Evsang[i]->GetYaxis()->SetRangeUser(0, new_max);
                }
                if (his_product_theta_E_lab[i]) {
                    his_product_theta_E_lab[i]->GetYaxis()->SetRangeUser(0, new_max);
                }
            }
        }
    }
    
    // Adjust decay product histograms
    if (decay_enabled) {
        for (int i = 0; i < his_decay_energy.size(); i++) {
            if (his_decay_energy[i] && his_decay_energy[i]->GetEntries() > 0) {
                // Find the actual maximum value in the data by scanning all bins
                double actual_max = 0.0;
                int n_bins = his_decay_energy[i]->GetNbinsX();
                
                for (int bin = 1; bin <= n_bins; bin++) {
                    if (his_decay_energy[i]->GetBinContent(bin) > 0) {
                        double bin_center = his_decay_energy[i]->GetXaxis()->GetBinCenter(bin);
                        actual_max = TMath::Max(actual_max, bin_center);
                    }
                }
                
                // Also check 2D histograms for maximum energy
                if (his_decay_theta_E_lab[i] && his_decay_theta_E_lab[i]->GetEntries() > 0) {
                    int n_y_bins = his_decay_theta_E_lab[i]->GetNbinsY();
                    for (int y_bin = 1; y_bin <= n_y_bins; y_bin++) {
                        for (int x_bin = 1; x_bin <= his_decay_theta_E_lab[i]->GetNbinsX(); x_bin++) {
                            if (his_decay_theta_E_lab[i]->GetBinContent(x_bin, y_bin) > 0) {
                                double y_center = his_decay_theta_E_lab[i]->GetYaxis()->GetBinCenter(y_bin);
                                actual_max = TMath::Max(actual_max, y_center);
                            }
                        }
                    }
                }
                
                if (actual_max > 0) {
                    // Set new range with generous margin (50% more)
                    double new_max = actual_max * 1.5;
                    double min_val = his_decay_energy[i]->GetXaxis()->GetXmin();
                    
                    cout << "Adjusting decay histogram " << i << " range: " << min_val << " - " << new_max 
                         << " MeV (actual max: " << actual_max << " MeV)" << endl;
                    
                    his_decay_energy[i]->GetXaxis()->SetRangeUser(min_val, new_max);
                    
                    // Also adjust 2D histograms
                    if (his_decay_Evsang[i]) {
                        his_decay_Evsang[i]->GetYaxis()->SetRangeUser(0, new_max);
                    }
                    if (his_decay_theta_E_lab[i]) {
                        his_decay_theta_E_lab[i]->GetYaxis()->SetRangeUser(0, new_max);
                    }
                }
            }
        }
    }
    
    // Adjust parent particle (heavy recoil) energy histograms
    if (decay_enabled && enable_energy_reconstruction) {
        // Adjust parent energy reconstructed histogram
        if (his_parent_energy_reconstructed && his_parent_energy_reconstructed->GetEntries() > 0) {
            double actual_max = 0.0;
            int n_bins = his_parent_energy_reconstructed->GetNbinsX();
            
            for (int bin = 1; bin <= n_bins; bin++) {
                if (his_parent_energy_reconstructed->GetBinContent(bin) > 0) {
                    double bin_center = his_parent_energy_reconstructed->GetXaxis()->GetBinCenter(bin);
                    actual_max = TMath::Max(actual_max, bin_center);
                }
            }
            
            if (actual_max > 0) {
                double new_max = actual_max * 1.5;
                double min_val = his_parent_energy_reconstructed->GetXaxis()->GetXmin();
                
                cout << "Adjusting parent energy reconstructed range: " << min_val << " - " << new_max 
                     << " MeV (actual max: " << actual_max << " MeV)" << endl;
                
                his_parent_energy_reconstructed->GetXaxis()->SetRangeUser(min_val, new_max);
            }
        }
        
        // Adjust parent energy actual histogram
        if (his_parent_energy_actual && his_parent_energy_actual->GetEntries() > 0) {
            double actual_max = 0.0;
            int n_bins = his_parent_energy_actual->GetNbinsX();
            
            for (int bin = 1; bin <= n_bins; bin++) {
                if (his_parent_energy_actual->GetBinContent(bin) > 0) {
                    double bin_center = his_parent_energy_actual->GetXaxis()->GetBinCenter(bin);
                    actual_max = TMath::Max(actual_max, bin_center);
                }
            }
            
            if (actual_max > 0) {
                double new_max = actual_max * 1.5;
                double min_val = his_parent_energy_actual->GetXaxis()->GetXmin();
                
                cout << "Adjusting parent energy actual range: " << min_val << " - " << new_max 
                     << " MeV (actual max: " << actual_max << " MeV)" << endl;
                
                his_parent_energy_actual->GetXaxis()->SetRangeUser(min_val, new_max);
            }
        }
    }
    
    // Adjust parent particle (heavy recoil) mass histograms
    if (decay_enabled && enable_mass_reconstruction) {
        // Adjust parent mass reconstructed histogram
        if (his_parent_mass_reconstructed && his_parent_mass_reconstructed->GetEntries() > 0) {
            double actual_max = 0.0;
            double actual_min = 1e6;
            int n_bins = his_parent_mass_reconstructed->GetNbinsX();
            
            for (int bin = 1; bin <= n_bins; bin++) {
                if (his_parent_mass_reconstructed->GetBinContent(bin) > 0) {
                    double bin_center = his_parent_mass_reconstructed->GetXaxis()->GetBinCenter(bin);
                    actual_max = TMath::Max(actual_max, bin_center);
                    actual_min = TMath::Min(actual_min, bin_center);
                }
            }
            
            if (actual_max > actual_min) {
                double range = actual_max - actual_min;
                double new_min = actual_min - range * 0.2;
                double new_max = actual_max + range * 0.2;
                
                cout << "Adjusting parent mass reconstructed range: " << new_min << " - " << new_max 
                     << " MeV (actual range: " << actual_min << " - " << actual_max << " MeV)" << endl;
                
                his_parent_mass_reconstructed->GetXaxis()->SetRangeUser(new_min, new_max);
            }
        }
        
        // Adjust parent mass actual histogram
        if (his_parent_mass_actual && his_parent_mass_actual->GetEntries() > 0) {
            double actual_max = 0.0;
            double actual_min = 1e6;
            int n_bins = his_parent_mass_actual->GetNbinsX();
            
            for (int bin = 1; bin <= n_bins; bin++) {
                if (his_parent_mass_actual->GetBinContent(bin) > 0) {
                    double bin_center = his_parent_mass_actual->GetXaxis()->GetBinCenter(bin);
                    actual_max = TMath::Max(actual_max, bin_center);
                    actual_min = TMath::Min(actual_min, bin_center);
                }
            }
            
            if (actual_max > actual_min) {
                double range = actual_max - actual_min;
                double new_min = actual_min - range * 0.2;
                double new_max = actual_max + range * 0.2;
                
                cout << "Adjusting parent mass actual range: " << new_min << " - " << new_max 
                     << " MeV (actual range: " << actual_min << " - " << actual_max << " MeV)" << endl;
                
                his_parent_mass_actual->GetXaxis()->SetRangeUser(new_min, new_max);
            }
        }
    }
}
