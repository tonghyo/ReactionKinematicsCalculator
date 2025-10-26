#include "FusionReaction.h"

// Check energy conservation
bool FusionReaction::CheckEnergyConservation() {
    // Use Invariant mass for both initial and final
    // This is frame-independent and should be conserved
    double E_beam_used = E_beam_current;
    
    // Calculate Initial Invariant mass (CM total energy)
    double E_total_lab = E_beam_used + M_beam + M_target;
    double p_beam_lab = sqrt(E_beam_used * (E_beam_used + 2 * M_beam));
    double s = E_total_lab * E_total_lab - p_beam_lab * p_beam_lab;
    double W_initial = sqrt(s);
    
    // Calculate Final Invariant mass (CM total energy)
    // Sum all 4-momenta in Lab frame, then calculate invariant mass
    double E_total_final = 0;
    double px_total_final = 0;
    double py_total_final = 0;
    double pz_total_final = 0;
    
    for (int i = 0; i < products.size(); i++) {
        double E_product = products[i].energy + products[i].mass; // Total energy in Lab frame (already calculated)
        
        // Use Lab frame momentum components (already calculated in CalculateProductKinematics)
        double px = products[i].px;
        double py = products[i].py;
        double pz = products[i].pz;
        
        E_total_final += E_product;
        px_total_final += px;
        py_total_final += py;
        pz_total_final += pz;
    }
    
    // Calculate invariant mass of final state
    double s_final = E_total_final * E_total_final - 
                     (px_total_final * px_total_final + 
                      py_total_final * py_total_final + 
                      pz_total_final * pz_total_final);
    double W_final = sqrt(s_final);
    
    double energy_diff = abs(W_initial - W_final);
    double tolerance = 5.0; // MeV (allow for numerical precision)
    
    if (energy_diff > tolerance) {
        cout << "WARNING: Energy conservation violated!" << endl;
        cout << "Initial W: " << W_initial << " MeV, Final W: " << W_final << " MeV" << endl;
        cout << "Difference: " << energy_diff << " MeV" << endl;
        return false;
    }
    return true;
}

// Reconstruct energy for analysis
void FusionReaction::ReconstructEnergy() {
    // Use Invariant mass (frame-independent)
    double E_beam_used = E_beam_current;
    
    // Calculate Invariant mass (CM total energy)
    double E_total_lab = E_beam_used + M_beam + M_target;
    double p_beam_lab = sqrt(E_beam_used * (E_beam_used + 2 * M_beam));
    double s = E_total_lab * E_total_lab - p_beam_lab * p_beam_lab;
    double E_total_initial = sqrt(s);  // Invariant mass
    
    // Final state: calculate invariant mass from all product 4-momenta
    double E_total_final = 0;
    double total_px = 0, total_py = 0, total_pz = 0;
    
    for (int i = 0; i < products.size(); i++) {
        double product_total_E = products[i].energy + products[i].mass;
        E_total_final += product_total_E;
        total_px += products[i].px;
        total_py += products[i].py;
        total_pz += products[i].pz;
    }
    
    // Calculate invariant mass of final state
    double s_final = E_total_final * E_total_final - 
                     (total_px * total_px + total_py * total_py + total_pz * total_pz);
    double E_total_final_invariant = sqrt(s_final);
    
    // Calculate total momentum magnitude
    double total_momentum = sqrt(total_px * total_px + total_py * total_py + total_pz * total_pz);
    
    // Calculate energy difference (invariant masses)
    double energy_diff = E_total_initial - E_total_final_invariant;
    
    // Fill histograms
    his_total_energy_initial->Fill(E_total_initial / 1000.0);  // Convert to GeV for display
    his_total_energy_final->Fill(E_total_final_invariant / 1000.0);  // Convert to GeV for display
    his_energy_difference->Fill(energy_diff);                   // Keep in MeV
    his_total_momentum_mag->Fill(total_momentum / 1000.0);     // Convert to GeV/c for display
}

void FusionReaction::ReconstructParentEnergy() {
    if (!decay_enabled || decay_momenta.empty()) return;
    
    // Reconstruct parent particle energy from decay products
    // Method: Sum all decay product 4-momenta to get parent 4-momentum
    
    double parent_px = 0.0, parent_py = 0.0, parent_pz = 0.0;
    double parent_E_total = 0.0;
    
    // Sum decay product momenta and energies
    for (int i = 0; i < decay_momenta.size(); i++) {
        // Get decay product momentum magnitude
        double p_decay = decay_momenta[i];
        
        // Get decay product energy (kinetic + rest mass)
        double E_decay_kinetic = decay_energies[i];
        double E_decay_total = E_decay_kinetic + decay_masses[i];
        
        // Get decay product angles (use CM frame angles for momentum direction)
        double theta_decay = decay_angles[i] * TMath::Pi() / 180.0;
        double phi_decay = 0.0; // Assume phi = 0 for simplicity
        
        // Calculate momentum components
        double px_decay = p_decay * sin(theta_decay) * cos(phi_decay);
        double py_decay = p_decay * sin(theta_decay) * sin(phi_decay);
        double pz_decay = p_decay * cos(theta_decay);
        
        // Sum to get parent momentum
        parent_px += px_decay;
        parent_py += py_decay;
        parent_pz += pz_decay;
        parent_E_total += E_decay_total;
    }
    
    // Calculate parent particle momentum magnitude
    double parent_p_mag = sqrt(parent_px*parent_px + parent_py*parent_py + parent_pz*parent_pz);
    
    // Calculate parent particle invariant mass (rest mass)
    double parent_invariant_mass = sqrt(parent_E_total*parent_E_total - parent_p_mag*parent_p_mag);
    
    // Calculate parent particle kinetic energy
    double parent_kinetic_energy = parent_E_total - parent_invariant_mass;
    
    // Get actual parent particle energy for comparison (original energy before decay)
    double actual_parent_energy = original_parent_energy;
    
    // Fill reconstruction histograms
    his_parent_energy_reconstructed->Fill(parent_kinetic_energy);
    his_parent_energy_actual->Fill(actual_parent_energy);
    his_parent_energy_difference->Fill(actual_parent_energy - parent_kinetic_energy);
}

void FusionReaction::ReconstructParentMass() {
    if (!decay_enabled || decay_momenta.empty()) return;
    
    // Reconstruct parent particle mass from decay products
    // Method: Calculate invariant mass from decay product 4-momenta
    
    double parent_px = 0.0, parent_py = 0.0, parent_pz = 0.0;
    double parent_E_total = 0.0;
    
    // Sum decay product momenta and energies
    for (int i = 0; i < decay_momenta.size(); i++) {
        // Get decay product momentum magnitude
        double p_decay = decay_momenta[i];
        
        // Get decay product energy (kinetic + rest mass)
        double E_decay_kinetic = decay_energies[i];
        double E_decay_total = E_decay_kinetic + decay_masses[i];
        
        // Get decay product angles (use CM frame angles for momentum direction)
        double theta_decay = decay_angles[i] * TMath::Pi() / 180.0;
        double phi_decay = 0.0; // Assume phi = 0 for simplicity
        
        // Calculate momentum components
        double px_decay = p_decay * sin(theta_decay) * cos(phi_decay);
        double py_decay = p_decay * sin(theta_decay) * sin(phi_decay);
        double pz_decay = p_decay * cos(theta_decay);
        
        // Sum to get parent momentum
        parent_px += px_decay;
        parent_py += py_decay;
        parent_pz += pz_decay;
        parent_E_total += E_decay_total;
    }
    
    // Calculate parent particle momentum magnitude
    double parent_p_mag = sqrt(parent_px*parent_px + parent_py*parent_py + parent_pz*parent_pz);
    
    // Calculate parent particle invariant mass (rest mass)
    double parent_mass_reconstructed = sqrt(parent_E_total*parent_E_total - parent_p_mag*parent_p_mag);
    
    // Get actual parent particle mass for comparison
    double actual_parent_mass = 0.0;
    if (decay_product_index < products.size()) {
        actual_parent_mass = products[decay_product_index].mass;
    }
    
    // Fill mass reconstruction histograms
    his_parent_mass_reconstructed->Fill(parent_mass_reconstructed);
    his_parent_mass_actual->Fill(actual_parent_mass);
    his_parent_mass_difference->Fill(actual_parent_mass - parent_mass_reconstructed);
}

void FusionReaction::ReconstructProductProperties() {
    if (!enable_product_reconstruction || selected_product1 < 0 || selected_product2 < 0) return;
    
    // Check if histograms are initialized
    if (!his_product1_mass_reconstructed || !his_product2_mass_reconstructed) return;
    
    // Reconstruct parent particle from two selected products (like decay reconstruction)
    // Method: Sum 4-momenta of two products to get parent particle properties
    
    if (selected_product1 >= products.size() || selected_product2 >= products.size()) return;
    
    Particle& p1 = products[selected_product1];
    Particle& p2 = products[selected_product2];
    
    // Get measured values (with experimental resolution) - use Lab frame
    double p1_measured = p1.momentum_lab;
    double E1_measured = p1.energy_lab;
    double p2_measured = p2.momentum_lab;
    double E2_measured = p2.energy_lab;
    
    // Get angles for momentum components - use Lab frame angles
    double theta1 = p1.theta_lab;
    double phi1 = p1.phi;
    double theta2 = p2.theta_lab;
    double phi2 = p2.phi;
    
    // Calculate momentum components for both products (Lab frame)
    double p1x = p1.px_lab;  // Use pre-calculated Lab frame momentum components
    double p1y = p1.py_lab;
    double p1z = p1.pz_lab;
    
    double p2x = p2.px_lab;
    double p2y = p2.py_lab;
    double p2z = p2.pz_lab;
    
    // Sum 4-momenta to get parent particle
    double parent_px = p1x + p2x;
    double parent_py = p1y + p2y;
    double parent_pz = p1z + p2z;
    double parent_p_mag = sqrt(parent_px*parent_px + parent_py*parent_py + parent_pz*parent_pz);
    
    // Total energy = kinetic energy + rest mass for each product
    double E1_total = E1_measured + p1.mass;
    double E2_total = E2_measured + p2.mass;
    double parent_E_total = E1_total + E2_total;
    
    // Calculate parent particle invariant mass (rest mass)
    double parent_mass_reconstructed = sqrt(parent_E_total*parent_E_total - parent_p_mag*parent_p_mag);
    
    // Calculate parent particle kinetic energy
    double parent_kinetic_energy_reconstructed = parent_E_total - parent_mass_reconstructed;
    
    // Get actual parent particle properties (from mass file)
    double parent_mass_actual = parent_mass;  // Actual parent mass from mass file
    double parent_kinetic_energy_actual = p1.energy_lab + p2.energy_lab;  // Sum of kinetic energies
    
    // Fill histograms
    his_product1_mass_reconstructed->Fill(parent_mass_reconstructed);
    his_product1_mass_actual->Fill(parent_mass_actual);
    his_product1_mass_difference->Fill(parent_mass_actual - parent_mass_reconstructed);
    
    his_product1_energy_reconstructed->Fill(parent_kinetic_energy_reconstructed);
    his_product1_energy_actual->Fill(parent_kinetic_energy_actual);
    his_product1_energy_difference->Fill(parent_kinetic_energy_actual - parent_kinetic_energy_reconstructed);
    
    // Product 2 histograms are not used in this approach, but keep them for compatibility
    his_product2_mass_reconstructed->Fill(0);  // Not used
    his_product2_mass_actual->Fill(0);         // Not used
    his_product2_mass_difference->Fill(0);     // Not used
    
    his_product2_energy_reconstructed->Fill(0);  // Not used
    his_product2_energy_actual->Fill(0);         // Not used
    his_product2_energy_difference->Fill(0);     // Not used
}

// Print event information
void FusionReaction::PrintEventInfo(int event_num) {
    cout << "\n========== Event " << event_num << " ==========" << endl;
    cout << "Q-value: " << fixed << setprecision(3) << CalculateQValue() << " MeV" << endl;
    cout << "Number of products: " << products.size() << endl;
    
    cout << "\nParticle Information (Lab frame):" << endl;
    cout << "Name\t\tA\tZ\tE_Lab(MeV)\tTheta_Lab(deg)\tP_Lab(MeV/c)" << endl;
    cout << "------------------------------------------------------------------------" << endl;
    
    for (int i = 0; i < products.size(); i++) {
        cout << products[i].name << "\t\t" 
             << products[i].A << "\t" 
             << products[i].Z << "\t"
             << fixed << setprecision(3) << products[i].energy << "\t\t"
             << fixed << setprecision(1) << products[i].theta * 180.0 / TMath::Pi() << "\t\t"
             << fixed << setprecision(3) << products[i].momentum << endl;
    }
}

// Print decay information
void FusionReaction::PrintDecayInfo(int event_num) {
    if (!decay_enabled) return;
    
    cout << "\n========== Decay Event " << event_num << " ==========" << endl;
    
    // Get the parent particle that will decay
    Particle& parent = products[decay_product_index];
    cout << "Parent particle: " << parent.name << " (A=" << parent.A << ", Z=" << parent.Z << ")" << endl;
    cout << "Parent energy (Lab): " << fixed << setprecision(3) << parent.energy << " MeV" << endl;
    cout << "Parent momentum (Lab): " << fixed << setprecision(3) << parent.momentum << " MeV/c" << endl;
    cout << "Parent angle (Lab): " << fixed << setprecision(1) << parent.theta * 180.0 / TMath::Pi() << " deg" << endl;
    
    // Calculate decay Q-value
    double Q_decay = parent.mass;
    for (int i = 0; i < decay_A.size(); i++) {
        Q_decay -= decay_masses[i];
    }
    cout << "Decay Q-value: " << fixed << setprecision(3) << Q_decay << " MeV" << endl;
    
    cout << "\nDecay products (Lab frame):" << endl;
    cout << "Name\t\tA\tZ\tMass(MeV)\tE_Lab(MeV)\tP_Lab(MeV/c)\tTheta_Lab(deg)" << endl;
    cout << "------------------------------------------------------------------------" << endl;
    
    for (int i = 0; i < decay_A.size(); i++) {
        cout << decay_names[i] << "\t\t" 
             << decay_A[i] << "\t" 
             << decay_Z[i] << "\t"
             << fixed << setprecision(1) << decay_masses[i] << "\t\t"
             << fixed << setprecision(3) << decay_energies[i] << "\t\t"
             << fixed << setprecision(3) << decay_momenta[i] << "\t\t"
             << fixed << setprecision(1) << decay_angles_lab[i] << endl;
    }
    
    // Calculate and show CM frame energies for comparison
    cout << "\nNote: Lab frame energies include parent's kinetic energy (27.877 MeV)" << endl;
    cout << "CM frame relative energies should be much smaller (~Q-value = " << fixed << setprecision(3) << Q_decay << " MeV)" << endl;
}

// Print product summary
void FusionReaction::PrintProductSummary() {
    cout << "\n========== Reaction Summary ==========" << endl;
    cout << "Beam: " << A_beam << " (Z=" << Z_beam << ") " << " (" << E_beam_initial << " MeV)" << endl;
    cout << "Target: " << A_target << " (Z=" << Z_target << ")" << endl;
    cout << "Q-value: " << fixed << setprecision(3) << CalculateQValue() << " MeV" << endl;
    cout << "\nProducts:" << endl;
    for (int i = 0; i < products.size(); i++) {
        cout << "  " << (i+1) << ". " << products[i].name << " (" 
             << products[i].A << ", Z=" << products[i].Z << ") - Mass: " 
             << fixed << setprecision(1) << products[i].mass << " MeV" << endl;
    }
}

// Run simulation
void FusionReaction::RunSimulation(int n_events, bool verbose) {
    cout << "Starting fusion reaction simulation..." << endl;
    PrintProductSummary();
    cout << "Number of events: " << n_events << endl;
    
    for (int event = 0; event < n_events; event++) {
        if (event % 10000 == 0) {
            cout << "Processing event " << event << endl;
        }
        
        CalculateProductKinematics();
        
        // Check energy conservation
        if (verbose && event < 3) {
            bool energy_ok = CheckEnergyConservation();
            if (!energy_ok) {
                cout << "Event " << event << " failed energy conservation!" << endl;
            }
        }
        
        // Lab frame already calculated in CalculateProductKinematics
        // TransformToLabFrame();  // No longer needed
        
        // Simulate decay if enabled
        if (decay_enabled) {
            SimulateDecay();
        }
        
        // Reconstruct total energy (if enabled)
        if (enable_total_energy_reconstruction) {
            ReconstructEnergy();
        }
        
        // Reconstruct parent energy from decay products (if enabled)
        if (enable_energy_reconstruction) {
            ReconstructParentEnergy();
        }
        
        // Reconstruct parent mass from decay products (if enabled)
        if (enable_mass_reconstruction) {
            ReconstructParentMass();
        }
        
        // Reconstruct product properties (if enabled)
        if (enable_product_reconstruction) {
            ReconstructProductProperties();
        }
        
        // Print detailed info for first 2 events if verbose
        // if (verbose && event < 2) {
        //     PrintEventInfo(event);
        //     PrintDecayInfo(event);
        // }
        
        // Fill beam histograms using the actual beam energy used in calculation
        his_beam_E->Fill(E_beam_current);
        
        double tar_x = fRandom->Gaus(0, tar_res);
        double tar_y = fRandom->Gaus(0, tar_res);
        his_beam_pos->Fill(tar_x, tar_y);
    }
    
    cout << "Simulation completed!" << endl;
}

// Save results to ROOT file
void FusionReaction::SaveResults(const char* filename) {
    TFile* file = new TFile(filename, "recreate");
    
    his_beam_E->Write();
    his_beam_pos->Write();
    his_multi_momentum->Write();
    
    // Write energy reconstruction histograms (if enabled)
    if (enable_total_energy_reconstruction) {
        his_total_energy_initial->Write();
        his_total_energy_final->Write();
        his_energy_difference->Write();
        his_total_momentum_mag->Write();
    }
    
    for (int i = 0; i < products.size(); i++) {
        his_product_angle[i]->Write();
        his_product_energy[i]->Write();
        his_product_Evsang[i]->Write();
        his_product_theta_E_lab[i]->Write();
    }
    
    // Write decay histograms (only if decay is enabled)
    if (decay_enabled) {
        for (int i = 0; i < his_decay_angle.size(); i++) {
            his_decay_angle[i]->Write();
            his_decay_energy[i]->Write();
            his_decay_Evsang[i]->Write();
            his_decay_theta_E_lab[i]->Write();
        }
    }
    
    // Write parent reconstruction histograms
    if (his_parent_energy_reconstructed) {
        his_parent_energy_reconstructed->Write();
        his_parent_energy_actual->Write();
        his_parent_energy_difference->Write();
    }
    
    // Write parent mass reconstruction histograms
    if (his_parent_mass_reconstructed) {
        his_parent_mass_reconstructed->Write();
        his_parent_mass_actual->Write();
        his_parent_mass_difference->Write();
    }
    
    // Write product reconstruction histograms
    if (enable_product_reconstruction && his_product1_mass_reconstructed) {
        his_product1_mass_reconstructed->Write();
        his_product1_mass_actual->Write();
        his_product1_mass_difference->Write();
        his_product1_energy_reconstructed->Write();
        his_product1_energy_actual->Write();
        his_product1_energy_difference->Write();
        
        his_product2_mass_reconstructed->Write();
        his_product2_mass_actual->Write();
        his_product2_mass_difference->Write();
        his_product2_energy_reconstructed->Write();
        his_product2_energy_actual->Write();
        his_product2_energy_difference->Write();
    }
    
    file->Close();
    cout << "Results saved to " << filename << endl;
}
