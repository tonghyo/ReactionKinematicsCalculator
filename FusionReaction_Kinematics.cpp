#include "FusionReaction.h"

// Calculate Q-value of the reaction
double FusionReaction::CalculateQValue() {
    double total_mass_initial = M_beam + M_target;
    double total_mass_final = 0;
    
    for (int i = 0; i < product_masses.size(); i++) {
        total_mass_final += product_masses[i];
    }
    
    return total_mass_initial - total_mass_final;
}

// Simplified N-body phase space generation
double FusionReaction::GeneratePhaseSpace() {
    double Q_val = CalculateQValue();
    double E_beam = fRandom->Gaus(E_beam_initial, E_beam_re) - E_loss * fRandom->Uniform();
    E_beam = fRandom->Gaus(E_beam, E_strag);
    
    return E_beam + Q_val; // Total available energy
}

// Calculate product kinematics using phase space
void FusionReaction::CalculateProductKinematics() {
    double E_beam = fRandom->Gaus(E_beam_initial, E_beam_re) - E_loss * fRandom->Uniform();
    E_beam = fRandom->Gaus(E_beam, E_strag);
    
    // Store current beam energy for reconstruction
    E_beam_current = E_beam;
    
    int n_products = products.size();
    if (n_products < 2) return;
    
    // Convert to GeV for TGenPhaseSpace
    double E_beam_GeV = E_beam / 1000.0;
    double M_beam_GeV = M_beam / 1000.0;
    double M_target_GeV = M_target / 1000.0;
    
    // Create initial state 4-vector in GeV (Lab frame)
    // Beam total energy = kinetic + mass
    double E_beam_total = E_beam_GeV + M_beam_GeV;
    double p_beam = sqrt(E_beam_total * E_beam_total - M_beam_GeV * M_beam_GeV);  // Relativistic momentum
    
    TLorentzVector beam(0.0, 0.0, p_beam, E_beam_total);
    TLorentzVector target(0.0, 0.0, 0.0, M_target_GeV);
    TLorentzVector W = beam + target;  // Total 4-momentum in Lab frame
    
    
    // Set masses for products in GeV (including excitation energy)
    Double_t masses[10];  // Maximum 10 products
    for (int i = 0; i < n_products; i++) {
        double excited_mass = products[i].mass;
        
        // If multiple excited states are enabled, randomly select one
        if (multiple_excited_states_enabled) {
            pair<int, int> nucleus_key = make_pair(products[i].A, products[i].Z);
            if (excited_states_energies.find(nucleus_key) != excited_states_energies.end()) {
                // Randomly select excited state based on branching ratios
                double random_value = fRandom->Uniform();
                double cumulative_probability = 0.0;
                
                for (int j = 0; j < excited_states_energies[nucleus_key].size(); j++) {
                    cumulative_probability += excited_states_ratios[nucleus_key][j];
                    if (random_value <= cumulative_probability) {
                        excited_mass += excited_states_energies[nucleus_key][j];
                        products[i].excitation_energy = excited_states_energies[nucleus_key][j];
                        break;
                    }
                }
            } else {
                // No excited states configured, use ground state
                excited_mass += products[i].excitation_energy;
            }
        } else {
            // Single excited state mode
            excited_mass += products[i].excitation_energy;
        }
        
        masses[i] = excited_mass / 1000.0;  // Convert to GeV
    }
    
    // Use TGenPhaseSpace for all reactions
    if (!fPhaseSpace->SetDecay(W, n_products, masses)) {
        cout << "ERROR: Phase space generation failed!" << endl;
        return;
    }
    
    Double_t weight = fPhaseSpace->Generate();
    
    // Get decay products (already in Lab frame from TGenPhaseSpace)
    for (int i = 0; i < n_products; i++) {
        TLorentzVector* p = fPhaseSpace->GetDecay(i);
        
        // Extract Lab frame kinematic variables (convert back to MeV)
        products[i].px = p->Px() * 1000.0;
        products[i].py = p->Py() * 1000.0;
        products[i].pz = p->Pz() * 1000.0;
        products[i].momentum = p->P() * 1000.0;
        
        // Set Lab frame momentum components
        products[i].px_lab = products[i].px;
        products[i].py_lab = products[i].py;
        products[i].pz_lab = products[i].pz;
        products[i].momentum_lab = products[i].momentum;
        
        double total_E_GeV = p->E();  // Total energy in GeV (Lab frame)
        products[i].energy = total_E_GeV * 1000.0 - products[i].mass;  // Kinetic energy in MeV (Lab frame)
        products[i].energy_lab = products[i].energy;  // Set Lab frame energy
        
        products[i].theta = p->Theta();
        products[i].phi = p->Phi();
        products[i].theta_lab = p->Theta();  // Lab frame theta (TGenPhaseSpace gives Lab frame results)
        
        // Add angular resolution (experimental uncertainty)
        double theta_with_resolution = products[i].theta + fRandom->Gaus(0, th_res);
        
        // Fill histograms with resolution (Lab frame)
        his_product_angle[i]->Fill(theta_with_resolution * 180.0 / TMath::Pi());
        his_product_energy[i]->Fill(products[i].energy);
        his_product_Evsang[i]->Fill(theta_with_resolution * 180.0 / TMath::Pi(), products[i].energy);
        his_product_theta_E_lab[i]->Fill(theta_with_resolution * 180.0 / TMath::Pi(), products[i].energy);
        his_multi_momentum->Fill(products[i].px, products[i].py);
        
    }
}

// Transform from CM frame to Lab frame
void FusionReaction::TransformToLabFrame() {
    int n_products = products.size();
    if (n_products == 0) return;
    
    // Use the same beam energy that was used in CalculateProductKinematics
    double E_beam = E_beam_current;
    
    // Calculate CM velocity (velocity of CM frame in Lab frame)
    double p_beam = sqrt(E_beam * (E_beam + 2 * M_beam));
    double E_beam_total = E_beam + M_beam;  // Beam total energy
    double E_target_total = M_target;  // Target at rest
    double E_cm_total = sqrt((E_beam_total + E_target_total) * (E_beam_total + E_target_total) - p_beam * p_beam);
    double beta_cm = p_beam / (E_beam_total + E_target_total);  // CM frame velocity
    double gamma_cm = 1.0 / sqrt(1.0 - beta_cm * beta_cm);
    
    // Transform each particle from CM to Lab frame
    for (int i = 0; i < n_products; i++) {
        // CM frame 4-momentum
        double E_cm = products[i].energy + products[i].mass;
        double p_cm = products[i].momentum;
        double pz_cm = products[i].pz;  // z-component in CM
        
        // Lorentz transformation to Lab frame
        double E_lab_total = gamma_cm * (E_cm + beta_cm * pz_cm);
        double pz_lab = gamma_cm * (pz_cm + beta_cm * E_cm);
        
        // Calculate lab frame momentum components
        // x,y components are unchanged in Lorentz boost along z-axis
        double px_lab = products[i].px;  
        double py_lab = products[i].py;
        double p_lab = sqrt(px_lab * px_lab + py_lab * py_lab + pz_lab * pz_lab);
        
        // Store lab frame momentum components
        products[i].px_lab = px_lab;
        products[i].py_lab = py_lab;
        products[i].pz_lab = pz_lab;
        products[i].momentum_lab = p_lab;
        
        // Lab frame kinetic energy
        products[i].energy_lab = E_lab_total - products[i].mass;
        
        // Lab frame angle
        products[i].theta_lab = acos(pz_lab / p_lab);
        
        // Add angular resolution (experimental uncertainty)
        double theta_lab_with_resolution = products[i].theta_lab + fRandom->Gaus(0, th_res);
        
        // Fill lab frame histogram with resolution
        his_product_theta_E_lab[i]->Fill(theta_lab_with_resolution * 180.0 / TMath::Pi(), 
                                        products[i].energy_lab);
    }
}

// Simulate decay of unbound state
void FusionReaction::SimulateDecay() {
    int n_decay_products = decay_A.size();
    if (n_decay_products < 2) return;
    
    // Get the parent particle that will decay
    Particle& parent = products[decay_product_index];
    
    // Only decay if the particle is in an excited state (excitation energy > 0)
    if (parent.excitation_energy <= 0.0) {
        // Ground state - no decay
        return;
    }
    
    // Store original parent energy before decay
    original_parent_energy = parent.energy_lab;
    
    // Calculate decay Q-value
    double Q_decay = parent.mass;
    for (int i = 0; i < n_decay_products; i++) {
        Q_decay -= decay_masses[i];
    }
    
    // Add excitation energy to Q-value (excited state has more energy available)
    Q_decay += parent.excitation_energy;
    
    if (Q_decay <= 0) {
        cout << "WARNING: Decay Q-value is negative or zero: " << Q_decay << " MeV" << endl;
        cout << "  DECAY: " << parent.name << " (excitation: " << parent.excitation_energy 
             << " MeV) -> Q-value: " << Q_decay << " MeV" << endl;
        cout << "  Parent mass: " << parent.mass << " MeV/c^2" << endl;
        cout << "  Decay products: ";
        for (int i = 0; i < n_decay_products; i++) {
            cout << decay_names[i] << " (" << decay_masses[i] << " MeV/c^2)";
            if (i < n_decay_products - 1) cout << " + ";
        }
        cout << endl;
        return;
    }
    
    // Create parent 4-momentum in Lab frame (convert to GeV for TLorentzVector)
    double E_parent_total = (parent.energy + parent.mass) / 1000.0;  // Convert to GeV
    double px_parent = parent.px / 1000.0;  // Convert to GeV/c
    double py_parent = parent.py / 1000.0;  // Convert to GeV/c
    double pz_parent = parent.pz / 1000.0;  // Convert to GeV/c
    
    TLorentzVector parent_4vec_lab(px_parent, py_parent, pz_parent, E_parent_total);
    
    // Boost to parent's CM frame for decay calculation
    TVector3 parent_boost = parent_4vec_lab.BoostVector();
    TLorentzVector parent_4vec_cm = parent_4vec_lab;
    parent_4vec_cm.Boost(-parent_boost);  // Boost to CM frame
    
    
    
    // Set masses for decay products in GeV
    Double_t decay_masses_GeV[10];
    for (int i = 0; i < n_decay_products; i++) {
        decay_masses_GeV[i] = decay_masses[i] / 1000.0;
    }
    
    // Use TGenPhaseSpace for decay in parent's CM frame
    TGenPhaseSpace decay_phase_space;
    if (!decay_phase_space.SetDecay(parent_4vec_cm, n_decay_products, decay_masses_GeV)) {
        cout << "ERROR: Decay phase space generation failed!" << endl;
        return;
    }
    
    Double_t decay_weight = decay_phase_space.Generate();
    
    // Get decay products in parent's CM frame, then boost to Lab frame
    for (int i = 0; i < n_decay_products; i++) {
        TLorentzVector* decay_p_cm = decay_phase_space.GetDecay(i);
        
        // Convert to Lab frame by boosting back
        TLorentzVector decay_p_lab = *decay_p_cm;
        decay_p_lab.Boost(parent_boost);  // Boost from CM to Lab frame
        
        // Extract decay product kinematic variables (convert back to MeV)
        double px_decay = decay_p_lab.Px() * 1000.0;
        double py_decay = decay_p_lab.Py() * 1000.0;
        double pz_decay = decay_p_lab.Pz() * 1000.0;
        double p_decay = decay_p_lab.P() * 1000.0;
        
        double E_decay_total = decay_p_lab.E() * 1000.0;
        double E_decay_kinetic = E_decay_total - decay_masses[i];
        
        double theta_decay = decay_p_lab.Theta();
        double phi_decay = decay_p_lab.Phi();
        
        // Add experimental resolution
        double theta_decay_with_resolution = theta_decay + fRandom->Gaus(0, th_res);
        double E_decay_kinetic_with_resolution = E_decay_kinetic + fRandom->Gaus(0, E_beam_re);
        
        // Calculate momentum with energy resolution effect
        // p = sqrt(E_kinetic * (E_kinetic + 2*mass))
        double E_total_with_resolution = E_decay_kinetic_with_resolution + decay_masses[i];
        double p_decay_with_resolution = sqrt(E_decay_kinetic_with_resolution * (E_decay_kinetic_with_resolution + 2 * decay_masses[i]));
        
        // Store decay product kinematics for display (Lab frame)
        // Use resolution-applied values for reconstruction analysis
        if (i < decay_energies.size()) {
            decay_energies[i] = E_decay_kinetic_with_resolution;  // Use resolution-applied energy
            decay_momenta[i] = p_decay_with_resolution;           // Use resolution-applied momentum
            decay_angles[i] = theta_decay * 180.0 / TMath::Pi();
            decay_angles_lab[i] = theta_decay_with_resolution * 180.0 / TMath::Pi();
        } else {
            decay_energies.push_back(E_decay_kinetic_with_resolution);  // Use resolution-applied energy
            decay_momenta.push_back(p_decay_with_resolution);           // Use resolution-applied momentum
            decay_angles.push_back(theta_decay * 180.0 / TMath::Pi());
            decay_angles_lab.push_back(theta_decay_with_resolution * 180.0 / TMath::Pi());
        }
        
        // Fill decay histograms with resolution
        his_decay_angle[i]->Fill(theta_decay_with_resolution * 180.0 / TMath::Pi());
        his_decay_energy[i]->Fill(E_decay_kinetic_with_resolution);
        his_decay_Evsang[i]->Fill(theta_decay_with_resolution * 180.0 / TMath::Pi(), E_decay_kinetic_with_resolution);
        his_decay_theta_E_lab[i]->Fill(theta_decay_with_resolution * 180.0 / TMath::Pi(), E_decay_kinetic_with_resolution);
    }
}
