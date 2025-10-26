#include "FusionReaction.h"

// Constructor
FusionReaction::FusionReaction() {
    fRandom = new TRandom3();
    fRandom->SetSeed(time(0));
    
    fPhaseSpace = new TGenPhaseSpace();
    
    // Default parameters (can be overridden by SetExperimentalParameters)
    E_loss = 1.0;
    E_strag = 0.05;
    E_beam_re = 0.05;
    tar_res = 0.5; // mm
    th_res = 0.1 * TMath::Pi() / 180.0; // radians
    
    products.clear();
    product_masses.clear();
    product_A.clear();
    product_Z.clear();
    
    // Initialize decay configuration
    decay_enabled = false;
    decay_product_index = -1;
    decay_A.clear();
    decay_Z.clear();
    decay_names.clear();
    decay_masses.clear();
    decay_energies.clear();
    decay_momenta.clear();
    decay_angles.clear();
    decay_angles_lab.clear();
    original_parent_energy = 0.0;
    
    // Initialize reconstruction control flags
    enable_energy_reconstruction = false;
    enable_mass_reconstruction = false;
    enable_total_energy_reconstruction = false;
    enable_product_reconstruction = false;
    
    // Initialize multiple excited states
    multiple_excited_states_enabled = false;
    excited_states_energies.clear();
    excited_states_ratios.clear();
    excited_states_product_index.clear();
    
    // Initialize product selection
    selected_product1 = -1;
    selected_product2 = -1;
    selected_product1_name = "";
    selected_product2_name = "";
    
    // Initialize parent particle info
    parent_A = 0;
    parent_Z = 0;
    parent_name = "";
    parent_mass = 0.0;
    
    // Initialize histogram pointers to nullptr
    his_parent_energy_reconstructed = nullptr;
    his_parent_energy_actual = nullptr;
    his_parent_energy_difference = nullptr;
    his_parent_mass_reconstructed = nullptr;
    his_parent_mass_actual = nullptr;
    his_parent_mass_difference = nullptr;
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

// Destructor
FusionReaction::~FusionReaction() {
    delete fRandom;
    delete fPhaseSpace;
}

// Set beam parameters (Energy, A, Z)
void FusionReaction::SetBeamParameters(double E_initial, int A, int Z) {
    E_beam_initial = E_initial;
    A_beam = A;
    Z_beam = Z;
}

// Set target parameters (A, Z)
void FusionReaction::SetTargetParameters(int A, int Z) {
    A_target = A;
    Z_target = Z;
}

// Add product (A, Z, name, excitation_energy)
void FusionReaction::AddProduct(int A, int Z, const string& name, double excitation_energy) {
    product_A.push_back(A);
    product_Z.push_back(Z);
    product_masses.push_back(0.0); // Will be filled from mass.dat
    product_names.push_back(name);
    
    Particle p;
    p.A = A;
    p.Z = Z;
    p.mass = 0.0; // Will be filled from mass.dat
    p.name = name;
    p.excitation_energy = excitation_energy; // Store excitation energy
    products.push_back(p);
    
    // Store product index for excited states if multiple excited states are enabled
    if (multiple_excited_states_enabled) {
        pair<int, int> nucleus_key = make_pair(A, Z);
        excited_states_product_index[nucleus_key] = products.size() - 1;
    }
    
    if (excitation_energy > 0.0) {
        cout << "Added product: " << name << " (A=" << A << ", Z=" << Z << ") with excitation energy: " 
             << excitation_energy << " MeV" << endl;
    } else {
        cout << "Added product: " << name << " (A=" << A << ", Z=" << Z << ") - ground state" << endl;
    }
}

// Set experimental parameters
void FusionReaction::SetExperimentalParameters(double E_loss, double E_strag, double E_beam_re, 
                                               double tar_res, double th_res) {
    this->E_loss = E_loss;
    this->E_strag = E_strag;
    this->E_beam_re = E_beam_re;
    this->tar_res = tar_res;
    this->th_res = th_res;
}

// Check conservation of A and Z numbers
bool FusionReaction::CheckConservation() {
    cout << "\n========== Conservation Check ==========" << endl;
    
    // Calculate initial A and Z
    int A_initial = A_beam + A_target;
    int Z_initial = Z_beam + Z_target;
    
    // Calculate final A and Z
    int A_final = 0;
    int Z_final = 0;
    
    for (int i = 0; i < product_A.size(); i++) {
        A_final += product_A[i];
        Z_final += product_Z[i];
    }
    
    // Check conservation
    bool A_conserved = (A_initial == A_final);
    bool Z_conserved = (Z_initial == Z_final);
    
    // Print results
    cout << "Initial: A = " << A_initial << ", Z = " << Z_initial << endl;
    cout << "Final:   A = " << A_final << ", Z = " << Z_final << endl;
    cout << "A conservation: " << (A_conserved ? "✓ PASS" : "✗ FAIL") << endl;
    cout << "Z conservation: " << (Z_conserved ? "✓ PASS" : "✗ FAIL") << endl;
    
    if (A_conserved && Z_conserved) {
        cout << "Overall: ✓ ALL CONSERVATION LAWS SATISFIED" << endl;
        return true;
    } else {
        cout << "Overall: ✗ CONSERVATION LAWS VIOLATED" << endl;
        cout << "ERROR: Program terminated due to conservation law violation!" << endl;
        exit(1);
    }
}

// Enable decay for a specific product
void FusionReaction::EnableDecay(int product_index) {
    if (product_index < 0 || product_index >= products.size()) {
        cout << "ERROR: Invalid product index for decay!" << endl;
        exit(1);
    }
    
    decay_enabled = true;
    decay_product_index = product_index;
    cout << "Decay enabled for product: " << product_names[product_index] << endl;
}

// Add decay product
void FusionReaction::AddDecayProduct(int A, int Z, const string& name) {
    if (!decay_enabled) {
        cout << "ERROR: Decay not enabled! Call EnableDecay() first." << endl;
        exit(1);
    }
    
    decay_A.push_back(A);
    decay_Z.push_back(Z);
    decay_names.push_back(name);
    decay_masses.push_back(0.0); // Will be filled from mass.dat
    
    cout << "Added decay product: " << name << " (A=" << A << ", Z=" << Z << ")" << endl;
}

// Disable decay
void FusionReaction::DisableDecay() {
    decay_enabled = false;
    decay_product_index = -1;
    decay_A.clear();
    decay_Z.clear();
    decay_names.clear();
    decay_masses.clear();
    cout << "Decay disabled." << endl;
}

// Enable/disable energy reconstruction
void FusionReaction::EnableEnergyReconstruction(bool enable) {
    enable_energy_reconstruction = enable;
    if (enable) {
        cout << "Energy reconstruction enabled" << endl;
    } else {
        cout << "Energy reconstruction disabled" << endl;
    }
}

// Enable/disable mass reconstruction
void FusionReaction::EnableMassReconstruction(bool enable) {
    enable_mass_reconstruction = enable;
    if (enable) {
        cout << "Mass reconstruction enabled" << endl;
    } else {
        cout << "Mass reconstruction disabled" << endl;
    }
}

// Enable/disable total energy reconstruction
void FusionReaction::EnableTotalEnergyReconstruction(bool enable) {
    enable_total_energy_reconstruction = enable;
    if (enable) {
        cout << "Total energy reconstruction enabled" << endl;
    } else {
        cout << "Total energy reconstruction disabled" << endl;
    }
}

// Enable/disable product reconstruction
void FusionReaction::EnableProductReconstruction(bool enable) {
    enable_product_reconstruction = enable;
    if (enable) {
        cout << "Product reconstruction enabled" << endl;
        
        // Auto-calculate parent particle info from selected products
        if (selected_product1 >= 0 && selected_product2 >= 0 && 
            selected_product1 < products.size() && selected_product2 < products.size()) {
            
            // Sum A and Z of selected products
            parent_A = product_A[selected_product1] + product_A[selected_product2];
            parent_Z = product_Z[selected_product1] + product_Z[selected_product2];
            parent_name = "Parent_" + selected_product1_name + "_" + selected_product2_name;
            parent_mass = 0.0;  // Will be loaded from mass file
            
            cout << "Auto-calculated parent particle info:" << endl;
            cout << "  Parent: " << parent_name << " (A=" << parent_A << ", Z=" << parent_Z << ")" << endl;
            cout << "  From: " << selected_product1_name << " (A=" << product_A[selected_product1] << ", Z=" << product_Z[selected_product1] << ")" << endl;
            cout << "       + " << selected_product2_name << " (A=" << product_A[selected_product2] << ", Z=" << product_Z[selected_product2] << ")" << endl;
        }
    } else {
        cout << "Product reconstruction disabled" << endl;
    }
}

// Select products for reconstruction
void FusionReaction::SelectProductsForReconstruction(int product1_index, int product2_index) {
    if (product1_index < 0 || product1_index >= products.size() || 
        product2_index < 0 || product2_index >= products.size() ||
        product1_index == product2_index) {
        cout << "ERROR: Invalid product indices for reconstruction!" << endl;
        cout << "Product indices must be different and within range [0, " << products.size()-1 << "]" << endl;
        return;
    }
    
    selected_product1 = product1_index;
    selected_product2 = product2_index;
    selected_product1_name = product_names[product1_index];
    selected_product2_name = product_names[product2_index];
    
    cout << "Selected products for reconstruction:" << endl;
    cout << "  Product 1: " << product_names[product1_index] << " (index " << product1_index << ")" << endl;
    cout << "  Product 2: " << product_names[product2_index] << " (index " << product2_index << ")" << endl;
}

// Select products for reconstruction by name
void FusionReaction::SelectProductsForReconstruction(const string& product1_name, const string& product2_name) {
    int product1_index = -1;
    int product2_index = -1;
    
    // Find indices by name
    for (int i = 0; i < product_names.size(); i++) {
        if (product_names[i] == product1_name) {
            product1_index = i;
        }
        if (product_names[i] == product2_name) {
            product2_index = i;
        }
    }
    
    // Check if both products were found
    if (product1_index == -1) {
        cout << "ERROR: Product '" << product1_name << "' not found!" << endl;
        cout << "Available products: ";
        for (int i = 0; i < product_names.size(); i++) {
            cout << product_names[i];
            if (i < product_names.size() - 1) cout << ", ";
        }
        cout << endl;
        return;
    }
    
    if (product2_index == -1) {
        cout << "ERROR: Product '" << product2_name << "' not found!" << endl;
        cout << "Available products: ";
        for (int i = 0; i < product_names.size(); i++) {
            cout << product_names[i];
            if (i < product_names.size() - 1) cout << ", ";
        }
        cout << endl;
        return;
    }
    
    if (product1_index == product2_index) {
        cout << "ERROR: Cannot select the same product twice!" << endl;
        return;
    }
    
    // Set selection
    selected_product1 = product1_index;
    selected_product2 = product2_index;
    selected_product1_name = product1_name;
    selected_product2_name = product2_name;
    
    cout << "Selected products for reconstruction:" << endl;
    cout << "  Product 1: " << product1_name << " (index " << product1_index << ")" << endl;
    cout << "  Product 2: " << product2_name << " (index " << product2_index << ")" << endl;
}

// Set parent particle info for reconstruction
void FusionReaction::SetParentParticleInfo(int A, int Z, const string& name) {
    parent_A = A;
    parent_Z = Z;
    parent_name = name;
    parent_mass = 0.0;  // Will be loaded from mass file
    
    cout << "Set parent particle info for reconstruction:" << endl;
    cout << "  Parent: " << name << " (A=" << A << ", Z=" << Z << ")" << endl;
}

// Enable/disable multiple excited states simulation
void FusionReaction::EnableMultipleExcitedStates(bool enable) {
    multiple_excited_states_enabled = enable;
    if (enable) {
        cout << "Multiple excited states simulation enabled" << endl;
        
        // Automatically set excited states for all products
        for (int i = 0; i < products.size(); i++) {
            // Use default excited states for common nuclei
            if (products[i].A == 26 && products[i].Z == 14) { // Si26
                vector<double> si26_energies = {0.0, 5.9, 6.3, 6.7, 8.0};
                vector<double> si26_ratios = {0.1, 0.2, 0.2, 0.3, 0.2};
                SetExcitedStates(26, 14, si26_energies, si26_ratios);
            }
            // Add more nuclei as needed
        }
    } else {
        cout << "Multiple excited states simulation disabled" << endl;
    }
}

// Set excited states for a specific nucleus
void FusionReaction::SetExcitedStates(int A, int Z, const vector<double>& excitation_energies, 
                                     const vector<double>& branching_ratios) {
    if (!multiple_excited_states_enabled) {
        cout << "ERROR: Multiple excited states not enabled! Call EnableMultipleExcitedStates(true) first." << endl;
        return;
    }
    
    if (excitation_energies.size() != branching_ratios.size()) {
        cout << "ERROR: Number of excitation energies must match number of branching ratios!" << endl;
        return;
    }
    
    // Normalize branching ratios
    double total_ratio = 0.0;
    for (double ratio : branching_ratios) {
        total_ratio += ratio;
    }
    
    vector<double> normalized_ratios;
    for (double ratio : branching_ratios) {
        normalized_ratios.push_back(ratio / total_ratio);
    }
    
    // Store excited states configuration
    pair<int, int> nucleus_key = make_pair(A, Z);
    excited_states_energies[nucleus_key] = excitation_energies;
    excited_states_ratios[nucleus_key] = normalized_ratios;
    
    cout << "Set excited states for nucleus A=" << A << ", Z=" << Z << ":" << endl;
    for (int i = 0; i < excitation_energies.size(); i++) {
        cout << "  State " << i << ": " << excitation_energies[i] << " MeV (ratio: " 
             << normalized_ratios[i] << ")" << endl;
    }
}
