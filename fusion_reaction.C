#include "FusionReaction.h"
#include "TApplication.h"

// Main simulation function
void run_fusion_simulation() {
    FusionReaction reaction;
    // Set beam (22Na) (Energy, A, Z)
    reaction.SetBeamParameters(142, 25, 13); // 143 MeV 22Na
    // Set target (24Mg) (A, Z)
    reaction.SetTargetParameters(2, 1);
    // Set experimental parameters (optional - uses defaults if not called)
    // Parameters: E_loss, E_strag, E_beam_re, tar_res, th_res
    reaction.SetExperimentalParameters(1.0, 0.05, 0.1, 0.5, 0.1*TMath::Pi()/180.0);
    
    
    // Add products first
    reaction.AddProduct(26, 14, "Si26");   // Si26 with random excited state
    reaction.AddProduct(1, 0, "n1");      // neutron 1
    
    // Heavy recoil (Si26) multiple excited states configuration
    // Define multiple excited states with their branching ratios
    vector<double> heavy_recoil_excitation_energies = {5.92, 5.9, 6.3, 6.7};  // MeV
    vector<double> heavy_recoil_branching_ratios = {1.0, 0.0, 0.0, 0.0};   // Relative probabilities
    
    // Enable multiple excited states simulation
    reaction.EnableMultipleExcitedStates(true);
    // Set custom excited states for Si26
    reaction.SetExcitedStates(26, 14, heavy_recoil_excitation_energies, heavy_recoil_branching_ratios);
    
    // Optional: Enable decay for 42V (index 0)
    reaction.EnableDecay(0);  // 42V will decay
    reaction.AddDecayProduct(25, 13, "25Al");  // 42V -> 41Ti + p
    reaction.AddDecayProduct(1, 1, "p");       // proton
    
    // Optional: Enable reconstruction analysis
    //reaction.EnableTotalEnergyReconstruction(true);  // Enable total energy reconstruction
    //reaction.EnableEnergyReconstruction(true);       // Enable parent energy reconstruction
    reaction.EnableMassReconstruction(true);         // Enable parent mass reconstruction
    
    // Optional: Enable product reconstruction
    //reaction.SelectProductsForReconstruction("Ti41", "p1");  // Select products for reconstruction
    //reaction.EnableProductReconstruction(true);      // Enable product reconstruction (must be called after SelectProductsForReconstruction)
    
    // Read mass file   
    reaction.ReadMassFile("mass.dat");
    
    // Initialize histograms (automatically initializes decay histograms if decay is enabled)
    reaction.InitializeHistograms();
    
    // Run simulation with verbose output for first 5 events
    reaction.RunSimulation(10000, true);
    // Save results
    reaction.SaveResults("fusion_results.root");
    // Draw results
    reaction.DrawResults();
}

// Auto-run when loading the macro
void fusion_reaction() {
    run_fusion_simulation();
}

// Main function for standalone execution
int main() {
    // Enable ROOT GUI
    TApplication app("FusionReaction", 0, 0);
    
    run_fusion_simulation();
    
    // Keep GUI alive
    app.Run();
    return 0;
}
