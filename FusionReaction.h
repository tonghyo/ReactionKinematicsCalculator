#ifndef FUSION_REACTION_H
#define FUSION_REACTION_H

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TFile.h>
#include <TGenPhaseSpace.h>
#include <TLorentzVector.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>

using namespace std;

// Particle structure for multi-particle reactions
struct Particle {
    int A;          // Mass number
    int Z;          // Atomic number
    double mass;    // Mass in MeV/c^2
    double energy;  // Kinetic energy
    double theta;   // Polar angle (radians) - CM frame
    double phi;     // Azimuthal angle (radians)
    double momentum; // Momentum magnitude
    double px, py, pz; // Momentum components
    double theta_lab; // Lab frame polar angle
    double energy_lab; // Lab frame energy
    double momentum_lab; // Lab frame momentum magnitude
    double px_lab, py_lab, pz_lab; // Lab frame momentum components
    double excitation_energy; // Excitation energy in MeV (0.0 for ground state)
    string name;    // Particle name
};

class FusionReaction {
private:
    // Beam parameters
    double E_beam_initial;
    double E_loss, E_strag, E_beam_re;
    double tar_res, th_res;
    int A_beam, Z_beam;
    double M_beam;
    double E_beam_current;  // Current event beam energy
    
    // Target parameters
    int A_target, Z_target;
    double M_target;
    
    // Reaction products
    vector<Particle> products;
    vector<double> product_masses;
    vector<int> product_A, product_Z;
    vector<string> product_names;
    
    // Multiple excited states configuration
    bool multiple_excited_states_enabled;
    map<pair<int, int>, vector<double>> excited_states_energies;  // (A,Z) -> excitation energies
    map<pair<int, int>, vector<double>> excited_states_ratios;   // (A,Z) -> branching ratios
    map<pair<int, int>, int> excited_states_product_index;      // (A,Z) -> product index
    
    // Decay configuration
    bool decay_enabled;
    int decay_product_index;  // Which product will decay
    vector<int> decay_A, decay_Z;  // Decay products A, Z
    vector<string> decay_names;    // Decay product names
    vector<double> decay_masses;   // Decay product masses
    
    // Decay product kinematics (for display)
    vector<double> decay_energies;
    vector<double> decay_momenta;
    vector<double> decay_angles;
    vector<double> decay_angles_lab;
    
    // Original parent particle energy (before decay)
    double original_parent_energy;
    
    // Reconstruction control flags
    bool enable_energy_reconstruction;
    bool enable_mass_reconstruction;
    bool enable_total_energy_reconstruction;
    bool enable_product_reconstruction;
    
    // Product selection for reconstruction
    int selected_product1;
    int selected_product2;
    string selected_product1_name;
    string selected_product2_name;
    
    // Original parent particle info (before separation)
    int parent_A;
    int parent_Z;
    string parent_name;
    double parent_mass;
    
    // Random number generator
    TRandom3* fRandom;
    
    // Phase space generator
    TGenPhaseSpace* fPhaseSpace;
    
public:
    // Histograms (public for drawing)
    TH1D* his_beam_E;
    TH2F* his_beam_pos;
    vector<TH1D*> his_product_angle;
    vector<TH1D*> his_product_energy;
    vector<TH2F*> his_product_Evsang;
    vector<TH2F*> his_product_theta_E_lab;  // Lab frame theta vs energy
    TH2F* his_multi_momentum;
    
    // Energy reconstruction histograms
    TH1D* his_total_energy_initial;
    TH1D* his_total_energy_final;
    TH1D* his_energy_difference;
    TH1D* his_total_momentum_mag;
    
    // Decay histograms
    vector<TH1D*> his_decay_angle;
    vector<TH1D*> his_decay_energy;
    vector<TH2F*> his_decay_Evsang;
    vector<TH2F*> his_decay_theta_E_lab;
    
    // Parent particle reconstruction histograms
    TH1D* his_parent_energy_reconstructed;
    TH1D* his_parent_energy_actual;
    TH1D* his_parent_energy_difference;
    
    // Parent particle mass reconstruction histograms
    TH1D* his_parent_mass_reconstructed;
    TH1D* his_parent_mass_actual;
    TH1D* his_parent_mass_difference;
    
    // Product reconstruction histograms
    TH1D* his_product1_mass_reconstructed;
    TH1D* his_product1_mass_actual;
    TH1D* his_product1_mass_difference;
    TH1D* his_product1_energy_reconstructed;
    TH1D* his_product1_energy_actual;
    TH1D* his_product1_energy_difference;
    
    TH1D* his_product2_mass_reconstructed;
    TH1D* his_product2_mass_actual;
    TH1D* his_product2_mass_difference;
    TH1D* his_product2_energy_reconstructed;
    TH1D* his_product2_energy_actual;
    TH1D* his_product2_energy_difference;
    
    // Constructor and Destructor
    FusionReaction();
    ~FusionReaction();
    
    // Setup functions
    void SetBeamParameters(double E_initial, int A, int Z);
    void SetTargetParameters(int A, int Z);
    void AddProduct(int A, int Z, const string& name, double excitation_energy = 0.0);
    void SetExperimentalParameters(double E_loss, double E_strag, double E_beam_re, 
                                  double tar_res, double th_res);
    
    // Multiple excited states functions
    void EnableMultipleExcitedStates(bool enable = true);
    void SetExcitedStates(int A, int Z, const vector<double>& excitation_energies, 
                         const vector<double>& branching_ratios);
    
    // Decay configuration
    void EnableDecay(int product_index);
    void AddDecayProduct(int A, int Z, const string& name);
    void DisableDecay();
    
    // Reconstruction control
    void EnableEnergyReconstruction(bool enable = true);
    void EnableMassReconstruction(bool enable = true);
    void EnableTotalEnergyReconstruction(bool enable = true);
    void EnableProductReconstruction(bool enable = true);
    
    // Product selection
    void SelectProductsForReconstruction(int product1_index, int product2_index);
    void SelectProductsForReconstruction(const string& product1_name, const string& product2_name);
    void SetParentParticleInfo(int A, int Z, const string& name);
    
    // Mass and Histogram functions
    void InitializeHistograms();
    void ReadMassFile(const char* filename);
    
    // Multi-body kinematics
    double CalculateQValue();
    double GeneratePhaseSpace();
    void CalculateProductKinematics();
    void TransformToLabFrame();
    
    // Particle information display
    void PrintEventInfo(int event_num);
    void PrintProductSummary();
    void PrintDecayInfo(int event_num);
    bool CheckEnergyConservation();
    void ReconstructEnergy();
    void ReconstructParentEnergy();
    void ReconstructParentMass();
    void ReconstructProductProperties();
    
    // Main simulation functions
    void RunSimulation(int n_events, bool verbose = false);
    void SaveResults(const char* filename);
    void DrawResults();
    bool CheckConservation();
    
    // Decay simulation functions
    void SimulateDecay();
    void InitializeDecayHistograms();
    void AutoAdjustHistogramRanges();
};

#endif // FUSION_REACTION_H
