#include "FusionReaction.h"
#include "TApplication.h"
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <algorithm>

// Simple helpers
static inline std::string Trim(const std::string &s) {
    auto start = s.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    auto end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

static inline std::vector<std::string> Split(const std::string &s, char delim) {
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) out.push_back(Trim(item));
    return out;
}

// Parse a comma separated list of doubles
static inline std::vector<double> ParseDoubles(const std::string &s) {
    std::vector<double> out;
    for (auto &tok : Split(s, ',')) {
        if (tok.empty()) continue;
        out.push_back(std::stod(tok));
    }
    return out;
}

// Parse product entries of form A,Z,label separated by semicolon or comma groups
static inline std::vector<std::tuple<int,int,std::string>> ParseProducts(const std::string &s) {
    std::vector<std::tuple<int,int,std::string>> out;
    // products may be separated by ';'
    for (auto &entry : Split(s, ';')) {
        if (entry.empty()) continue;
        // allow either A,Z,label or A,Z,label with commas
        auto parts = Split(entry, ',');
        if (parts.size() >= 3) {
            int A = std::stoi(parts[0]);
            int Z = std::stoi(parts[1]);
            std::string label = parts[2];
            out.emplace_back(A, Z, label);
        }
    }
    return out;
}

// Read key=value parameter file; lines starting with # are comments
static inline std::map<std::string,std::string> ReadParamFile(const std::string &path) {
    std::map<std::string,std::string> params;
    std::ifstream ifs(path);
    if (!ifs) return params;
    std::string line;
    while (std::getline(ifs, line)) {
        auto pos = line.find('#');
        if (pos != std::string::npos) line = line.substr(0, pos);
        auto eq = line.find('=');
        if (eq == std::string::npos) continue;
        std::string key = Trim(line.substr(0, eq));
        std::string val = Trim(line.substr(eq+1));
        if (!key.empty()) params[key] = val;
    }
    return params;
}

// Main simulation function: optional param file path
void run_fusion_simulation(const char *paramFilePath = "params.txt") {
    FusionReaction reaction;

    // Read parameters from file (if exists)
    auto params = ReadParamFile(paramFilePath ? paramFilePath : "");

    // 1) Beam parameters: beam = Energy,A,Z
    if (params.count("beam")) {
        auto parts = Split(params["beam"], ',');
        if (parts.size() >= 3) {
            double E = std::stod(parts[0]);
            int A = std::stoi(parts[1]);
            int Z = std::stoi(parts[2]);
            reaction.SetBeamParameters(E, A, Z);
        }
    } else {
        cerr << "No beam parameters specified in param file." << endl;
        return;
    }

    // 2) Target parameters: target = A,Z
    if (params.count("target")) {
        auto parts = Split(params["target"], ',');
        if (parts.size() >= 2) {
            int A = std::stoi(parts[0]);
            int Z = std::stoi(parts[1]);
            reaction.SetTargetParameters(A, Z);
        }
    } else {
        cerr << "No target parameters specified in param file." << endl;
        return;
    }

    // 3) Experimental parameters: experimental = E_loss,E_strag,E_beam_re,tar_res,th_res_deg
    if (params.count("experimental")) {
        auto parts = Split(params["experimental"], ',');
        if (parts.size() >= 5) {
            double E_loss = std::stod(parts[0]);
            double E_strag = std::stod(parts[1]);
            double E_beam_re = std::stod(parts[2]);
            double tar_res = std::stod(parts[3]);
            double th_res_deg = std::stod(parts[4]);
            reaction.SetExperimentalParameters(E_loss, E_strag, E_beam_re, tar_res, th_res_deg*TMath::Pi()/180.0);
        }
    } else {
        cerr << "No experimental parameters specified in param file." << endl;
        return;
    }

    // 4) Products: products = A,Z,label;A,Z,label;...
    if (params.count("products")) {
        auto prods = ParseProducts(params["products"]);
        for (auto &t : prods) {
            reaction.AddProduct(std::get<0>(t), std::get<1>(t), std::get<2>(t).c_str());
        }
    } else {
        cerr << "No products specified in param file." << endl;
        return;
    }

    // 5) Multiple excited states
    if (params.count("multiple_excited_states")) {
        std::string val = params["multiple_excited_states"];
        std::transform(val.begin(), val.end(), val.begin(), ::tolower);
        reaction.EnableMultipleExcitedStates(val == "1" || val == "true" || val == "yes");
    } else {
        cerr << "No multiple excited states parameter specified in param file." << endl;
        return;
    }

    // 6) Excited energies and branching
    if (params.count("excited_energies") && params.count("excited_branching")) {
        auto energies = ParseDoubles(params["excited_energies"]);
        auto branches = ParseDoubles(params["excited_branching"]);
        // If sizes match and non-empty, set for the first heavy product by A,Z from params (optional keys heavy_A heavy_Z)
        int heavyA = 26, heavyZ = 14;
        if (params.count("heavy_A")) heavyA = std::stoi(params["heavy_A"]);
        if (params.count("heavy_Z")) heavyZ = std::stoi(params["heavy_Z"]);
        if (!energies.empty() && energies.size() == branches.size()) {
            reaction.SetExcitedStates(heavyA, heavyZ, energies, branches);
        }
    }

    // 7) Decay: enable_decay = index (int) or enable_decay=true with decay_products key
    if (params.count("enable_decay")) {
        std::string v = params["enable_decay"];
        std::transform(v.begin(), v.end(), v.begin(), ::tolower);
        if (v == "true" || v == "1" || v == "yes") {
            // expect decay_products = A,Z,label;A,Z,label
            reaction.EnableDecay(0);
            if (params.count("decay_products")) {
                auto dprods = ParseProducts(params["decay_products"]);
                for (auto &t : dprods) {
                    reaction.AddDecayProduct(std::get<0>(t), std::get<1>(t), std::get<2>(t).c_str());
                }
            }
        }
    } else {
        reaction.DisableDecay();
    }

    // 8) Reconstruction flags
    if (params.count("enable_mass_reconstruction")) {
        std::string v = params["enable_mass_reconstruction"]; std::transform(v.begin(), v.end(), v.begin(), ::tolower);
        reaction.EnableMassReconstruction(v == "1" || v == "true" || v == "yes");
    } else {
        reaction.EnableMassReconstruction(false);
    }
    if (params.count("enable_total_energy_reconstruction")) {
        std::string v = params["enable_total_energy_reconstruction"]; std::transform(v.begin(), v.end(), v.begin(), ::tolower);
        reaction.EnableTotalEnergyReconstruction(v == "1" || v == "true" || v == "yes");
    } else {
        reaction.EnableTotalEnergyReconstruction(false);
    }
    if (params.count("enable_energy_reconstruction")) {
        std::string v = params["enable_energy_reconstruction"]; std::transform(v.begin(), v.end(), v.begin(), ::tolower);
        reaction.EnableEnergyReconstruction(v == "1" || v == "true" || v == "yes");
    } else {
        reaction.EnableEnergyReconstruction(false);
    }

    // 9) Product reconstruction flags
    if (params.count("select_product")) {
        auto parts = Split(params["select_product"], ',');
        if (parts.size() >= 2) {
            reaction.SelectProductsForReconstruction(parts[0], parts[1]);
            reaction.EnableProductReconstruction(true);
        } else {
            cerr << "Invalid select_product parameter format." << endl;
        }
    } else {
        reaction.EnableProductReconstruction(false);
    }

    // 10) Mass file
    if (params.count("mass_file")) reaction.ReadMassFile(params["mass_file"].c_str());
    else reaction.ReadMassFile("mass.dat");

    // 11) Initialize histograms
    reaction.InitializeHistograms();

    // 12) Run simulation: n_events (default 10000), verbose_events (bool)
    int n_events = 10000;
    bool verbose = true;
    if (params.count("n_events")) n_events = std::stoi(params["n_events"]);
    if (params.count("verbose_events")) {
        std::string v = params["verbose_events"]; std::transform(v.begin(), v.end(), v.begin(), ::tolower);
        verbose = (v == "1" || v == "true" || v == "yes");
    }
    reaction.RunSimulation(n_events, verbose);

    // 13) Results file
    if (params.count("output_file")) reaction.SaveResults(params["output_file"].c_str());
    else reaction.SaveResults("fusion_results.root");

    // Draw results if requested
    if (!params.count("no_draw")) reaction.DrawResults();
}

// Auto-run when loading the macro (ROOT will call this)
void fusion_reaction() {
    run_fusion_simulation();
}

// Main function for standalone execution
int main(int argc, char **argv) {
    const char *paramFile = "params.txt";
    if (argc >= 2 && argv[1] && argv[1][0] != '\0') paramFile = argv[1];
    cout << "Using parameter file: " << paramFile << endl;

    // Enable ROOT GUI - create TApplication after we've captured user args
    TApplication app("FusionReaction", &argc, argv);

    run_fusion_simulation(paramFile);

    // Keep GUI alive
    app.Run();
    return 0;
}
