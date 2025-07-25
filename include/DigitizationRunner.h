/*
 * Copyright (C) 2025 CYGNO Collaboration
 *
 *
 * Author: Stefano Piacentini
 * Created in 2025
 *
 */

#pragma once

#include "ConfigManager.h"
#include "Globals.h"
#include "cygnolib.h"
#include <string>
#include <vector>
#include <map>
#include <memory>
#include <utility>
#include "TH2F.h"

// Forward declarations
class TFile;

/**
 * @class DigitizationRunner
 * @author Stefano Piacentini
 * @brief Manages the digitization process from MC simulation input to ROOT image output.
 *
 * @details This class initializes the configuration file, sets up camera parameters, processes ROOT files,
 * and handles digitization including background addition, and application of vignetting maps.
 */
class DigitizationRunner {
public:
    /**
     * @brief Constructor.
     * @param[in] configFile Path to the configuration file.
     * @param[in] inputDir Path to the directory containing input ROOT files.
     * @param[in] outputDir Path to the directory where output will be written.
     */
    DigitizationRunner(const std::string& configFile, const std::string& inputDir, const std::string& outputDir);

    /**
     * @brief Executes the full digitization pipeline.
     */
    void run();

private:

    /**
     * @brief Initializes global variables like GEM gains and extraction efficiencies.
     */
     void initializeGlobals();

     /**
      * @brief Prepares camera parameters such as pixel size, optics, and sensor dimensions.
      */
     void prepareCameraSettings();
 
     /**
      * @brief Validates consistency between simulation dimensions and camera dimensions.
      */
     void checkDimensionConsistency();
 
     /**
      * @brief Loops through ROOT files in the input directory and processes them.
      */
     void processRootFiles();
 
     /**
      * @brief Checks whether a given file is a valid input ROOT file.
      * @param[in] filename Name of the file to check.
      * @return True if valid, false otherwise.
      */
     bool isValidInputFile(const std::string& filename) const;
 
     /**
      * @brief Sets the random seed for reproducibility.
      * @param[in] seed Seed value (default is 10).
      */
     void setSeed(int seed = 10);
 
     /**
      * @brief Loads ion list for SRIM format events.
      * @param[in] config ConfigManager object.
      * @param[in] SRIM_events Container to store parsed events.
      * @param[in] filename File name of the ion list.
      * @param[in] infolder Path to the input folder.
      */
     void loadIonList4SRIM(ConfigManager& config,
                           std::vector<std::vector<std::string>>& SRIM_events,
                           const std::string& filename,
                           const std::string& infolder);
 
     /**
      * @brief Initializes the source directory based on configFile variable.
      */
     void initSourceDir();
 
     /**
      * @brief Saves global and run parameters into the output ROOT file.
      * @param[in] outfile Shared pointer to the output TFile.
      */
     void SaveValues(std::shared_ptr<TFile>& outfile);
 
     /**
      * @brief Fills the vignetting map histogram.
      * @param[in] VignMap Reference to the 2D histogram containing the vignetting map.
      */
     void fillVigmap(TH2F& VignMap);
 
     /**
      * @brief Determines whether an event is a nuclear recoil (NR).
      * @param[in] pdgID_hits List of PDG codes for all hits.
      * @param[in] pdg PDG code to compare against.
      * @return True if NR, false otherwise.
      */
     bool is_NR(std::vector<int> pdgID_hits, int pdg);
 
     /**
      * @brief Adds pedestal to an image, according the configfile options.
      * @param[in] background Reference to 2D image matrix that will contain the pedestal.
      */
     void AddBckg(std::vector<std::vector<int>>& background);

     /**
     * @brief Parses a MC axis label string and returns the corresponding digititazion axis and sign.
     *
     * This function interprets a string like "x", "-y", or "z" to determine which
     * digitization axis corresponds to a MC axis, and whether that axis is inverted.
     *
     * For example:
     * - "x"   → returns ('x', 1.0)
     * - "-z"  → returns ('z', -1.0)
     *
     * @param axis_label A string representing the axis mapping from the config file.
     *                   Expected format is one of: "x", "y", "z", "-x", "-y", or "-z".
     *
     * @return A pair where:
     *         - first: the digitization axis as a char ('x', 'y', or 'z')
     *         - second: a double representing the sign (1.0 or -1.0)
     */
     std::pair<char, double> getAxisMapping(const std::string& axis_label);
 
     ConfigManager config;           ///< Configuration manager
     std::string configFile;         ///< Path to the config file
     std::string inputDir;           ///< Path to input directory
     std::string outputDir;          ///< Path to output directory
     int runCount = 0;               ///< Current run number
     int NMAX_EVENTS = 200;          ///< Max number of events to process per output file
     std::string SOURCE_DIR = "";    ///< Path to event source directory
 
};
