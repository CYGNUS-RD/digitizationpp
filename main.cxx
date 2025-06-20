/*
 * Copyright (C) 2025 CYGNO Collaboration
 *
 *
 * Author: Stefano Piacentini
 * Created in 2025
 *
 */

#include "ConfigManager.h"
#include "DigitizationRunner.h"

#include <iostream>
#include <unistd.h>
#include <limits.h>
#include <filesystem>

using namespace std;

int main(int argc, char** argv) {
    if (argc < 2) {
        cerr << "No Configfile given!\nUsage: ./mc_sim Configfile -I inputDir -O outputDir" << endl;
        return EXIT_FAILURE;
    }

    string configFile = argv[1];
    string inputDir, outputDir;

    // Handle optional input/output dirs
    if (argc == 6) {
        string arg2 = argv[2];
        string val2 = argv[3];
        string arg3 = argv[4];
        string val3 = argv[5];

        if (arg2 == "-I" && arg3 == "-O") {
            inputDir = val2;
            outputDir = val3;
        } else if (arg2 == "-O" && arg3 == "-I") {
            outputDir = val2;
            inputDir = val3;
        } else {
            cerr << "Invalid argument structure. Use -I inputDir -O outputDir" << endl;
            return EXIT_FAILURE;
        }
    } else {
        // default fallback: current working directory + relative paths
        char cwd[PATH_MAX];
        char* ret = getcwd(cwd, sizeof(cwd));
        if(!ret) {
            cerr << "Failed to get current directory... exiting" << endl;
            return EXIT_FAILURE;
        }
        string base = std::filesystem::path(configFile).parent_path().string();
        inputDir = string(cwd) + "/" + base + "/../input/";
        outputDir = string(cwd) + "/" + base + "/../OutDir/";
    }

    
    
    try {
        DigitizationRunner runner(configFile, inputDir, outputDir);
        runner.run();
    } catch (const std::exception& e) {
        cerr << "Error during digitization: " << e.what() << endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
