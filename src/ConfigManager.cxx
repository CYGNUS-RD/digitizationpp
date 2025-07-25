/*
 * Copyright (C) 2025 CYGNO Collaboration
 *
 *
 * Author: Stefano Piacentini
 * Created in 2025
 *
 */

#include "ConfigManager.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <filesystem>
#include <set>
#include <stdexcept>
#include "TFile.h"

bool ConfigManager::loadConfig(const std::string& configFile) {
    std::ifstream config(configFile);
    if (!config.is_open()) return false;

    std::string line;
    while (getline(config, line)) {
        line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
        line.erase(remove_if(line.begin(), line.end(), [](char c) { return c == '\''; }), line.end());

        if (line.empty() || line[0] == '#' || line[0] == '{' || line[0] == '}') continue;

        auto delim1 = line.find(":");
        auto delim2 = line.find(",");
        if (delim2 == std::string::npos) delim2 = line.size();

        std::string key = line.substr(0, delim1);
        std::string val = line.substr(delim1 + 1,  delim2 - delim1 - 1);
        
        options[key] = val;
    }
    return true;
}

void ConfigManager::loadIsotopes(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) return;

    std::string line;
    while (getline(file, line)) {
        line.erase(remove_if(line.begin(),line.end(),[] (char c){return isspace(c);}),line.end());
        line.erase(remove_if(line.begin(),line.end(),[] (char c){return c=='\'';}),line.end());
        if(line[0] == '#' || line.empty() || line[0] == '{' || line[0] == '}') continue;

        std::stringstream to_split(line.c_str());
        std::string element;
        
        int counter = 0;
        std::string key;
        std::string val;

        while(getline(to_split, element, ',')) {
            if(counter == 0) key = element;
            else if (counter == 1) val = element;
            else if (counter == 2) val += Form("%03d", std::stoi(element));
            counter++;
        }
        
        isotopes[key]=val;
        
    }
    
    return;
}


void ConfigManager::loadIonList(const std::string& filename, std::vector<std::vector<std::string>>& ionlist) {
    std::ifstream file(filename);
    if (!file.is_open()) return;

    std::string line;
    while (getline(file, line)) {
        line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
        if (line.find("ionlist") != std::string::npos) continue;

        std::stringstream outer(line);
        std::string block;
        while (getline(outer, block, '[')) {
            if(!block.empty()){
                block.erase(remove_if(block.begin(), block.end(), [](char c) { return c == ']'; }), block.end());
    
                std::stringstream inner(block);
                std::string item;
                std::vector<std::string> row;
                while (getline(inner, item, ',')) {
                    if (!item.empty()) row.push_back(item);
                }
                
                ionlist.push_back(row);
            }
        }
        break;
    }
    return;
}

void ConfigManager::printConfig() const {
    std::cout << "Configuration Options:\n";
    for (const auto& [key, val] : options) {
        std::cout << "  " << key << " = " << val << "\n";
    }
}


const std::map<std::string, std::string>& ConfigManager::getOptions() const {
    return options;
}

std::string ConfigManager::getIsotope(const std::string& key) const {
    auto it = isotopes.find(key);
    return (it != isotopes.end()) ? it->second : "";
}

std::string ConfigManager::get(const std::string& key) const {
    auto it = options.find(key);
    return (it != options.end()) ? it->second : "";
}

bool ConfigManager::getBool(const std::string& key) const {
    std::string val = get(key);
    return val == "True" || val == "true" || val == "1";
}

double ConfigManager::getDouble(const std::string& key) const {
    return std::stod(get(key));
}

int ConfigManager::getInt(const std::string& key) const {
    return std::stoi(get(key));
}

void ConfigManager::validateAxisMappings() {
    std::set<char> axes_found;

    // List of keys to check
    const std::vector<std::string> keys = {"MC_xaxis", "MC_yaxis", "MC_zaxis"};

    for (const std::string& key : keys) {
        if (options.find(key) == options.end()) {
            throw std::runtime_error("Missing required option: " + key);
        }

        const std::string& val = options.at(key);

        if (val.empty() || (val[0] != '-' && val.size() != 1) || (val[0] == '-' && val.size() != 2)) {
            throw std::runtime_error("Invalid axis format for option: " + key + " = '" + val + "'");
        }

        char axis_char = (val[0] == '-') ? val[1] : val[0];

        if (axis_char != 'x' && axis_char != 'y' && axis_char != 'z') {
            throw std::runtime_error("Invalid axis character in option: " + key + " = '" + val + "'");
        }

        axes_found.insert(axis_char);
    }

    if (axes_found.size() != 3) {
        throw std::runtime_error("Invalid axis mapping: Must include one 'x', one 'y', and one 'z'.");
    }

    return;
}