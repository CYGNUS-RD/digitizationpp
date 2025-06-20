/*
 * Copyright (C) 2025 CYGNO Collaboration
 *
 *
 * Author: Stefano Piacentini
 * Created in 2025
 *
 */

#pragma once

#include <string>
#include <map>
#include <vector>

/**
 * @class ConfigManager
 * @author Stefano Piacentini
 * @brief Handles loading and accessing configuration file for the digitization.
 *
 * @details Supports parsing of general settings, isotope mappings, and ion lists.
 */
class ConfigManager {
public:
    /**
     * @brief Loads configuration key-value pairs from a file.
     * @param configFile Path to the configuration file.
     * @return True if loaded successfully, false otherwise.
     */
    bool loadConfig(const std::string& configFile);

    /**
     * @brief Loads isotope conversion data (e.g. Z to ZAI mapping) from a file.
     * @param isotopeFile Path to the isotope definition file.
     */
    void loadIsotopes(const std::string& isotopeFile);

    /**
     * @brief Loads an ion list file into a 2D vector structure.
     * @param ionFile Path to the ion list file.
     * @param ionlist Reference to the container to populate.
     */
    void loadIonList(const std::string& ionFile, std::vector<std::vector<std::string>>& ionlist);

    /**
     * @brief Prints the current configuration options to stdout.
     */
    void printConfig() const;

    /**
     * @brief Retrieves the full map of configuration options.
     * @return Constant reference to the internal key-value map.
     */
    const std::map<std::string, std::string>& getOptions() const;

    /**
     * @brief Gets the isotope code (e.g. ZAI) for a given atomic number key.
     * @param key The atomic number or element symbol key.
     * @return The isotope string code, or empty string if not found.
     */
    std::string  getIsotope(const std::string& key) const;

    /**
     * @brief Retrieves a configuration value as a string.
     * @param key Configuration key.
     * @return String value associated with the key.
     */
    std::string  get(const std::string& key) const;

    /**
     * @brief Retrieves a configuration value and interprets it as a boolean.
     * @param key Configuration key.
     * @return True if the value is \"True\", \"true\", or \"1\".
     */
    bool     getBool(const std::string& key) const;

    /**
     * @brief Retrieves a configuration value and converts it to a double.
     * @param key Configuration key.
     * @return Double-precision value.
     */
    double getDouble(const std::string& key) const;

    /**
     * @brief Retrieves a configuration value and converts it to an integer.
     * @param key Configuration key.
     * @return Integer value.
     */
    int       getInt(const std::string& key) const;

private:
    std::map<std::string, std::string> options;  ///< Configuration key-value pairs
    std::map<std::string, std::string> isotopes; ///< Isotope lookup table
};
