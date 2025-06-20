/*
 * Copyright (C) 2025 CYGNO Collaboration
 *
 *
 * Author: Stefano Piacentini
 * Created in 2025
 *
 */

#pragma once

#include <vector>
#include <string>
#include <filesystem>
#include <algorithm>


/**
 * @namespace Utils
 * @brief A collection of general-purpose utility functions for math, string processing, and file handling.
 * @author Stefano Piacentini
 */
namespace Utils {

    /**
     * @brief Converts a relative file path to an absolute path using std::filesystem.
     * @param relativePath Relative path to resolve.
     * @return Absolute file path as a string.
     */
    std::string resolvePath(const std::string& relativePath);

    
    /**
     * @brief Splits a string by a given delimiter into a vector of substrings.
     * @param input The input string.
     * @param delimiter The character used to split.
     * @return A vector of substrings.
     */
    std::vector<std::string> splitString(const std::string& input, char delimiter);

    /**
     * @brief Computes the angle in radians between two 3D vectors.
     * @param v1 First 3D vector.
     * @param v2 Second 3D vector.
     * @return Angle in radians.
     */
    double angleBetween(const std::vector<double>& v1, const std::vector<double>& v2);

    /**
     * @brief Computes the cross product of two 3D vectors.
     * @param a First 3D vector.
     * @param b Second 3D vector.
     * @return The cross product vector.
     */
    std::vector<double> crossProduct(const std::vector<double>& a, const std::vector<double>& b);

    /**
     * @brief Rotates a vector around a specified axis by a given angle (in radians).
     * @param vec The vector to rotate.
     * @param angle The rotation angle in radians.
     * @param axis The axis to rotate around (should be normalized).
     * @return The rotated vector.
     */
    std::vector<double> rotateByAngleAndAxis(const std::vector<double>& vec, double angle, const std::vector<double>& axis);

    /**
     * @brief Rounds a floating-point number up to the next even integer.
     * @param value The input value.
     * @return The nearest even integer greater than or equal to value.
     */
    double roundUpToEven(double value);

    /**
     * @brief Generates a std::vector of double numbers similar to Python's `numpy.arange`.
     * @param start Starting value of the sequence.
     * @param stop End value (not included in result).
     * @param step Increment between values.
     * @return A vector containing the sequence.
     */
    std::vector<double> arange(double start, double stop, double step);

}
