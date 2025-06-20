/*
 * Copyright (C) 2025 CYGNO Collaboration
 *
 *
 * Author: Stefano Piacentini
 * Created in 2025
 *
 */

#pragma once

/**
 * @file Globals.h
 * @author Stefano Piacentini
 * @brief Declares global variables used throughout the digitization and simulation pipeline.
 */

/**
 * @brief Gain of the first GEM (Gas Electron Multiplier) stage.
 */
extern double GEM1_gain;

/**
 * @brief Gain of the second GEM stage.
 */
extern double GEM2_gain;

/**
 * @brief Gain of the third GEM stage.
 */
extern double GEM3_gain;

/**
 * @brief Extraction efficiency of the first GEM stage.
 */
extern double extraction_eff_GEM1;

/**
 * @brief Extraction efficiency of the second GEM stage.
 */
extern double extraction_eff_GEM2;

/**
 * @brief Extraction efficiency of the third GEM stage.
 */
extern double extraction_eff_GEM3;

/**
 * @brief Optical transfer efficiency factor (depends on aperture, demagnification, etc.).
 */
extern double omega;

/**
 * @brief Number of pixels along the X (horizontal) axis of the camera sensor.
 */
extern int x_pix;

/**
 * @brief Number of pixels along the Y (vertical) axis of the camera sensor.
 */
extern int y_pix;

/**
 * @brief Optical counts (ADU or electrons) per detected photon.
 */
extern double optcounts_per_photon;

/**
 * @brief Physical size of the camera sensor along the Y axis (in mm).
 */
extern double y_sensor_size;

/**
 * @brief Time to expose the whole camera sensor (ms).
 */
extern double readout_time;

/**
 * @brief Number of threads used in parallel processing (TBB-based).
 */
extern int num_threads;
