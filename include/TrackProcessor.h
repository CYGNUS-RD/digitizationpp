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
#include "TH2F.h"
#include <vector>
#include <string>
#include <map>

/**
 * @class TrackProcessor
 * @author Stefano Piacentini
 * @brief Executes the digitization of Monte Carlo hits including charge smearing, gain application, and vignetting.
 *
 * @details This class applies physical processes such as 2D/3D smearing, saturation, and vignetting
 * to simulate how an event would appear on a detector image.
 */
class TrackProcessor {
public:
    /**
     * @brief Constructor.
     * @param[in] config Reference to a ConfigManager containing digitization parameters.
     */
    TrackProcessor(ConfigManager& config);

    /**
     * @brief Applies full charge diffusion and saturation model to energy deposits.
     *
     * @param[in] x_hits_tr Vector of x coordinates of energy deposits.
     * @param[in] y_hits_tr Vector of y coordinates of energy deposits.
     * @param[in] z_hits_tr Vector of z coordinates of energy deposits.
     * @param[in] energy_hits_tr Energy deposited at each hit.
     * @param[in] energy Total event energy.
     * @param[in] NR_flag Flag indicating if the event is a nuclear recoil.
     * @param[in] image 2D image to be filled with the simulated event.
     */
    void computeWithSaturation(const std::vector<double>& x_hits_tr,
                               const std::vector<double>& y_hits_tr,
                               const std::vector<double>& z_hits_tr,
                               const std::vector<double>& energy_hits_tr,
                               float energy,
                               bool NR_flag,
                               std::vector<std::vector<double>>& image);

    /**
     * @brief Computes the image without applying saturation.
     *
     * @param[in] x_hits_tr Vector of x coordinates of energy deposits.
     * @param[in] y_hits_tr Vector of y coordinates of energy deposits.
     * @param[in] z_hits_tr Vector of z coordinates of energy deposits.
     * @param[in] energy_hits_tr Energy deposited at each hit.
     * @param[in] image 2D image to be filled with the simulated event.
     */
    void computeWithoutSaturation(const std::vector<double>& x_hits_tr,
                                  const std::vector<double>& y_hits_tr,
                                  const std::vector<double>& z_hits_tr,
                                  const std::vector<double>& energy_hits_tr,
                                  std::vector<std::vector<double>>& image);

    
    /**
     * @brief Applies vignetting correction to a 2D image using a vignette map.
     *
     * @param[in] image 2D pixel matrix to apply correction to.
     * @param[in] xpix Image width in pixels.
     * @param[in] ypix Image height in pixels.
     * @param[in] VignMap Vignetting map.
     */
    void TrackVignetting(std::vector<std::vector<double>>& image, int xpix, int ypix, const TH2F& VignMap);

private:
    ConfigManager& config; ///< Reference to configuration

    /**
     * @brief Performs 3D smearing of the ionization cloud in space.
     * @param[in] x_hits_tr Vector of x coordinates of energy deposits.
     * @param[in] y_hits_tr Vector of y coordinates of energy deposits.
     * @param[in] z_hits_tr Vector of z coordinates of energy deposits.
     * @param[in] energy_hits_tr Energy deposited at each hit.
     * @param[in] S3D_x Vector of x coordinates of smeared electrons.
     * @param[in] S3D_y Vector of y coordinates of smeared electrons.
     * @param[in] S3D_z Vector of z coordinates of smeared electrons.
     */
    void cloud_smearing3D(  const std::vector<double>& x_hits_tr,
                            const std::vector<double>& y_hits_tr,
                            const std::vector<double>& z_hits_tr,
                            const std::vector<double>& energy_hits_tr,
                            std::vector<float>& S3D_x,
                            std::vector<float>& S3D_y,
                            std::vector<float>& S3D_z);

     /**
     * @brief Applies Gaussian smearing to a single axis using local sigma and charge.
     * @param[in] axis_hit Vector of position along the considered axis of energy deposits.
     * @param[in] axis_sigma Vector of total diffusion sigma along the considered axis to apply to each energy deposits (it depends on z!)
     * @param[in] nel Number of electrons deposited at each hit.
     * @return Vector of smeared positions
     */
    std::vector<float> smear(const std::vector<double>& axis_hit,
                             const std::vector<double>& axis_sigma,
                             const std::vector<double>& nel);

    /**
     * @brief Parallel version of 3D smearing using axis-specific sigmas and charges.
     * @param[in] x_axis_hit Vector of x coordinates of energy deposits.
     * @param[in] y_axis_hit Vector of y coordinates of energy deposits.
     * @param[in] z_axis_hit Vector of z coordinates of energy deposits.
     * @param[in] x_axis_sigma Vector of total diffusion sigma along the x axis to apply to each energy deposits (it depends on z!)
     * @param[in] y_axis_sigma Vector of total diffusion sigma along the y axis to apply to each energy deposits (it depends on z!)
     * @param[in] z_axis_sigma Vector of total diffusion sigma along the z axis to apply to each energy deposits (it depends on z!)
     * @param[in] nel Number of electrons deposited at each hit.
     * @param[in] X Vector of x coordinates of smeared electrons.
     * @param[in] Y Vector of y coordinates of smeared electrons.
     * @param[in] Z Vector of z coordinates of smeared electrons.
     */
    void smear_parallel(const std::vector<double>& x_axis_hit,
                        const std::vector<double>& y_axis_hit,
                        const std::vector<double>& z_axis_hit,
                        const std::vector<double>& x_axis_sigma,
                        const std::vector<double>& y_axis_sigma,
                        const std::vector<double>& z_axis_sigma,
                        const std::vector<double>& nel,
                        std::vector<float>& X,
                        std::vector<float>& Y,
                        std::vector<float>& Z);

    /**
     * @brief Performs 2D smearing in case of no saturation (and therefore no voxelization).
     * @param[in] x_hits_tr Vector of x coordinates of energy deposits.
     * @param[in] y_hits_tr Vector of y coordinates of energy deposits.
     * @param[in] z_hits_tr Vector of z coordinates of energy deposits.
     * @param[in] energy_hits_tr Energy deposited at each hit.
     * @param[in] S2D_x Vector of x coordinates of smeared electrons.
     * @param[in] S2D_y Vector of y coordinates of smeared electrons.
     */
    void ph_smearing2D( const std::vector<double>& x_hits_tr,
                        const std::vector<double>& y_hits_tr,
                        const std::vector<double>& z_hits_tr,
                        const std::vector<double>& energy_hits_tr,
                        std::vector<float>& S2D_x,
                        std::vector<float>& S2D_y);

    /**
     * @brief Computes diffusion sigma values for given path lengths.
     */
    std::vector<double> compute_sigma(const double diff_const, const double diff_coeff, const std::vector<double>& dz);
    
    /**
     * @brief Applies GEM1 gain to primary electrons.
     */
    std::vector<double> NelGEM1(const std::vector<double>& N_ioniz_el);

    /**
     * @brief Applies GEM2 gain and charge losses to electrons based on z position.
     */
    std::vector<double> NelGEM2(const std::vector<double>& energyDep,const std::vector<double>& drift_l);

    /**
     * @brief Applies saturation function to a certain voxel.
     *
     * @param[in] nel Number of electrons in the voxel.
     * @param[in] A Saturation parameter A.
     * @param[in] beta Saturation parameter beta.
     * @return Saturated output value.
     */
    double Nph_saturation(int nel, double A, double beta);



};
