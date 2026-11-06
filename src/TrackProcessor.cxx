/*
 * Copyright (C) 2025 CYGNO Collaboration
 *
 *
 * Author: Stefano Piacentini
 * Created in 2025
 *
 */
 
#include "TrackProcessor.h"
#include "ConfigManager.h"
#include "Utils.h"
#include "Globals.h"
#include <cmath>
#include <chrono>
#include <numeric>
#include <random>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include "TRandom3.h"
#include <oneapi/tbb.h>
#include <oneapi/tbb/info.h>
#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/task_arena.h>
#include <oneapi/tbb/mutex.h>
#include <oneapi/tbb/enumerable_thread_specific.h>
#include <thread>

using namespace std;

TrackProcessor::TrackProcessor(ConfigManager& configMgr)
    : config(configMgr) {
        N_primaries_reaching_GEMs = -1;
    }

bool TrackProcessor::computeWithSaturation(const vector<double>& x_hits_tr,
                                           const vector<double>& y_hits_tr,
                                           const vector<double>& z_hits_tr,
                                           const vector<double>& energy_hits,
                                           const TH2F& VignMap,
                                           float energy,
                                           bool NR_flag,
                                           vector<vector<double>>& image,
                                           vector<PMTVoxel>& pmt_voxels)
{
    // vectorized smearing
    vector<float> S3D_x;
    vector<float> S3D_y;
    vector<float> S3D_z;
    
    // if there are no electrons on GEM3, just use empty image
    if (x_hits_tr.size() == 0) {
        N_primaries_reaching_GEMs = 0;
        return true;
    }
    // if there are electrons on GEM3, apply saturation effect
    else {
        // Sample primaries and GEM gains once for all input hits of this call.
        // nel[h] is the post-GEM2 population, not the primary count.
        // If the caller filtered hits, the primary metadata covers that subset.
        const vector<double> nel = NelGEM2(energy_hits, z_hits_tr);
        if (nel.empty()) return false;

        double OFF = 15;
        double OFFz = 10.;

        double xmin = (*min_element(x_hits_tr.begin(), x_hits_tr.end()))-OFF;
        double xmax = (*max_element(x_hits_tr.begin(), x_hits_tr.end()))+OFF;
        double ymin = (*min_element(y_hits_tr.begin(), y_hits_tr.end()))-OFF;
        double ymax = (*max_element(y_hits_tr.begin(), y_hits_tr.end()))+OFF;
        double zmin = (*min_element(z_hits_tr.begin(), z_hits_tr.end()))-OFFz;
        double zmax = (*max_element(z_hits_tr.begin(), z_hits_tr.end()))+OFFz;
        
        // FIXME: create a function for the saturation loop
        // This limits dense histogram cells, not total process memory.
        long int max_3Dhisto_volume=(long int)5e+7;  // 5e7 doubles: about 400 MB, excluding vector overhead.

        int xy_vox_scale = config.getInt("xy_vox_scale");
        
        double xbin_dim = config.getDouble("x_dim") / x_pix / xy_vox_scale;
        double ybin_dim = config.getDouble("y_dim") / y_pix / xy_vox_scale;
        double zbin_dim = config.getDouble("z_vox_dim");
        
        
        long int x_n_bin = Utils::roundUpToEven((xmax-xmin)/xbin_dim);
        long int y_n_bin = Utils::roundUpToEven((ymax-ymin)/ybin_dim);
        long int z_n_bin_MAX = Utils::roundUpToEven((zmax-zmin)/zbin_dim);
        
        // Preserve the existing energy-per-voxel heuristic: sparse map at low
        // occupancy, dense slabs above the threshold.
        bool map_algorithm = true;
        if(energy/(x_n_bin * y_n_bin * z_n_bin_MAX)>1e-6) map_algorithm = false;

        long int z_n_bin;
        
        // DEBUG
        //cout<<"size = "<<S3D_z.size()<<endl;
        //std::chrono::duration<double> dur=endsmear-startsmear;
        //std::cout << "Time smear " << dur.count() << " seconds" <<std::endl;
        
        double optA                 = config.getDouble("A");
        double optbeta              = config.getDouble("beta");
        double optphotons_per_el    = config.getDouble("photons_per_el");

        //Timing and debug variables
        long int size_tot=0;
        std::chrono::duration<double> dur_smear;
        std::chrono::duration<double> dur_criti;
        std::chrono::duration<double> dur_ampli;
        
        vector<vector<int>> hout(x_n_bin, vector<int>(y_n_bin, 0.0));
        
        if(map_algorithm) {
            // The sparse map combines all hit batches before saturation.
            // WID bounds hits per batch, not amplified electrons or map size.
            z_n_bin = z_n_bin_MAX;
            map<long int, int> hcmap;
            //long int L = x_n_bin-1;
            long int M = y_n_bin;
            long int N = z_n_bin;
            
            if((NR_flag==true && energy>100)||(NR_flag==false && energy>500)) {
                
                int WID = config.getInt("WID");
                if (WID <= 0) return false;
                int nparts = (x_hits_tr.size() + WID - 1)/WID;
                for(int part = 0; part < nparts; part++) {
                    cout<<"   part "<<part<<"/"<<nparts<<"..."<<endl;
                    int hit_tr_idx = part*WID;
                    int hit_tr_idx_up = min((part+1)*WID,(int)x_hits_tr.size());
                    S3D_x.clear();
                    S3D_y.clear();
                    S3D_z.clear();
                    
                    vector<double> x_hits_tr_i(x_hits_tr.begin() + hit_tr_idx, x_hits_tr.begin() + hit_tr_idx_up);
                    vector<double> y_hits_tr_i(y_hits_tr.begin() + hit_tr_idx, y_hits_tr.begin() + hit_tr_idx_up);
                    vector<double> z_hits_tr_i(z_hits_tr.begin() + hit_tr_idx, z_hits_tr.begin() + hit_tr_idx_up);
                    vector<double> nel_i(nel.begin() + hit_tr_idx, nel.begin() + hit_tr_idx_up);
                
                    auto startsmear = std::chrono::steady_clock::now();
                    cloud_smearing3D_from_nel(x_hits_tr_i, y_hits_tr_i, z_hits_tr_i, nel_i, S3D_x, S3D_y, S3D_z);
                    auto endsmear = std::chrono::steady_clock::now();
                    dur_smear=dur_smear+endsmear-startsmear;
                    size_tot+=S3D_z.size();
                    
                    vector<size_t> indices(S3D_z.size());
                    // Fill indices with 0, 1, 2, ..., numbers.size() - 1
                    iota(indices.begin(), indices.end(), 0);
                    
                    auto startcriti = std::chrono::steady_clock::now();
                    // THIS IS THE COMPUTATIONALLY EXPENSIVE PART
                    //Parallel version
                    //Types of Mutex: spin -> unfair, unscalable (better for locking things that require little time)
                    //queue ->> fair, scalable
                    //tbb::queuing_mutex myMutex;
                    tbb::spin_mutex myMutex;
                    // Run the default parallelism
                    tbb::task_arena arena(num_threads);
                    arena.execute([&]{
                        parallel_for(tbb::blocked_range<size_t>(0, indices.size(),max<size_t>(1, indices.size()/max(1, num_threads))),
                            [&] (auto &range) {
                            map<long int, int> paralmap;
                            //cout << "Thread ID: " << this_thread::get_id() << endl;
                            for(auto ihit=range.begin();ihit<range.end();++ihit)
                            {
                                long int xx = floor((S3D_x[ihit]- xmin)/ xbin_dim);
                                long int yy = floor((S3D_y[ihit]- ymin)/ ybin_dim);
                                long int zz = floor((S3D_z[ihit]- zmin)/ zbin_dim);
                                long int map_index = xx * (M * N) + yy * N + zz;
                                
                                if (paralmap.find(map_index) == paralmap.end()) paralmap[map_index] = 1;
                                else paralmap[map_index] += 1;
                            }
                            //cout<<"Mutex\n";
                            //cout << "Thread ID: " << this_thread::get_id() << endl;
                            //tbb::queuing_mutex::scoped_lock myLock(myMutex);
                            tbb::spin_mutex::scoped_lock myLock(myMutex);	//maybe better to use oneapi::tbb::spin_mutex to use spi_mutex since the critical operation is small
                            //merge of paralmap in hcmap
                            hcmap.merge(paralmap);
                            for(auto it=paralmap.begin();it!=paralmap.end();++it) hcmap[it->first] = hcmap[it->first] + paralmap[it->first];               //Here auto would be std:map<long int,int>::iterator
                        });
                    });
                    
                    //End parallel version

                    //Sequencial version
                    //for_each(indices.begin(), indices.end(), [&](int ihit) {
                    //    long int xx = floor((S3D_x[ihit]- xmin)/ xbin_dim);
                    //    long int yy = floor((S3D_y[ihit]- ymin)/ ybin_dim);
                    //    long int zz = floor((S3D_z[ihit]- zmin)/ zbin_dim);
                    //    long int map_index = xx * (M * N) + yy * N + zz;
                    //    
                    //    if (hcmap.find(map_index) == hcmap.end()) hcmap[map_index] = 1;
                    //    else hcmap[map_index] += 1;
                    //    
                    //});
                    //End sequencial version
                    auto endcriti = std::chrono::steady_clock::now();
                    dur_criti=dur_criti+endcriti-startcriti;
                }
                auto startampli = std::chrono::steady_clock::now();
                for(auto const& [key, val] : hcmap) {
                    int zz = key % N ;
                    int yy = ((key - zz) / N) % M;
                    int xx = (key - yy * N - zz) / N / M;
                    // optbeta is multiplied by factors to normalize on volume chosen
                    double n_ph_val = Nph_saturation(val, optA, optbeta*xy_vox_scale*xy_vox_scale*0.1/zbin_dim);
                    hout[xx][yy]+=n_ph_val;

                    // --- PMT Simulation Integration ---
                    if (n_ph_val > 0) {
                        PMTVoxel voxel;
                        // Convert bin indices back to physical coordinates (mm)
                        voxel.x = (xx * xbin_dim) + xmin;
                        voxel.y = (yy * ybin_dim) + ymin;
                        voxel.z = (zz * zbin_dim) + zmin; 
                        
                        voxel.n_photons = n_ph_val;
                        pmt_voxels.push_back(voxel);
                    }
                }
                auto endampli = std::chrono::steady_clock::now();
                dur_ampli=dur_ampli+endampli-startampli;

            } else {

                auto startsmear = std::chrono::steady_clock::now();
                cloud_smearing3D_from_nel(x_hits_tr, y_hits_tr, z_hits_tr, nel, S3D_x, S3D_y, S3D_z);
                auto endsmear = std::chrono::steady_clock::now();
                dur_smear=dur_smear+endsmear-startsmear;
                size_tot+=S3D_z.size();

                vector<size_t> indices(S3D_z.size());
                // Fill indices with 0, 1, 2, ..., numbers.size() - 1
                iota(indices.begin(), indices.end(), 0);
                
                auto startcriti = std::chrono::steady_clock::now();
                //Parallel version
                //Types of Mutex: spin -> unfair, unscalable (better for locking things that require little time)
                //queue ->> fair, scalable
                //tbb::queuing_mutex myMutex;
                tbb::spin_mutex myMutex;
                // Run the default parallelism
                tbb::task_arena arena(num_threads);
                arena.execute([&]{
                    parallel_for(tbb::blocked_range<size_t>(0, indices.size(),max<size_t>(1, indices.size()/max(1, num_threads))),
                        [&] (auto &range) {
                        map<long int, int> paralmap;
                        //cout << "Thread ID: " << this_thread::get_id() << endl;
                        for(auto ihit=range.begin();ihit<range.end();++ihit)
                        {
                            long int xx = floor((S3D_x[ihit]- xmin)/ xbin_dim);
                            long int yy = floor((S3D_y[ihit]- ymin)/ ybin_dim);
                            long int zz = floor((S3D_z[ihit]- zmin)/ zbin_dim);
                            long int map_index = xx * (M * N) + yy * N + zz;
                            
                            if (paralmap.find(map_index) == paralmap.end()) paralmap[map_index] = 1;
                            else paralmap[map_index] += 1;
                        }
                        //cout<<"Mutex\n";
                        //cout << "Thread ID: " << this_thread::get_id() << endl;
                        //tbb::queuing_mutex::scoped_lock myLock(myMutex);
                        tbb::spin_mutex::scoped_lock myLock(myMutex);	//maybe better to use oneapi::tbb::spin_mutex to use spi_mutex since the critical operation is small
                        //merge of paralmap in hcmap
                        hcmap.merge(paralmap);
                        for(auto it=paralmap.begin();it!=paralmap.end();++it) hcmap[it->first] = hcmap[it->first] + paralmap[it->first];               //Here auto would be std:map<long int,int>::iterator
                    });
                });
                
                //End parallel version

                //Sequencial version
                //for_each(indices.begin(), indices.end(), [&](int ihit) {
                //    long int xx = floor((S3D_x[ihit]- xmin)/ xbin_dim);
                //    long int yy = floor((S3D_y[ihit]- ymin)/ ybin_dim);
                //    long int zz = floor((S3D_z[ihit]- zmin)/ zbin_dim);
                //    long int map_index = xx * (M * N) + yy * N + zz;
                    
                //    if (hcmap.find(map_index) == hcmap.end()) hcmap[map_index] = 1;
                //    else hcmap[map_index] += 1;
                    
                //});
                //End sequencial version
                auto endcriti = std::chrono::steady_clock::now();
                dur_criti=dur_criti+endcriti-startcriti;

                auto startampli = std::chrono::steady_clock::now();
                for(auto const& [key, val] : hcmap) {
                    int zz = key % N ;
                    int yy = ((key - zz) / N) % M;
                    int xx = (key - yy * N - zz) / N / M;
                    // optbeta is multiplied by factors to normalize on volume chosen
                    double n_ph_val = Nph_saturation(val, optA, optbeta*xy_vox_scale*xy_vox_scale*0.1/zbin_dim);
                    hout[xx][yy]+=n_ph_val;
                    
                    // --- PMT Simulation Integration ---
                    if (n_ph_val > 0) {
                        PMTVoxel voxel;
                        // Convert bin indices back to physical coordinates (mm)
                        voxel.x = (xx * xbin_dim) + xmin;
                        voxel.y = (yy * ybin_dim) + ymin;
                        voxel.z = (zz * zbin_dim) + zmin; 
                        
                        voxel.n_photons = n_ph_val;
                        pmt_voxels.push_back(voxel);
                    }
                }
                auto endampli = std::chrono::steady_clock::now();
                dur_ampli=dur_ampli+endampli-startampli;
                
            }

        } else {

            // Slabs are groups of bins of one global z grid.
            // hc has exactly these dimensions (no unused +1 borders).
            const long int xy_volume = x_n_bin * y_n_bin;
            if (xy_volume <= 0 || xy_volume > max_3Dhisto_volume) {
                cerr << "A single xy plane exceeds max_3Dhisto_volume." << endl;
                return false;
            }
            const long int slab_bins = max_3Dhisto_volume / xy_volume;
            const auto sigma_x = compute_sigma(config.getDouble("diff_const_sigma0T"), config.getDouble("diff_coeff_T"), z_hits_tr);
            const auto sigma_z = compute_sigma(config.getDouble("diff_const_sigma0L"), config.getDouble("diff_coeff_L"), z_hits_tr);
            // One integer counter per hit; no full-event coordinate cloud.
            // NelGEM2 already rounded these populations to integer values.
            vector<long long> nel_remaining(nel.size());
            for (size_t h = 0; h < nel.size(); ++h) {
                if (!isfinite(nel[h]) || nel[h] < 0 || nel[h] >= static_cast<double>(numeric_limits<long long>::max()) ||
                    !isfinite(sigma_x[h]) || !isfinite(sigma_z[h])) return false;
                nel_remaining[h] = static_cast<long long>(nel[h]);
            }

            // Adapter to use an efficient binomial with the existing ROOT RNG.
            // This adapter owns no RNG state: each call advances gRandom.
            // The standard distribution expects an integer-producing callable.
            // 4294967296 is 2^32; scaling Rndm() supplies 32-bit integers.
            // min/max declare the adapter range; operator() makes random() callable.
            struct RootEngine {
                using result_type = uint32_t;
                static constexpr result_type min() { return 0; }
                static constexpr result_type max() { return numeric_limits<result_type>::max(); }
                result_type operator()() { return static_cast<result_type>(gRandom->Rndm() * 4294967296.0); }
            } random;
            // Draw a count, not coordinates: Binomial(n, probability).
            // Handle certain outcomes without consuming random numbers.
            auto draw = [&](long long n, double probability) -> long long {
                if (n == 0 || probability <= 0.) return 0;
                if (probability >= 1.) return n;
                return binomial_distribution<long long>(n, probability)(random);
            };
            // P(Z >= edge) for a Gaussian with the given mean and sigma.
            // erfc evaluates the upper tail directly instead of subtracting
            // a CDF close to one. A zero-width Gaussian is a point mass.
            auto survival = [](double edge, double mean, double sigma) {
                if (sigma == 0.) return mean >= edge ? 1. : 0.;
                return 0.5 * erfc((edge - mean) / (sqrt(2.) * sigma));
            };

            // Remove the lower Gaussian tail once; the upper tail is
            // whatever remains after the final bin. Neither is forced inside.
            for (size_t h = 0; h < nel.size(); ++h)
                nel_remaining[h] = draw(nel_remaining[h], survival(zmin, z_hits_tr[h], sigma_z[h]));

            // Persistent streams for logical workers: initialized once, never
            // reset between bins/blocks/slabs. This avoids repeating clouds.
            // Each logical worker owns one stream, regardless of the OS thread
            // TBB schedules it on. Seeding is serial; worker draws are disjoint.
            const int workers = max(1, num_threads);
            vector<TRandom3> paralrandom;
            paralrandom.reserve(workers);
            for (int w = 0; w < workers; ++w)
                paralrandom.emplace_back(1u + gRandom->Integer(4294967295u));
            tbb::task_arena arena(workers);

            // Two float arrays of this length use about 512 KiB together.
            // Unlike WID (previous implementation), this also bounds a single highly amplified hit.
            const long long max_smear_electrons = 65536;

            for (long int first_z = 0; first_z < z_n_bin_MAX; first_z += slab_bins) {
                // first_z is global; zz below is local to this slab.
                // The last slab may contain fewer bins.
                z_n_bin = min(slab_bins, z_n_bin_MAX - first_z);
                vector<vector<vector<double>>> hc(x_n_bin,
                    vector<vector<double>>(y_n_bin, vector<double>(z_n_bin, 0.)));

                for (long int zz = 0; zz < z_n_bin; ++zz) {
                    const double low = zmin + (first_z + zz) * zbin_dim;
                    const double high = zmin + (first_z + zz + 1) * zbin_dim;
                    for (size_t h = 0; h < nel.size(); ++h) {
                        // Conditional binomial: assign each electron to at most
                        // one z bin. This is Gaussian smearing directly in bins.
                        const double qlow = survival(low, z_hits_tr[h], sigma_z[h]);
                        const double qhigh = survival(high, z_hits_tr[h], sigma_z[h]);
                        // Remaining electrons are conditioned on Z >= low:
                        // P(low <= Z < high | Z >= low) = 1 - qhigh/qlow.
                        // Clamp floating-point roundoff. If qlow underflows to
                        // zero, the fallback assigns any numerical residue here.
                        const double probability = qlow > 0. ? clamp(1. - qhigh / qlow, 0., 1.) : 1.;
                        long long count = draw(nel_remaining[h], probability);
                        // Remove the assigned population before proceeding:
                        // the same electron cannot be assigned to another bin.
                        nel_remaining[h] -= count;
                        while (count > 0) {
                            const long long block = min(count, max_smear_electrons);
                            auto startsmear = chrono::steady_clock::now();
                            // Bounded x/y cloud; z is already a bin.
                            S3D_x.resize(block); S3D_y.resize(block);
                            // Every element is overwritten, even when resize
                            // reuses capacity. Workers write disjoint intervals.
                            arena.execute([&] {
                                tbb::parallel_for(tbb::blocked_range<int>(0, workers, 1), [&](auto& range) {
                                    for (int w = range.begin(); w < range.end(); ++w) {
                                        const long long begin = block * w / workers;
                                        const long long end = block * (w + 1) / workers;
                                        for (long long j = begin; j < end; ++j) {
                                            S3D_x[j] = paralrandom[w].Gaus(x_hits_tr[h], sigma_x[h]);
                                            S3D_y[j] = paralrandom[w].Gaus(y_hits_tr[h], sigma_x[h]);
                                        }
                                    }
                                });
                            });
                            dur_smear += chrono::steady_clock::now() - startsmear;
                            size_tot += block;
                            // arena.execute has completed: accumulate serially
                            // so histogram increments do not race. Gaussian x/y
                            // tails outside the grid are discarded explicitly.
                            auto startcriti = chrono::steady_clock::now();
                            for (size_t j = 0; j < S3D_x.size(); ++j) {
                                const long int xx = floor((S3D_x[j] - xmin) / xbin_dim);
                                const long int yy = floor((S3D_y[j] - ymin) / ybin_dim);
                                if (xx >= 0 && xx < x_n_bin && yy >= 0 && yy < y_n_bin)
                                    hc[xx][yy][zz] += 1.;
                            }
                            dur_criti += chrono::steady_clock::now() - startcriti;
                            count -= block;
                        }
                    }
                }
                // Apply the nonlinear response to complete voxel populations:
                // saturating separate hit/block contributions would be incorrect.
                auto startampli = chrono::steady_clock::now();
                for (long int xx = 0; xx < x_n_bin; ++xx)
                    for (long int yy = 0; yy < y_n_bin; ++yy)
                        for (long int zz = 0; zz < z_n_bin; ++zz) {
                            if (hc[xx][yy][zz] == 0.) continue;
                            // Nph_saturation currently accepts int. If this single voxel contains more electrons than the maximum that can be represented by a double, stop the digitization
                            if (hc[xx][yy][zz] > numeric_limits<int>::max()) return false;

                            const double n_ph_val = Nph_saturation(static_cast<int>(hc[xx][yy][zz]), optA, optbeta * xy_vox_scale * xy_vox_scale * 0.1 / zbin_dim);
                            hout[xx][yy] += n_ph_val;
                            
                            if (n_ph_val > 0.) {
                                PMTVoxel voxel;
                                voxel.x = xx * xbin_dim + xmin;
                                voxel.y = yy * ybin_dim + ymin;
                                // Keep the existing lower-bin-edge convention,
                                // using the global z index for PMT coordinates.
                                voxel.z = (first_z + zz) * zbin_dim + zmin;
                                voxel.n_photons = n_ph_val;
                                pmt_voxels.push_back(voxel);
                            }
                        }
                dur_ampli += chrono::steady_clock::now() - startampli;
            }
        }

        cout<<"size = "<<size_tot<<endl;
        std::cout << "Time smear " << dur_smear.count() << " seconds" <<std::endl;
        std::cout << "Time Critical " << dur_criti.count() << " seconds" <<std::endl;
        std::cout << "Time ampli " << dur_ampli.count() << " seconds" <<std::endl;
        
        // Padding [Now before camera digitization because needed for VignMap]
        // FIXME: Write a function padding()
        //Define a translation vector
        
        // DEBUG
        // int sommatotale =0;
        
        // sommatotale = accumulate(hout.cbegin(), hout.cend(), 0, [](auto sum, const auto& row) {
        //             return accumulate(row.cbegin(), row.cend(), sum);
        //         });
        
        // cout<<"DEBUG: INT = "<<sommatotale<<endl;

        int x_center_cloud=(int)round(((xmax+xmin)/2.)/xbin_dim);
        int y_center_cloud=(int)round(((ymax+ymin)/2.)/ybin_dim);

        //DEBUG
        //cout<<"cc   = ("<<x_center_cloud<<","<<y_center_cloud<<")"<<endl;
        vector<int> translation = {x_center_cloud, y_center_cloud};
        // Calculate the center position of the original array in the padded array
        vector<int> center = {(int)(static_cast<double>(x_pix)*xy_vox_scale/2.)+translation[0],
                              (int)(static_cast<double>(y_pix)*xy_vox_scale/2.)+translation[1]
                             };
                             
        // DEBUG
        //cout<<"c   = ("<<center[0]<<","<<center[1]<<")"<<endl;
        
        // cout<<"Center: "<<center[0]<<", "<<center[1]<<endl;
        const int x_origin = center[0] - (int)hout.size()/2;
        const int y_origin = center[1] - (int)hout[0].size()/2;
        int x_start = max(0, x_origin);
        int y_start = max(0, y_origin);
        int x_end   = min(x_pix*xy_vox_scale, x_origin + (int)hout.size());
        int y_end   = min(y_pix*xy_vox_scale, y_origin + (int)hout[0].size());
        // cout<<"PADDING ["<<x_start<<":"<<x_end<<","<<y_start<<":"<<y_end<<"]"<<endl;

        // Adding voxels to the same pixel before applying camera response
        for(int xx=x_start; xx<x_end; xx++){
            for(int yy=y_start; yy<y_end; yy++){
                image[xx/xy_vox_scale][yy/xy_vox_scale]+=hout[xx-x_origin][yy-y_origin];
            }
        }
        
        // Applying Poisson smearing + camera response
        for(unsigned int ii = 0; ii<image.size(); ii++) {
            for(unsigned int jj = 0; jj<image[0].size(); jj++) {
                if (image[ii][jj] != 0) {
                    // Vignetting
                    double _vigtmp = 1.0;
                    if(config.getBool("Vignetting")) {
                        _vigtmp = VignMap.GetBinContent(VignMap.GetXaxis()->FindBin(ii),
                                                        VignMap.GetYaxis()->FindBin(jj)
                                                          );
                    }
                    
                    // Generate number of photons reaching camera
                    const int _number_of_photons = gRandom->Poisson(static_cast<double>(image[ii][jj])
                                                                   * _vigtmp
                                                                   * omega
                                                                   * camera_quantum_efficiency
                                                                   * optphotons_per_el);
                    
                    // Average number of electrons
                    const double _mu = static_cast<double>(_number_of_photons) * optcounts_per_photon;

                    image[ii][jj] = static_cast<int>(_mu);
                    // Fluctuations of electrons // not needed if we sum random pedestal!!!
                    // if(_mu!=0) {
                    //     const double _smeared = gRandom-> Gaus(_mu, camera_electron_rms*optcounts_per_photon);

                    //     // If negative -> zero, but should happen very rarely
                    //     const double _nonneg  = (_smeared < 0.0) ? 0.0 : _smeared;
    
                    //     // Return integer number
                    //     image[ii][jj] =  static_cast<int>(std::lround(_nonneg));
                    // } else {
                    //     image[ii][jj] = 0;
                    // }
                }
            }
        }
        
    }
    return true;

}

// Sparse-path helper: smear an already sampled post-GEM2 population.
// The caller supplies one nel entry per hit. Unlike cloud_smearing3D, this
// method never samples charge. Clear outputs because smear_parallel appends.
void TrackProcessor::cloud_smearing3D_from_nel(const vector<double>& x_hits_tr,
                                              const vector<double>& y_hits_tr,
                                              const vector<double>& z_hits_tr,
                                              const vector<double>& nel,
                                              vector<float>& S3D_x,
                                              vector<float>& S3D_y,
                                              vector<float>& S3D_z) {
    S3D_x.clear();
    S3D_y.clear();
    S3D_z.clear();

    const auto sigma_x = compute_sigma(config.getDouble("diff_const_sigma0T"), config.getDouble("diff_coeff_T"), z_hits_tr);
    const auto sigma_y = sigma_x;
    const auto sigma_z = compute_sigma(config.getDouble("diff_const_sigma0L"), config.getDouble("diff_coeff_L"), z_hits_tr);
    smear_parallel(x_hits_tr, y_hits_tr, z_hits_tr, sigma_x, sigma_y, sigma_z,
                   nel, S3D_x, S3D_y, S3D_z);
}

// Legacy energy-based entry point: this still samples charge internally.
// computeWithSaturation deliberately uses the from_nel helper instead.
bool TrackProcessor::cloud_smearing3D(const vector<double>& x_hits_tr,
                                      const vector<double>& y_hits_tr,
                                      const vector<double>& z_hits_tr,
                                      const vector<double>& energy_hits,
                                      vector<float>& S3D_x,
                                      vector<float>& S3D_y,
                                      vector<float>& S3D_z) {

    vector<double> nel = NelGEM2(energy_hits, z_hits_tr);
    if (nel.empty()) {
        std::cerr << "Warning: TrackProcessor::cloud_smearing3D: skipping this track because of empty result from TrackProcessor::NelGEM2."<<std::endl;
        return false;
    }
    
    //DEBUG
    //for(unsigned int i=0; i<nel.size(); i++) {
    //    cout<<nel[i]<<"\n";
    //}

    vector<double> sigma_x = compute_sigma(config.getDouble("diff_const_sigma0T"), config.getDouble("diff_coeff_T"), z_hits_tr);
    vector<double> sigma_y = compute_sigma(config.getDouble("diff_const_sigma0T"), config.getDouble("diff_coeff_T"), z_hits_tr);
    vector<double> sigma_z = compute_sigma(config.getDouble("diff_const_sigma0L"), config.getDouble("diff_coeff_L"), z_hits_tr);

    //Here this is the slowest part
    //Sequential
    /*S3D_x = smear(x_hits_tr, sigma_x, nel);
    S3D_y = smear(y_hits_tr, sigma_y, nel);
    S3D_z = smear(z_hits_tr, sigma_z, nel);*/
    //Parallel
    smear_parallel(x_hits_tr,y_hits_tr,z_hits_tr,sigma_x,sigma_y,sigma_z,nel,S3D_x,S3D_y,S3D_z);
    // DEBUG
    //for(unsigned int i=0; i<S3D_x.size(); i++) {
    //    cout<<S3D_x[i]<<endl;
    //}

    return true;
}


vector<float> TrackProcessor::smear(const vector<double>& axis_hit, const vector<double>& axis_sigma, const vector<double>& nel) {

    long int nelsum = accumulate(nel.begin(), nel.end(), (long int)0);
    
    vector<float> X;
    X.reserve(nelsum);
    
    // Create a vector of indices where each index i is repeated nel[i] times
    vector<long int> indices(nelsum);
    vector<long int> positions(axis_hit.size() + 1, 0);
    
    // Compute cumulative sum of nel to determine positions
    partial_sum(nel.begin(), nel.end(), positions.begin() + 1);

    // Fill the indices vector
    for_each(positions.begin(), positions.end() - 1, [&, i = 0](long int pos) mutable {
        fill(indices.begin() + pos, indices.begin() + positions[i + 1], i);
        ++i;
    });
    
    // Fill X with Gaussian-distributed values based on axis_hit and axis_sigma
    //This here is the slowest part
    transform(indices.begin(), indices.end(), back_inserter(X), [&](int i) {
        return (float)gRandom->Gaus(axis_hit[i], axis_sigma[i]);
    });
    
    return X;
}

void TrackProcessor::smear_parallel(const vector<double>& x_axis_hit,const vector<double>& y_axis_hit,const vector<double>& z_axis_hit,const vector<double>& x_axis_sigma,const vector<double>& y_axis_sigma,const vector<double>& z_axis_sigma,const vector<double>& nel,vector<float>& X,vector<float>& Y,vector<float>& Z)
{
    long int nelsum = accumulate(nel.begin(), nel.end(), (long int)0);
    
    X.reserve(nelsum);
    Y.reserve(nelsum);
    Z.reserve(nelsum);

    // Fixed random seed
    bool fixedSeed = config.getBool("fixed_seed");
    
    // Create a vector of indices where each index i is repeated nel[i] times
    vector<long int> indices(nelsum);
    vector<long int> positions(x_axis_hit.size() + 1, 0);
    
    // Compute cumulative sum of nel to determine positions
    partial_sum(nel.begin(), nel.end(), positions.begin() + 1);

    // Fill the indices vector
    for_each(positions.begin(), positions.end() - 1, [&, i = 0](long int pos) mutable {
        fill(indices.begin() + pos, indices.begin() + positions[i + 1], i);
        ++i;
    });
    
    // Fill X with Gaussian-distributed values based on axis_hit and axis_sigma
    //This here is the slowest part
    //auto startstep4 = std::chrono::steady_clock::now();

    tbb::spin_mutex myMutex;
    // Run the default parallelism
    tbb::task_arena arena(num_threads);
    arena.execute([&]{
        parallel_for(tbb::blocked_range<size_t>(0, indices.size(),max<size_t>(1, indices.size()/max(1, num_threads))),
            [&] (auto &range) {
            vector<float> x_paralvec,y_paralvec,z_paralvec;
            TRandom3 paralrandom;
            myMutex.lock();
            //deterministic seed from chunk start
            if(fixedSeed) paralrandom.SetSeed(3+static_cast<UInt_t>(range.begin()));
            else paralrandom.SetSeed(floor(gRandom->Rndm()*10000));
            myMutex.unlock();
            for(auto iterator=range.begin();iterator<range.end();++iterator)
            {
                int index = indices[iterator];
                x_paralvec.push_back((float)paralrandom.Gaus(x_axis_hit[index], x_axis_sigma[index]));
                y_paralvec.push_back((float)paralrandom.Gaus(y_axis_hit[index], y_axis_sigma[index]));
                z_paralvec.push_back((float)paralrandom.Gaus(z_axis_hit[index], z_axis_sigma[index]));
            }
            //cout<<"Mutex\n";
            //cout << "Thread ID: " << this_thread::get_id() << endl;
            //tbb::queuing_mutex::scoped_lock myLock(myMutex);
            tbb::spin_mutex::scoped_lock myLock(myMutex);	//maybe better to use oneapi::tbb::spin_mutex to use spi_mutex since the critical operation is small
            //merge of paralvec in X
            X.insert(X.end(),x_paralvec.begin(),x_paralvec.end());
            Y.insert(Y.end(),y_paralvec.begin(),y_paralvec.end());
            Z.insert(Z.end(),z_paralvec.begin(),z_paralvec.end());
        });
    });    

    //auto endstep4 = std::chrono::steady_clock::now();
    //cout<<X.size()<<endl;
    //cout<<Y.size()<<endl;
    //cout<<Z.size()<<endl;

    //std::chrono::duration<double> dur=endstep4-startstep4;
    //std::cout << "Slowest Time smear part " << dur.count() << " seconds" <<std::endl;

    return ;
}

vector<double> TrackProcessor::compute_sigma(const double diff_const, const double diff_coeff, const vector<double>& dz) {
    vector<double> sigmas;
    transform(dz.begin(), dz.end(), back_inserter(sigmas), [&] (double a) {
        return sqrt(diff_const + diff_coeff * a /10.0);
    });
    return sigmas;
}


vector<double> TrackProcessor::NelGEM2(const vector<double>& energyDep, const vector<double>& drift_l) {
    vector<double> n_ioniz_el_ini;
    double opt_pot=config.getDouble("ion_pot");
    transform(energyDep.begin(),energyDep.end(),back_inserter(n_ioniz_el_ini), [&] (double a) { return a/opt_pot;});

    // Old version: z_extra was added here
    // vector<double> drift_l;
    // int opt_gem=config.getDouble("z_extra");
    // transform(z_hit.begin(),z_hit.end(),back_inserter(drift_l), [&] (double a) { return a+opt_gem;});

    // Check if hits with negative drift length
    vector<bool> negatives(drift_l.size());
    transform(drift_l.begin(), drift_l.end(), negatives.begin(),
              [](double x) { return x < 0; });
    if (any_of(negatives.begin(), negatives.end(), [](bool b) { return b; })) {
        std::cerr << "Warning: TrackProcessor::NelGEM2: track contains hits with negative drift length."<<std::endl;
        return {};  // empty vector signals "skip"
    }
    
    vector<double> n_ioniz_el_mean(n_ioniz_el_ini.size(), 0.0);
    
    double optabsorption_l=config.getDouble("absorption_l");
    for(unsigned int i=0;i<n_ioniz_el_mean.size();i++) n_ioniz_el_mean[i]=abs(n_ioniz_el_ini[i]*exp(-drift_l[i]/optabsorption_l));
    
    vector<double> n_ioniz_el(n_ioniz_el_ini.size(), 0);
    transform(n_ioniz_el_mean.begin(), n_ioniz_el_mean.end(), n_ioniz_el.begin(), [&] (double a) {
        return gRandom->Poisson(a);
    });

    // Primary metadata is recorded before GEM multiplication.
    // The saturation path calls this once for its complete input hit set.
    N_primaries_reaching_GEMs = std::accumulate(n_ioniz_el.begin(), n_ioniz_el.end(), 0.0);
    std::cout << "Number of primaries reaching GEMs: " << N_primaries_reaching_GEMs << std::endl;
    
    
    // total number of secondary electrons considering the gain in the 2nd GEM foil
    vector<double> n_tot_el = NelGEM1(n_ioniz_el);
    transform(n_tot_el.begin(), n_tot_el.end(), n_tot_el.begin(), [&] (double a) {
        return round(a * GEM2_gain * extraction_eff_GEM2);
    });
    
    return n_tot_el;
    
}


vector<double> TrackProcessor::NelGEM1(const vector<double>& N_ioniz_el) {
    
    vector<double> n_tot_el(N_ioniz_el.size(), 0);
    
    for(unsigned int i = 0; i<N_ioniz_el.size(); i++) {
        for(int j = 0; j<(int)round(N_ioniz_el[i]); j++) {
            double nsec = gRandom->Exp(GEM1_gain) * extraction_eff_GEM1;
            n_tot_el[i] += nsec;
        }
    }

    //DEBUG
    //for(unsigned int i=0; i<n_tot_el.size(); i++) {
    //    cout<<n_tot_el[i]<<"\n";
    //}
    
    return n_tot_el;
}


double TrackProcessor::Nph_saturation(int nel, double A, double beta) {
    return nel * A * GEM3_gain / (1.0 + beta * GEM3_gain * nel);
}


bool TrackProcessor::ph_smearing2D( const vector<double>& x_hits_tr,                        //Function still to be corrected
                                    const vector<double>& y_hits_tr,
                                    const vector<double>& z_hits_tr,
                                    const vector<double>& energy_hits,
                                    vector<float>& S2D_x,
                                    vector<float>& S2D_y) {

    // Electrons in GEM2
    vector<double> nel = NelGEM2(energy_hits, z_hits_tr);
    if (nel.empty()) {
        std::cerr << "Warning: TrackProcessor::ph_smearing2D: skipping this track because of empty result from TrackProcessor::NelGEM2."<<std::endl;
        return false;
    }

    double optphotons_per_el    = config.getDouble("photons_per_el");
    double optA                 = config.getDouble("A");
    // Photons in GEM3 (the factor A is added to be able to compare saturated and non-saturated results)
    vector<double> nph;
    transform(nel.begin(), nel.end(), back_inserter(nph), [&](double nel_i) {
        return nel_i * optA * GEM3_gain * omega * optphotons_per_el * optcounts_per_photon;
    });

    vector<double> sigma_xy = compute_sigma(config.getDouble("diff_const_sigma0T"), config.getDouble("diff_coeff_T"), z_hits_tr);

    S2D_x = smear(x_hits_tr, sigma_xy, nph);
    S2D_y = smear(y_hits_tr, sigma_xy, nph);

    return true;
}

void TrackProcessor::TrackVignetting(vector<vector<double>>& image, int xpix, int ypix, const TH2F & VignMap) {
    
    for(int xx = 0; xx < xpix; xx++) {
        for(int yy = 0; yy < ypix; yy++) {
            if(image[xx][yy] != 0) {
                image[xx][yy]=round(image[xx][yy] *
                                    VignMap.GetBinContent(VignMap.GetXaxis()->FindBin(xx),
                                                          VignMap.GetYaxis()->FindBin(yy)
                                                          )
                                    );
            }
        }
    }
    
    return;
}


int TrackProcessor::ApplyExposureCut(vector<vector<double>>& image,
                                      vector<double>& x_hits_tr,
                                      vector<double>& y_hits_tr,
                                      vector<double>& z_hits_tr,
                                      vector<double>& energy_hits) {

    
    int row_cut = -1;
    
    //CUT TRACKS due to exposure of camera
    double randcut = gRandom->Uniform(config.getDouble("exposure_time")+readout_time);
    // DEBUG
    //cout<<"DEBUG: randcut = "<<randcut<<endl;

    if (randcut<readout_time) {
                            
        double y_cut_tmp = config.getDouble("y_dim") * (0.5 - randcut/readout_time);
        
        // Removing elements from x_hits_tr
        x_hits_tr.erase(std::remove_if(x_hits_tr.begin(), x_hits_tr.end(), [&](const double& x) {
            return y_hits_tr[&x-&*x_hits_tr.begin()] > y_cut_tmp;
        }), x_hits_tr.end());
        // Removing elements from z_hits_tr
        z_hits_tr.erase(std::remove_if(z_hits_tr.begin(), z_hits_tr.end(), [&](const double& z) {
            return y_hits_tr[&z-&*z_hits_tr.begin()] > y_cut_tmp;
        }), z_hits_tr.end());
        // Removing elements from energy_hits
        energy_hits.erase(std::remove_if(energy_hits.begin(), energy_hits.end(), [&](const double& e) {
            return y_hits_tr[&e-&*energy_hits.begin()] > y_cut_tmp;
        }), energy_hits.end());
        // Removing elements from y_hits_tr [must be done after the previous ones]
        y_hits_tr.erase(std::remove_if(y_hits_tr.begin(), y_hits_tr.end(),[&](const double& y) {
            return y > y_cut_tmp;
        }), y_hits_tr.end());
            
        row_cut = (int)(randcut * static_cast<double>(y_pix) / readout_time);

        for(int ix = 0; ix < (int)image.size(); ix++) {
            for(int iy = 0; iy < (int)image[0].size(); iy++) {
                if(iy > row_cut) {
                    image[ix][iy] = 0.0;
                }
            }
        }
 
    } else if (randcut>config.getDouble("exposure_time")) {
        double y_cut_tmp = config.getDouble("y_dim") * (0.5 - (randcut - config.getDouble("exposure_time")) / readout_time);
        
        // Removing elements from x_hits_tr
        x_hits_tr.erase(std::remove_if(x_hits_tr.begin(), x_hits_tr.end(), [&](const double& x) {
            return y_hits_tr[&x-&*x_hits_tr.begin()] < y_cut_tmp;
        }), x_hits_tr.end());
        // Removing elements from z_hits_tr
        z_hits_tr.erase(std::remove_if(z_hits_tr.begin(), z_hits_tr.end(), [&](const double& z) {
            return y_hits_tr[&z-&*z_hits_tr.begin()] < y_cut_tmp;
        }), z_hits_tr.end());
        // Removing elements from energy_hits
        energy_hits.erase(std::remove_if(energy_hits.begin(), energy_hits.end(), [&](const double& e) {
            return y_hits_tr[&e-&*energy_hits.begin()] < y_cut_tmp;
        }), energy_hits.end());
            
        // Removing elements from y_hits_tr [must be done after the previous ones]
        y_hits_tr.erase(std::remove_if(y_hits_tr.begin(), y_hits_tr.end(),[&](const double& y) {
            return y < y_cut_tmp;
        }), y_hits_tr.end());
        
        row_cut = (int)((randcut-config.getDouble("exposure_time")) * static_cast<double>(y_pix) / readout_time);

        for(int ix = 0; ix < (int)image.size(); ix++) {
            for(int iy = 0; iy < (int)image[0].size(); iy++) {
                if(iy < row_cut) {
                    image[ix][iy] = 0.0;
                }
            }
        }

    }
    
    // DEBUG
    //cout<<"DEBUG: row_cut  = "<<row_cut<<endl;
    return row_cut;
}



double TrackProcessor::GetTrackVariable(const std::string& varname,
                                        std::vector<double>& x_hits_tr,
                                        std::vector<double>& y_hits_tr,
                                        std::vector<double>& z_hits_tr) {
    double ret = -1.0;
    int numhits = (int)x_hits_tr.size();

    int varcode = varname_map[varname];

    switch (varcode) {
        case 0: { // proj_track_2D
            double proj_track_2D = 0;
            for (int i = 0; i < numhits - 1; ++i) {
                proj_track_2D += std::hypot(x_hits_tr[i + 1] - x_hits_tr[i],
                                            y_hits_tr[i + 1] - y_hits_tr[i]);
            }
            ret = proj_track_2D;
            break;
        }
        case 1: // x_vertex
            ret = (x_hits_tr[0] + 0.5 * config.getDouble("x_dim")) * static_cast<double>(x_pix) / config.getDouble("x_dim");
            break;
        case 2: // y_vertex
            ret = (y_hits_tr[0] + 0.5 * config.getDouble("y_dim")) * static_cast<double>(y_pix) / config.getDouble("y_dim");
            break;
        case 3: // z_vertex
            ret = z_hits_tr[0];
            break;
        case 4: // x_vertex_end
            ret = (x_hits_tr[numhits - 1] + 0.5 * config.getDouble("x_dim")) * static_cast<double>(x_pix) / config.getDouble("x_dim");
            break;
        case 5: // y_vertex_end
            ret = (y_hits_tr[numhits - 1] + 0.5 * config.getDouble("y_dim")) * static_cast<double>(y_pix) / config.getDouble("y_dim");
            break;
        case 6: // z_vertex_end
            ret = z_hits_tr[numhits - 1];
            break;
        case 7: // x_min
            ret = (*std::min_element(x_hits_tr.begin(), x_hits_tr.end()) + 0.5 * config.getDouble("x_dim")) * static_cast<double>(x_pix) / config.getDouble("x_dim");
            break;
        case 8: // x_max
            ret = (*std::max_element(x_hits_tr.begin(), x_hits_tr.end()) + 0.5 * config.getDouble("x_dim")) * static_cast<double>(x_pix) / config.getDouble("x_dim");
            break;
        case 9: // y_min
            ret = (*std::min_element(y_hits_tr.begin(), y_hits_tr.end()) + 0.5 * config.getDouble("y_dim")) * static_cast<double>(y_pix) / config.getDouble("y_dim");
            break;
        case 10: // y_max
            ret = (*std::max_element(y_hits_tr.begin(), y_hits_tr.end()) + 0.5 * config.getDouble("y_dim")) * static_cast<double>(y_pix) / config.getDouble("y_dim");
            break;
        case 11: // z_min
            ret = (*std::min_element(z_hits_tr.begin(), z_hits_tr.end()));
            break;
        case 12: // z_max
            ret = (*std::max_element(z_hits_tr.begin(), z_hits_tr.end()));
            break;
        default:
            throw std::logic_error("Unhandled varname_map code: " + std::to_string(varcode));
    }

    return ret;
}
