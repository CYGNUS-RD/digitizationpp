# From GEANT4 energy deposits to camera images and PMT waveforms

Source review: [piacent/digitizationpp](https://github.com/piacent/digitizationpp), commit `a33fb5439406ecec4858ed94d3f722ed985f1ce3`, inspected 17 November 2026. This document describes the implemented forward model, including differences between configuration comments and executable behavior. It is a source-based explanation, not a numerical validation against detector data.

## 1. Introductory scheme: physics and decision flow

GEANT4 supplies the positions and energies of deposits in the gas. This program does not rerun particle transport: it converts those deposits into surviving ionization electrons, amplified and diffused charge, and finally measured camera counts and PMT voltages. The camera integrates a two-dimensional projection; the PMTs use the longitudinal distribution to construct a time profile.

The shared physical chain is:

**Energy deposition → ionization and attachment → drift diffusion → GEM1/GEM2 multiplication → local three-dimensional charge density → saturated GEM3 response.**

From that shared distribution, two detector branches follow:

- **Camera:** integrate along drift → map onto pixels → optional optical blur → vignetting, collection efficiency and quantum efficiency → Poisson detection → camera counts → optional exposure cut → pedestal addition.
- **PMTs:** retain voxel positions → effective detection yield → relative electron-arrival times → geometrical collection and Poisson detection at each PMT → single-photoelectron pulses → noise → sampling and quantization.

### 1.1 Ordered logical scheme, including branches

```text
Read configuration; establish GEM gains, camera geometry and collection factor
Read ROOT nTuple events (start_event, events)
  Classify event using energyDep_NR > 0
    NR event and NR=False: skip event
    NR event and NR=True: use supplied energyDep_hits_NRQF
    Otherwise: use energyDep_hits
  If total energyDep < ion_pot: save background-only event
  Transform axes, add offsets/extras, apply one rigid random XYZ shift
  Select deposits inside [min_drift, max_drift]
    No active deposits: zero camera signal and no PMT simulation
  Sample primary electrons surviving attachment
    primary_granularity=True: drift each primary, then spread its avalanche
    False: spread avalanche electrons with combined drift/intrinsic width
  Build post-GEM2 electron occupancy in voxels
    Low estimated occupancy: sparse map
      High event energy: process input populations in WID batches
      Otherwise: process in one batch
    High estimated occupancy: dense z slabs with bounded electron blocks
  Apply GEM3 saturation once to each complete voxel population
  Store voxel response for PMTs; project response for camera
  Camera: optional sigma0_light; optional Vignetting; Poisson and counts
  PMT branch if pmt_mode=True and retained effective light exists
    Fast / Slow / Both selects which pulse arrays receive signal
    Add noise and quantize arrays
  Save signal integral and sparse pixels BEFORE exposure cut
  If exposure_time_effect=True: cut camera rows and corresponding hits
  Add camera pedestal if bckg=True
  redpix_output=True: save sparse signal and event_info, no image histogram
  redpix_output=False: also save camera image histogram
```

Saturation is always used by the active event loop. Setting `beta=0` removes its nonlinear suppression; the separate `ph_smearing2D` helper is not the selected alternative pipeline.

### 1.2 Flags and subcases at a glance

| Control | Implemented consequence |
|---|---|
| `NR` | Allows NR-containing events when true; false skips them, rather than converting them as ER. |
| `primary_granularity` | Chooses correlated, per-primary drift clouds versus independently smeared avalanche electrons. Missing boolean reads false. |
| `sigma0_light` | Positive value adds camera light blur after saturation; absent/zero disables it. |
| `Vignetting`, `Vig_Map` | Enable and select the camera spatial response map. |
| `Camera_type` | `Fusion`, `Quest1`, `Quest2` select pixel dimensions, conversion, QE and readout time. |
| `exposure_time_effect`, `exposure_time` | Enable the random rolling-exposure row cut and set exposure duration. |
| `bckg`, `random_ped` | Add measured camera pedestal; `-1` selects a random frame, other values a fixed frame. |
| `redpix_output` | Select sparse-only versus image output. The implementation still calls pedestal loading when `bckg=True`, even in sparse mode. |
| `pmt_mode` | Enable the PMT branch; it is also conditional on effective voxel light being nonempty. |
| `digitizers` | `Fast`, `Slow`, `Both` gate pulse generation; both arrays still receive noise and are returned. |
| `randX/Y/Z_min`, `randX/Y/Z_max` | Uniform rigid event shifts; equal bounds give fixed shifts. Both bounds are required. |
| Legacy `randX/Y/Z_range` | Symmetric shift interval about zero, used only without explicit bounds. |
| `min_drift`, `max_drift` | Inclusive response-deposit selection. Missing values actually default to 0 and 510 mm. |
| `fixed_seed` | Sets ROOT seed to 3 and affects electron smearing streams; does not seed the independent PMT random generators. |
| `WID`, `Parallel_threads` | Computational batching/concurrency, not new physical processes. |
| `start_event`, `events`, `start_run_number`, `queue` | Event selection, output numbering and output location; `events=-1` requests remaining events, `queue` forces local output. |
| `save_mc_first3_hits` | Save up to three original GEANT4 hit positions as truth metadata. |
| `GEANT4isotopes`, `rootfiles`, `donotremove`, `tag`, `zcloud` | Present in the supplied configuration but not active physical switches in the inspected pipeline. In particular, `run()` calls ROOT processing unconditionally. |

The remaining parameters are physical coefficients, geometry, calibration or file locations; they are explained below. Boolean parsing accepts `True`, `true` or `1`.

## 2. GEANT4 inputs and recoil-energy interpretation

**Physics.** Deposited energy becomes ionization. Nuclear recoils can deposit energy in channels that do not produce the same ionization yield as electronic recoils. The digitizer relies on an already supplied quenched energy vector; it does not calculate a quenching-factor law itself.

**Implementation.** `DigitizationRunner::processRootFiles` reads `nTuple`, including `energyDep`, `energyDep_NR`, `energyDep_hits`, `energyDep_hits_NRQF`, positions and PDG identifiers. If `energyDep_NR>0`, the entire event uses `energyDep_hits_NRQF`; otherwise it uses `energyDep_hits`. The total `energy` metadata remains `energyDep`. Thus response energy, total deposited energy and active-volume energy are not interchangeable.

The preliminary background-only threshold uses total `energyDep < ion_pot`, before active-volume filtering. This is an explicit event-level shortcut, not the later Poisson ionization process. Positions are interpreted in mm and energies in keV.

## 3. Detector placement and active volume

**Physics.** Drift distance sets attachment, diffusion and electron-arrival time. In the digitization frame the GEM plane is at z=0 and positive z points into the drift region.

**Implementation.** `MC_xaxis`, `MC_yaxis`, `MC_zaxis` specify where each input axis goes, with an optional minus sign. Each output coordinate is the signed input coordinate plus its `*_offset` and `*_extra`. One additional shift per axis is sampled for the whole event; individual hits are not relocated independently. The sampled shifts are saved as `placement_dx_mm`, `placement_dy_mm`, `placement_dz_mm`.

Only hits satisfying `min_drift <= z <= max_drift` contribute to detector response. Complete-track geometry metadata is computed before that selection. Camera edges clip projected light; placement does not enforce full-track containment. The comment suggesting that omission of the upper bound removes the cut is inconsistent with the implemented 510 mm fallback.

## 4. Primary ionization and attachment during drift

For deposit energy E at drift distance z, the mean number of primary electrons reaching the GEMs is

`mu_primary = (E / ion_pot) × exp(-z / absorption_l)`.

A Poisson draw supplies the surviving count. The supplied configuration uses `ion_pot=0.035 keV` and `absorption_l=1350 mm`.

**Physical interpretation.** The energy-per-electron parameter is an effective ionization work scale. Attachment removes electrons exponentially with path length. Ionization and survival fluctuations are combined in a single Poisson sample; there is no separately sampled initial population followed by binomial attachment, and no explicit Fano-factor model.

**Implementation.** `prepare_primaries` or `NelGEM2` performs this calculation. `N_primaries_reaching_GEMs` records the count before multiplication, for the active deposits passed to the processor. Negative input drift distances cause these routines to fail; normal event processing filters them out beforehand.

## 5. Diffusion and GEM1/GEM2 amplification

### 5.1 Gain model

Each GEM has a mean gain

`G_i = c_G × exp(alpha_G × GEMi_HV)`

and an extraction factor

`epsilon_i = 0.87319885 × exp(-0.002 × GEMi_HV)`.

GEM1 avalanche size is sampled from an exponential distribution with mean G1. It is multiplied by epsilon1, G2 and epsilon2. GEM2 gain is deterministic in this model. Although epsilon3 is initialized, it does not appear in the active GEM3 saturation expression.

### 5.2 Per-primary diffusion: primary_granularity=True

Each surviving primary gets one Gaussian displacement with widths

`sd_T = sqrt(diff_coeff_T × z/10)` and `sd_L = sqrt(diff_coeff_L × z/10)`.

The factor 10 converts mm to cm for the diffusion coefficients. The primary then receives a rounded post-GEM2 avalanche population, sampled as `round(Exp(G1) × epsilon1 × G2 × epsilon2)`. Its individual avalanche electrons are spread around that displaced primary with widths `sqrt(diff_const_sigma0T)` and `sqrt(diff_const_sigma0L)`.

**Consequence.** Electrons from one primary share its drift displacement. This preserves local avalanche correlations, which matter because saturation depends nonlinearly on local density.

### 5.3 Legacy diffusion: primary_granularity=False

GEM1 exponential contributions are summed per GEANT4 hit, then multiplied by G2 and epsilon2 and rounded. Each resulting electron is independently smeared around the original hit using

`sigma_T,L = sqrt(diff_const_sigma0T,L + diff_coeff_T,L × z/10)`.

The marginal widths combine the same two variance terms, but the correlation structure differs. Equal overall cloud widths therefore need not imply equal saturation or image texture. The intrinsic constants are variances, despite their sigma-like names.

Changing `drift_field` does not automatically update these diffusion coefficients: it affects the PMT drift-velocity calculation separately.

## 6. Voxel creation and filling

**Physics.** Charge arriving in the same small volume competes for the GEM response. Voxelization supplies a local occupancy on which the saturation model acts; it is not a microscopic tracking of individual GEM holes.

**Grid.** The active-hit bounding box is expanded by 15 mm on each transverse side and 10 mm on each longitudinal side. Bin counts are rounded upward to even integers. Dimensions are

- `dx = x_dim / x_pix / xy_vox_scale`;
- `dy = y_dim / y_pix / xy_vox_scale`;
- `dz = z_vox_dim`.

The grid is local to the event. `xy_vox_scale` must be an integer at least one. A position is assigned using `floor((coordinate - grid_min)/bin_size)`.

### 6.1 Sparse-map route

The code chooses a sparse map unless `energy/(Nx×Ny×Nz) > 1e-6`. This is a computational heuristic using total event energy, not an additional interaction law.

Gaussian electron coordinates are generated explicitly. Occupancy is stored with flattened index `ix×Ny×Nz + iy×Nz + iz`; thread-local maps are merged by summing counts. If the event is NR above 100 keV or ER above 500 keV, populations are processed in batches of `WID`. All batches contribute to the same map before saturation. With per-primary mode, the batched entries are primaries rather than original hits. WID does not bound the electron count of a large avalanche or total map memory.

### 6.2 Dense-slab route

The dense grid is split into z slabs capped at 50 million cells, approximately 400 MB of doubles before other allocations. A transverse plane too large for this cap causes failure.

For each input population, sequential conditional binomial draws allocate electrons to Gaussian z intervals. Allocated electrons are subtracted from the remaining population, so no electron can be reused in another z bin. Transverse coordinates are then sampled in blocks of at most 65,536 electrons and accumulated into the slab. Gaussian tails outside the grid are discarded. Saturation is applied only after all populations have contributed to the voxel.

This route avoids storing all three coordinates for every amplified electron simultaneously. The two algorithms target the same occupancy model but use different sampling sequences; identical random seeds do not imply identical event realizations across routes. Sparse-mode out-of-grid checks are explicitly conditional on `primary_granularity`, unlike the dense route.

## 7. GEM3 saturation and the shared camera–PMT quantity

For a voxel containing n post-GEM2 electrons, the response is

`S(n) = n × A × G3 / (1 + beta_eff × G3 × n)`,

where `beta_eff = beta × xy_vox_scale² × 0.1/z_vox_dim`.

At low occupancy the response is linear. At high occupancy it tends toward `A/beta_eff`. The voxel-size rescaling compensates for the configured transverse subdivision and longitudinal bin size relative to the reference convention. It is not a general recalibration for arbitrary detector pixel pitch.

**Implementation.** `Nph_saturation` implements this formula. Contributions must be combined before applying it: saturating two charge populations separately would generally give a different answer from saturating their sum.

Despite names such as `n_ph_val` and `PMTVoxel::n_photons`, S is the saturated avalanche response **before the camera's `photons_per_el` factor**. The camera sums S into the transverse integer array `hout`, so fractional contributions are truncated during accumulation. PMT voxels retain the unprojected response and physical coordinates. Coordinates are reconstructed at voxel lower edges, not centers.

This distinction is essential: the camera and PMT branches share S, but do not implement the same subsequent photon-yield normalization.

## 8. Camera projection, optics and photon detection

**Projection.** The z contributions are summed, the local transverse grid is positioned within the camera field, and subpixels are merged into camera pixels. Light outside the camera is discarded.

**Optical blur.** Positive `sigma0_light` applies a separable, pixel-integrated Gaussian with approximately four-sigma support, expressed in pixels using `x_dim/x_pix`. This happens after saturation. It redistributes light without reducing charge density in the saturation calculation; light lost at the image boundary is not restored. It does not blur PMT voxels.

**Collection.** With demagnification `m = y_dim/y_sensor_size` and aperture f, the collection factor is

`omega = 1/[4(m+1)f]²`.

For projected response S_pixel, the detected count is sampled as

`Ndet ~ Poisson(S_pixel × vignette × omega × camera_QE × photons_per_el)`.

The vignette factor is one when disabled, otherwise read from the map at the pixel position. Although the variable is named number of photons, inclusion of QE makes this an effective detected-photoelectron count.

**Conversion.** The camera stores `int(Ndet × counts_per_photon)`.

| Camera | Pixels | QE | Conversion | Readout time |
|---|---:|---:|---:|---:|
| Fusion | 2304×2304 | 0.80 | 1/0.21 counts per detected quantum | 184.4 ms |
| Quest1 | 4096×2304 | 0.85 | 1/0.107 | 199.0 ms |
| Quest2 | 4096×2304 | 0.85 | 1/0.107 | 39.0 ms |

The configured camera-electron RMS values are not sampled in the active response: the Gaussian read-noise code is commented out in favor of adding measured pedestal frames. There is no explicit camera full-well saturation or ADC clipping in this path; GEM saturation is a separate effect.

## 9. PMT light selection and timing

The PMT branch runs when `pmt_mode=True` and the voxel list is nonempty. For each voxel it computes

`Neff = round(S × 0.0136)`.

Only positive Neff values survive. The code labels 0.0136 as an average PMT QE. It does **not** additionally multiply by `photons_per_el`. Consequently this factor should be treated as the implemented effective normalization, whose physical calibration must be checked before interpreting absolute PMT photoelectron yields. Sub-threshold voxel contributions are discarded individually rather than pooled.

The earliest retained voxel defines z_ref. Arrival times are

`t_arrival = abs(z - z_ref)/v_drift`.

`Utils::compute_drift_velocity` linearly interpolates a hard-coded He/CF4 60/40 field–velocity table, with endpoint clamping. Dividing its cm/us result by 100 gives mm/ns. At the configured 0.8 kV/cm, the interpolation gives 5.296 cm/us, or 0.05296 mm/ns.

The z distribution already contains longitudinal diffusion. Timing preserves relative drift spread but removes absolute event drift time. No GEANT4 hit-time input, explicit optical flight time, or prompt-scintillation component is used here. The voxel z coordinate is a drift-time proxy; optical emission is treated as occurring at the GEM plane.

## 10. PMT geometrical collection and photoelectron statistics

For PMT radius r, GEM-to-PMT distance d, and lateral displacement `(dx,dy)`, define `R²=dx²+dy²+d²`. The implemented mean is

`mu_PMT = Neff × r² × d² / (4 R⁴)`.

Each PMT independently receives a Poisson draw with that mean. This is the implemented approximate collection law; it is not a ray-tracing calculation or an exact finite-disk solid angle. It combines inverse-distance and angular dependence, with no explicit reflections or optical attenuation.

**Implementation.** `PhotonPropagation::pmt_hits` loops over retained voxels and configured PMTs. Only positive draws are stored, paired with the voxel arrival time. Geometry uses `pmt_number`, `pmt_radius`, `dist_gem_pmt` and `pmt_i_x/y`.

## 11. Single-photoelectron pulse formation

For every detected count, `SignalSimulation::spe_signal` draws a Gaussian transit delay and an exponential dispersion delay, then adds an exponentially modified Gaussian voltage pulse.

The transit mean is `transit_time` in ns. **The implementation treats `transit_time_spread` as FWHM and divides by 2 sqrt(2 ln 2)**, although the example configuration calls it a standard deviation. The exponential dispersion parameter `exp_dispersion_scale` is a rate in s⁻¹, so its mean delay is its reciprocal.

The pulse has charge `Q = -1.6e-19 × pmt_gain` coulombs, load resistance 50 ohms, Gaussian width `pmt_sigma` in seconds and exponential-tail rate `pmt_lambda` in s⁻¹. Its negative polarity follows from Q. Pulses are summed linearly; no per-photoelectron gain distribution or PMT saturation model is added.

Fast and slow windows add fixed offsets of 200 and 1466 samples respectively. When both are enabled, the same photoelectron delay feeds both samplings. The code evaluates pulses at sample times; it does not integrate them over ADC sample intervals. Finite windows can truncate pulse area.

## 12. PMT noise, digitizers and quantization

Time samples are `t_i=i/Fs`. The example configuration requests 1024 samples at 750 MHz and 4000 samples at 250 MHz. `digitizers` controls pulse insertion into these arrays.

Noise templates are selected by `detector`, channel and Fast/Slow label, loaded from the local `pmt` directory or downloaded using `pmt_noise_cloud_dir`. The routine treats the input values as PSD values, draws random Fourier phases, sets magnitudes proportional to `sqrt(PSD×Fs/N)`, constructs a conjugate spectrum and applies an inverse Fourier transform. The resulting noise is added to the waveform.

Finally, each voltage sample is rounded as `round(V×4096)/4096`. This is a step of 1/4096 in waveform units, not a complete ADC model with a configured voltage range or rail clipping.

Implementation qualifications:

- Both fast and slow arrays receive noise even if only one receives signal.
- If no voxel survives effective-light selection, waveform generation is skipped entirely; the event does not automatically receive a noise-only waveform.
- Propagation and noise setup use `pmt_number`, but waveform generation and output assignments explicitly enumerate PMTs 1–4. Arbitrary PMT counts therefore are not consistently supported.
- PMT generators use `random_device`, independently of `fixed_seed`.

## 13. Camera exposure and measured backgrounds

**Exposure.** With `exposure_time_effect=True`, a time is sampled uniformly over exposure plus readout duration. During the initial readout interval, one side of a row boundary is removed; during the final readout interval, the other side is removed. In the central interval there is no cut. The same routine filters hit vectors by y and zeros corresponding image rows. This approximates partial sensor exposure, rather than tracking individual electron arrival times through a shutter simulation.

**Ordering matters.** PMT waveforms, the ordinary signal integral and sparse pixels have already been generated before this cut. PMTs are therefore unaffected by the camera exposure mask. Separate cut quantities describe the modified camera image and hits.

**Background.** `AddBckg` loads a pedestal image from `bckg_path` plus `bckg_name`; if absent it assembles the configured `noiserun` entries from `ped_cloud_dir`. The pedestal is added pixel by pixel to the integer signal image. `bckg=False` leaves the background zero. Pedestals represent measured detector/electronic backgrounds; they are not additional GEANT4 interactions.

## 14. Outputs, interpretation and special cases

Image mode writes `histograms_RunXXXXX.root`, including `TH2I` images and `event_info`. Sparse mode writes `digi_RunXXXXX.root` with signal pixels and metadata, omitting image histograms. Sparse coordinate and content vectors use `uint16_t`; they should not be assumed to preserve arbitrary large pixel values without range limits.

`N_photons` is the sum of digitized camera signal counts, despite its name. It is not the total number of emitted optical photons. `N_primaries_reaching_GEMs` is a different quantity: surviving primary electrons before GEM multiplication.

`joinPeds` reconstructs images from saved sparse pixels and adds pedestal frames. Because sparse pixels are recorded before the exposure cut, this reconstruction does not reproduce an exposure-cut image automatically.

The low-energy shortcut writes a background-only event with reported energy zero. An explicit response failure writes a background-only event marked with energy -1. An event with no active deposits returns an empty signal successfully. NR events rejected by `NR=False` are skipped instead of written as background-only events.

Truth positions saved by `save_mc_first3_hits` refer to original GEANT4 coordinates. Direction inferred from the first two deposited-hit positions is not necessarily the generated primary-particle direction; that requires joining the MC truth.

## 15. Source map and model boundaries

The principal implementation references, pinned to the reviewed revision, are:

- [DigitizationRunner.cxx](https://github.com/piacent/digitizationpp/blob/a33fb5439406ecec4858ed94d3f722ed985f1ce3/src/DigitizationRunner.cxx): input selection, coordinate transforms, calibration, branch ordering, camera background, PMT preparation and output.
- [TrackProcessor.cxx](https://github.com/piacent/digitizationpp/blob/a33fb5439406ecec4858ed94d3f722ed985f1ce3/src/TrackProcessor.cxx): ionization, attachment, correlated diffusion, amplification, voxel occupancy, saturation, camera response and exposure cuts.
- [PMTHits.cxx](https://github.com/piacent/digitizationpp/blob/a33fb5439406ecec4858ed94d3f722ed985f1ce3/src/PMTHits.cxx): geometrical collection and Poisson hit sampling.
- [PMTSignal.cxx](https://github.com/piacent/digitizationpp/blob/a33fb5439406ecec4858ed94d3f722ed985f1ce3/src/PMTSignal.cxx): single-photoelectron pulse, timing, PSD noise and quantization.
- [Utils.cxx](https://github.com/piacent/digitizationpp/blob/a33fb5439406ecec4858ed94d3f722ed985f1ce3/src/Utils.cxx): drift-velocity interpolation and grid rounding.
- [Example configuration](https://github.com/piacent/digitizationpp/blob/a33fb5439406ecec4858ed94d3f722ed985f1ce3/config/ConfigFile_new.txt): supplied numerical settings and comments.

This is a parameterized detector-response model. It does not explicitly solve electric fields, trace electron trajectories through GEM holes, track optical photons, or simulate a detailed electronic circuit. The most consequential interpretation points are the externally supplied NR quenching, per-primary charge correlations, voxel-dependent saturation, different camera/PMT yield normalizations, relative PMT timing, and the distinction between pre-exposure sparse signal and final pedestal-added images.
