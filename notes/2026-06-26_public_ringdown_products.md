# Public GW250114 Ringdown Products

This note records the addition of the larger public GW250114 spectroscopy data
release.

## Data

- File: `data/raw/GW250114_data_release.tar.gz`
- Checked compressed size: 108.52 MB
- Extracted directory: `data/raw/GW250114_data_release/`

The archive contains ringdown-specific products that are more directly useful
than the full IMR posterior for the beyond-Kerr branch.

## Implemented

- Added `config/gw250114_public_ringdown_products.wl`.
- Added `scripts/wolfram/gw250114_public_ringdown_products.wl`.
- The script summarizes:
  - `pyring_220_posterior.dat`;
  - `posterior_with_qnm_frequencies.dat`;
  - `220+221+df221+dg221_6M_f220meas_f221meas_df221meas_120Ksamps.hdf5`.
- It verifies from the public Figure 4 script that pyRing stores
  `domega_221`, while the plotted frequency deviation is
  `df_221 = log(1 + domega_221)`.
- It creates a first 1D proxy mapping public `df_221` samples to EFT coupling
  samples using the already imported higher-derivative QNM sensitivities.

## Caveat

The 1D proxy is not a full EFT constraint. It ignores damping time, mode
correlations, amplitudes, polarization mixing, and first-order EFT validity
limits.

## Next Step

Use the public ringdown variables `f_220`, `f_221`, and `df_221` to build a
proper likelihood-level projection for each EFT fingerprint.
