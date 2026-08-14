# External data sources

Large public inputs are not committed to this repository. This avoids
duplicating authoritative archives and keeps the code repository small.

## GW250114 public data release

Source record:

- Zenodo record: `16877102`
- DOI family: `10.5281/zenodo.16877101`
- release file: `GW250114_data_release.tar.gz`

After downloading and extracting the release, the scripts expect:

```text
data/raw/GW250114_data_release/data/pre_-40M_with_lnL_cut.dat
data/raw/GW250114_data_release/data/pyring_220_221_delta_posterior_without_frequencies.dat
data/raw/GW250114_data_release/data/220+221+df221+dg221_6M_f220meas_f221meas_df221meas_120Ksamps.hdf5
data/raw/GW250114_data_release/data/220+221_amps_fs_gammas_merged_timescans.hdf5
```

Running `scripts/python/article1_independent_prior.py` converts the public
pre-peak posterior into the compact remnant-prior samples used downstream.

## Full-IMR remnant posterior

The full-IMR robustness branch uses:

```text
data/raw/posterior_samples_NRSur7dq4.h5
```

This is the public `posterior_samples_NRSur7dq4.h5` product from the GW250114
release. Running `scripts/wolfram/gw250114_posterior_calibration.wl` creates the
compact selected columns used by the final analysis.

## Numerical-relativity calibration

The evolving-Kerr control was calibrated using `SXS:BBH:3617`, level 3:

```text
data/raw/SXS_BBH_3617/Lev3_metadata.json
data/raw/SXS_BBH_3617/Lev3_Horizons.h5
data/raw/SXS_BBH_3617/Lev3_Strain_N2.h5
data/raw/SXS_BBH_3617/Lev3_Strain_N2.json
```

The calibration parameters used by the final control are frozen in
`config/hairy_gw250114_publication.json`; the large SXS files are needed only
to reconstruct that calibration from scratch.

## Published rotating-hairy QNM table

`data/ovalle_rotating_hairy_qnm_zhen_li_2022.csv` is a compact transcription
used for implementation validation. Its scientific source is Zhen Li,
arXiv:2212.08112. It is included as numerical reference data, not as manuscript
content.

## Integrity and licensing

Users should obtain external data from the original repositories and cite the
corresponding data releases and papers. Those files are governed by the terms
and licenses of their original providers, not by this repository's MIT
license.
