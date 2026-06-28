# Data sources and large-file policy

## Public GW250114 inputs

Primary public event products:

- Zenodo record `16877102`, DOI family `10.5281/zenodo.16877101`
- `GW250114_data_release.tar.gz`
- `posterior_samples_NRSur7dq4.h5`

These files are not redistributed in this repository because they are large and
already citable at their source. Place local copies under:

```text
data/raw/GW250114_data_release.tar.gz
data/raw/GW250114_data_release/
data/raw/posterior_samples_NRSur7dq4.h5
```

The `.gitignore` intentionally excludes `data/raw/`, `*.h5`, `*.hdf5`, and
`*.tar.gz`.

## Theory inputs

The higher-derivative Kerr QNM fingerprints are imported from the public
BeyondKerrQNM/Cano et al. fit material and stored in compact derived tables
under:

```text
data/
results/higher_derivative_qnm_complete/
```

## Static QNM references

The static supplied-potential examples use published reference tables for
Schwarzschild, braneworld tidal charge, Bardeen, and Hayward axial
gravitational QNM benchmarks. The repository stores compact validation outputs
and scripts, not third-party PDFs.
