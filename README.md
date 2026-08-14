# GW250114 rotating-hairy QNM analysis code

This repository contains the public analysis code accompanying the article
"Testing black-hole quasinormal-mode grids with GW250114: a rotating-hairy
case study".

The repository is intentionally **code only**. It does not duplicate the
manuscript, bibliography, submission files, or publication figures maintained
separately in Overleaf.

## Publication snapshot policy

This repository is intended as the frozen code snapshot for one article. After
publication, scientific extensions should be released in a new repository
linked to the new article. The history of this repository should remain
available so readers can recover the code that accompanied the published work.

## Scientific scope

The code implements:

- a Dudley--Finley continued-fraction calculation for the rotating hairy
  black-hole spectrum studied by Zhen Li;
- construction and interpolation validation of the production QNM grid;
- projection of the grid onto public pyRing and RINGDOWN posterior products
  for GW250114;
- static supplied-potential validation examples;
- a numerical-relativity-calibrated evolving-Kerr false-hair control;
- positive-injection and robustness calculations.

The event-level calculation uses public marginalized posterior products. It is
not a new detector-strain likelihood. The Dudley--Finley equation is an
approximate spectral prescription and does not supply the full coupled
perturbation sector of every theory sharing the background metric.

## Repository contents

- `scripts/python/`: numerical solvers, inference, controls, and audits.
- `scripts/wolfram/`: public full-IMR posterior calibration.
- `config/`: fixed analysis settings used for the article.
- `data/`: two compact numerical input tables needed by the scripts; no large
  detector or numerical-relativity data are redistributed.
- `DATA_SOURCES.md`: where to obtain the external public data.
- `REPRODUCIBILITY.md`: run order and expected outputs.
- `MANIFEST.csv`: SHA-256 inventory of the public snapshot.

## Installation

The frozen environment was tested with Python 3.12. Install the core packages:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

The optional reconstruction of the pre-peak remnant prior and the SXS
calibration requires additional packages:

```powershell
python -m pip install -r requirements-optional.txt
```

On Linux or macOS, activate the environment with `source .venv/bin/activate`.

## Quick validation

Check the public folder structure and Python syntax without installing the
scientific dependencies:

```powershell
python scripts/validate_repository.py
```

The shortest independent check reproduces the published QNM table and Kerr
limit:

```powershell
python scripts/python/hairy_continued_fraction.py
```

A smaller development grid can be generated with:

```powershell
python scripts/python/build_hairy_qnm_grid.py
```

The full production grid contains 35,343 direct complex roots and is
computationally more expensive:

```powershell
python scripts/python/build_hairy_qnm_production_grid.py
```

See `REPRODUCIBILITY.md` for the full sequence and for the copy commands that
place the supplied compact derived inputs in their expected runtime paths.

## External data

Large public GW250114 and SXS files are deliberately excluded. Download
instructions, expected paths, and source identifiers are in `DATA_SOURCES.md`.
The repository's `.gitignore` prevents those files from being committed by
accident.

## Citation

If you use this code, please cite the associated article. Bibliographic details
can be added to `CITATION.cff` after the article receives its final DOI.

## License

The original code in this repository is released under the MIT License. Public
data, third-party software, and published reference values remain subject to
their own licenses and citation requirements.
