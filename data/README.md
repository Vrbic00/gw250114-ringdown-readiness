# Data

Keep raw downloaded GWOSC data in `data/raw/` and derived/cache products in
`data/cache/`.

Large data products should not be committed unless explicitly needed.

`raw/posterior_samples_NRSur7dq4.h5` is the preferred public GW250114
NRSur7dq4 posterior sample file from the GWOSC event API / Zenodo record
`10.5281/zenodo.16877101`. Its checked size is 25.71 MB. It is used to
calibrate remnant mass/spin posterior widths and Kerr QNM pushforwards, not as
a ringdown-only EFT posterior.

`raw/GW250114_data_release.tar.gz` is the larger public GW250114 spectroscopy
data release from the same Zenodo record. Its checked size is 108.52 MB. The
archive has been extracted into `raw/GW250114_data_release/`; it contains the
ringdown-specific RINGDOWN and pyRing products needed for the project’s
event-level branch.

`beyond_kerr_qnm_selected_fits.csv` is a compact generated import of the
public GPL-3.0 BeyondKerrQNM fit files for modes `220`, `221`, `222`, `330`,
and `440`. Its rows record the exact upstream commit, source file, repository,
paper, and license. Regenerate it with:

```powershell
python scripts/python/import_beyond_kerr_qnm_fits.py <Fits-directory> data/beyond_kerr_qnm_selected_fits.csv --commit <commit-hash>
```
