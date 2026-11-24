# Phone Spectrometer for Public Outreach — Advanced Version (v2.1.0)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24506272.svg)](https://doi.org/10.5281/zenodo.24506272)

Open-source code and validation data for turning ordinary smartphone photographs of a projected spectrum
into absorption spectra, with gamma (sRGB) linearization, a measured camera response, and validation against
a laboratory UV–vis spectrophotometer.

This is the **advanced version**. The **basic version** (v1.0.0, the code used at European Researchers' Night;
https://doi.org/10.5281/zenodo.23400260) remains available unchanged. Everything in the basic version is
included here, plus the extensions and data listed below.



## Contents

| Folder | What it contains |
|---|---|
| `matlab_code/` | MATLAB code: `v1_as_used/` (event scripts), `v2_user_friendly/` (point-and-click app, batch script, UV–vis tool, `core/`, `experiments/` for the validation protocol). No toolboxes needed (tested in R2026b). See `matlab_code/v2_user_friendly/HOW_TO_USE.txt`. |
| `analysis/` | Python analysis that produced every number and figure in the paper. `analyse_paperB.py` → `results.json`; `response.py` (camera response, Debevec–Malik); `make_figs.py`, `fig1_setup.py`, `cuvette_panel.py`, `graphical_abstract.py` (figures). |
| `photos/raw/ern_photos/` | All validation photographs (JPG), Samsung Galaxy S25 FE, 19–20 Nov 2026, named by time stamp. `photos/photo_sets.md` and `photos/inventory.tsv` say which photo is which (dark, ethanol references, DCM series, 2026 mixed-dye set, lasers, shutter series). |
| `validation_photos_DNG.zip` (separate file in the Zenodo record) | RAW (DNG) versions of the same 69 photographs (350 MB; 1.7 GB unzipped). Unzip into `photos/raw/ern_photos/` to run `raw_pipeline.py`. |
| `photos/setup/` | Top view of the spectrometer (Figure 1b); the 2025 mixed-dye set used at the event (Figure S10). |
| `uvvis/` | Laboratory UV–vis spectra (CSV, converted from the instrument's UTF-16 export to UTF-8): `a1` ethanol; DCM `b1` 3 µL, `b2` 4 µL, `a2` 5 µL, `a3` 10 µL, `a4` 15 µL, `a5` 20 µL of 5 mM stock in 2.00 mL ethanol; `real1`–`real5` the 2026 mixed-dye set. |
| `docs/` | **`PYTHON_GUIDE.md`: step-by-step guide to the Python analysis, with a picture for every step** (`docs/guide/`, made by `analysis/guide_figures.py`). |
| `figs/` | All figures of the paper and supplement (PNG 300 dpi and PDF), the 3D schematic (SVG), and the graphical abstract. |

## Reproduce the paper's numbers

New to the code? Start with **[docs/PYTHON_GUIDE.md](docs/PYTHON_GUIDE.md)**, which explains each step with pictures.

```
pip install -r requirements.txt
cd analysis
python analyse_paperB.py      # ~1.5 min, writes results.json and arrays.npz
python make_figs.py           # figures into ../figs
python raw_pipeline.py        # same analysis from the RAW (DNG) files (needs validation_photos_DNG.zip from the Zenodo record, unzipped into photos/raw/ern_photos)
python blank_baseline.py      # baseline correction from the repeated ethanol blanks
python calib_compare.py       # calibration options for the mixed-dye set
python guide_figures.py       # pictures for docs/PYTHON_GUIDE.md (one per analysis step)
python gamma_2026.py          # encoded vs linearized absorbance, 2026 set and DCM (Figure 3b)
python fig_basic2025.py       # Figure S10: 2025 set with the basic version (needs the basic record's example_data in basic_example_data/)
python fig3_merged.py         # Figure 3: gamma encoding, camera response, wavelength calibration
python fig5_panels.py         # Figure 4: DCM cuvettes, spectra, linearity
python fig8_improved.py       # Figure 5 (mixed-dye cuvettes + corrected spectra), Figures S13 (normalized) and S14 (absolute); chooses the red-end axis adjustment
```

Figure 5b uses an empirical adjustment of the wavelength axis above 587 nm (the last calibration band): nm' = nm + c·max(0, px − px_S1)², with c chosen by the highest mean correlation over the five mixed-dye samples (c = −4×10⁻⁵ nm/px², red end 663 nm instead of 670 nm). It is fitted to the samples, not an independent calibration.

`results.json` shipped in this archive is the output used in the paper; re-running reproduces it exactly.

## Key results (DCM in ethanol, 7.5–49.5 µM)

| Processing | Slope vs UV–vis | R² |
|---|---|---|
| sRGB-linearized | 1.09 | 0.996 |
| Encoded (no correction) | 0.60 | 0.9999 |
| Measured camera response | 0.96 | 0.998 |

Detection limit (cuvette removed and replaced): A = 0.11 (3 µM DCM). Reference-based wavelength calibration:
cross-validated error 8.5 nm.

## What has been tested

| Part | Tested with | Result |
|---|---|---|
| Python analysis (`analysis/`) | the 2026 validation photos in this repository | reproduces every number in the paper; rerunning on a clean copy gives an identical `results.json` |
| MATLAB `batch_example.m` | `matlab_code/example_data` (2026 photos), MATLAB R2026b | same spectra and band positions as the Python analysis |
| MATLAB app `PhoneSpectrometer.m`: measuring, peak table, wavelength axis from known peaks, lab overlay, normalizing, saving and loading settings | same data, MATLAB R2026b | same band positions as the Python analysis (within 5 nm); axis 434–658 nm |
| MATLAB `uvvis_process.m` | `example_data/uvvis`, MATLAB R2026b | same lamp-artefact correction as the Python code (identical values) |
| Live mode (watching a folder), drawing the strip by hand, calibration from laser or fluorescent-lamp photos, `experiments/e1_calibration.m` and `e3_dilution.m` | not used for the paper | provided for new measurements; please report problems through GitHub issues |

The analysis was validated with one phone (Samsung Galaxy S25 FE). Other phones should work with the same
protocol (Pro mode, fixed settings; SI section S12 of the paper), but their photos have not yet been analyzed here.

## How to cite

Raj, A. *Phone spectrometer for public outreach: advanced version with validation data* (v2.1.0). Zenodo, 2026. https://doi.org/10.5281/zenodo.24506272

Please also cite the companion paper (details will be added on publication). `CITATION.cff` gives the same information in machine-readable form.

## License

MIT (see `matlab_code/LICENSE`). Photographs and spectra are released under the same terms.

## Figure files and their numbers in the JCE manuscript

| Figure | File (figs/) |
|---|---|
| 1 | Fig1_setup |
| 2 | Fig3_photo_to_absorbance |
| 3 | Fig3_gamma_response_calibration |
| 4 | Fig5_linearity |
| 5 | Fig5_event_samples |
| S10 | S_basic2025 (2025 set at the event, basic version) |
| S13, S14 | S_event_samples_normalized, S_event_samples_absolute |
