# Publishing the advanced version on Zenodo (separate record)

The basic version (10.5281/zenodo.23400260) stays as it is. The advanced version is a NEW record, so it gets its own DOI.

1. zenodo.org → log in → **New upload** (do NOT use "New version" on the basic record — that would make it a version of the same concept DOI).
2. Upload:
   - `phone-spectrometer_advanced_v2.1.0.zip` (code, Python analysis, JPG photos, UV–vis, figures)
   - `validation_photos_DNG.zip` — zip the 69 `.dng` files from your phone/computer (folder of 19–20 Nov 2026), ~1.7 GB
3. Fill the form from `.zenodo.json` (title, version 2.1.0, MIT, keywords, description).
4. Related works: add `10.5281/zenodo.23400260` as **"Is derived from"** (or "Continues"). After acceptance, add the paper DOI as **"Is supplement to"**.
5. Click **Reserve DOI** before publishing, and send me the DOI: I will put it in the paper (ref. 13, Data availability, SI S10) and in `CITATION.cff`.
6. Optionally add a line to the basic record's description: "Advanced version with validation data: <new DOI>".
7. Publish.
