# Lab day checklist: Paper B validation experiments

Print this and tick boxes as you go. The file names below are the ones the scripts expect by default. If you use other names, change the SETTINGS block of the script.

## 0. Before the lab day

- [ ] **Dye for E3: DCM in ethanol** (absorption peak about 470 nm, broad band without sharp structure; ε ≈ 4 × 10⁴ M⁻¹ cm⁻¹; molar mass 303.4 g/mol).
  - Stock: about 3.0 mg DCM in 10 mL ethanol = 1 mM. Dilute 1:10 to 100 µM as a working stock.
- [ ] **Dilution series** (volumetric, from the one working stock, made the same day, kept capped and in the dark): about **2.5, 6, 12, 18, 24 and 36 µM**, plus an ethanol blank. This gives A(470) ≈ 0.1–1.5.
  - The exact values do not matter: the UV-Vis absorbance of the same cuvette is the reference. Just write down what you made.
- [ ] **E2 samples.** Dilute the five event samples to A_max of 0.3–1.0.
- [ ] Use **glass or quartz cuvettes** (ethanol attacks plastic), the same ones for the phone and the UV-Vis. Label on the frosted side only.
- [ ] Phone (Samsung Galaxy S25) charged, with free storage.

## 1. Record the setup (needed for the paper; fill in on the day)

| Item | Value |
|---|---|
| Lamp type / power / supply setting | |
| Condenser lens f, diameter | |
| Iris: max aperture; working opening (mm) | |
| Focusing lens f | |
| Long lens f = 1000 mm, diameter | |
| Prism glass / side length | |
| Grating used? (lines/mm) | |
| Prism → screen distance (cm); spectrum length on screen (cm) | |
| Screen paper | |
| Phone model(s) | |
| Pro mode: ISO / shutter / white balance (K) / focus | |
| Image size (px); JPG / DNG | |
| UV-Vis instrument model; slit / bandwidth | |
| Dyes: names, supplier, lot; solvent; concentrations made | |
| Room temperature | |

Take two or three overview photos of the setup, with the lamp on and off, for Figure 1 and SI S1.

## 2. Camera settings: Samsung Galaxy S25 (all validation photos)

- [ ] **Lens:** use only the **1× main camera**. Do not use 0.6× (ultrawide), 3× or any zoom.
- [ ] **Camera settings (gear icon):** *High efficiency pictures* **off**, so photos are saved as JPG, not HEIC (the scripts read JPG). *Pro mode RAW copies* **on**, so a DNG is saved next to each JPG. *Auto HDR / Super HDR*, *Scene optimiser*, *Shot suggestions* and *Watermark* **off**. Keep one picture size (12 MP) all day. Menu names can differ slightly between One UI versions.
- [ ] **Pro mode:** ISO **50 or 100**; WB a fixed Kelvin value (e.g. **5500 K**), never AUTO; FOCUS **MF** with focus peaking on, focused on the edge of the spectrum; leave any colour/tone adjustments at zero. EV is greyed out once ISO and SPEED are both manual, which is correct.
- [ ] **SPEED (shutter):** turn the histogram on. Choose the speed so that the histogram stays clear of the right edge for the reference photo (water/ethanol cuvette): no pure-white area anywhere in the spectrum. Then **keep this speed unchanged** for E1–E3.
- [ ] Clamp the phone. Use the **2 s timer** or a remote, so touching the screen does not move the phone. Do not touch the phone or the screen during a series.
- [ ] Darken the room, and keep monitor and phone light off the screen.
- [ ] Write ISO, SPEED, WB and picture size in the setup table.

## 3. E1: Wavelength calibration and resolution → `E1_photos/`

Take these at the **end of the E3 session, without moving the phone**. E3 reuses the calibration.

- [ ] `ref_white.jpg`: lamp on, cuvette with solvent.
- [ ] Block the lamp. For each laser, direct it along the beam path through the sample position: `laser405.jpg`, `laser532.jpg`, `laser650.jpg`.
  - Use a **shorter shutter time** if a spot looks pure white. The script warns about saturation.
  - *or* put a CFL at the spectrometer entrance and take `cfl.jpg`, then set `lineFiles = {'cfl.jpg'}` and `lineNm = {[435.8 487 546.1 611]}`.
- [ ] **Resolution.** Put tracing paper or ground glass just before the iris, so the laser fills the iris the way the lamp does. Then photograph the 532 nm line:
  - `res532_iris<working>mm.jpg` at the working opening;
  - `res532_iris<small>mm.jpg` with the iris closed further.
  - Note both openings in mm.
- [ ] Run `e1_calibration.m`. Check panel (a): every dotted line must sit on the right peak.

## 4. E3: Dilution series → `E3_photos/`, `E3_uvvis/`

Take the photos in this order, keeping the same settings throughout:

- [ ] `dark.jpg`: lamp blocked.
- [ ] `ref_01.jpg`: ethanol cuvette.
- [ ] `blank_1.jpg` … `blank_5.jpg`: the ethanol cuvette taken out and put back each time.
- [ ] `c1_1.jpg`, `c1_2.jpg`, `c1_3.jpg`: lowest concentration, three photos each. Take the cuvette out and put it back between photos.
- [ ] `c2_…`, then `ref_02.jpg`; `c3_…`, `c4_…`, then `ref_03.jpg`; `c5_…`, `c6_…`, then `ref_04.jpg`. This puts a reference photo after every two concentrations.
- [ ] **The same cuvettes** go into the UV-Vis (380–800 nm at least, ethanol baseline). Save `c1.csv` … `c6.csv` in `E3_uvvis/`.
- [ ] Then do E1, without moving anything.
- [ ] Run `e3_dilution.m`, and fill in `conc` (the concentrations you actually made) first.
  - **Expected:** linearised slope ≈ 1, encoded slope ≈ 0.5, and curvature above A ≈ 1–1.5 (stray light).

## 5. E2: Peak positions of the five samples (same session if possible)

- [ ] Take `ref.jpg`, then `s1_1.jpg` … `s5_3.jpg` (three photos per sample), with UV-Vis of the same cuvettes.
- [ ] Process with `batch_example.m`, setting `laserFile` to the E1 photo or using the E1 calibration, and `offsetMode = 'none'`. Set `alignOnBlue = false` there too, because several of the samples absorb blue light.

## 6. E4: Exposure (recommended)

- [ ] Use one sample (A ≈ 0.7) at 3–4 shutter speeds. At each speed, take a matched reference (same speed) and the sample: `exp1_ref.jpg`, `exp1_s.jpg`, and so on.
- [ ] Take one extra sample photo at a different shutter speed from its reference (the deliberate mismatch).

## 7. E5: Repeatability and other phones (recommended)

- [ ] With one phone, take the cuvette out and put it back five times, one photo each.
- [ ] With 2–3 other phones (different manufacturers), repeat E2 or three E3 concentrations. Note each model.

## 8. Optional

- [ ] **E6:** a transmission grating setup with three E3 concentrations and a linear calibration.
- [ ] **E7/E8:** as planned in the paper draft.
- [ ] If the phone saves DNG, keep the DNG files of the E3 series (RAW vs JPG comparison).

## After the lab

- [ ] Copy all photos (originals, not screenshots) and UV-Vis files into the folders above.
- [ ] Run E1 → E3 and send `E1_results.txt`, `E3_results.txt` and the PNG figures. I fill in the paper and the SI from these.
- [ ] Back up the raw data. It goes to Zenodo with the paper.
