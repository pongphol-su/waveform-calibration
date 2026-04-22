# waveform-calibration
This repository contains a python code to calibrate spectrum waveform with a reference.

## waveformCalibration.py
- Prompts the user to select an uncalibrated spectrum and another reference spectrum.
- Calibrates each spectrum of input file with the reference waveform, individually.
- Save the results.

## config.json
- Contains adjustable parameters and system configuration.
- Default band range 1 and 2 are B-band and A-band of telluric oxygen absorption lines, respectively.
- If the computation is too heavy for your system or it takes too long, you can adjust the shift resolution.
- *Adjustments should be made with care.*

## Data preparation
- The input and reference spectrum should be normalized before they can be used with this pipeline.
- Both spectrum **MUST** cover the range of band range 1 and 2 as specified in the 'config.json' file.
- Acceptable file formats: *.csv .xlsx .pkl .p (pickle)*

## Reference spectrum
- Must contain exactly two columns, one header row.
- The first column stores wavelength while the second column stores intensity.

## Input (uncalibrated) spectrum
- Contains two or more columns, one header row.
- The first column stores wavelength while other columns store intensities of each different spectrum measurement.

## Results
The results will be saved as a folder inside a selected directory.
1. corrected_intensities.csv
2. corrected_wavelengths.csv
3. correction_coefficients.csv
4. avg_corrected_waveform.csv