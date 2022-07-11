# Physics Potentials of Supernova Detection

Study the physics potentials like absolute neutrino mass, neutrino mass ordering with supernova neutrino detection in JUNO. The main chain includes:

1. visible event spectra PDF generation;
2. ToyMC datasets generation;
3. Fitting and sensitivity analysis.

## PDF Generation:
### 1. Code Tools
A local numerical simulation package developed by HuiLing Li is modified and used as the PDF generation tool. It can provide information like:

- Initial neutrino luminosity, event rate, average energy;
- neutrino rate at the detector considering flavour conversion and distance scale;
- visible event spectra for different channels in the detector (only proton quenching added)...

The production/nuFlux/snewpy_models.py can generate nuFlux data from SNEWPY website.

### 2. Datasets
1. Garching models: provided directly in the code package (total 32 models);
2. Burrows 2D: from SNEWPY website, different data format, requiring modifications in BurrowsClass;
3. ...

### 3. ROOT File 
Currently the data from generator are in txt files (requiring improvements). We use genPDF.py in the production/PDFs/10kpc path to generate TH1Ds and save them into root files. 

## Data Generation:
Based on the TH1D of PDFs, toyMC datasets can be generated by sampling with genData.py in the production/PDFs/10kpc path.

## Fitting
### 1. Time Scanning
The first try is to scan time and calcualte the likelihood value without any true fitting. The analysis codes locate in analysis/MH/. The scan_NLL_combine.py is the NLL scanner with figures and values outputs. And the sens.py can compare the NLL values of different models.