# Fitting of frerromagnetic resonance (FMR) data by Lorentzian model

**Input Data**: FMR Spectrum as a function of the applied magnetic field, microwave frequency and temperature.

**Output Data**: Resonance field and linewidth (full width at half maximum with 95% confidence intervals).
**Fitting Model**: Derivative of Lorentzian with real + imaginary parts. Fitting is done by least squared method, with confidence intervals of 95%.

![Image](https://github.com/YCHEN-NYU/FMR/blob/master/data/8GHz_060K.png?raw=true)

**LabVIEW UI Interface**: Hardware programming interface with LabVIEW in data acquisition in the Quantum Design PPMS setup. Communications are performed by IEEE GPIB interfaces and S-parameters from the Vector Network Analyzer Model E8364B.
![Image](https://github.com/YCHEN-NYU/FMR/blob/master/ui-LabVIEW/LabVIEW_UI.jpeg)

**PPMS-VNA setups**: Hardware setup for PPMS-FMR probes for thin film FMR experiments.
![Image](https://github.com/YCHEN-NYU/FMR/blob/master/setup-hardware/PPMS-VNA%20setup.jpeg);
