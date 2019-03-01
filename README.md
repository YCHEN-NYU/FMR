# FMR
Fitting of rerromagnetic resonance (FMR) data with Lorentzian model
Input Data: FMR Spectrum as a function of the applied magnetic field, microwave frequency and temperature.
Output Data: Resonance field and linewidth (full width at half maximum with 95% confidence intervals).
Fitting Model: Derivative of Lorentzian with real + imaginary parts. Fitting is done by least squared method, with confidence intervals of 95%.
LabVIEW UI interface: Hardware programming interface with LabVIEW in data acquisition in the Quantum Design PPMS setup. Communications are performed by IEEE GPIB interfaces and S-parameters from the Vector Network Analyzer Model E8364B.
