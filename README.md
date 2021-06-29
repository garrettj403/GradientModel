Gradient Model
==============

*Calculate the surface impedance of a rough metal surface using the Gradient Model (GM)*

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4694878.svg)](https://doi.org/10.5281/zenodo.4694878)

Example: Surface Impedance of Gold
----------------------------------

Surface impedance of gold (conductivity: 4.1e7 S/m) between 280 GHz and 360 GHz:

<p align="center">
<img src="https://raw.githubusercontent.com/garrettj403/GradientModel/main/examples/results/wr3p0-surface-impdance-300k.png" width="800">
</p>

To use these results in HFSS, select `Assign boundary > Impedance...` and then copy/paste these values into the dialog box:

```
50 nm surface roughness:
   Real:       -1.6644e-25 * Freq^2 + 5.1199e-13 * Freq^1 + 6.1107e-02
   Imaginary:  -3.4813e-25 * Freq^2 + 2.2603e-12 * Freq^1 + 7.2291e-02
   ur:         -1.3448e-23 * Freq^2 + 5.0514e-11 * Freq^1 + 3.9464e+00

75 nm surface roughness:
   Real:       -1.5186e-25 * Freq^2 + 6.0379e-13 * Freq^1 + 5.8012e-02
   Imaginary:  -4.5447e-25 * Freq^2 + 3.0960e-12 * Freq^1 + 7.9819e-02
   ur:         -2.5105e-23 * Freq^2 + 9.6127e-11 * Freq^1 + 5.8256e+00

100 nm surface roughness:
   Real:       -1.4427e-25 * Freq^2 + 7.0489e-13 * Freq^1 + 5.5407e-02
   Imaginary:  -5.6445e-25 * Freq^2 + 3.8871e-12 * Freq^1 + 8.8834e-02
   ur:         -3.9840e-23 * Freq^2 + 1.5236e-10 * Freq^1 + 8.0508e+00
```

References
----------

Gradient model:

   - G. Gold and K. Helmreich, “A Physical Surface Roughness Model and Its Applications,” IEEE Trans. Microw. Theory Tech., vol. 65, no. 10, pp. 3720–3732, Oct. 2017, doi: [10.1109/TMTT.2017.2695192](https://doi.org/10.1109/TMTT.2017.2695192).

   - K. Lomakin, G. Gold, and K. Helmreich, “Analytical Waveguide Model Precisely Predicting Loss and Delay Including Surface Roughness,” IEEE Trans. Microw. Theory Tech., vol. 66, no. 6, pp. 2649–2662, Jun. 2018, doi: [10.1109/TMTT.2018.2827383](https://doi.org/10.1109/TMTT.2018.2827383).

Closed-form solution:

   - D. N. Grujic, “Closed-Form Solution of Rough Conductor Surface Impedance,” IEEE Trans. Microw. Theory Tech., vol. 66, no. 11, pp. 4677–4683, 2018, doi: [10.1109/TMTT.2018.2864586](https://doi.org/10.1109/TMTT.2018.2864586).
