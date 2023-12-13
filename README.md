# Extrapolate_ADG_VIStoUV_Distribution
---

The ADG extrapolation model extrapolates the chromophoric dissolved organic matter (CDOM) absorption coefficient, a<sub>g</sub>(λ), depigmented particulate absorption coefficient, a<sub>d</sub>(λ), and the non-phytoplankton absorption coefficient, a<sub>dg</sub>(λ), from the visible (VIS) to the ultraviolet (UV) electromagnetic spectrum. The model utilizes inputs of a<sub>g</sub>(λ) and associated a<sub>d</sub>(λ) spectra in the VIS to extrapolate to the UV. The complete development and validation of the ADG extrapolation model is described in [Kehrli et al., 2023](https://opg.optica.org/oe/fulltext.cfm?uri=oe-31-11-17450&id=530517). The presented version of the ADG extrapolation model source code is in MATLAB file format. In this version of the code, we introduced some changes with the primary purpose of streamlining the structure of the code and to facilitate its application by users.

This README document provides information about the files within the Extrapolate_ADG_VIStoUV_Distribution repository.

---

## Extrapolate_ADG_VIStoUV.m:
ADG extrapolation model to estimate a<sub>g</sub>(λ), a<sub>d</sub>(λ), and a<sub>dg</sub>(λ) in the near-UV from 350 to 399 nm for a given pair of a<sub>g</sub>(λ) and a<sub>d</sub>(λ) spectra in the VIS. Returns output wavelength (lambda_out), a<sub>g</sub>(λ) (agout), a<sub>d</sub>(λ) (adout), and a<sub>dg</sub>(λ) (adgout). See supporting documentation for further details.

## ADG_cor_LUT.mat:
Look-up tables (LUTs) necessary to run and perform corrections in ADG_ext.m. The structure contains two fields to perform the short-wavelength low signal a<sub>g</sub>(λ) correction, and the a<sub>d</sub>(λ) bias corrections described in Tables 3 and 4 of the Kehrli et al., 2023. See Extrapolate_ADG_VIStoUV.m.m function documentation for further details about the .mat file.

## Extrapolate_ADG_VIStoUV_Test.m:
MATLAB script that executes Extrapolate_ADG_VIStoUV.m on 1 sample input of spectral a<sub>g</sub>(λ) and a<sub>d</sub>(λ) defined from 400 to 700 nm. The code extrapolates a<sub>g</sub>(λ), a<sub>d</sub>(λ), and a<sub>dg</sub>(λ) from 350 to 399 nm in 1 nm increments.

## Extrapolate_ADG_VIStoUV_Test_Run.xls:
Microsoft Excel spreadsheet containing the input and resulting output data obtained from the application of the Extrapolate_ADG_VIStoUV.m on 1 sample input to determine a<sub>g</sub>(λ). a<sub>d</sub>(λ) and a<sub>dg</sub>(λ) in the near-UV. The file is the original output file generated by Extrapolate_ADG_VIStoUV_Test.m. 

---

Contributors: Matthew Kehrli, Dariusz Stramski, Rick A. Reynolds, and Ishan Joshi\
Contacts: Matthew Kehrli (mdkehrli@ucsd.edu | mdkehrli@gmail.com), Dariusz Stramski (dstramski@ucsd.edu), Rick Reynolds(rreynolds@ucsd.edu), Ishan Joshi(dstramski@ucsd.edu)\
Ocean Optics Research Laboratory, Scripps Institution of Oceanography, University of California San Diego

---
