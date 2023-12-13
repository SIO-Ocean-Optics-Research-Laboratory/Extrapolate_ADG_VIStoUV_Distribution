function [lambda_out,agout,adout,adgout] = Extrapolate_ADG_VIStoUV(lambda,ag,ad,ADG_cor_LUT)
%Implements the extrapolation model to extend the chromophoric dissolved
%organic matter (CDOM) absorption coefficient(ag), depigmented (also
%referred to as non-algal) particulate absorption coefficient (ad), and
%non-phytoplankton absorption coefficient (adg) to the near-UV spectral
%range of 350 to 400 nm in 1 nm increments from input of spectral ag and ad
%(and hence adg as well) that spans 400 to 450 nm in the visible spectral
%region.
%
%Reference:
%
%Kehrli, M. D., Stramski, D., Reynolds, R. A., & Joshi, I. D. (2023).
%Estimation of chromophoric dissolved organic matter and non-algal
%particulate absorption coefficients of seawater in the ultraviolet by
%extrapolation from the visible spectral region. Optics Express, 31(11),
%17450. https://doi.org/10.1364/OE.486354
%
%Required function inputs: lambda, ag, ad, adg_cor_LUT
%   lambda [m-by-1 numeric]: Input values of light wavelength [nm]
%   corresponding to the input values of spectral absorption coefficients. 
%
%   ag [m-by-1 numeric]: Input values of spectral absorption coefficient of
%   CDOM [m^-1]. 
%
%   ad [m-by-1 numeric]: Input values of spectral absorption coefficient of
%   depigmented particulate matter [m^-1]. 
%
%   ADG_cor_LUT [1-by-1 structure]: Structure containing the required
%   look-up tables (LUTs) to perform corrections described in Step 3 of
%   Figure 9 from Kehrli et al., 2023; loaded via load('adg_cor_LUT.mat').
%
%	    ADG_cor_LUT.agcor: LUT with coefficients and weights from Table 3 of
%	    Kehrli et al., 2023.
%
%	    ADG_cor_LUT.adcor: LUT with MdR (median ratio) values from Table 4 of
%	    Kehrli et al., 2023.
%
%Outputs: lambda_out, agout, adout, adgout
%   lambda_out [50-by-1 double]: Light wavelengths [nm] corresponding to
%   output values of spectral absorption coefficients ag, ad, and adg from
%   the partitioning model. Output light wavelengths span from 350 to 399
%   nm with a 1 nm spectral interval.
%   
%   agout [50-by-1 double]: An array containing the estimated spectral
%   values of the CDOM absorption coefficient, ag, in the UV. The output
%   spectrum of ag spans from 350 to 399 nm with a 1 nm spectral
%   interval.
%
%   adout [50-by-1 double]: An array containing the estimated spectral
%   values of the depigmented particulate absorption coefficient, ag, in
%   the UV. The output spectrum of ad spans from 350 to 399 nm with a 1 nm
%   spectral interval.
%
%   adgout [50-by-1 double]: An array containing the estimated spectral
%   values of the non-phytoplankton absorption coefficient, adg, in the UV.
%   The output spectrum of adg spans from 350 to 399 nm with a 1 nm
%   spectral interval.
%
%Version 1.0 (v1.0)
%
%Version history: 
%2023-11-09: Revised ADG extrapolation model and Matlab version, M. D.
%Kehrli, D. Stramski, R. A. Reynolds, I. D. Joshi.
%202X-XX-XX: Final revised MATLAB version (v1.0), M. Kehrli, D. Stramski,
%R. A. Reynolds, and I. D. Joshi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check input arguments
    arguments
        lambda (:,1) {mustBeNumeric}
        ag (:,1) {mustBeNumeric}
        ad (:,1) {mustBeNumeric}
        ADG_cor_LUT (1,1) struct
    end

%input adg
adg = ag + ad;

%spectral range in the visible (blue) region for fitting the exponential
%function
range = [400 450];
%spectral range with wavelengths at 1 nm intervals for UV extrapolation
UVLambda = (350:399)';

%perform initial ag and adg extrapolations described in Step 2 of Figure 9
%from Kehrli et al., 2023

%ag extrapolation from visible to 350 nm (1 nm intervals)
[ag_linextrap] = linearextrap(lambda, log(ag), range, UVLambda);
%adg extrapolation from visible to 350 nm (1 nm intervals)
[adg_linextrap] = linearextrap(lambda, log(adg), range, UVLambda);

%convert initial output of ag and adg from log-transformed space
ag_init = exp(ag_linextrap);
adg_init = exp(adg_linextrap);

%calculate initial output of ad (before application of bias correction)
ad_init = adg_init-ag_init;

%load LUTs with coefficients and weights for ag (Table 3 in Kehrli et al.,
%2023) and the values of MdR for ad (Table 4 in Kehrli et al., 2023) which
%are used for correction of the initial extrapolation results
agcor = ADG_cor_LUT.agcor;
adcor = ADG_cor_LUT.adcor;

%store the coefficients for ag correction (Table 3 in Kehrli et al., 2023)
%to workspace
alpha = agcor.alpha;
beta = agcor.beta;
w1 = agcor.w1;
w2 = agcor.w2;

%apply ag correction Equation (2) in Kehrli et al., 2023
y = 1 + alpha.^(log10(ag_init)-beta);
%equation (3) in Kehrli et al., 2023
CF = w1.*y+w2;
ag_cor = ag_init.*CF;

%store the median ratio MdR for correction of ad (Table 4 in Kehrli et al.,
%2023)
%to workspace
MdRadUV = adcor.MdR;

%apply ad correction; see Fig. 9 in Kehrli et al., 2023
ad_cor = ad_init.*(1./MdRadUV);

%add corrected ad and ag to estimate adg in the UV
adg_cor = ad_cor + ag_cor;

%function output
lambda_out = UVLambda;
agout = ag_cor;
adout = ad_cor;
adgout = adg_cor;
end
%end of main code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Ancillary subfunction
function [extrapout] = linearextrap(lambda,loga,lamrange,lamextrap)
%Compute a linear regression over a specified wavelength range to
%extrapolate to a desired output wavelength range.
%
%Inputs: lambda, loga, lamrange, lamextrap
%   lambda [m-by-1 numeric]: Input values of light wavelength [nm]
%   corresponding to the input values of spectral absorption coefficient.
%
%   loga [m-by-1 numeric]: Input values of log-e transformed spectral
%   absorption coefficient.
%
%   lamrange [1-by-2 numeric]: Wavelength range over which to calculate
%   linear regression of loga vs. lambda (e.g. [400 450] spectral values
%   which range from 400 to 450 nm).
%
%   lamextrap [n-by-1 numeric]: Wavelengths in the spectral extrapolation
%   range (e.g. [350:1:400]).
%
%Outputs:
%   extrapout [n-by-1 numeric]: Extrapolated values of the absorption
%   coefficient in log-e space using a 'poly1' matlab regression applied to
%   log-e transformed spectral values of the input absorption coefficient
%   corresponding to lamrange.
%
%Created: July 28, 2023
%Completed: July 28, 2023
%Updates: October 27, 2023 - Removed uncertainty and output coefficients
%produced by the linear extrapolation to simplify the function.
%Additionally set to match boundary wavelength of 400 nm only if values of
%the log-e transformed absorption coefficient are available at 400 nm.
%
%Matthew Kehrli Ocean Optics Research Laboratory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%indices of input wavelengths bound within lamrange
lineidx = find(lambda >= lamrange(1) & lambda <= lamrange(2));

%check if input absorption spectrum is nan; set function output to nan if
%true, otherwise perform regression
if isnan(loga(lambda==400))
    extrapout = nan(numel(lambda),1);
else
    %perform regression routine with 1-degree polynomial
    fitag = fit(lambda(lineidx),loga(lineidx),'poly1');
    
    %store the extrapolation output before potentially matching values at
    %400 nm if available
    extrap_pre_shift = fitag(lamextrap);

    %set up shift of extrapolation to match input absorption coefficient at
    %400 nm if available
    if lambda(lineidx(1)) == 400
        shiftidx = lamrange(1);
        diff = fitag(shiftidx) - loga(lineidx(1));
    
        %shift all values of extrapolation
        extrapout = extrap_pre_shift - diff;
    end
end
end