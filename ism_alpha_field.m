function [alpha] = ism_alpha_field(acoeff,vv, pp, gg, oo)
%% Calculate Basal Slipperiness from alpha coefficients
% Inputs:
%   vv      solution variables
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs:
%   vv       solution variables

%reconstruct alpha parameter       
alpha = exp(acoeff);                   



end