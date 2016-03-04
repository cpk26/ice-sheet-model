function [Cslip] = ism_cslip_field(vv, pp, gg, oo)
%% Calculate Basal Slipperiness from alpha coefficients
% Inputs:
%   vv      solution variables
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs:
%   vv       solution variables

%reconstruct basal slipperiness
if isequal(oo.Cdisc, 'dct2')
Cslip = exp(idct2(vv.acoeff));                   

elseif isequal(oo.Cdisc, 'grid')
Cslip = exp(vv.acoeff);                   
   
end



end