function [acoeff] = ism_cslip_acoeff(vv, pp, gg, oo)
%% Calculate Basal Slipperiness alpha coefficients
% Inputs:
%   vv      solution variables
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs:
%   vv       solution variables

if ~isfield(pp,'acoeff_nx'), pp.acoeff_nx = 45; end;    %Number of terms in x,y directions to keep in DCT2
if ~isfield(pp,'acoeff_ny'), pp.acoeff_ny = 45; end;    


aField = reshape(log(vv.C),gg.nJ,gg.nI);                %Set slipperiness, C = exp(a(:)*F(:))
aField(vv.C == 0) = mean(gg.S_h*log(vv.C(:)));          %Remove INF values outside of mask

if isequal(oo.Cdisc, 'dct2')
%% Discretize Basal slipperiness using discrete cosine series
n_x = pp.acoeff_nx;   %Number of terms in x,y directions to keep in DCT2
n_y = pp.acoeff_ny; 

AA = zeros(gg.nJ,gg.nI); AA(1:n_y,1:n_x) = 1; %Create Mask

acoeff = dct2(aField,gg.nJ,gg.nI);                      %2D DCT, removing high frequency terms
acoeff(~AA) = 0;                                         

elseif isequal(oo.Cdisc, 'grid')
acoeff = aField;    
end



end