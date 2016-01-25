function [vv] = ism_cslip_idct(vv,pp,gg,oo )
%% Calculate basal slipperiness from alpha coeff and 2d sin/cos basis
% Inputs:
%   vv      solution variables
%   gg      grid and operators
%   oo      options
% Outputs: 
%   vv      updated basal slipperiness in solution struct

n_x = vv.n_x;       %Number of coefficients (x/y dirs)
n_y = vv.n_y;
Lx = gg.Lx;     %Length of Domain (x/y dimensions)
Ly = gg.Ly;


B = zeros(gg.nJ*gg.nI,1);

%% Construct basal slipperiness from DCT coefficients
for j = 0:n_x
for k = 0:n_y
aInd = (j+1) + k*(n_x+1);
B = B + vv.acoeff(aInd)*cos(pi*j*gg.xx(:)/Lx).*cos(pi*k*gg.yy(:)/Ly);
end, end

vv.C = B.^2;
end
