function [vv2] = ism_inversion(vv,aa,pp,gg,oo )
%% Inversion using Shallow Stream Model 
% Inputs:
%   vv      struct containing initial solution variables
%   aa      prescribed fields, including inputs and boundary conditions
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs:
%   vv2     struct containing new solution variables

vv2 = struct();


%% Discretize Basal slipperiness using fourier series
n_x = ceil(log2(gg.Lx/pp.c5)); vv2.n_x = n_x; %length of series in x-dir
n_y = ceil(log2(gg.Ly/pp.c5)); vv2.n_y = n_y; %length of series in y-dir
vv2.C = vv.C;

[vv2] = ism_cslip_dct(vv2,pp,gg,oo );

%% Initial Basal Slipperiness Field   
[vv2] = ism_cslip_idct(vv2,pp,gg,oo );

%% Solve forward problem for initial guess
[ii] = ism_sia(aa.s,aa.h,vv2.C, pp,gg,oo);  %SIA 
vv2.u = ii.u; vv2.v = ii.v;
[vv2] = ism_sstream(vv2,aa,pp,gg,oo );      %SSA 


%% Iterate
for j = 1:15

%% solve for mu, lambda
[vv2] = ism_adjoint(vv2,aa,pp,gg,oo );

%% solve for gradient
[vv2] = ism_cslip_grad(vv2, pp, gg, oo);

%% update coefficients
[vv2] = ism_acoeff_linesearch(vv2,aa, pp, gg, oo);

%% Update Basal Slipperiness Field   
[vv2] = ism_cslip_basis(vv2, gg, oo);

disp(j)
%% End Loop
end




end
