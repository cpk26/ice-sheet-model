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

if ~isfield(vv,'inv_iter'), vv.inv_iter = 50; end;        %Initial step coefficient

vv2 = struct();

%% Discretize Basal slipperiness using discrete cosine series
%n_x = ceil(log2(gg.Lx/pp.c5)); vv2.n_x = n_x;   %Number of terms in x,y directions to keep
%n_y = ceil(log2(gg.Ly/pp.c5)); vv2.n_y = n_y; 

vv2.n_x = 20;   %TEST
vv2.n_y = 20;
n_x = 20;
n_y = 20;
AA = zeros(gg.nJ,gg.nI); AA(1:n_x,1:n_y) = 1;

acoeff = dct2(reshape(vv.C,gg.nJ,gg.nI));          %Calculate 2d DCT
acoeff(~AA) = 0;                                %Keep a subset of the terms
vv2.acoeff = acoeff;


%% Initial Basal Slipperiness Field   
vv2.C = idct2(acoeff);                      %reconstruct basal slipperiness

%% Solve forward problem for initial guess
[ii] = ism_sia(aa.s,aa.h,vv2.C, pp,gg,oo);  %SIA 
vv2.u = ii.u; vv2.v = ii.v;
[vv2] = ism_sstream(vv2,aa,pp,gg,oo );      %SSA 

inv_iternorm = zeros(vv.inv_iter,1);

%% Iterate
for j = 1:vv.inv_iter
fprintf('Inversion iteration: %i of %i \n',[j,vv.inv_iter])
%% solve for mu, lambda
[vv2] = ism_adjoint(vv2,aa,pp,gg,oo );

%% solve for gradient
[vv2] = ism_cslip_grad(vv2, pp, gg, oo);

%% newton-raphson step 
[vv2] = ism_acoeff_linesearch(vv2,aa, pp, gg, oo);

inv_iternorm(j) = ism_inversion_misfit(vv2.u,vv2.v,aa,pp,gg, oo);
end


pause()


end
