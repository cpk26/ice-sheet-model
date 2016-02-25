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

if ~isfield(oo,'inv_iter'), oo.inv_iter = 15; end;    %Number of iterations for inversion

vv2 = struct();

%% Discretize Basal slipperiness using discrete cosine series
vv2.n_x = ceil(log2(gg.Lx/pp.c5)); vv2.n_x = 45;   %Number of terms in x,y directions to keep
vv2.n_y = ceil(log2(gg.Ly/pp.c5)); vv2.n_y = 45; 
AA = zeros(gg.nJ,gg.nI); AA(1:vv2.n_y,1:vv2.n_x) = 1;


acoeff = dct2(reshape(sqrt(vv.C),gg.nJ,gg.nI));  %2D DCT, removing high frequency terms
acoeff(~AA) = 0;                                         
vv2.acoeff = acoeff;


%% Initial Basal Slipperiness Field   
vv2.C = idct2(acoeff).^2;                   %reconstruct basal slipperiness

%% Solve forward problem for initial guess
[vv2] = ism_sia(aa.s,aa.h,vv2.C,vv2,pp,gg,oo);  %SIA 
[vv2] = ism_sstream(vv2,aa,pp,gg,oo );      %SSA 


inv_norm = zeros(oo.inv_iter+1,1);            %Store velocity misfit
inv_norm(1) = ism_vel_misfit(vv2.u,vv2.v,aa,pp,gg, oo);

%% Iterate
for j = 1:oo.inv_iter
fprintf('Inversion iteration: %i of %i \n',[j,oo.inv_iter])
%% solve for mu, lambda
[vv2] = ism_adjoint(vv2,aa,pp,gg,oo );

%% solve for gradient
[vv2] = ism_cslip_grad(vv2, pp, gg, oo);

%% newton-raphson step 
[vv2] = ism_acoeff_nrs(vv2,aa, pp, gg, oo);
if isfield(vv2,'armflag'); break; end

inv_norm(j+1) = ism_vel_misfit(vv2.u,vv2.v,aa,pp,gg, oo);
end

vv2.inv_norm = inv_norm;
fprintf('Inversion Complete \n \n')

end
