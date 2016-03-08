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

if ~isfield(oo,'inv_iter'), oo.inv_iter = 40; end;    %Number of iterations for inversion

vv2 = struct();

%% Discretize Basal slipperiness into coefficients                                    
vv2.acoeff = ism_cslip_acoeff(vv, pp, gg, oo);

%% Initial Basal Slipperiness Field   
vv2.C = ism_cslip_field(vv2, pp, gg, oo);    %reconstruct basal slipperiness

%% Solve forward problem for initial guess
[vv2] = ism_sia(aa.s,aa.h,vv2.C,vv2,pp,gg,oo);  %SIA 
[vv2] = ism_sstream(vv2,aa,pp,gg,oo );      %SSA 


inv_norm = zeros(oo.inv_iter+1,1);            %Store velocity misfit
inv_norm(1) = ism_inv_cost(vv2,aa,pp,gg, oo);

%% Iterate
for j = 1:oo.inv_iter
fprintf('Inversion iteration: %i of %i \n',[j,oo.inv_iter])

%% Line search step 
[vv2] = ism_acoeff_ls(vv2,aa, pp, gg, oo);
if isfield(vv2,'armflag'); break; end

inv_norm(j+1) = ism_inv_cost(vv2,aa,pp,gg, oo);
end

vv2.inv_norm = inv_norm;
fprintf('Inversion Complete \n \n')

end
