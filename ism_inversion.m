function [vv2] = ism_inversion(vv,aa,pp,gg,oo )
%% Inversion using Depth Integrated Ice Sheet Model (DEISM)
% Inputs:
%   vv      struct containing initial solution variables
%   aa      prescribed fields, including inputs and boundary conditions
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs:
%   vv2     struct containing new solution variables

%% Initialize
vv2 = vv;

%% 

%% Discretize Basal slipperiness into coefficients                                    
vv2.acoeff = ism_alpha_acoeff(vv.alpha, vv, pp, gg, oo);

% if oo.hybrid;                                           %Rediscretize [Unnecessary if removing fourier discretization]
% Cb = ism_cslip_field(vv2, pp, gg, oo); 
% C = Cb(:)./(1 + (pp.c13*Cb(:)).*(gg.S_h'*vv.F2));
% vv2.Cb = Cb; vv2.C = C; 
%     
% else vv2.C = ism_cslip_field(vv2, pp, gg, oo);   end

 
%% Solve Initial Forward problem
%[vv2] = ism_sia(aa.s,aa.h,vv2.C,vv2,pp,gg,oo);  %SIA [can remove]
[vv2] = ism_deism(vv2,aa,pp,gg,oo );          %SSA 

%% Initial Cost
if oo.hybrid
F1 = ism_falpha(1,vv2.uv,vv2.nEff_lyrs,vv2,aa,pp,gg,oo );          %Calculate F alpha factors [Hybrid]
F2 = ism_falpha(2,vv2.uv,vv2.nEff_lyrs,vv2,aa,pp,gg,oo );
cst = ism_inv_cost(vv2.uv,vv2.Cb,vv2.alpha,F1,F2,vv2,aa,pp,gg, oo);  %Current misfit

else
cst = ism_inv_cost(vv2.uv,vv2.Cb,vv2.alpha,[],[],vv2,aa,pp,gg, oo);  %Current misfit
end


%% Optimization Options

%% MINFUNC [https://www.cs.ubc.ca/~schmidtm/Software/minFunc.html]
options = struct();
options.Method = 'lbfgs';
options.Display = 'full';
options.LS_init = 3; %Double previous stepsize for intitial step guess
options.LS_type = 1; %Custom Armijo backtracking
options.Corr = 35;
options.MaxIter = oo.inv_iter;
options.MaxFunEvals = oo.inv_funcEval;
options.progTol = 1e-7;
 
%% Optimization


if strcmp(oo.inv_meth, 'AD')
ism_adjAD_generate( vv2,aa, pp, gg, oo );

if strcmp(oo.inv_opt,'gd')
[vv2.acoeff,cst,exitflag,output] = ism_steepDesc(@(x)ism_adjAD_optWrapper(x,vv2,aa, pp, gg, oo),vv2.acoeff(:));
elseif strcmp(oo.inv_opt,'lbfgs')
[vv2.acoeff,cst,exitflag,output] = minFunc(@(x)ism_adjAD_optWrapper(x,vv2,aa, pp, gg, oo),vv2.acoeff(:),options);
else
error('Optimization method not specified')
end  

else
error('Inversion Method not specified')
end
     
%% Finish up 
vv2.output = output;
% if oo.hybrid, vv2.Cb = ism_cslip_field(vv2, pp, gg, oo); 
% else vv2.C = ism_cslip_field(vv2, pp, gg, oo); end;

vv2.alpha = ism_alpha_field(vv2.acoeff,vv2, pp, gg, oo);
if oo.hybrid
vv2.Cb = ism_slidinglaw(vv2.alpha,vv2.uv,vv2.Cb,vv2.F2,vv2,aa,pp,gg,oo);
else vv2.Cb = ism_slidinglaw(vv2.alpha,vv2.uv,vv2.Cb,[],vv2,aa,pp,gg,oo); end;

[vv2] = ism_deism(vv2,aa,pp,gg,oo );   %Optimized Velocities and reconstructed C       



fprintf('Inversion Complete \n \n')

end

