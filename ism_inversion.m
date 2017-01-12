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
[vv2] = ism_deism(vv2,aa,pp,gg,oo );          %SSA 

%% Initial Cost
[cst] = ism_inv_cost(vv2.uv,vv2.Cb,vv2.alpha,vv2.F1,vv2.F2,vv2,aa,pp,gg, oo);

%% Optimization Options

%% MINFUNC [https://www.cs.ubc.ca/~schmidtm/Software/minFunc.html]
options = struct();
options.Method = 'lbfgs';
options.Display = 'full';
options.LS_type = 0; %Armijo Backtracking
options.Corr = 45;
options.MaxIter = oo.inv_iter;
options.MaxFunEvals = oo.inv_funcEval;
options.progTolFunc = 0.000001;
 
%% Optimization


if strcmp(oo.inv_meth, 'AD')
%ism_adjAD_generate( vv2,aa, pp, gg, oo );
ism_adjAD_generate2( vv2,aa, pp, gg, oo );

disp(['Initial Cost: ', num2str(cst,'%10.2e\n')])
if strcmp(oo.inv_opt,'gd')
[vv2.acoeff,cst,exitflag,output] = ism_steepDesc(@(x)ism_adjAD_optWrapper(x,vv2,aa, pp, gg, oo),vv2.acoeff(:));
elseif strcmp(oo.inv_opt,'lbfgs')
[vv2.acoeff,cst,exitflag,output] = minFunc_cpk(@(x)ism_adjAD_optWrapper(x,vv2,aa, pp, gg, oo),vv2.acoeff(:),options);
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
vv2.uv = zeros(gg.nua + gg.nva,1);
if oo.hybrid
vv2.Cb = ism_slidinglaw(vv2.alpha,vv2.uv,vv2.Cb,vv2.F2,vv2,aa,pp,gg,oo);
else vv2.Cb = ism_slidinglaw(vv2.alpha,vv2.uv,vv2.Cb,[],vv2,aa,pp,gg,oo); end;

[vv2] = ism_deism(vv2,aa,pp,gg,oo );   %Optimized Velocities and reconstructed C       



fprintf('Inversion Complete \n \n')

end

