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

%% Discretize Basal slipperiness into coefficients                                    
vv2.acoeff = ism_cslip_acoeff(vv, pp, gg, oo);
vv2.C = ism_cslip_field(vv2, pp, gg, oo);    
 
%% Solve Initial Forward problem
[vv2] = ism_sia(aa.s,aa.h,vv2.C,vv2,pp,gg,oo);  %SIA 
[vv2] = ism_deism(vv2,aa,pp,gg,oo );          %SSA 

%% Optimization Options

%% FMINUNC
  options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton',...
  'HessUpdate', 'bfgs','GradObj','on','TolFun', 1e-10);

%% MINFUNC [https://www.cs.ubc.ca/~schmidtm/Software/minFunc.html]
options = struct();
options.Method = 'lbfgs';
options.Display = 'full';
options.LS_init = 2;
options.Corr = 35;
options.MaxIter = 100;
options.optTol = 5e-4;
 
%% Optimization

if strcmp(oo.inv_meth, 'LM')
    
if strcmp(oo.inv_opt,'gd')
[vv2.acoeff,cst,exitflag,output] = ism_steepDesc(@(x)ism_adjLM_optWrapper(x,{},aa, pp, gg, oo),vv2.acoeff(:));
elseif strcmp(oo.inv_opt,'lbfgs')
[vv2.acoeff,cst,exitflag,output] = fminunc(@(x)ism_adjLM_optWrapper(x,vv2,aa, pp, gg, oo),vv2.acoeff(:),options);
else
error('Optimization method not specified')
end  

elseif strcmp(oo.inv_meth, 'AD')
ism_adjAD_generate( vv2,aa, pp, gg, oo );

if strcmp(oo.inv_opt,'gd')
[vv2.acoeff,cst,exitflag,output] = ism_steepDesc(@(x)ism_adjAD_optWrapper(x,{},aa, pp, gg, oo),vv2.acoeff(:));
elseif strcmp(oo.inv_opt,'lbfgs')
[vv2.acoeff,cst,exitflag,output] = minFunc(@(x)ism_adjAD_optWrapper(x,vv2,aa, pp, gg, oo),vv2.acoeff(:),options);
else
error('Optimization method not specified')
end  

%[vv2.acoeff,cst,exitflag,output] = ism_steepDesc(@(x)ism_adjAD_optWrapper(x,{},aa, pp, gg, oo),vv2.acoeff(:));
%[vv2.acoeff,cst,exitflag,output] = fminunc(@(x)ism_adjAD_optWrapper(x,vv2,aa, pp, gg, oo),vv2.acoeff(:),options);
%[vv2.acoeff,cst,exitflag,output] = minFunc(@(x)ism_adjAD_optWrapper(x,vv2,aa, pp, gg, oo),vv2.acoeff(:),options);
else
error('Inversion Method not specified')
end
     
%% Finish up 
vv2.output = output;
vv2.C = ism_cslip_field(vv2, pp, gg, oo);  
[vv2] = ism_deism(vv2,aa,pp,gg,oo );   %Optimized Velocities and reconstructed C       



fprintf('Inversion Complete \n \n')

end

