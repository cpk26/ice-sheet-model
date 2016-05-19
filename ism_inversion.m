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
vv2.C = ism_cslip_field(vv2, pp, gg, oo);    %Reconstruct basal slipperiness

%% Solve Initial Forward problem
[vv2] = ism_sia(aa.s,aa.h,vv2.C,vv2,pp,gg,oo);  %SIA 
[vv2] = ism_sstream(vv2,aa,pp,gg,oo );          %SSA 


%% Optimization
 %options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton',...
 %    'HessUpdate', 'bfgs','GradObj','on');
 options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton',...
     'HessUpdate', 'steepdesc','GradObj','on', 'TolFun', 1e-15);
 
 if strcmp(oo.inv_meth, 'LM')
 [vv2.acoeff,cst,exitflag,output] = fminunc(@(x)ism_adjLM_optWrapper(x,vv2,aa, pp, gg, oo),vv2.acoeff(:),options);
 elseif strcmp(oo.inv_meth, 'AD')
 ism_adjAD_generate( vv2,aa, pp, gg, oo );
 [vv2.acoeff,cst,exitflag,output] = fminunc(@(x)ism_adjAD_optWrapper(x,vv2,aa, pp, gg, oo),vv2.acoeff(:),options);
 else
 error('Inversion Method not specified')
 end
     
     
%% cs.ubc.ca function minimizer

% options = struct();
% %options.Method = 'cg';
% options.Display = 'iter';
% %options.LS_init = 0;
% 
% options.LS_type = 2;
% [vv2.acoeff,cst,exitflag,output] = minFunc(@(x)ism_invoptWrapper(x,vv2,aa, pp, gg, oo),vv2.acoeff(:),options);

%% Finish up
vv2.C = ism_cslip_field(vv2, pp, gg, oo);    %reconstruct basal slipperiness


fprintf('Inversion Complete \n \n')

end

