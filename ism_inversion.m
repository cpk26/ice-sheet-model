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

vv2 = vv;

%% Discretize Basal slipperiness into coefficients                                    
vv2.acoeff = ism_cslip_acoeff(vv, pp, gg, oo);


% %% Initial Basal Slipperiness Field   
vv2.C = ism_cslip_field(vv2, pp, gg, oo);    %Reconstruct basal slipperiness
 
% %% Solve Initial Forward problem
[vv2] = ism_sia(aa.s,aa.h,vv2.C,vv2,pp,gg,oo);  %SIA 
[vv2] = ism_deism(vv2,aa,pp,gg,oo );          %SSA 

% %Initial Cost
% if oo.hybrid
% F1 = ism_falpha(1,vv2.nEff,vv2,aa,pp,gg,oo );          %Calculate F alpha factors 
% F2 = ism_falpha(2,vv2.nEff,vv2,aa,pp,gg,oo );
% else
% F1 =[]; F2 = [];
% end
% 
% cst1 = ism_inv_cost(vv2.U,gg.S_h*vv2.C(:),F1,F2,vv2,aa,pp,gg, oo);



%% Optimization
  options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton',...
  'HessUpdate', 'bfgs','GradObj','on','TolFun', 1e-10);
% options = optimoptions(@fminunc, 'Display','iter','Algorithm','quasi-newton',...
%     'HessUpdate', 'steepdesc','GradObj','on','TolFun', 1e-10);
%  
 if strcmp(oo.inv_meth, 'LM')
 %[vv2.acoeff,cst,exitflag,output] = fminunc(@(x)ism_adjLM_optWrapper(x,vv2,aa, pp, gg, oo),vv2.acoeff(:),options);
 [vv2.acoeff,cst,exitflag,output] = ism_steepDesc(@(x)ism_adjLM_optWrapper(x,{},aa, pp, gg, oo),vv2.acoeff(:));

 elseif strcmp(oo.inv_meth, 'AD')
 ism_adjAD_generate( vv2,aa, pp, gg, oo );
 [vv2.acoeff,cst,exitflag,output] = ism_steepDesc(@(x)ism_adjAD_optWrapper(x,{},aa, pp, gg, oo),vv2.acoeff(:));
  %[vv2.acoeff,cst,exitflag,output] = fminunc(@(x)ism_adjAD_optWrapper(x,vv2,aa, pp, gg, oo),vv2.acoeff(:),options);
 else
 error('Inversion Method not specified')
 end
     
     
% %% cs.ubc.ca function minimizer
% 
% options = struct();
% options.Method = 'lbgs';
% options.Display = 'excessive';
% options.optTol = 1e-10;
% options.progTol = 1e-15;
% options.optTol = 1e-10;
% 
% % [vv2.acoeff,cst,exitflag,output] = minFunc(@(x)ism_invoptWrapper(x,vv2,aa, pp, gg, oo),vv2.acoeff(:),options);
% 
% ism_adjAD_generate( vv2,aa, pp, gg, oo );
% [vv2.acoeff,cst,exitflag,output] = minFunc(@(x)ism_adjAD_optWrapper(x,vv2,aa, pp, gg, oo),vv2.acoeff(:),options);

%% Finish up
vv2.output = output;
vv2.C = ism_cslip_field(vv2, pp, gg, oo);    %reconstruct basal slipperiness


fprintf('Inversion Complete \n \n')

end

