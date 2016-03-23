
function [cst,grad] = ism_invoptWrapper(acoeff,vv,aa, pp, gg, oo)
% Inputs:
%   vv      struct containing initial solution variables
%   aa      prescribed fields, including inputs and boundary conditions
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs:
%   vv2     updated struct with new alpha coefficients

vv.acoeff = acoeff;

[vv] = ism_sia(aa.s,aa.h,vv.C,vv, pp,gg,oo);   %SIA                                        
[vv] = ism_sstream(vv,aa,pp,gg,oo );          %SSA 
[vv] = ism_cost_jac(vv, pp, gg, oo);

grad = vv.cJac;                      %Gradient
cst = ism_inv_cost(vv,aa,pp,gg, oo); %Current misfit


end   

