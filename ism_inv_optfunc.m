
function [cst,grad] = ism_inv_optfunc(acoeff,vv,aa, pp, gg, oo)
% Inputs:
%   vv      struct containing initial solution variables
%   aa      prescribed fields, including inputs and boundary conditions
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs:
%   vv2     updated struct with new alpha coefficients

vv.acoeff = reshape(acoeff,gg.nJ,gg.nI);   %Array=>matrix
vv.C = ism_cslip_field(vv, pp, gg, oo);    %reconstruct basal slipperiness

[vv] = ism_sia(aa.s,aa.h,vv.C,vv, pp,gg,oo);   %SIA                                        
[vv] = ism_sstream(vv,aa,pp,gg,oo );          %SSA 
[vv] = ism_adjoint(vv,aa,pp,gg,oo );
[vv] = ism_cost_jac(vv, pp, gg, oo);

cst = ism_inv_cost(vv,aa,pp,gg, oo); %Current misfit

if nargout > 1 % gradient required
    grad = vv.cJac(:);                      %Gradient (as array)
end

end   

