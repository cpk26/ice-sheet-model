
function [cst,gradN] = ism_adjLM_optWrapper(acoeff,vv,aa, pp, gg, oo)
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

[vv] = ism_sia(aa.s,aa.h,vv.C,vv, pp,gg,oo);    %SIA                                        
[vv] = ism_deism(vv,aa,pp,gg,oo );            %SSA 
cst = ism_inv_cost(vv.U,[],[],[],vv,aa,pp,gg, oo);            %Current misfit

if nargout > 1 % gradient required
    [vv] = ism_adjLM_main(vv,aa,pp,gg,oo );
    [vv] = ism_adjLM_costJac(vv, pp, gg, oo);
    grad = vv.cJac(:);                      %Gradient (as array)
    gradN = grad./max(abs(grad));
end

end   

