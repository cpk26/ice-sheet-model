
function [tau,sdir] = ism_acoeff_minalg(cst,vv,aa, pp, gg, oo)
%% Function which calculates step size and direction to minimize acoeff.
% Inputs:
%   vv      struct containing initial solution variables
%   aa      prescribed fields, including inputs and boundary conditions
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs:
%   vv2     updated struct with new alpha coefficients



if strcmp(oo.inv_opt,'dg')
tau = cst/norm(vv.cJac(:),2);  
sdir = -vv.cJac / norm(vv.cJac(:),2);
    
    
elseif strcmp(oo.inv_opt,'lbfgs')
    
    
    
end


end   

