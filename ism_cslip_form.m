function [ C2 ] = ism_cslip_form(flag,F2,C,aa,pp,gg,oo )
%ism_cslip_form Switch between Cslip effective, and Cslip basal; Needed for
%hybrid model.
% Inputs:
%   flag    Set to 1 to return B_eff, 0 for B_basal
%   F2      Falpha integral
%   vv      struct containing initial solution variables
%   aa      prescribed fields, including inputs and boundary conditions
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs:
%   C       Cslip parameter

if flag
    C2 = C(:)./(1 + (pp.c13*C(:)).*F2);             %Cslip basal -> Cslip effective
else
    C2 = C(:)./(1 - (pp.c13*C(:)).*F2);             %Cslip effective -> Cslip basal
end

end

