function [D_diag] = ism_dim_Ddiag(C,nEff,aa,pp,gg,oo)
%% Field Equations for SSA velocities 
% Inputs:
%   u,v     velocities in the x,y directions
%   aa      prescribed fields, including inputs and boundary conditions
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs: Ax = b
%   LHS     Left hand side. 
%   RHS     Right hand side.


%% Determine Basal Slipperines on U/V grids 

Cslip_u = (gg.c_hu*C)./(gg.c_hu*(C > 0));            %Slipperiness on u,v grids
Cslip_v = (gg.c_hv*C)./(gg.c_hv*(C > 0));

Cslip_u(logical(gg.S_u*gg.nmgn_ugrid(:))) = 0;                           %Slipperiness is zero at margin for BC
Cslip_v(logical(gg.S_v*gg.nmgn_vgrid(:))) = 0;

%% Build the diagonal of the D matrix
h = gg.S_h*aa.h(:);
D_diag = [3*nEff.*h; nEff.*h; nEff.*h; -pp.c3*Cslip_u; -pp.c3*Cslip_v];



end

