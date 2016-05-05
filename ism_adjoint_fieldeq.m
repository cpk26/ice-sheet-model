function [LHS, RHS] = ism_adjoint_fieldeq(vv,aa,pp,gg,oo)
%% Field Equations for Adjoint  
% Inputs:
%   u,v     velocities in the x,y directions
%   aa      prescribed fields, including inputs and boundary conditions
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs: Ax = b
%   LHS     Left hand side. 
%   RHS     Right hand side.

if ~isfield(oo,'inv_msft'), oo.inv_msft = 'abs'; end 

n = pp.n_Glen;                          

%% Variables (Non-Dimensionalized)
nha = sum(gg.S_h(:));                               %number of active h-grid nodes
nua = sum(gg.S_u(:));
nva = sum(gg.S_v(:));

Cslip_u = (gg.c_hu*gg.S_h*vv.C(:))./(gg.c_hu*gg.S_h*(vv.C(:) > 0));            %Slipperiness on u,v grids
Cslip_v = (gg.c_hv*gg.S_h*vv.C(:))./(gg.c_hv*gg.S_h*(vv.C(:) > 0));

Cslip_u(logical(gg.S_u*gg.nmgn_ugrid(:))) = 0;                           %Slipperiness is zero at margin for BC
Cslip_v(logical(gg.S_v*gg.nmgn_vgrid(:))) = 0;


h_diag = spdiags(gg.S_h*aa.h(:),0,nha,nha);         %Diagonalize  
Cslip_udiag = spdiags(Cslip_u(:),0,nua,nua);
Cslip_vdiag = spdiags(Cslip_v(:),0,nva,nva);       

exx = gg.du_x*vv.u;                                %Strain Rates
eyy = gg.dv_y*vv.v;
exy = 0.5*(gg.dhu_y*vv.u + gg.dhv_x*vv.v);
edeff = sqrt(exx.^2 + eyy.^2 + exx.*eyy + exy.^2 + pp.n_rp.^2);

nEff =  edeff.^((1-n)/n);                       %Effective Viscosity [dimensionless]
nEff_diag = spdiags(nEff(:),0,nha,nha);                                  


%% Field equations for lambda and mu
X = [gg.du_x gg.dv_y; gg.du_x -gg.dv_y; gg.dhu_y gg.dhv_x; speye(nua,nua) sparse(nua,nva); sparse(nva,nua) speye(nva,nva)];
X2 = [gg.dh_x gg.dh_x gg.duh_y speye(nua,nua) sparse(nva,nua)'; gg.dh_y -gg.dh_y gg.dvh_x sparse(nua,nva)' speye(nva,nva)];
D = blkdiag(3*pp.c6*nEff_diag*h_diag, pp.c6*nEff_diag*h_diag, pp.c6*nEff_diag*h_diag, -Cslip_udiag, -Cslip_vdiag);

LHS = X2*D*X;

if strcmp(oo.inv_msft,'abs')                        %RHS Adjoint equations
E1 = pp.L_vel*pp.c7*(aa.u - vv.u);                                
E2 = pp.L_vel*pp.c7*(aa.v - vv.v);
elseif strcmp(oo.inv_msft,'rel')
E1 = (pp.L_vel/pp.c7)*(aa.u - vv.u)./(aa.u.^2);                               
E2 = (pp.L_vel/pp.c7)*(aa.v - vv.v)./(aa.v.^2);
end

RHS = [E1; E2];


end