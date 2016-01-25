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


n = pp.n_Glen;       

h = aa.h; h_diag = spdiags(h(:),0,gg.nIJ,gg.nIJ);
Cslip_diag = spdiags(vv.C(:),0,gg.nIJ,gg.nIJ);
u = vv.u; v = vv.v;

exx = gg.du_x*u;                                %Strain Rates
eyy = gg.dv_y*v;
exy = 0.5*(gg.du_y*u + gg.dv_x*v);
edeff = sqrt(exx.^2 + eyy.^2 + exx.*eyy + exy.^2 + pp.n_rp.^2);

nEff =  edeff.^((1-n)/n);        %Effective Viscosity [dimensionless]
nEff_diag = spdiags(nEff(:),0,gg.nIJ,gg.nIJ);     


A1 = gg.S_u*gg.dh_x * (pp.c6 * 4 * nEff_diag*h_diag) * gg.du_x*gg.S_u';   %LHS Adjoint equations
A2 = gg.S_u*gg.dhu_y * (pp.c6 * nEff_diag*h_diag) * gg.du_y*gg.S_u';
A3 = gg.S_u*gg.c_hu * (Cslip_diag) * gg.c_uh*gg.S_u';
AA = A1 + A2 - A3;

B1 = gg.S_u*gg.dh_x * (pp.c6 * 2 * nEff_diag*h_diag) * gg.dv_y*gg.S_v';
B2 = gg.S_u*gg.dhu_y * (pp.c6 * nEff_diag*h_diag) * gg.dv_x*gg.S_v';
BB = B1 + B2;

C1 = gg.S_v*gg.dh_y * (pp.c6 * 2 * nEff_diag * h_diag) * gg.du_x*gg.S_u';
C2 = gg.S_v*gg.dhv_x * (pp.c6 * nEff_diag * h_diag) * gg.du_y*gg.S_u';
CC = C1 + C2;

D1 = gg.S_v*gg.dh_y * (pp.c6 * 4 * nEff_diag * h_diag) * gg.dv_y*gg.S_v';
D2 = gg.S_v*gg.dhv_x * (pp.c6 * nEff_diag * h_diag) *gg.dv_x*gg.S_v';
D3 = gg.S_v*gg.c_hv * (Cslip_diag) * gg.c_vh*gg.S_v';
DD = D1 + D2 - D3;

LHS = [AA BB; CC DD];   

E1 = pp.c7*(aa.u - vv.u);                                %RHS Adjoint equations
E2 = pp.c7*(aa.v - vv.v);
RHS = [E1; E2];


end