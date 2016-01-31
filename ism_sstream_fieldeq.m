function [LHS, RHS] = ism_sstream_fieldeq(u,v,C,aa,pp,gg,oo)
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


n = pp.n_Glen;                          

%% Variables (Non-Dimensionalized)
s = aa.s; s_diag = spdiags(s(:),0,gg.nIJ,gg.nIJ);      %Topography
h = aa.h; h_diag = spdiags(h(:),0,gg.nIJ,gg.nIJ);
Cslip_diag = spdiags(C(:),0,gg.nIJ,gg.nIJ);

%Use gradient instead of gg.nddx/y 
%since periodic BC conditions do not apply
[Sx,Sy] = gradient(s, gg.dx, gg.dy);        
Sx = Sx(:); Sy = Sy(:); 


exx = gg.du_x*(gg.S_u'*u);                                %Strain Rates
eyy = gg.dv_y*(gg.S_v'*v);
exy = 0.5*(gg.du_y*(gg.S_u'*u) + gg.dv_x*(gg.S_v'*v));
edeff = sqrt(exx.^2 + eyy.^2 + exx.*eyy + exy.^2 + pp.n_rp.^2);

nEff =  edeff.^((1-n)/n);        %Effective Viscosity [dimensionless]
nEff_diag = spdiags(nEff(:),0,gg.nIJ,gg.nIJ);                                  


%% Field equations for velocities
A1 = gg.S_u*gg.dh_x*4*nEff_diag*h_diag*gg.du_x*gg.S_u';     %LHS SSA 
A2 = gg.S_u*gg.dhu_y*nEff_diag*h_diag*gg.du_y*gg.S_u';
A3 = pp.c3*gg.S_u*gg.c_hu*Cslip_diag*gg.c_uh*gg.S_u';
AA = A1 + A2 - A3;


B1 = gg.S_u*gg.dh_x*2*nEff_diag*h_diag*gg.dv_y*gg.S_v';
B2 = gg.S_u*gg.dhu_y*nEff_diag*h_diag*gg.dv_x*gg.S_v';
BB =  B1+B2;

C1 = gg.S_v*gg.dh_y*2*nEff_diag*h_diag*gg.du_x*gg.S_u';
C2 = gg.S_v*gg.dhv_x*nEff_diag*h_diag*gg.du_y*gg.S_u';
CC = C1 + C2;


D1 = gg.S_v*gg.dh_y*4*nEff_diag*h_diag*gg.dv_y*gg.S_v';
D2 = gg.S_v*gg.dhv_x*nEff_diag*h_diag*gg.dv_x*gg.S_v';
D3 = pp.c3*gg.S_v*gg.c_hv*Cslip_diag*gg.c_vh*gg.S_v';
DD = D1 + D2 - D3;

LHS = [AA BB; CC DD];

f1a = gg.S_u*pp.c4*gg.c_hu*h_diag*Sx;                              %RHS SSA
f1b = gg.S_v*pp.c4*gg.c_hv*h_diag*Sy;

RHS = [f1a;f1b];



end

