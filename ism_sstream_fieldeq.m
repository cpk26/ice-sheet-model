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
nha = sum(gg.S_h(:));                               %number of active h-grid nodes
nua = sum(gg.S_u(:));
nva = sum(gg.S_v(:));

h_diag = spdiags(gg.S_h*aa.h(:),0,nha,nha);         %Diagonalize     
Cslip_diag = spdiags(gg.S_h*C(:),0,nha,nha);        

[Sx,Sy] = gradient(aa.s, gg.dx, gg.dy); %Use gradient instead of gg.nddx/y since periodic BC conditions do not apply      
Sx = Sx(:); Sy = -Sy(:); 


exx = gg.du_x*u;                                %Strain Rates
eyy = gg.dv_y*v;
exy = 0.5*(gg.dhu_y*u + gg.dhv_x*v);
edeff = sqrt(exx.^2 + eyy.^2 + exx.*eyy + exy.^2 + pp.n_rp.^2);

nEff =  edeff.^((1-n)/n);                       %Effective Viscosity [dimensionless]
nEff_diag = spdiags(nEff(:),0,nha,nha);                                  


%% Field equations for velocities
X = [gg.du_x gg.dv_y; gg.du_x -gg.dv_y; gg.dhu_y gg.dhv_x; gg.c_uh zeros(nha,nva); zeros(nha,nua) gg.c_vh];
X2 = [gg.dh_x gg.dh_x gg.duh_y gg.c_hu zeros(nha,nua)'; gg.dh_y -gg.dh_y gg.dvh_x zeros(nha,nva)' gg.c_hv];
D = blkdiag(3*nEff_diag*h_diag, nEff_diag*h_diag, nEff_diag*h_diag, -pp.c3*Cslip_diag, -pp.c3*Cslip_diag);

LHS = X2*D*X;

A1 = gg.c_hu*h_diag*gg.S_h*Sx;             %RHS SSA
A2 = (gg.c_hu*gg.S_h*(aa.h(:) > 0));       %Interpolate within mask, extrap at edges             
f1a = pp.c4*(A1./A2);   

A1 = gg.c_hv*h_diag*gg.S_h*Sy;           
A2 = (gg.c_hv*gg.S_h*(aa.h(:) > 0));
f1b = pp.c4*(A1./A2);   

RHS = [f1a;f1b];



end

