function [LHS, RHS] = ism_deism_fieldeq(C,nEff,aa,pp,gg,oo)
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


% nha = full(gg.nha);                               %number of active h/u/v grid nodes
% nua = full(gg.nua);
% nva = full(gg.nva);


nha = gg.nha;
nua = gg.nua;
nva = gg.nva;
NU = nua + nva;

%% Surface Gradient
Sx = aa.Sx;
Sy = aa.Sy;

Cslip_u = (gg.c_hu*C)./(gg.c_hu*gg.S_h*(gg.m(:) ==2));            %Slipperiness on u,v grids
Cslip_v = (gg.c_hv*C)./(gg.c_hv*gg.S_h*(gg.m(:) ==2));

Cslip_u(logical(gg.S_u*gg.nmgn_ugrid(:))) = 0;                           %Slipperiness is zero at margin for BC
Cslip_v(logical(gg.S_v*gg.nmgn_vgrid(:))) = 0;

h = gg.S_h*aa.h(:);




%% Field equations for velocities
Dvec = ([3*nEff(:).*h(:); nEff(:).*h(:); nEff(:).*h(:); -pp.c3*Cslip_u(:); -pp.c3*Cslip_v(:)]);
L = numel(Dvec);
D = sparse([1:L],[1:L],Dvec,L,L);

LHS = X2*D*X;                              %LHS

A1 = gg.c_hu*h_diag*gg.S_h*Sx;             %Driving Stress
A2 = (gg.c_hu*gg.S_h*(aa.h(:) > 0));       %Interpolate within mask, extrap at edges             
f1a = pp.c4*(A1./A2);   

A1 = gg.c_hv*h_diag*gg.S_h*Sy;           
A2 = (gg.c_hv*gg.S_h*(aa.h(:) > 0));
f1b = pp.c4*(A1./A2);   

RHS = [f1a;f1b];



end

