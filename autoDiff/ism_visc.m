function [nEff_z] = ism_visc(z,U,nEff,C,aa,pp,gg,oo)
%% Calculate Ice Viscosity 
% Inputs:
%   z       Depth at each grid cell at which to evaluate viscosity
%   U       [u;v] ice velocities in the x,y directions
%   visEff  Current effective viscosity
%   C       Basal Slip
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs: Ax = b
%   LHS     Left hand side. 
%   RHS     Right hand side.

u = U(1:gg.nua); u_h = gg.c_uh*u;
v = U(gg.nua+1:end); v_h = gg.c_vh*v;
s = gg.S_h*aa.s(:);
h = gg.S_h*aa.h(:);
n = pp.n_Glen;


exx = gg.du_x*u;                                %Strain Rates
eyy = gg.dv_y*v;
exy = 0.5*(gg.dhu_y*u + gg.dhv_x*v);
exz = pp.c11*C.*u_h.*(s-z)./(nEff.*h);
eyz = pp.c11*C.*v_h.*(s-z)./(nEff.*h);
edeff = sqrt(exx.^2 + eyy.^2 + exx.*eyy + exy.^2 + (1/4)*exz.^2 + (1/4)*eyz.^2 + pp.n_rp.^2);

nEff_z =  edeff.^((1-n)/n);                       %Effective Viscosity [dimensionless]
                             
end


