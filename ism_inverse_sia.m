function [ vv ] = ism_inverse_sia(s,h,u_obs,v_obs,vv,pp,gg,oo )
%% Inverse Shallow ice approximation model
% Inputs:
%   s       surface elevation
%   h       ice sheet depth
%   u_obs   observed ice sheet velocity, x-dir, on h-grid, dimensionless
%   v_obs   observed ice sheet velocity, y-dir, on h-grid, dimensionless
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs:
%   vv     struct containing new solution variables

%Use gradient instead of gg operators as periodic BC do not apply to topography

[Sx,Sy] = gradient(s, gg.dx, gg.dy);            %Topography   
Sx = Sx(:); 
Sy = -Sy(:); 
h = h(:); 

U_obs = sqrt((gg.S_h * u_obs).^2 + (gg.S_h * v_obs).^2); %Magnitude of observed velocities (h-grid)
U_b = U_obs;

%% Determine Basal Slip

A1 = -pp.c1 * gg.S_h * (h .* Sx);                   %Basal Slip multiplier
A2 = -pp.c1 * gg.S_h * (h .* Sy);
AA = sqrt(A1.^2 + A2.^2);
 
C = gg.S_h' * (AA./U_b + pp.C_rp);    %Basal slipperiness

vv.C = C;
end

