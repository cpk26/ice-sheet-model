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

%% Calculate Deformational Velocities
if oo.hybrid,                                                 %Deformational velocity
u_d = -gg.S_h * (pp.c2*(Sx.^2 + Sy.^2).*(h.^4).*Sx);
v_d = -gg.S_h * (pp.c2*(Sx.^2 + Sy.^2).*(h.^4).*Sy);

u_d(logical(gg.S_h*gg.nmgn(:))) = 0;                          %Except at ice margin
v_d(logical(gg.S_h*gg.nmgn(:))) = 0;

else                                                          %Assumed Plug Flow
u_d = 0;
v_d = 0;
end

U_d = sqrt((u_d).^2 + (v_d).^2);

%% Determine Basal Velocities
F = 0.80;
U_d(U_d > F*U_obs) = F*U_obs(U_d > F*U_obs);        %Limit deformational velocity to fraction of observed velocity. 
U_b = U_obs - U_d;                                  %This acts a smoothing on C. Note: allow a minimum slip to prevent C = Inf 


%% Determine Basal Slip

A1 = -pp.c1 * gg.S_h * (h .* Sx);                   %Basal Slip multiplier
A2 = -pp.c1 * gg.S_h * (h .* Sy);
AA = sqrt(A1.^2 + A2.^2);
 
C = gg.S_h' * (AA./U_b);    %Basal slipperiness
C = C + pp.C_rp;            %Add regularization parameter

vv.C = C;
end

