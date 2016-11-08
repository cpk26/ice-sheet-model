function [ vv ] = ism_sia(s,h,C,vv,pp,gg,oo )
%% Shallow ice approximation model
% Inputs:
%   s       surface elevation
%   h       ice sheet depth
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs:
%   vv     struct containing new solution variables

%Use gradient instead of gg operators as periodic BC do not apply to the
%topography of the the ISMIP-C experiment

[Sx,Sy] = gradient(s, gg.dx, gg.dy);            %Topography   
Sx = Sx(:); 
Sy = -Sy(:); 
h = h(:); 
C = C(:); 
C(C < pp.C_rp) = pp.C_rp;             %Regularization to avoid C=0 -> InF

%% Basal Velocities
u_b = -pp.c1 * (gg.S_h * C.^-1) .* (gg.S_h * (h .* Sx));      
v_b = -pp.c1 * (gg.S_h * C.^-1) .* (gg.S_h * (h .* Sy));

%% Deformational Velocity
if oo.hybrid,                                                 %Deformational velocity
u_d = -gg.S_h * (pp.c2*(Sx.^2 + Sy.^2).*(h.^4).*Sx);
v_d = -gg.S_h * (pp.c2*(Sx.^2 + Sy.^2).*(h.^4).*Sy);

u_d(logical(gg.S_h*gg.nmgn(:))) = 0;                          %Except at ice margin
v_d(logical(gg.S_h*gg.nmgn(:))) = 0;

else                                                          %Plug Flow
u_d = 0;
v_d = 0;
  
end

%% Depth Integrated Velocity

u_di = u_b + pp.mdR*u_d;                                              
v_di = v_b + pp.mdR*v_d;

u = gg.c_hu*u_di(:)./(gg.c_hu*(gg.S_h*gg.m(:) == 2));  %Transfer onto u/v grids
v = gg.c_hv*v_di(:)./(gg.c_hv*(gg.S_h*gg.m(:) == 2));  

vv.u = u;                              %Non-dimensional velocity
vv.v = v;
vv.uv = [u;v];

U = sqrt( (gg.c_uh*u).^2 + (gg.c_vh*v).^2 );
vv.U = U;

% vv.u_d = 0.5*gg.c_hu*u_d(:)./(gg.c_hu*(gg.S_h*gg.m(:) == 2));  %Transfer onto u/v grids
% vv.v_d = 0.5*gg.c_hv*v_d(:)./(gg.c_hv*(gg.S_h*gg.m(:) == 2));  
% vv.U_d = [u_d;v_d];



end

