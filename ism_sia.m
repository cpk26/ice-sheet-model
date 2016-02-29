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

%Because our main forward/inverse model uses the SSA approximation, we only
%calculate the surface velocity using the basal sliding component of the SIA
%approximation, ignoring the contribution to velocity from internal deformation.
%See the Greve + Blatter book, 'Dynamics of Ice Sheets and Glaciers'.


%Use gradient instead of gg operators as periodic BC do not apply to the
%topography of the the ISMIP-C experiment

[Sx,Sy] = gradient(s, gg.dx, gg.dy);            %Topography   
Sx = Sx(:); 
Sy = -Sy(:); 
h = h(:); 
C = C(:) + pp.C_rp;                             %Regularziation to avoid C=0 -> InF

n = pp.n_Glen;

u_s = -pp.c1 * (gg.S_h * C.^-1) .* (gg.S_h * (h .* Sx));      %Sliding Velocities (h-grid)
v_s = -pp.c1 * (gg.S_h * C.^-1) .* (gg.S_h * (h .* Sy));

%u_g = -pp.c2*(h.^(n+1)) * sqrt(Sx + xy)


u = gg.c_hu*u_s(:);                      %Transfer onto u/v grids
v = gg.c_hv*v_s(:);

vv.u = u;                              %Non-dimensional velocity
vv.v = v;



end

