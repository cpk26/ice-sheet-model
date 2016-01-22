function [ vv ] = ism_sia(s,h,C,pp,gg,oo )
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
Sy = Sy(:); 
Sg = sqrt(Sx.^2 + Sy.^2); 
h = h(:); 
C = C(:) + pp.C_rp;                             %Regularziation to avoid C=0 -> InF

n = pp.n_Glen;                                  %Ice Flow Parameters

ubas = -pp.c1.* C.^-1 .* h.* Sx;                  %Velocities (h-grid)
udef = -pp.c2 .* h.^(n+1) .* Sg.^(n-1) .* Sx; 
u = ubas + udef;

vbas = -pp.c1.* C.^-1 .* h .* Sy;
vdef = -pp.c2 .* h.^(n+1) .* Sg.^(n-1) .* Sy; 
v = vbas(:) + vdef;        

u = gg.S_u * gg.c_hu*u(:);                      %Transfer onto u/v grids
v = gg.S_v * gg.c_hv*v(:);

vv.u = u;                                       %Non-dimensional velocity
vv.v = v;



end

