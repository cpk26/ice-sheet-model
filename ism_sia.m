function [ vv ] = ism_sia(s,h,pp,gg,oo )
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

n = pp.n_Glen;                                  %Ice Flow Parameters
           
u = pp.c3 * h.^(n+1) .* Sg.^(n-1) .* Sx;        %Velocities (h-grid)
v = pp.c3 * h.^(n+1) .* Sg.^(n-1) .* Sy;        

u = gg.S_u * gg.c_hu*u(:);                            %Transfer onto u/v grids
v = gg.S_v * gg.c_hv*v(:);

vv.U = [u;v];
vv.u = u;
vv.v = v;



end

