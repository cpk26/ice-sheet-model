function [ vv ] = ism_sia(s,h,pp,gg,oo )
%% Shallow ice approximation model
% Inputs:
%   vv      struct containing initial solution variables
%   aa      prescribed fields, including inputs and boundary conditions
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs:
%   vv     struct containing new solution variables

                                     
%Use gradient instead of gg.nddx/gg.nddy 
%since periodic BC conditions do not apply
[Sx,Sy] = gradient(s, gg.dx, gg.dy);               %Topography   
Sx = Sx(:); 
Sy = Sy(:); 
Sg = sqrt(Sx.^2 + Sy.^2); 
h = h(:); 

n = pp.n_Glen;                                     %Ice Flow Parameters
           
u = pp.c3 * h.^(n+1) .* Sg.^(n-1) .* Sx;           %Velocities
v = pp.c3 * h.^(n+1) .* Sg.^(n-1) .* Sy;

vv.u = u;
vv.v = v;



end

