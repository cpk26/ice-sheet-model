function [ msft ] = ism_vel_misfit(u,v,aa,pp,gg, oo)
%% Inversion cost function
% Inputs:
%   vv      struct containing initial solution variables
%   aa      prescribed fields, including inputs and boundary conditions
%   oo      options
% Outputs:
%   misfit  Least squares misfit

msft = pp.c9*0.5*sum( gg.c_uh*(u-aa.u).^2 + gg.c_vh*(v-aa.v).^2 )*gg.dx*gg.dy;


end

