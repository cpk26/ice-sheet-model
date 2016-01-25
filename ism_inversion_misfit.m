function [ cost ] = ism_inversion_misfit(u,v,aa,pp,gg, oo)
%% Inversion cost function
% Inputs:
%   vv      struct containing initial solution variables
%   aa      prescribed fields, including inputs and boundary conditions
%   oo      options
% Outputs:
%   cost     value of the cost

cost = pp.c9*sum( (aa.u-u).^2 + (aa.v-v).^2)*gg.dx*gg.dy;


end

