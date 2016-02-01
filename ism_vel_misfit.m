function [ msft ] = ism_inversion_misfit(u,v,aa,pp,gg, oo)
%% Inversion cost function
% Inputs:
%   vv      struct containing initial solution variables
%   aa      prescribed fields, including inputs and boundary conditions
%   oo      options
% Outputs:
%   misfit  Least squares misfit

msft = pp.c9*0.5*sum( gg.c_uh*gg.S_u'*(aa.u-u).^2 + gg.c_vh*gg.S_v'*(aa.v-v).^2 )*gg.dx*gg.dy;


end

