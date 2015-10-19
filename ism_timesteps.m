function [tt,vv,info] = ism_timesteps(t_span,vv,aa,pp,gg,oo)
% timestep nevis variables over time t_span
% Inputs:
%   t_span  start and end times and times for saving output
%   vv      struct containing initial solution variables
%   aa      prescribed fields, including inputs and boundary conditions
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs:
%   tt      struct containing time series of average quantities
%   vv      struct containing final solution variables
%   info    information about last computation [optional]
%

tInit = t_span(0);
tFinal = t_span(end);
dt = pp.ts;

t = tInit;
while t < tFinal

    [vv1,vv2,info] = ism_timestep(dt,vv,aa,pp,gg,oo);
    

    t = t + dt;
end
