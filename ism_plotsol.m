function [ vv2 ] = ism_plotsol(vv, dd, ps, gg, oo )
%% Plot Solution
% Inputs:
%   vv      struct containing solution variables
%   aa      prescribed fields, including inputs and boundary conditions
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs:
%   vv2     struct containing dimensionalized solution variables


u = vv.u * ps.u; u = reshape(u, gg.nJ, gg.nI);          %Velocities
v = vv.v * ps.u; v = reshape(v, gg.nJ, gg.nI);
U = sqrt(u.^2 + v.^2);

h = dd.h;                                               %Topography
b = dd.b;
s = dd.s;

x = dd.x;                                               %Grid
y = dd.y;

vv2.u = u;                                              %Dimensionalized solutions
vv2.v = v;

figure()                                                %Plot Velocities
subplot(3,1,1)
imagesc(x,y,U);
title('Velocity');
colorbar()

subplot(3,1,2)
imagesc(x,y,u);
title('X component of velocity');
colorbar()

subplot(3,1,3)
imagesc(x,y,v);
title('Y component of velocity');
colorbar()

figure()                                                %Plot Topography
subplot(3,1,1)
imagesc(x,y,s)
title('Surface Elevation');
colorbar()

subplot(3,1,2)
imagesc(x,y,b)
title('Bed Elevation');
colorbar()

subplot(3,1,3)
imagesc(x,y,h)
title('Ice Thickness');
colorbar()



end

