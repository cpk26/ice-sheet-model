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


u_h = gg.c_uh*gg.S_u'*vv.u; u_h = reshape(u_h, gg.nJ, gg.nI);      %Velocities
v_h = gg.c_vh*gg.S_v'*vv.v; v_h = reshape(v_h, gg.nJ, gg.nI);      %u,v grids onto h-grid 
U = sqrt(u_h.^2 + v_h.^2);

u_h = u_h*ps.u;                                         %Dimensionalize
v_h = v_h*ps.u;
U = U*ps.u;

h = dd.h;                                               %Topography
b = dd.b;
s = dd.s;

x = dd.x;                                               %Grid
y = dd.y;

vv2.u = u_h;                                            %Dimensionalized solutions
vv2.v = v_h;
vv2.U = U;

figure()                                                %Plot Velocities
subplot(3,1,1)
imagesc(x,y,U);
title('Velocity');
colorbar()

subplot(3,1,2)
imagesc(x,y,u_h);
title('X component of velocity');
colorbar()

subplot(3,1,3)
imagesc(x,y,v_h);
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

