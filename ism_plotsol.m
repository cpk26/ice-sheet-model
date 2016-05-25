function [ vv2 ] = ism_plotsol(vv, aa, dd, pp, pd, gg, oo )
%% Plot Solution
% Inputs:
%   vv      struct containing solution variables
%   aa      prescribed fields, including inputs and boundary conditions
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs:
%   vv2     struct containing dimensionalized solution variables


if strcmp(oo.pT, 'forward'); C = aa.C; end; %Problem Type
if strcmp(oo.pT, 'inverse'); C = vv.C; end;



u = vv.u;
v = vv.v;
U = vv.U;

u_h = gg.c_uh*u;                    %Velocities
v_h = gg.c_vh*v;                    %u,v grids onto h-grid 

if oo.hybrid                        %For hybrid model, convert u_Eff to u_surface
F1 = ism_falpha(1,vv,aa,pp,gg,oo );
F2 = ism_falpha(2,vv,aa,pp,gg,oo );
tmpa = (1 + C(:).*F1)./(1 + C(:).*F2);

u_h = u_h.*tmpa;
v_h = v_h.*tmpa;
end


u_h = reshape(u_h, gg.nJ, gg.nI);      %Velocities
v_h = reshape(v_h, gg.nJ, gg.nI);      %Reshape for plotting
U_h = sqrt(u_h.^2 + v_h.^2);

u_h = u_h*pp.u*pd.ty;                                  %Dimensionalize [m/yr]
v_h = v_h*pp.u*pd.ty;
U_h = U_h*pp.u*pd.ty;


h = dd.h;                                          %Topography
b = dd.b;
s = dd.s;

x = dd.x;                                          %Grid
y = dd.y;

vv2.u = u_h;                                            %Dimensionalized solutions
vv2.v = v_h;
vv2.U = U;

figure()                                                %Plot Solution Velocities
subplot(3,1,1)
imagesc(U_h);
title('Velocity');
colorbar()


subplot(3,1,2)
imagesc(u_h);
title('X component of velocity');
colorbar()

subplot(3,1,3)
imagesc(v_h);
title('Y component of velocity');
colorbar()

figure()                                                %Plot Topography
subplot(3,1,1)
imagesc(s)
title('Surface Elevation');
colorbar()

subplot(3,1,2)
imagesc(b)
title('Bed Elevation');
colorbar()

subplot(3,1,3)
imagesc(h)
title('Ice Thickness');
colorbar()



end

