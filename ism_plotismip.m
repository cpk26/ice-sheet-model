function [ ] = ism_plotsol(vv, aa, dd, pp, pd, gg, oo )
%% Plot Solution
% Inputs:
%   vv      struct containing solution variables
%   aa      prescribed fields, including inputs and boundary conditions
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs:
%   vv2     struct containing dimensionalized solution variables


Cb = vv.Cb;

u = vv.u;
v = vv.v;
U = vv.U;

u_h = gg.S_h'*gg.c_uh*vv.u;       %Solved Velocities
v_h = gg.S_h'*gg.c_vh*vv.v;       %u,v grids onto h-grid 

if oo.hybrid                        %For hybrid model, convert u_Eff to u_surface; Need to d
F1 = ism_falpha(1,vv.uv,vv.nEff_lyrs,vv,aa,pp,gg,oo ); F1 = gg.S_h'*F1;
F2 = ism_falpha(2,vv.uv,vv.nEff_lyrs,vv,aa,pp,gg,oo ); F2 = gg.S_h'*F2;

tmpa = (1 + pp.c13*Cb(:).*F1)./(1 + pp.c13*Cb(:).*F2);

u_h = u_h.*tmpa;
v_h = v_h.*tmpa;
end


u_h = reshape(u_h, gg.nJ, gg.nI);      %Velocities
v_h = reshape(v_h, gg.nJ, gg.nI);      %Reshape for plotting
U_h = sqrt(u_h.^2 + v_h.^2);

u_h = u_h*pp.ty*pp.u;                                  %Dimensionalize [m/yr]
v_h = v_h*pp.ty*pp.u;  
U_h = U_h*pp.ty*pp.u; 

figure()                                         %Plot Solution Line
rX = dd.x(:)/dd.L; 
rY = dd.y(:)/dd.L;
rUx = reshape(u_h,gg.nJ,gg.nI);
X = linspace(0,1,gg.nI);
Y = 0.25;
Ux = griddata(rX,rY,rUx,X,Y);
plot(X,Ux,'k.','LineWidth',3);

figure()                                        %Plot Solution Field             
imagesc(U_h);
title('Velocity');
colorbar()





end

