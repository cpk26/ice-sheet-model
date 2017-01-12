function [ ] = ism_plotismip_inv( vv, aa, dd, pp, pd, gg, oo )
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

effSurfFac = (1 + pp.c13*Cb(:).*F1)./(1 + pp.c13*Cb(:).*F2);

u_h = u_h.*effSurfFac;
v_h = v_h.*effSurfFac;
end


u_h = reshape(u_h, gg.nJ, gg.nI);      %Velocities
v_h = reshape(v_h, gg.nJ, gg.nI);      %Reshape for plotting
U_h = sqrt(u_h.^2 + v_h.^2);

u_h = u_h*pp.ty*pp.u;                                  %Dimensionalize [m/yr]
v_h = v_h*pp.ty*pp.u;  
U_h = U_h*pp.ty*pp.u; 

figure()                                         %Plot Solution vs known
subplot(1,4,1)
rX = dd.x(:)/dd.L; 
rY = dd.y(:)/dd.L;
rUx = reshape(u_h,gg.nJ,gg.nI);
X = linspace(0,1,gg.nI);
Y = 0.25;
Ux = griddata(rX,rY,rUx,rX,Y);
plot(X,Ux,'b+','LineWidth',1); hold on;
rUx_ismip = reshape(dd.vx/pp.u,gg.nJ,gg.nI);
Ux_ismip = griddata(rX,rY,rUx_ismip,rX,Y);


plot(X,Ux_ismip,'k.-','LineWidth',1);
ylabel('Surface Velocity (ma^{-1})')
xlabel('Normalized x')  

subplot(1,4,2)                                      %error 
plot(Ux_ismip-Ux)
ylabel('Surface Velocity Difference (ma^{-1})')  
xlabel('Normalized x') 


subplot(1,4,3)                                      %Convergence
loglog(vv.output.trace.fval)
ylabel('Cost Function')
xlabel('Iteration')
grid on

subplot(1,4,4)                                      %Convergence
rC = reshape(vv.alpha,gg.nJ,gg.nI);
C = griddata(rX,rY,rC,rX,Y);
plot(X,C,'b+','LineWidth',1); hold on;

rC_ismip = reshape(dd.alpha,gg.nJ,gg.nI);
C_ismip = griddata(rX,rY,rC_ismip,rX,Y);


plot(X,C_ismip,'k.-','LineWidth',1);


end

