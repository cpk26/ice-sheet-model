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
effSurfFac = (1 + pp.c13*Cb(:).*vv.F1)./(1 + pp.c13*Cb(:).*vv.F2);

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
subplot(2,2,1)
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
ylabel('U_s (m yr^{-1})') 
xlabel('x (normalized)') 
text(0.05,.9,'a','Units','normalized', 'FontWeight','bold', 'FontSize', 12)




subplot(2,2,2)                                      %error 
plot(X,Ux_ismip-Ux)
ylabel('U_s (m yr^{-1})')  
xlabel('x (normalized)') 
text(0.05,.9,'b','Units','normalized', 'FontWeight','bold', 'FontSize', 12)


subplot(2,2,3)                                      %Convergence
loglog(vv.output.trace.fval)
ylabel('J')
xlabel('Iteration')
grid on
text(0.05,.9,'c','Units','normalized', 'FontWeight','bold', 'FontSize', 12)

subplot(2,2,4)                                      %Basal Drag
rC = reshape(vv.alpha,gg.nJ,gg.nI);
C = griddata(rX,rY,rC,rX,Y);
plot(X,C,'b+','LineWidth',1); hold on;

rC_ismip = reshape(dd.alpha,gg.nJ,gg.nI);
C_ismip = griddata(rX,rY,rC_ismip,rX,Y);


plot(X,C_ismip,'k.-','LineWidth',1);
xlabel('x (normalized)') 
text(0.05,.9,'d','Units','normalized', 'FontWeight','bold', 'FontSize', 12)



set(gcf,'color','w');
set(gcf,'PaperUnits','inches','PaperPosition',[0 0.25 7 5])

end

