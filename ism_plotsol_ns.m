function ism_plotsol_ns(vv, aa, dd, pp, gg, oo )
%% Plot Solution
% Inputs:
%   vv      struct containing solution variables
%   aa      prescribed fields, including inputs and boundary conditions
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs:
%   vv2     struct containing dimensionalized solution variables


if ~isfield(oo,'plot_vel'), oo.plot_vel = 1; end        % Plotting Options                 
if ~isfield(oo,'plot_obs'), oo.plot_obs = 1; end
if ~isfield(oo,'plot_err'), oo.plot_err = 1; end
if ~isfield(oo,'plot_C'), oo.plot_C = 1; end  
if ~isfield(oo,'plot_topo'), oo.plot_topo = 1; end                        


blnk = ones(size(gg.m));                                       %Mask for background of images
blnkImage = double(cat(3,blnk, blnk, blnk));
imMask = ~(gg.m ==2);


Cb = vv.Cb;


u_h = gg.c_uh*vv.u;       %Solved Velocities
v_h = gg.c_vh*vv.v;       %u,v grids onto h-grid 

if oo.hybrid                        %For hybrid model, convert u_Eff to u_surface; Need to d
effSurfFac = (1 + pp.c13*Cb(:).*vv.F1)./(1 + pp.c13*Cb(:).*vv.F2);

u_h = u_h.*effSurfFac;
v_h = v_h.*effSurfFac;
end

u_h = reshape(gg.S_h'*u_h, gg.nJ, gg.nI);      %Velocities
v_h = reshape(gg.S_h'*v_h, gg.nJ, gg.nI);      %Reshape for plotting
U = sqrt(u_h.^2 + v_h.^2);

u_h = u_h*pp.u*pp.ty;                                              %Dimensionalize [m/yr]
v_h = v_h*pp.u*pp.ty;
U = U*pp.u*pp.ty;


obs = struct();                                                    %Observed Velocities [m/yr] 
obs.u = dd.vx *pp.ty; obs.u(imMask) = 0;
obs.v = dd.vy *pp.ty; obs.v(imMask) = 0;
obs.U = sqrt(obs.u.^2 + obs.v.^2);


umin = min([u_h(:); obs.u(:)]); umax = max([u_h(:); obs.u(:)]);     %Solved Velocity extrema
vmin = min([v_h(:); obs.v(:)]); vmax = max([v_h(:); obs.v(:)]);
Umin = min([U(:); obs.U(:)]); Umax = max([U(:); obs.U(:)]);
Umax = 500;


h = dd.h;                                                           %Topography
b = dd.b;
s = dd.s;


if oo.plot_vel
figure()                                                %Plot Solution Velocities
imagesc(U);
title('Velocity [Solved]');
colorbar()
caxis([Umin Umax])
hold on
% hh = imagesc(blnkImage);
% set(hh, 'AlphaData', imMask)


if oo.hybrid                                    %Plot Ratio of Basal to Surface Velocities

figure
velRatio = (gg.S_h'*(1+pp.c13*Cb.*vv.F1).^-1);
imagesc(reshape(velRatio,gg.nJ,gg.nI));
title('Ratio of basal to sliding velocity')
colorbar()  
caxis([.50 1])
axis equal
axis tight
    
    
end

end

if oo.plot_obs
figure()                                                %Plot Observed velocities
imagesc(obs.U);
title('Velocity [observed]');
colorbar()
caxis([Umin Umax])
hold on
% hh = imagesc(blnkImage);
% set(hh, 'AlphaData', imMask)

end

if oo.plot_topo
figure()                                                %Plot Topography
imagesc(s)
axis equal
axis tight
title('Surface Elevation');
colorbar()
caxis([0 2500])

hold on
% hh = imagesc(blnkImage);
% set(hh, 'AlphaData', imMask)

figure;
imagesc(b)
axis equal
axis tight
title('Bed Elevation');
colorbar()
caxis([-500 500])
hold on

axis equal
axis tight
set(gca, 'FontSize', 9)
set(gca, 'FontName', 'Arial')

XTicks = [1, 176];
XLabels = {'-214075', '-99505'};
YTicks = [1,140];
YLabels = {'-2196035';'-2271905'};

%set(gca,'YDir','normal')
axis equal
axis([XTicks YTicks])

set(gca,'XTick',XTicks)
set(gca,'XTickLabel',XLabels)

set(gca,'YTick',YTicks)
set(gca,'YTickLabel',YLabels)

xlabel('Easting (m)')
ylabel('Northing (m)')

set(gcf,'color','w');
set(gca,'TickLength',[0 0]);

figure
imagesc(h)
title('Ice Thickness');
colorbar()
hold on
% hh = imagesc(blnkImage);
% set(hh, 'AlphaData', imMask)

end

if oo.plot_C                                    %Basal Drag Figure
figure()
if oo.hybrid, imagesc(reshape(gg.S_h'*Cb,gg.nJ,gg.nI))
else, imagesc(reshape(C,gg.nJ,gg.nI)); end
axis equal
axis tight
title('Basal Slipperiness');
c = colorbar();
hold on
hh = imagesc(blnkImage);
set(hh, 'AlphaData', imMask) 
end

if oo.plot_err

figure()                                        %Absolute Error Figure
imagesc(abs((U-obs.U)))
axis equal
axis tight
set(gca, 'FontSize', 9)
set(gca, 'FontName', 'Arial')
XTicks = [1, 176];
XLabels = {'-214075', '-99505'};
YTicks = [1,140];
YLabels = {'-2196035';'-2271905'};

axis equal
axis([XTicks YTicks])

set(gca,'XTick',XTicks)
set(gca,'XTickLabel',XLabels)

set(gca,'YTick',YTicks)
set(gca,'YTickLabel',YLabels)

xlabel('Easting (m)')
ylabel('Northing (m)')

set(gcf,'color','w');
set(gca,'TickLength',[0 0]);

title('|u_{obs} - u_{s}|');
c = colorbar();
caxis([0 50])
ylabel(c,'Error (m/yr)') 
hold on
%hh = imagesc(blnkImage);
%set(hh, 'AlphaData', imMask)   


figure()                                    %Relative Error Figure
imagesc(abs((100*(U-obs.U)./obs.U)))

axis equal
axis tight
set(gca, 'FontSize', 9)
set(gca, 'FontName', 'Arial')

XTicks = [1, 176];
XLabels = {'-214075', '-99505'};
YTicks = [1,140];
YLabels = {'-2196035';'-2271905'};

%set(gca,'YDir','normal')
axis equal
axis([XTicks YTicks])

set(gca,'XTick',XTicks)
set(gca,'XTickLabel',XLabels)

set(gca,'YTick',YTicks)
set(gca,'YTickLabel',YLabels)

xlabel('Easting (m)')
ylabel('Northing (m)')

set(gcf,'color','w');
set(gca,'TickLength',[0 0]);

caxis([0 25])
c = colorbar();
%colormap((brewermap(10,'Blues')))
colormap(jet(5))

ylabel(c,'Error (%)')
title('|u_{obs} - u_{s}|/|u_{obs}|');

figure()                                            %Error in Measurements
velErr = sqrt(...
    ((aa.erru + 0.03*abs(aa.u))).^2 +...
    ((aa.errv + 0.03*abs(aa.v))).^2 ...
    );
relErr = 100*velErr.*((gg.S_h*U(:)).^-1);
imagesc(reshape(gg.S_h'*relErr,gg.nJ,gg.nI))
colorbar()
caxis([0 25])
colormap(jet(5))
title('Observed measurement error')

figure()
histogram(abs(gg.S_h*(U(:)-obs.U(:)))./velErr,[0:14],'Normalization', 'Probability')
xlabel('Standard deviation of velocity observations')
ylabel('Normalized frequency')

figure()
histogram(abs(gg.S_h*(U(:)-obs.U(:))), [0:5:115],'Normalization', 'Probability')
xlabel('|u_{obs}-u_{s}| (ma^-1)')
ylabel('Normalized frequency')


end


% hold on
% hh = imagesc(blnkImage);
% set(hh, 'AlphaData', imMask)   
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0.25 7.0 5.0])

end

