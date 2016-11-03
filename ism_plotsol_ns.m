function [ vv2 ] = ism_plotsol_ns(vv, aa, dd, pp, gg, oo )
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

if strcmp(oo.pT, 'forward');                            %Problem Type
if oo.hybrid, Cb = aa.Cb; C = vv.C;             
else, C = aa.C; end; end;

if strcmp(oo.pT, 'inverse'); C = vv.C; 
if oo.hybrid, Cb = vv.Cb; end;
end;

blnk = ones(size(gg.m));                                       %Mask for background of images
blnkImage = double(cat(3,blnk, blnk, blnk));
imMask = ~(gg.m ==2);

u_h = gg.S_h'*gg.c_uh*vv.u;       %Solved Velocities
v_h = gg.S_h'*gg.c_vh*vv.v;       %u,v grids onto h-grid 
U = sqrt(u_h.^2 + v_h.^2);

if oo.hybrid                        %For hybrid model, convert u_Eff to u_surface; Need to d
F1 = ism_falpha(1,vv.U,vv.nEff_lyrs,vv,aa,pp,gg,oo ); F1 = gg.S_h'*F1;
F2 = ism_falpha(2,vv.U,vv.nEff_lyrs,vv,aa,pp,gg,oo ); F2 = gg.S_h'*F2;

tmpa = (1 + pp.c13*Cb(:).*F1)./(1 + pp.c13*Cb(:).*F2);

u_h = u_h.*tmpa;
v_h = v_h.*tmpa;
end

u_h = reshape(u_h, gg.nJ, gg.nI);      %Velocities
v_h = reshape(v_h, gg.nJ, gg.nI);      %Reshape for plotting
U = sqrt(u_h.^2 + v_h.^2);

u_h = u_h*pp.u*pp.ty;                                              %Dimensionalize [m/yr]
v_h = v_h*pp.u*pp.ty;
U = U*pp.u*pp.ty;

vv2.u = u_h;                                            
vv2.v = v_h;
vv2.U = U;

obs = struct();                                                    %Observed Velocities [m/yr] 
obs.u = dd.vx *pp.ty; obs.u(imMask) = 0;
obs.v = dd.vy *pp.ty; obs.v(imMask) = 0;
obs.U = sqrt(obs.u.^2 + obs.v.^2);


umin = min([u_h(:); obs.u(:)]); umax = max([u_h(:); obs.u(:)]);     %Solved Velocity extrema
vmin = min([v_h(:); obs.v(:)]); vmax = max([v_h(:); obs.v(:)]);
Umin = min([U(:); obs.U(:)]); Umax = max([U(:); obs.U(:)]);
Umax = 500;

if strcmp(oo.pT, 'inverse');                                         %Basal Slipperiness
vv2.C = reshape(vv.C, gg.nJ, gg.nI); 
if oo.hybrid, 
vv2.Cb = reshape(vv.Cb, gg.nJ, gg.nI); 
vv2.F1 = F1; vv2.F2 = F2;
end;  
end;


h = dd.h;                                                           %Topography
b = dd.b;
s = dd.s;


if oo.plot_vel
figure()                                                %Plot Solution Velocities
subplot(3,1,1)
imagesc(U);
title('Velocity [Solved]');
colorbar()
caxis([Umin Umax])
hold on
% hh = imagesc(blnkImage);
% set(hh, 'AlphaData', imMask)


subplot(3,1,2)
imagesc(u_h);
title('X component of velocity');
colorbar()
caxis([umin umax])
hold on
hh = imagesc(blnkImage);
set(hh, 'AlphaData', imMask)

subplot(3,1,3)
imagesc(v_h);
title('Y component of velocity');
colorbar()
caxis([vmin vmax])
hold on
hh = imagesc(blnkImage);
set(hh, 'AlphaData', imMask)

if oo.hybrid                                    %Plot Ratio of Basal to Surface Velocities

figure
velRatio = gg.S_h'*(gg.S_h*(1+pp.c13*Cb.*F1).^-1);
imagesc(reshape(velRatio,gg.nJ,gg.nI));
title('Ratio of basal to sliding velocity')
colorbar()  
axis equal
axis tight
    
    
end

end

if oo.plot_obs
figure()                                                %Plot Observed velocities
subplot(3,1,1)
imagesc(obs.U);
title('Velocity [observed]');
colorbar()
caxis([Umin Umax])
hold on
% hh = imagesc(blnkImage);
% set(hh, 'AlphaData', imMask)


subplot(3,1,2)
imagesc(obs.u);
title('X component of velocity');
colorbar()
caxis([umin umax])
hold on
% hh = imagesc(blnkImage);
% set(hh, 'AlphaData', imMask)

subplot(3,1,3)
imagesc(obs.v);
title('Y component of velocity');
colorbar()
caxis([vmin vmax])
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
if oo.hybrid, imagesc(reshape(Cb,gg.nJ,gg.nI))
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
imagesc(((U-obs.U)))
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

title('Velocity Error');
c = colorbar();
caxis([-50 50])
ylabel(c,'Error (m/yr)') 
hold on
%hh = imagesc(blnkImage);
%set(hh, 'AlphaData', imMask)   


figure()                                    %Relative Error Figure
imagesc(((100*(U-obs.U)./obs.U)))

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

caxis([-25 25])
c = colorbar();
%colormap((brewermap(10,'Blues')))
colormap(jet(5))

ylabel(c,'Error (%)')


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

figure()
histogram(abs(gg.S_h*(U(:)-obs.U(:)))./velErr,[0:14],'Normalization', 'Probability')
xlabel('Standard Deviation of Measurement Error')
ylabel('Probability')

figure()
histogram(abs(gg.S_h*(U(:)-obs.U(:))), [0:5:115],'Normalization', 'Probability')
xlabel('Difference between predicted vs measured velocities (ma^-1)')
ylabel('Probability')


end


% hold on
% hh = imagesc(blnkImage);
% set(hh, 'AlphaData', imMask)   
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0.25 7.0 5.0])

end

