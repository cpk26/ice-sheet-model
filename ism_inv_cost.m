function [ cst ] = ism_inv_cost(U,C,F1,F2,vv,aa,pp,gg, oo)
%% Inversion cost function
% Inputs:
%   vv      struct containing initial solution variables
%   aa      prescribed fields, including inputs and boundary conditions
%   oo      options
% Outputs:
%   cst     Inversion cost

%if ~isfield(pp,'L_smooth'), pp.L_smooth = 0; end 


u = U(1:gg.nua); 
v = U(gg.nua+1:end);

Cb = C;                                           %Basal Slip

if oo.hybrid                                      %Hybrid Model; Compute surface velocities

tmpa = (1 + (pp.c13*Cb(:)).*F1)./(1 + (pp.c13*Cb(:)).*F2);          %U_eff to U_surface factor on h-grid
tmpa_u = (gg.c_hu*tmpa)./(gg.c_hu*(gg.S_h*gg.m(:) ==2));    %Interpolate onto u/v grids
tmpa_v = (gg.c_hv*tmpa)./(gg.c_hv*(gg.S_h*gg.m(:) ==2));    %Extrapolate at edges

if any(gg.nmgn(:))                                  %Ice Margin Nodes
tmpa_u(gg.nmgn_ugrid) = 1;                          %No vertical variation in velocity at ice margin                                       
tmpa_v(gg.nmgn_vgrid) = 1;
end

u = u.*tmpa_u;                                      %Surface velocities                    
v = v.*tmpa_v;                    
end

%Velocity Misfit
if strcmp(oo.inv_cst,'abs') 
    cu = (gg.c_uh*u-aa.u).^2;
    cv = (gg.c_vh*v-aa.v).^2;   
    cst = pp.c9*0.5*sum( [cu; cv] )*gg.dx*gg.dy;
elseif strcmp(oo.inv_cst,'rel') 
    N = 1e-10;
    cu = ((gg.c_uh*u-aa.u)./aa.u).^2;
    cv = ((gg.c_vh*v-aa.v)./aa.v).^2;
    cst = N*pp.c10*0.5*sum( [cu; cv] )*gg.dx*gg.dy; 
elseif strcmp(oo.inv_cst,'log') 
    N = 1e-10;
    vn = (10/pp.ty)/pp.u;               %Velocity normalization 
    cm = sqrt( (gg.c_uh*u).^2 + (gg.c_vh*v).^2 ) + vn;
    co = sqrt( (aa.u).^2 + (aa.v).^2 ) + vn;
    cst = N*pp.c10*0.5*sum( (log(cm./co)).^2 );
elseif strcmp(oo.inv_cst,'wls')
    N = 1e-10;
    cu = ( (gg.c_uh*u-aa.u)./(0.03*abs(aa.u)+aa.erru) ).^2;
    cv = ( (gg.c_vh*v-aa.v)./(0.03*abs(aa.v)+aa.errv) ).^2;  
    cst = N*pp.c10*0.5*sum( [cu; cv] )*gg.dx*gg.dy;
end

cst = pp.L_vel*cst;

%% Tikanov Regularization
if pp.L_smooth > 0

in_u = ((gg.c_hu*gg.S_h*ones(gg.nIJ,1)) == 1);      %Interior nodes on ugrid
in_u(gg.S_u*gg.nperbc_ugrid(:) > 0) = 0;            %Handle BC
in_u(gg.S_u*gg.nfxd_ugrid(:) > 0) = 0;         
C1 = (1/pp.x)*in_u.*(gg.dh_x*Cb);                   %gradient of alpha coefficients, x-dir/ugrid


in_v = ((gg.c_hv*gg.S_h*ones(gg.nIJ,1)) == 1);      %Interior nodes on vgrid 
in_v(gg.S_v*gg.nperbc_vgrid(:) > 0) = 0;            %Handle BC
in_v(gg.S_v*gg.nfxd_vgrid(:) > 0) = 0;              %Handle BC         
C2 = (1/pp.x)*in_v.*(gg.dh_y*Cb);                    %gradient of alpha coefficients, y-dir/vgrid

tik = 0.5*pp.c10*(sum(C1(:).^2) + sum(C2(:).^2))*gg.dx*gg.dy;
cst = cst + pp.L_smooth*tik;
end


end

