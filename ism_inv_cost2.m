function [ cst,tik ] = ism_inv_cost2(uv,C,alpha,F1,F2,vv,aa,pp,gg, oo)
%% Inversion cost function
% Inputs:
%   vv      struct containing initial solution variables
%   aa      prescribed fields, including inputs and boundary conditions
%   oo      options
% Outputs:
%   cst     Inversion cost

%if ~isfield(pp,'L_smooth'), pp.L_smooth = 0; end 


u = uv(1:gg.nua); 
v = uv(gg.nua+1:end);

Cb = C;                                           %Basal Slip

esf = 1;
if oo.hybrid                                     %Hybrid Model; Compute surface velocities
esf = (1 + (pp.c13*Cb(:)).*F1)./(1 + (pp.c13*Cb(:)).*F2);          %U_eff to U_surface factor on h-grid
end

%Velocity Misfit
if isequal(oo.inv_cst,1)                            %Absolute
    cu = ((gg.c_uh*u.*esf)-aa.u).^2;
    cv = ((gg.c_vh*v.*esf)-aa.v).^2;     
    cst = pp.c9*0.5*sum( [cu; cv] )*gg.dx*gg.dy;
elseif isequal(oo.inv_cst,2)                        %Relative
    cu = (((gg.c_uh*u.*esf)-aa.u)./aa.u).^2;
    cv = (((gg.c_vh*v.*esf)-aa.v)./aa.v).^2;
    cst = pp.c10*0.5*sum( [cu; cv] )*gg.dx*gg.dy; 
elseif isequal(oo.inv_cst,3)                        %Log 
    vn = (10/pp.ty)/pp.u;               %Velocity normalization 
    cm = sqrt( (gg.c_uh*u.*esf).^2 + (gg.c_vh*v.*esf).^2 ) + vn;
    co = sqrt( (aa.u).^2 + (aa.v).^2 ) + vn;
    cst = pp.c10*0.5*sum( (log(cm./co)).^2 );
elseif isequal(oo.inv_cst,4)                        %Weighted least squares 
    cu = ( ((gg.c_uh*u.*esf)-aa.u)./(0.03*abs(aa.u)+aa.erru) ).^2;
    cv = ( ((gg.c_vh*v.*esf)-aa.v)./(0.03*abs(aa.v)+aa.errv) ).^2;  
    cst = pp.c10*0.5*sum( [cu; cv] )*gg.dx*gg.dy;
end

cst = pp.L_vel*cst;

%% Tikanov Regularization
tik = 0;

if pp.L_smooth > 0
h = gg.S_h*ones(gg.nIJ,1);
umask = 1 - abs(gg.dh_x*h)*gg.dx;                              %Remove mgn nodes
vmask = 1 - abs(gg.dh_y*h)*gg.dy;

C1 = (umask.*(gg.dh_x*(pp.bd*alpha)))./pp.x;                   %gradient of alpha coefficients, x-dir/ugrid
C2 = (vmask.*(gg.dh_y*(pp.bd*alpha)))./pp.x;                   %gradient of alpha coefficients, y-dir/vgrid

tik = 0.5*pp.c10*(sum(C1(:).^2) + sum(C2(:).^2))*gg.dx*gg.dy;
cst = cst + pp.L_smooth*tik;
end



end

