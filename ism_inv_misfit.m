function [ msft ] = ism_inv_misfit(vv,aa,pp,gg, oo)
%% Inversion cost function
% Inputs:
%   vv      struct containing initial solution variables
%   aa      prescribed fields, including inputs and boundary conditions
%   oo      options
% Outputs:
%   misfit  Least squares misfit

if ~isfield(oo,'inv_msft'), oo.inv_msft = 'abs'; end  

u = vv.u;
v = vv.v;

%Velocity Misfit
if strcmp(oo.inv_msft,'abs') 
msft = pp.c9*0.5*sum( gg.c_uh*(u-aa.u).^2 + gg.c_vh*(v-aa.v).^2 )*gg.dx*gg.dy;
elseif strcmp(oo.inv_msft,'rel')    
msft = pp.c10*0.5*sum( gg.c_uh*((u-aa.u)./aa.u).^2 + gg.c_vh*((v-aa.v)./aa.v).^2 )*gg.dx*gg.dy;  
end

msft = pp.L_vel*msft;

%Add Cost function enforcing smoothness for 'grid' discretization
if strcmp(oo.Cdisc,'grid') 
C1 = (1/pp.x)*gg.dh_x*gg.S_h*vv.acoeff(:).*((gg.c_hu*gg.S_h*ones(gg.nIJ,1)) == 1); %gradient of alpha coefficients, x,y directions
C2 = (1/pp.x)*gg.dh_y*gg.S_h*vv.acoeff(:).*((gg.c_hv*gg.S_h*ones(gg.nIJ,1)) == 1);
cst = 0.5*pp.c10*(sum(C1(:).^2) + sum(C2(:).^2))*gg.dx*gg.dy;
msft = msft + pp.L_smooth*cst;
end


end

