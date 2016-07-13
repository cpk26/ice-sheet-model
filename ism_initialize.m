function [aa,vv2] = ism_initialize(b,s, u, v, C,vv, dd,gg,pp,oo)
% Assign default prescribed fields and boundary conditions and initial
% conditions [ probably need to change initial hs and phi ]
% Inputs
%   b bed elevation on nodes
%   s surface elevation on nodes
%   u velocity in x-direction
%   v velocity in y-direction
%   C basal slipperiness
%   gg grid structure
%   pp parameters
%   oo options [optional]
% Outputs
%   aa structure containing prescribed fields and boundary conditions
%   vv2 structure containing initial variables
%

if nargin<5, oo = struct; end
if ~isfield(dd,'nfxd'), dd.nfxd = []; end
if ~isfield(dd,'errvx_u'), dd.errvx_u = ones(size(u))*pp.u; end 
if ~isfield(dd,'errvy_v'), dd.errvy_v = ones(size(v))*pp.u; end
if ~isfield(dd,'errv'), dd.errv = ones(size(C))*pp.u; end


%% put fields and variables in structs

%Topography
aa.s = s;
aa.b = b;
aa.h = s-b;

%Fixed Boundary Conditions. 

if any(dd.nfxd(:))
aa.nfxd_uval = dd.vx_u.*gg.nfxd_ugrid/pp.u;
aa.nfxd_vval = dd.vy_v.*gg.nfxd_vgrid/pp.u;
else
aa.nfxd_uval = zeros(gg.nJ,gg.nI+1); aa.nfxd_vval = zeros(gg.nJ+1,gg.nI);
end

if strcmp(oo.pT, 'forward')
aa.C = C(:);                   %Forward Problem
vv2.U = [u(:); v(:)];
vv2.u = u(:);
vv2.v = v(:);

elseif strcmp(oo.pT, 'inverse')
vv2.C = C(:);                   %Inverse Problem
vv2.u = vv.u(:);
vv2.v = vv.v(:);

aa.U = [u(:); v(:)];
aa.u = u(:); 
aa.v = v(:);

aa.errvx_u = gg.S_u*dd.errvx_u(:)/pp.u; 
aa.errvy_v = gg.S_v*dd.errvy_v(:)/pp.u; 
aa.errv = gg.S_h*dd.errv(:)/pp.u;

end

if oo.hybrid,                               %Determine initial viscosity for Hybrid Approximation                           
Cb = C(:);                                  %Basal Slipperiness
U = [vv2.u(:);vv2.v(:)];
nEff = ism_visc(U,vv2,aa,pp,gg,oo);          %Initial viscosity;

for j=[1:10]                                 %Self consistent viscosity
F2 = ism_falpha(2,nEff,vv2,aa,pp,gg,oo );    %Effective Basal Slipperiness
C = Cb(:)./(1 + Cb(:).*(gg.S_h'*F2));
nEff = ism_visc_di(U,nEff,gg.S_h*C(:),aa,pp,gg,oo); %Updated Viscosity

end   

vv2.nEff = nEff;
end



% aa.u(aa.u <= 0) = aa.u(aa.u <= 0) - pp.U_rp; %Regularize
% aa.u(aa.u > 0) = aa.u(aa.u > 0) + pp.U_rp;
% aa.v(aa.v <= 0) = aa.v(aa.v <= 0) - pp.U_rp; 
% aa.v(aa.v > 0) = aa.v(aa.v > 0) + pp.U_rp;



end