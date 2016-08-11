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

if ~isfield(oo,'nl'), oo.nl = 50; end                   %%Number of layers, must even so vl+1 is odd.
nl = oo.nl;      

if nargin<5, oo = struct; end
if ~isfield(dd,'nfxd'), dd.nfxd = []; end
if ~isfield(dd,'errvx'), dd.errvx = ones(size(u))*pp.u; end 
if ~isfield(dd,'errvy'), dd.errvy = ones(size(v))*pp.u; end
if ~isfield(dd,'errv'), dd.errv = ones(size(C))*pp.u; end


%% put fields and variables in structs

%Topography
aa.s = s;
aa.b = b;
aa.h = s-b;
[Bx,By] = gradient(aa.b, gg.dx, gg.dy);             %Bed Gradient
aa.prj = sqrt(1+ Bx.^2 + By.^2);


%Fixed Boundary Conditions. 

if any(dd.nfxd(:))
aa.nfxd_uval = dd.vx_u.*gg.nfxd_ugrid/pp.u;
aa.nfxd_vval = dd.vy_v.*gg.nfxd_vgrid/pp.u;
else
aa.nfxd_uval = zeros(gg.nJ,gg.nI+1); aa.nfxd_vval = zeros(gg.nJ+1,gg.nI);
end

if strcmp(oo.pT, 'forward')                     %Forward Problem
if oo.hybrid, Cb = C(:); aa.Cb = Cb;   
else aa.C = C(:); end                   
vv2.U = [u(:); v(:)];
vv2.u = u(:);
vv2.v = v(:);

elseif strcmp(oo.pT, 'inverse')                 %Inverse Problem
if oo.hybrid, Cb = C(:); vv2.Cb = Cb; 
else vv2.C = C(:); end;         

vv2.U = [u(:); v(:)];
vv2.u = u(:);
vv2.v = v(:);


aa.u = gg.S_h*dd.vx(:)/pp.u; 
aa.v = gg.S_h*dd.vy(:)/pp.u; 


aa.err = gg.S_h*dd.errv(:)/pp.u;
aa.erru = gg.S_h*dd.errvx(:)/pp.u; 
aa.errv = gg.S_h*dd.errvy(:)/pp.u; 

end

if oo.hybrid,                               %Determine initial viscosity for Hybrid Approximation                           
U = [vv2.u(:);vv2.v(:)];
nEff = ism_visc(U,vv2,aa,pp,gg,oo);          %Initial viscosity (SSA);
nEff_lyrs = repmat(nEff,1,nl+1);

for j=[1:50]                                 %Self consistent viscosity
F2 = ism_falpha(2,U,nEff_lyrs,vv2,aa,pp,gg,oo );    %Effective Basal Slipperiness
C = Cb(:)./(1 + (pp.c13*Cb(:)).*(gg.S_h'*F2)); 
[nEff, nEff_lyrs] = ism_visc_di(U,nEff_lyrs,gg.S_h*C(:),aa,pp,gg,oo); %Updated Viscosity

end   


vv2.C = C;
vv2.nEff = nEff;
vv2.nEff_lyrs = nEff_lyrs;
end



% aa.u(aa.u <= 0) = aa.u(aa.u <= 0) - pp.U_rp; %Regularize
% aa.u(aa.u > 0) = aa.u(aa.u > 0) + pp.U_rp;
% aa.v(aa.v <= 0) = aa.v(aa.v <= 0) - pp.U_rp; 
% aa.v(aa.v > 0) = aa.v(aa.v > 0) + pp.U_rp;



end