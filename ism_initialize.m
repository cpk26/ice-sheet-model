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

if ~isfield(oo,'nl'), oo.nl = 50; end                   %Number of layers, must even so vl+1 is odd.
nl = oo.nl;      

if nargin<5, oo = struct; end
if ~isfield(dd,'nfxd'), dd.nfxd = []; end
if ~isfield(dd,'errvx'), dd.errvx = ones(gg.nIJ,1)*pp.u; end 
if ~isfield(dd,'errvy'), dd.errvy = ones(gg.nIJ,1)*pp.u; end


%% Put fields and variables in structs

%% Topography
%Use gradient instead of gg.nddx/y since periodic BC conditions do not apply
aa.s = s;
aa.b = b;
aa.h = s-b;

%% Bed Gradient
[Bx,By] = gradient(aa.b, gg.dx, gg.dy);             
aa.prj = sqrt(1+ Bx.^2 + By.^2);

%% Surface Gradient
 
%For interior of ice Sheet
[Sx,Sy] = gradient(aa.s, gg.dx, gg.dy);                 

%Correction for ice margin
su = (gg.c_hu*gg.S_h*aa.s(:))./(gg.c_hu*gg.S_h*(aa.h(:) > 0));  %Thickness on u,v grids,
sv = (gg.c_hv*gg.S_h*aa.s(:))./(gg.c_hv*gg.S_h*(aa.h(:) > 0));  %linear extrapolation at the edges

ddx = gg.S_h'*gg.du_x*su;   %Surface gradient on h-grid for margin
ddy = gg.S_h'*gg.dv_y*sv;

Sx = Sx(:);                                   %Vectorize, flip the sign in y-direction due to convention
Sx(logical(gg.nmgn)) = ddx(logical(gg.nmgn)); %Update gradient at ice margin
Sy(logical(gg.nmgn)) = -ddy(logical(gg.nmgn));
Sy = -Sy(:); 

aa.Sx = Sx;     %Save
aa.Sy = Sy;



%% Fixed Boundary Conditions. 
if any(dd.nfxd(:))
aa.nfxd_uval = dd.vx_u.*gg.nfxd_ugrid/pp.u;
aa.nfxd_vval = dd.vy_v.*gg.nfxd_vgrid/pp.u;
else
aa.nfxd_uval = zeros(gg.nJ,gg.nI+1); aa.nfxd_vval = zeros(gg.nJ+1,gg.nI);
end

%% Forward/Inverse Problem Settings
vv2.U = [u(:); v(:)];
vv2.u = u(:);
vv2.v = v(:);


if strcmp(oo.pT, 'forward')                     %Forward Problem
if oo.hybrid, Cb = C(:); aa.Cb = Cb;   
else aa.C = C(:); end   


elseif strcmp(oo.pT, 'inverse')                 %Inverse Problem
if oo.hybrid, Cb = C(:); vv2.Cb = Cb; 
else vv2.C = C(:); end;  

aa.u = gg.S_h*dd.vx(:)/pp.u;                    %h-grid
aa.v = gg.S_h*dd.vy(:)/pp.u; 
                
aa.erru = gg.S_h*dd.errvx(:)/pp.u;              %h-grid 
aa.errv = gg.S_h*dd.errvy(:)/pp.u; 

end

%% Initialize Hybrid or SSA variables
U = vv2.U;

if oo.hybrid,                                       %Hybrid                                                      
nEff = ism_visc(U,vv2,aa,pp,gg,oo);         
nEff_lyrs = repmat(nEff,1,nl+1);
F2 = ism_falpha(2,U,nEff_lyrs,vv2,aa,pp,gg,oo );    %Effective Basal Slipperiness
C = Cb(:)./(1 + (pp.c13*Cb(:)).*(gg.S_h'*F2)); 

vv2.C = C;
vv2.F2 = F2;
vv2.nEff = nEff;
vv2.nEff_lyrs = nEff_lyrs;

else                                                %SSA
nEff = ism_visc(U,vv2,aa,pp,gg,oo);         
vv2.nEff = nEff;   
end


end