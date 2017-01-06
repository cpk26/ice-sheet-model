function [aa,vv2] = ism_initialize(b,s, u, v, C,vv, dd,gg,pp,oo)
% Assign default prescribed fields and boundary conditions and initial
% conditions 
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
if ~isfield(dd,'N'), dd.N = 0.1*(s-b)*pp.g*pp.rho_i/pp.phi; end


%% Put fields and variables in structs

%% Topography
%Use gradient instead of gg.nddx/y since periodic BC conditions do not apply
aa.s = s;
aa.b = b;
aa.h = s-b;

%% Bed Gradient
[Bx,By] = gradient(aa.b, gg.dx, gg.dy);             

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

A1 = gg.c_hu*gg.S_h*(aa.h(:).*Sx);             %Driving Stress
A2 = (gg.c_hu*gg.S_h*(aa.h(:) > 0));       %Interpolate within mask, extrap at edges             
DSa = pp.c4*(A1./A2);   

A1 = gg.c_hv*gg.S_h*(aa.h(:).*Sy);           
A2 = (gg.c_hv*gg.S_h*(aa.h(:) > 0));
DSb = pp.c4*(A1./A2);   

aa.DRIVSTRESS = [DSa;DSb];





%% Fixed Boundary Conditions. 
if any(dd.nfxd(:))
aa.nfxd_uval = dd.vx_u.*gg.nfxd_ugrid/pp.u;
aa.nfxd_vval = dd.vy_v.*gg.nfxd_vgrid/pp.u;
else
aa.nfxd_uval = zeros(gg.nJ,gg.nI+1); aa.nfxd_vval = zeros(gg.nJ+1,gg.nI);
end

%% Forward/Inverse Problem Settings
vv2.uv = [u(:); v(:)];
vv2.u = u(:);
vv2.v = v(:);



Cb = gg.S_h*C(:); 
vv2.Cb = Cb;

if isequal(oo.slidinglaw, 'linear')  
alpha = gg.S_h*C(:);
else
alpha = dd.alpha;
aa.N =gg.S_h*dd.N(:)/pp.phi; end

if strcmp(oo.pT, 'forward'), aa.alpha = alpha; 
else
vv2.alpha = alpha; 

aa.u = gg.S_h*dd.vx(:)/pp.u;                    %h-grid
aa.v = gg.S_h*dd.vy(:)/pp.u; 
                
aa.erru = gg.S_h*dd.errvx(:)/pp.u;              %h-grid 
aa.errv = gg.S_h*dd.errvy(:)/pp.u; 

end


% if strcmp(oo.pT, 'forward')                     %Forward Problem
% if oo.hybrid, Cb = gg.S_h*C(:); vv2.Cb = Cb; vv2.C = NaN(size(Cb));  
% else aa.Cb = gg.S_h*C(:); vv2.C = vv2.Cb; end   
% 
% if isequal(oo.slidinglaw, 'linear')  
% aa.alpha = gg.S_h*C(:);
% else aa.alpha = dd.alpha; 
% aa.N =gg.S_h*dd.N(:)/pp.phi; end
% 
% 
% elseif strcmp(oo.pT, 'inverse')                 %Inverse Problem
% if oo.hybrid, Cb = gg.S_h*C(:); vv2.Cb = Cb; 
% else vv2.Cb =gg.S_h*C(:); vv2.C = Cb; end;  
% 
% if isequal(oo.slidinglaw, 'linear')   
% vv2.alpha = gg.S_h*C(:);
% else vv2.alpha = dd.alpha;
% aa.N =gg.S_h*dd.N(:)/pp.phi; end
% 
% aa.u = gg.S_h*dd.vx(:)/pp.u;                    %h-grid
% aa.v = gg.S_h*dd.vy(:)/pp.u; 
%                 
% aa.erru = gg.S_h*dd.errvx(:)/pp.u;              %h-grid 
% aa.errv = gg.S_h*dd.errvy(:)/pp.u; 
% 
% end


%% Initialize Hybrid or SSA variables
uv = vv2.uv;

if oo.hybrid,                                       %Hybrid                                                      
nEff = ism_visc(uv,vv2,aa,pp,gg,oo);         
nEff_lyrs = repmat(nEff,1,nl+1);
F2 = ism_falpha(2,uv,nEff_lyrs,vv2,aa,pp,gg,oo );    %Effective Basal Slipperiness
C = Cb./(1 + (pp.c13*Cb).*(F2)); 

vv2.C = C;
vv2.F2 = F2;
vv2.nEff = nEff;
vv2.nEff_lyrs = nEff_lyrs;

else                                                %SSA
vv2.C = vv2.Cb;
nEff = ism_visc(uv,vv2,aa,pp,gg,oo);         
vv2.nEff = nEff;   
vv2.F2 =[];
end


end