function [aa,vv] = ism_initialize(b,s, u, v, C, dd,gg,pp,oo)
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
%   vv structure containing initial variables
%

if nargin<5, oo = struct; end
if ~isfield(dd,'nfxd'), dd.nfxd = []; end


%% put fields and variables in structs

%Topography
aa.s = s;
aa.b = b;
aa.h = max(s-b,0);

%Fixed Boundary Conditions. 

if ~isempty(dd.nfxd)
aa.nfxd_uval = dd.vx_u(gg.nfxd_uind)/pp.u;
aa.nfxd_vval = dd.vy_v(gg.nfxd_vind)/pp.u;
else
aa.nfxd_uval = []; aa.nfxd_vval = [];
end

if strcmp(oo.pT, 'forward')
aa.C = C(:);                   %Forward Problem
vv.u = u(:);
vv.v = v(:);
elseif strcmp(oo.pT, 'inverse')
vv.C = C(:);                   %Inverse Problem
aa.u = u(:); 
aa.v = v(:); 

aa.u(aa.u <= 0) = aa.u(aa.u <= 0) - pp.U_rp; %Regularize
aa.u(aa.u > 0) = aa.u(aa.u > 0) + pp.U_rp;
aa.v(aa.v <= 0) = aa.v(aa.v <= 0) - pp.U_rp; 
aa.v(aa.v > 0) = aa.v(aa.v > 0) + pp.U_rp;
end


end