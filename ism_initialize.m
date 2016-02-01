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


% if ~isempty(dd.nfxd)
% 
% i = 1; aa.nfxd_uval = NaN(numel(gg.nfxd_uind),1);
% for p = gg.nfxd_uind'
% [r,c] = ind2sub(size(gg.xx_u),p);
% aa.nfxd_uval(i) = nanmean(dd.vx(r, max(c-1,1):min(c+1,gg.nI)));  
% i=i+1;
% end;
%  
% offset = (gg.nI+1)*gg.nJ;
% i = 1; aa.nfxd_vval = NaN(numel(gg.nfxd_vind),1);
% for p = gg.nfxd_vind'
% [r,c] = ind2sub(size(gg.xx_v),p);
% aa.nfxd_vval(i) = nanmean(dd.vy(min(r-1,1):max(r+1,gg.nJ), c));
% i=i+1;
% end;
% 
% end


if strcmp(oo.pT, 'forward')
aa.C = C(:);                   %Forward Problem
vv.u = u(:);
vv.v = v(:);
elseif strcmp(oo.pT, 'inverse')
vv.C = C(:);                   %Inverse Problem
aa.u = u(:);
aa.v = v(:); 
end


end