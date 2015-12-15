function [aa,vv] = ism_initialize(b,s, u, v, C, gg,pp,oo)
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
if ~isfield(pp,'l_c'), pp.l_c = 1; end
if ~isfield(pp,'phi_s'), pp.phi_s = -inf; end
if ~isfield(oo,'prob'), oo.prob = 1; end

%% put fields and variables in structs

%Topography
aa.s = s;
aa.b = b;
aa.h = max(s-b,0);

if oo.pT == 1
aa.C = C(:);                   %Forward Problem
vv.u = u(:);
vv.v = v(:);
else
vv.C = C(:);                   %Inverse Problem
aa.u = u(:);
aa.v = v(:); 
end


end