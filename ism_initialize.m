function [aa,vv] = ism(b,s,gg,pp,oo)
% Assign default prescribed fields and boundary conditions and initial
% conditions [ probably need to change initial hs and phi ]
% Inputs
%   b bed elevation on nodes
%   s surface elevation on nodes
%   gg grid structure
%   pp parameters
%   oo options [optional]
% Outputs
%   aa structure containing prescribed fields and boundary conditions
%   vv structure containing initial variables
%

if nargin<5, oo = struct; end
if ~isfield(pp,'u_b'), pp.u_b = 1; end
if ~isfield(pp,'l_c'), pp.l_c = 1; end
if ~isfield(pp,'phi_s'), pp.phi_s = -inf; end


%% basal velocity [ value taken from pp ]
Ub = pp.u_b*ones(gg.nJ,gg.nI); 
Ub(gg.next) = NaN;
Ub(gg.nmgn) = 0;

%% Strain
exx = zeros(gg.nJ,gg.nI);
exy = zeros(gg.nJ,gg.nI);
eyx = zeros(gg.nJ,gg.nI);
eyy = zeros(gg.nJ,gg.nI);


%% boundary conditions


%% put fields and variables in structs
aa.s = s;
aa.b = b;
aa.H = max(s-b,0);

aa.Ub = Ub; 
aa.exx = exx; 
aa.exy = exy; 
aa.eyx = eyx; 
aa.eyy = eyy; 


vv.Ub = Ub; 
vv.exx = exx; 
vv.exy = exy; 
vv.eyx = eyx; 
vv.eyy = eyy; 


end