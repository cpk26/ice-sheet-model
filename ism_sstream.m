function [vv2] = ism_sstream(vv,aa,pp,gg,oo )
%% Shallow Stream Model 
% Inputs:
%   vv      struct containing initial solution variables
%   aa      prescribed fields, including inputs and boundary conditions
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs:
%   vv2     struct containing new solution variables
%   J       Jacobian matrix

numIter = 10;                               %Solver parameters
sstream_norm = zeros(numIter,1);

if strcmp(oo.pT, 'forward'); C = aa.C; end; %Problem Type
if strcmp(oo.pT, 'inverse'); C = vv.C; end;

u = vv.u;                                   %Initial iterate velocity 
v = vv.v;

%% Remap indices [from whole region, to masked area]

A = sum(gg.S_u); A2 = cumsum(A);            %u-grid
nfxd_uind = A2(gg.nfxd_uind);               %Fixed nodes                   
nperbc_u1ind  = A2(gg.nperbc_u1ind);        %Periodic BC nodes
nperbc_u2ind  = A2(gg.nperbc_u2ind);
nmgn_uind = A2(gg.nmgn_uind);               %Margin Nodes

vOff = sum(A);                              %offset to v values [number of u values]

A = sum(gg.S_v); A2 = cumsum(A);            %v-grid
nfxd_vind = A2(gg.nfxd_vind) + vOff;        %Fixed nodes   
nperbc_v1ind = A2(gg.nperbc_v1ind) + vOff;  %Periodic BC nodes
nperbc_v2ind = A2(gg.nperbc_v2ind) + vOff;
nmgn_vind = A2(gg.nmgn_vind)+ vOff;         %Margin Nodes


%% Picard Iterations
for j = 1:numIter


[LHS, RHS] = ism_sstream_fieldeq(u,v,C,aa,pp,gg,oo);      %Field Equations
U = Inf(size(LHS,2),1);                                   %Velocity vector, unmodified length   


%% Boundary Conditions

DEL = [];                   %Columns to delete
DEL2 = [];                  %Rows to delete

if any(gg.nmgn(:))  %Ice Margin Nodes
hu = (gg.c_hu*gg.S_h*aa.h(:))./(gg.c_hu*gg.S_h*(aa.h(:) > 0));            %Thickness on u,v grids
hv = (gg.c_hv*gg.S_h*aa.h(:))./(gg.c_hv*gg.S_h*(aa.h(:) > 0));

B = zeros(numel(RHS),1); 
B(nmgn_uind) = 1/(pp.x*gg.dx); B(nmgn_vind) = 1/(pp.x*gg.dy);

f2a = 0.5*pp.c4*hu.^2;
f2b = 0.5*pp.c4*hv.^2;
F = [f2a; f2b];

RHS(B>0) = 0;           %Replace forcing at the edge 
RHS = RHS + B.*F;
end

if any(gg.nfxd(:))  %Dirichlet BC Nodes
RHS = RHS - LHS(:,nfxd_uind)*aa.nfxd_uval;
RHS = RHS - LHS(:,nfxd_vind)*aa.nfxd_vval;    
DEL = union(nfxd_uind, nfxd_vind);
DEL2 = DEL;
end

if any(gg.nperbc(:))    %Periodic BC nodes
LHS(:, nperbc_u1ind) = LHS(:, nperbc_u1ind) + LHS(:, nperbc_u2ind); % Apply Periodic BC
LHS(:, nperbc_v1ind) = LHS(:, nperbc_v1ind) + LHS(:, nperbc_v2ind); 
DEL = union(DEL, [nperbc_u2ind; nperbc_v2ind]);
end


LHS(:,DEL) = [];
LHS(DEL2,:) = [];
RHS(DEL2,:) = [];

%Solve 
Um = LHS\RHS;               %Solve modified field equations

%Return to original velocity vector
U(DEL) = NaN;
U(~isnan(U)) = Um;

if any(gg.nfxd(:))
    U(nfxd_uind) = aa.nfxd_uval;
    U(nfxd_vind) = aa.nfxd_vval;
end

if any(gg.nperbc)
    U(nperbc_u2ind) = U(nperbc_u1ind);
    U(nperbc_v2ind) = U(nperbc_v1ind);
end

u = U(1:vOff);    %u,v velocity fields
v = U(vOff+1:end);

sstream_norm(j) = norm(RHS-LHS*Um,oo.norm); %iteration norm (using Um)

end

vv.u = u;
vv.v = v;
vv.sstream_norm = sstream_norm;

vv2=vv;

end



