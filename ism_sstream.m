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

numIter = 10;                               %Solver parameters
sstream_norm = zeros(numIter,1);

if strcmp(oo.pT, 'forward'); C = aa.C; end; %Problem Type
if strcmp(oo.pT, 'inverse'); C = vv.C; end;

u = vv.u;                                   %Initial iterate velocity 
v = vv.v;

%% Picard Iterations
for j = 1:numIter


[LHS, RHS] = ism_sstream_fieldeq(u,v,C,aa,pp,gg,oo);      %Field Equations
U = Inf(size(LHS,2),1);                                   %Velocity vector, unmodified length   


%% Apply Boundary Conditions

DEL = zeros(gg.nua+gg.nva,1);                   %Columns to delete
DEL2 = DEL;                                     %Rows to delete

if any(gg.nmgn(:))  %Ice Margin Nodes
hu = (gg.c_hu*gg.S_h*aa.h(:))./(gg.c_hu*gg.S_h*(aa.h(:) > 0));            %Thickness on u,v grids, linear extrapolation at the edges
hv = (gg.c_hv*gg.S_h*aa.h(:))./(gg.c_hv*gg.S_h*(aa.h(:) > 0));

t_mgn = [0.5*pp.c4*(0.5*hu).^2; 0.5*pp.c4*(0.5*hv).^2];                   %Boundary Condition at Ice Margin
mgn_mask = [gg.S_u*gg.nmgn_ugrid(:); gg.S_v*gg.nmgn_vgrid(:)];            %Ice Margin Node Mask

RHS(logical(mgn_mask)) = 0;                                               %Insert Ice Margin BC
RHS = RHS + mgn_mask.*t_mgn;
end

if any(gg.nfxd(:))  %Dirichlet BC Nodes
tmp_a = [gg.S_u*aa.nfxd_uval(:); gg.S_v*aa.nfxd_vval(:)];                 %Vector of fixed values
tmp_b = [gg.S_u*gg.nfxd_ugrid(:); gg.S_v*gg.nfxd_vgrid(:)];               %Location of fixed values

RHS = RHS - LHS*tmp_a;
DEL = DEL + tmp_b;
DEL2 = DEL2 + DEL;

clear tmp_a tmp_b;
end

if any(gg.nperbc(:))    %Periodic BC nodes
tmp_a = gg.S_u*(gg.nperbc_ugrid(:)>0); tmp_a = logical(tmp_a);
tmp_b = gg.S_u*(gg.nperbc_ugrid(:)<0); tmp_b = logical(tmp_b);
tmp_c = gg.S_v*(gg.nperbc_vgrid(:)>0); tmp_c = logical(tmp_c);
tmp_d = gg.S_v*(gg.nperbc_vgrid(:)<0); tmp_d = logical(tmp_d);

LHS(:, tmp_a) = LHS(:, tmp_a) + LHS(:, tmp_b); % Apply Periodic BC
LHS(:, tmp_c) = LHS(:, tmp_c) + LHS(:, tmp_d); 
DEL = DEL + [tmp_b; tmp_d];

clear tmp_i tmp_j tmp_k tmp_l

end

DEL = logical(DEL);
DEL2 = logical(DEL2);

LHS(:,DEL) = [];
LHS(DEL2,:) = [];
RHS(DEL2,:) = [];

%Solve 
Um = LHS\RHS;               %Solve modified field equations

%Return to original velocity vector
U(DEL) = NaN;
U(~isnan(U)) = Um;

if any(gg.nfxd(:))
tmp_a = [gg.S_u*gg.nfxd_ugrid(:); gg.S_v*gg.nfxd_vgrid(:)]; tmp_a = logical(tmp_a);
tmp_b = [gg.S_u*aa.nfxd_uval(:); gg.S_v*aa.nfxd_vval(:)];              

U(tmp_a) = tmp_b(tmp_a);

clear tmp_a tmp_b;
end

if any(gg.nperbc)
tmp_a = [gg.S_u*(gg.nperBC_ugrid(:) < 0); gg.S_v*(gg.nperBC_vgrid(:) < 0)]; tmp_a = logical(tmp_a);
tmp_b = [gg.S_u*(gg.nperBC_ugrid(:) > 0); gg.S_v*(gg.nperBC_vgrid(:) > 0)]; tmp_b = logical(tmp_b); 

U(tmp_a) = U(tmp_a) + tmp_b(tmp_a);

clear tmp_a tmp_b;
end

u = U(1:gg.nua);    %u,v velocity fields
v = U(gg.nua+1:end);

sstream_norm(j) = norm(RHS-LHS*Um,oo.norm); %iteration norm (using Um)

end

vv.u = u;
vv.v = v;
vv.sstream_norm = sstream_norm;

vv2=vv;

end



