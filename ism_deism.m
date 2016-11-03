function [vv2, rr] = ism_deism(vv,aa,pp,gg,oo )
%% Shallow Stream Model 
% Inputs:
%   vv      struct containing initial solution variables
%   aa      prescribed fields, including inputs and boundary conditions
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs:
%   vv2     struct containing new solution variables

numIter = oo.pic_iter;                      %Solver parameters
sstream_norm = zeros(numIter,1);
eq_norm = zeros(numIter,1);


if strcmp(oo.pT, 'forward'); 
    if oo.hybrid, Cb = aa.Cb; C = vv.C; nEff = vv.nEff; nEff_lyrs = vv.nEff_lyrs;
    else C = aa.C; nEff = vv.nEff; end; end
if strcmp(oo.pT, 'inverse'); 
    if oo.hybrid, Cb = vv.Cb; C = vv.C; nEff = vv.nEff; nEff_lyrs = vv.nEff_lyrs;
    else C = vv.C; nEff = vv.nEff; end; end


U = vv.U;                                   %Initial iterate velocity

rr = struct();                              %Preallocate arrays if we are saving picard iterations
if isequal(oo.savePicIter,1)
    rr.An = cell(1,numIter);
    rr.Un = zeros(gg.nua+gg.nva,numIter+1);
    rr.nEffn = zeros(gg.nha,numIter);
    if oo.hybrid; rr.F2n = zeros(gg.nha,numIter+1);end
end

if isequal(oo.savePicIter,1),   %Save initial velocity/viscosity; 
rr.Un(:,1) = U; 
rr.nEffn(:,1) = nEff(:);
if oo.hybrid
rr.F2n(:,1) = vv.F2(:); end;
end                      

        
%% Picard Iterations
for j = 1:numIter



%% Field Equations
[LHS, RHS] = ism_deism_fieldeq(C,nEff, aa,pp,gg,oo);              %Field Equations
LHSf = LHS;
RHSf = RHS;
sstream_norm(j) = norm(RHS-LHS*U)./norm(RHS); %iteration norm (using last iterations U)

U = Inf(size(LHS,2),1);                                                 %Velocity vector, full length   

if isequal(oo.savePicIter,1),                               %Save
rr.An{j} = LHS;  
end                        

%% Apply Boundary Conditions

DEL = zeros(gg.nua+gg.nva,1);                   %Columns to delete
DEL2 = DEL;                                     %Rows to delete

if any(gg.nmgn(:))  %Ice Margin Nodes
hu = (gg.c_hu*gg.S_h*aa.h(:))./(gg.c_hu*gg.S_h*(aa.h(:) > 0));            %Thickness on u,v grids, linear extrapolation at the edges
hv = (gg.c_hv*gg.S_h*aa.h(:))./(gg.c_hv*gg.S_h*(aa.h(:) > 0));

du = (gg.c_hu*gg.S_h*abs(min(0,aa.b(:))))./(gg.c_hu*gg.S_h*(aa.h(:) > 0));   %Draft on u,v grids, linear extrapolation at the edges
dv = (gg.c_hv*gg.S_h*abs(min(0,aa.b(:))))./(gg.c_hv*gg.S_h*(aa.h(:) > 0));


t_mgn = [0.5*pp.c4*(hu.^2 - pp.c5*du.^2);...                            %Boundary Condition at Ice Margin
        0.5*pp.c4*(hv.^2 - pp.c5*dv.^2)];                            
  
mgn_mask = [gg.S_u*gg.nmgn_ugrid(:); gg.S_v*gg.nmgn_vgrid(:)];            %Ice Margin Node Mask

RHS(logical(mgn_mask)) = 0;                                               %Insert Ice Margin BC
RHS = RHS + mgn_mask.*t_mgn;
end

if any(gg.nfxd(:))  %Dirichlet BC Nodes
tmp_a = [gg.S_u*gg.nfxd_ugrid(:); gg.S_v*gg.nfxd_vgrid(:)];               %Location of fixed values
tmp_b = [gg.S_u*aa.nfxd_uval(:); gg.S_v*aa.nfxd_vval(:)];                 %Vector of fixed values

if oo.hybrid,                                     %Convert surface vel -> effective vel
F1 = ism_falpha(1,U,nEff_lyrs,vv,aa,pp,gg,oo );          %Calculate F alpha factors 
F2 = ism_falpha(2,U,nEff_lyrs,vv,aa,pp,gg,oo );
tmp_c = (1 + pp.c13*(gg.S_h*Cb(:)).*F1)./(1 + pp.c13*(gg.S_h*Cb(:)).*F2);          %Un -> Us factor
tmpc_u = (gg.c_hu*tmp_c)./(gg.c_hu*(gg.S_h*gg.m(:)==2));    %Interpolate onto u/v grids
tmpc_v = (gg.c_hv*tmp_c)./(gg.c_hv*(gg.S_h*gg.m(:)==2));    %Extrapolate at edges
effVelFac = [tmpc_u; tmpc_v];            

tmp_b = tmp_b./effVelFac;                                %effective velocities
end


RHS = RHS - LHS*tmp_b;
DEL = DEL + tmp_a;
DEL2 = DEL2 + DEL;

clear tmp_a tmp_b;
end

if any(gg.nperbc(:))    %Periodic BC nodes
    
tmp_a = [gg.S_u*(gg.nperbc_ugrid(:) < 0); gg.S_v*(gg.nperbc_vgrid(:) < 0)]; tmp_a = logical(tmp_a);
tmp_b = [gg.S_u*(gg.nperbc_ugrid(:) > 0); gg.S_v*(gg.nperbc_vgrid(:) > 0)]; tmp_b = logical(tmp_b); 

LHS(:, tmp_b) = LHS(:, tmp_b) + LHS(:, tmp_a);
DEL = DEL + [tmp_a];

%DEL2 = DEL2 + [tmp_a];
clear tmp_a tmp_b tmp_c tmp_d

end

DEL = logical(DEL);
DEL2 = logical(DEL2);

LHS(:,DEL) = [];
LHS(DEL2,:) = [];
RHS(DEL2,:) = [];

%% Solve 

Um = LHS\RHS;               %Solve modified field equations

%% Return to original velocity vector

U(DEL) = NaN;
U(~isnan(U)) = Um;

if any(gg.nfxd(:))
tmp_a = [gg.S_u*gg.nfxd_ugrid(:); gg.S_v*gg.nfxd_vgrid(:)]; tmp_a = logical(tmp_a);
tmp_b = [gg.S_u*aa.nfxd_uval(:); gg.S_v*aa.nfxd_vval(:)];              

if oo.hybrid,                                     %Convert surface vel -> effective vel
tmp_b = tmp_b./effVelFac;                         %Use previous factor for efficiency
end

U(tmp_a) = tmp_b(tmp_a);
clear tmp_a tmp_b;
end

if any(gg.nperbc(:))

tmp_a = [gg.S_u*(gg.nperbc_ugrid(:) < 0); gg.S_v*(gg.nperbc_vgrid(:) < 0)]; tmp_a = logical(tmp_a);
tmp_b = [gg.S_u*(gg.nperbc_ugrid(:) > 0); gg.S_v*(gg.nperbc_vgrid(:) > 0)]; tmp_b = logical(tmp_b); 

U(tmp_a) = U(tmp_b);

clear tmp_a tmp_b;
end


u = U(1:gg.nua);    %u,v velocity fields
v = U(gg.nua+1:end);

%% Viscosity, effective basal drag    
    
if oo.hybrid                                            %Determine viscosity, basal slipperiness appropriately

[nEff, nEff_lyrs] = ism_visc_di(U,nEff_lyrs,gg.S_h*C(:),aa,pp,gg,oo); %Updated Viscosity
F2 = ism_falpha(2,U,nEff_lyrs,vv,aa,pp,gg,oo );
[Cb] = ism_slidinglaw(vv,aa,pp,gg,oo);

C = Cb(:)./(1 + (pp.c13*Cb(:)).*(gg.S_h'*F2));                    %Effective Basal Slipperiness

else nEff = ism_visc(U,vv,aa,pp,gg,oo); end             %SSA Viscosity



if isequal(oo.savePicIter,1),
    rr.Un(:,j+1) = U; 
    rr.nEffn(:,j+1) = nEff(:);
    if oo.hybrid; rr.F2n(:,j+1) = F2; end
end                %Save Intermediate velocity array


%% Plot Velocities [For manual testing of picard iteration]
% u_h = gg.S_h'*gg.c_uh*u; %Solution Velocities
% v_h = gg.S_h'*gg.c_vh*v;
% 
% if oo.hybrid                        %Hybrid Model specific
% F1 = ism_falpha(1,U,nEff_lyrs,vv,aa,pp,gg,oo );
% F2 = ism_falpha(2,U,nEff_lyrs,vv,aa,pp,gg,oo );
% 
% 
% tmpa = (1 + pp.c13*Cb(:).*F1)./(1 + pp.c13*Cb(:).*F2);          %Use to Surface Velocities
% tmpa = F1./F2;
% 
% u_h = u_h.*tmpa;                    
% v_h = v_h.*tmpa;
% end
% 
% imagesc(reshape(u_h,gg.nJ,gg.nI))
% %%


eq_norm(j) = norm(RHSf-LHSf*U)./norm(RHS);


vv.U = U;
vv.u = u;
vv.v = v;
vv.nEff = nEff;
vv.sstream_norm = sstream_norm;
vv.eq_norm = eq_norm;

if oo.hybrid, vv.C = C; vv.nEff_lyrs = nEff_lyrs; vv.F2 = F2; end;

end


vv2=vv;

end

      

