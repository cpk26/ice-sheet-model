function [uv] = adi_deism_backbone(uvprev, C,nEff,vv,aa,pp,gg,oo )

%% Field Equations
[LHS, RHS] = ism_deism_fieldeq(C,nEff, aa,pp,gg,oo);              %Field Equations

uv = Inf(size(LHS,2),1);                                          %Velocity vector, full length   
                       

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
F1 = adi_falpha(1,uvprev,C,vv,aa,pp,gg,oo );          %Calculate F alpha factors 
F2 = adi_falpha(2,uvprev,C,vv,aa,pp,gg,oo );
tmp_c = (1 + (pp.c13*C).*(F2-F1));          %uvn -> uvs factor
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

uvm = LHS\RHS;               %Solve modified field equations

%% Return to original velocity vector

uv(DEL) = NaN;
uv(~isnan(uv)) = uvm;

if any(gg.nfxd(:))
tmp_a = [gg.S_u*gg.nfxd_ugrid(:); gg.S_v*gg.nfxd_vgrid(:)]; tmp_a = logical(tmp_a);
tmp_b = [gg.S_u*aa.nfxd_uval(:); gg.S_v*aa.nfxd_vval(:)];              

if oo.hybrid,                                     %Convert surface vel -> effective vel
tmp_b = tmp_b./effVelFac;                         %Use previous factor for efficiency
end

uv(tmp_a) = tmp_b(tmp_a);
clear tmp_a tmp_b;
end

if any(gg.nperbc(:))

tmp_a = [gg.S_u*(gg.nperbc_ugrid(:) < 0); gg.S_v*(gg.nperbc_vgrid(:) < 0)]; tmp_a = logical(tmp_a);
tmp_b = [gg.S_u*(gg.nperbc_ugrid(:) > 0); gg.S_v*(gg.nperbc_vgrid(:) > 0)]; tmp_b = logical(tmp_b); 

uv(tmp_a) = uv(tmp_b);

clear tmp_a tmp_b;
end

end

