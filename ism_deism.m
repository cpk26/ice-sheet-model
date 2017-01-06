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

numPicIter = oo.pic_iter;                      %Solver parameters
numAdjIter = oo.adj_iter;

adjFlag = (numAdjIter > 0);
adjIterThresh = max(numPicIter - numAdjIter - 1, 0);
adjIterPtr = 1;

sstream_norm = zeros(numPicIter,1);
eq_norm = zeros(numPicIter,1);



uv = vv.uv;                                   %Initial Values
Cb = vv.Cb; 
C = vv.C; 
nEff = vv.nEff; 
if oo.hybrid
nEff_lyrs = vv.nEff_lyrs; 
F2 = vv.F2;
end

if strcmp(oo.pT, 'forward'), alpha = aa.alpha;  
else alpha = vv.alpha; end;

rr = struct();                              %Preallocate arrays if we are saving picard iterations
if numAdjIter > 0
    if isequal(numAdjIter,numPicIter), k = 1; %Save previous iteration if applicable
    else, k = 2; end;
    rr.An = cell(1,numAdjIter);
    rr.uvn = zeros(gg.nua+gg.nva,numAdjIter+k);
    rr.nEffn = zeros(gg.nha,numAdjIter);
    %if oo.hybrid; 
        rr.nEff_lyrsn = cell(1,numAdjIter);
        rr.F2n = zeros(gg.nha,numAdjIter+k);
        rr.Cbn = zeros(gg.nha,numAdjIter+k); %end
end

if isequal(adjIterThresh, 0)   %Save initial velocity/viscosity; 
rr.alpha = alpha;
rr.uvn(:,1) = uv; 
rr.nEffn(:,1) = nEff(:);
rr.Cbn(:,1) = Cb(:); 
%rr.Cn(:,1) = C(:); 
if oo.hybrid
rr.nEff_lyrsn{1} = nEff_lyrs;
rr.F2n(:,1) = vv.F2(:); end;
end                      

        
%% Picard Iterations
for j = 1:numPicIter

    

%% Field Equations
[LHS, RHS] = ism_deism_fieldeq(C,nEff, aa,pp,gg,oo);              %Field Equations
LHSf = LHS;
RHSf = RHS;
sstream_norm(j) = norm(RHS-LHS*uv)./norm(RHS); %iteration norm (using last iterations uv)

uv = Inf(size(LHS,2),1);                                                 %Velocity vector, full length   

if adjFlag && j > adjIterThresh,   %Save
rr.An{adjIterPtr} = LHS;  
adjIterPtr = adjIterPtr + 1;    %Increment
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
F1 = ism_falpha(1,uv,nEff_lyrs,vv,aa,pp,gg,oo );          %Calculate F alpha factors 
F2 = ism_falpha(2,uv,nEff_lyrs,vv,aa,pp,gg,oo );
tmp_c = (1 + (pp.c13*Cb).*F1)./(1 + (pp.c13*Cb).*F2);          %uvn -> uvs factor
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


u = uv(1:gg.nua);    %u,v velocity fields
v = uv(gg.nua+1:end);

%% Viscosity, effective basal drag    
    
if oo.hybrid                                            %Determine viscosity, basal slipperiness appropriately

[nEff, nEff_lyrs] = ism_visc_di(uv,nEff_lyrs,C,aa,pp,gg,oo); %Updated Viscosity
F2 = ism_falpha(2,uv,nEff_lyrs,vv,aa,pp,gg,oo );
[Cb] = ism_slidinglaw(alpha,vv.uv,Cb,F2,vv,aa,pp,gg,oo);

C = Cb(:)./(1 + (pp.c13*Cb).*(F2));                    %Effective Basal Slipperiness

else
nEff = ism_visc(uv,vv,aa,pp,gg,oo);
[C] = ism_slidinglaw(alpha,vv.uv,Cb,[],vv,aa,pp,gg,oo);
end             %SSA Viscosity



if adjFlag && j >= adjIterThresh
    rr.uvn(:,adjIterPtr) = uv; 
    rr.nEffn(:,adjIterPtr) = nEff(:);
    rr.Cbn(:,adjIterPtr) = Cb;
    %rr.Cn(:,adjIterPtr) = C;
    if oo.hybrid; 
        rr.F2n(:,adjIterPtr) = F2; 
    rr.nEff_lyrsn{adjIterPtr} = nEff_lyrs;
    end
end                %Save Intermediate velocity array


%% Plot Velocities [For manual testing of picard iteration]
% u_h = gg.S_h'*gg.c_uh*u; %Solution Velocities
% v_h = gg.S_h'*gg.c_vh*v;
% 
% if oo.hybrid                        %Hybrid Model specific
% F1 = ism_falpha(1,uv,nEff_lyrs,vv,aa,pp,gg,oo );
% F2 = ism_falpha(2,uv,nEff_lyrs,vv,aa,pp,gg,oo );
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


eq_norm(j) = norm(RHSf-LHSf*uv)./norm(RHS);


vv.uv = uv;          %Store variables
vv.u = u;
vv.v = v;
vv.nEff = nEff;            
vv.C = C;
if oo.hybrid
vv.Cb = Cb;
vv.nEff_lyrs = nEff_lyrs; 
vv.F2 = F2; 
else
vv.Cb = C;
end;      



end

%% Append solution norms
vv.sstream_norm = sstream_norm;
vv.eq_norm = eq_norm;

%% Calculate effective velocity and basal velocity
U = sqrt( (gg.c_uh*u).^2 + (gg.c_vh*v).^2 );
if oo.hybrid
tmpa = (1 + (pp.c13*vv.Cb).*(F2)); %use corresponding Cb to u/v
Ub = U./tmpa;

vv.U = U;
vv.Ub = Ub;
    
else
vv.U = U;
vv.Ub = vv.U;
end

%% Save + return
vv2=vv;

end

      

