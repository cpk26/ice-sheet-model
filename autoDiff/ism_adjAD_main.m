function [rr] = ism_adjAD_main(vv,rr ,aa,pp,gg,oo )
%% Shallow Stream Model 
% Inputs:
%   vv      struct containing initial solution variables
%   aa      prescribed fields, including inputs and boundary conditions
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs:
%   vv2     struct containing new solution variables

numIter = oo.pic_iter;                               %Number of Picard Iterations

if isfield(rr,'runC'), runC = rr.runC;               %Variable accumalating the adjoint of C
else runC = zeros(gg.nha,1); end;                    %initialize to vv.runC if provided

                                             
if oo.hybrid, Cb = vv.Cb;                           %Initialize Basal Drag
else C = vv.C; end

X = [gg.du_x gg.dv_y; gg.du_x -gg.dv_y; gg.dhu_y gg.dhv_x; speye(gg.nua,gg.nua) sparse(gg.nua,gg.nva); sparse(gg.nva,gg.nua) speye(gg.nva,gg.nva)];
X2 = [gg.dh_x gg.dh_x gg.duh_y speye(gg.nua,gg.nua) sparse(gg.nva,gg.nua)'; gg.dh_y -gg.dh_y gg.dvh_x sparse(gg.nua,gg.nva)' speye(gg.nva,gg.nva)];


DEL = zeros(gg.nua+gg.nva,1);                   %Columns to delete
DEL2 = DEL;                                     %Rows to delete                 

firstpass = 1;

disp('Adjoint Calculation')


for j = numIter:-1:1
disp(['Inverse Picard Iteration: ', num2str(j)])

%% Initiate Variables
A_r = rr.An{j};                             %A matrix from current iteration (r indicates that it will be reduced during application of BC)
A_sp = (rr.An{j} ~= 0);
Uf = rr.Un(:,j+1);                          %Velocity from iteration n+1 (forward iteration)
b = zeros(numel(Uf),1);                     %Preallocate  b array

if oo.hybrid, 
F2 = rr.F2n(:,j);
C = Cb(:)./(1 + (pp.c13*Cb(:)).*(gg.S_h'*F2)); 

end;

if firstpass
    adjU = rr.adjU; adjU_r = adjU;               %Adjoint of U from cost function for first iteration
    firstpass = 0;    
    
else
    U_adi = struct('f', Uf, 'dU',ones(gg.nua+gg.nva,1));             %Otherwise calculate from adjoint of viscosity
    
    if oo.hybrid
    UAD = ism_visc_diSAD(U_adi,rr.nEffn(:,j),gg.S_h*C(:),aa,pp,gg,oo);    
    U_visc = sparse(UAD.dU_location(:,1),UAD.dU_location(:,2), UAD.dU, UAD.dU_size(1), UAD.dU_size(2));        
    else
    UAD = ism_visc_AD(U_adi,vv,aa,pp,gg,oo);    
    U_visc = sparse(UAD.dU_location(:,1),UAD.dU_location(:,2), UAD.dU, UAD.dU_size(1), UAD.dU_size(2));
    end
    
    adjU = U_visc'*adjnEff; adjU_r = adjU;
    clear U_adi UAD U_visc


end

%% Handle BC for A matrix and Velocity Array 

if any(gg.nmgn(:)); %ice margin nodes
tmp_a = [gg.S_u*gg.nmgn_ugrid(:); gg.S_v*gg.nmgn_vgrid(:)];

DEL = DEL + tmp_a;
DEL2 = DEL2 + tmp_a;

clear tmp_a 
end; 

if any(gg.nfxd(:))  %Dirichlet BC Nodes
tmp_a = [gg.S_u*gg.nfxd_ugrid(:); gg.S_v*gg.nfxd_vgrid(:)];               %Location of fixed values

DEL = DEL + tmp_a;
DEL2 = DEL2 + tmp_a;

clear tmp_a 
end

if any(gg.nperbc)   %Periodic BC

tmp_a = [gg.S_u*(gg.nperbc_ugrid(:) < 0); gg.S_v*(gg.nperbc_vgrid(:) < 0)]; tmp_a = logical(tmp_a);
tmp_b = [gg.S_u*(gg.nperbc_ugrid(:) > 0); gg.S_v*(gg.nperbc_vgrid(:) > 0)]; tmp_b = logical(tmp_b);

A_r(tmp_b,:) = A_r(tmp_b,:) + A_r(tmp_a,:);

DEL2 = DEL2 + tmp_a;
clear tmp_a tmp_b tmp_c;

end

DEL = logical(DEL);
DEL2 = logical(DEL2);


A_r(:,DEL) = [];
A_r(DEL2,:) = [];


adjU_r(DEL) = [];   

%% Intermediate step to determining the adjoint of A matrix
b_r = A_r'\adjU_r;

%% Return to original velocity vector
b(~DEL2) = b_r;

%% Apply BC to b array

if any(gg.nperbc)

tmp_a = [gg.S_u*(gg.nperbc_ugrid(:) < 0); gg.S_v*(gg.nperbc_vgrid(:) < 0)]; tmp_a = logical(tmp_a);
tmp_b = [gg.S_u*(gg.nperbc_ugrid(:) > 0); gg.S_v*(gg.nperbc_vgrid(:) > 0)]; tmp_b = logical(tmp_b);

b(tmp_a) = b(tmp_b);
clear tmp_a tmp_b;

end   
  


%% Determine adjoint of A matrix

adjA = -spdiags(b,0,gg.nua+gg.nva,gg.nua+gg.nva)*A_sp*spdiags(Uf,0,gg.nua+gg.nva,gg.nua+gg.nva);

clear b tmp_a;

%% Determine adjoint of Viscosity                                   
nEff_adi = struct('f', rr.nEffn(:,j), 'dnEff',ones(gg.nha,1));   %Calculate main diagonal of D matrix (A = X2*D*X)
nEffAD = ism_dim_Ddiag_ADnEff(gg.S_h*C(:),nEff_adi,aa,pp,gg,oo);
nEff_Ddiag = sparse(nEffAD.dnEff_location(:,1),nEffAD.dnEff_location(:,2), nEffAD.dnEff, nEffAD.dnEff_size(1), nEffAD.dnEff_size(2));


adjnEff = zeros(gg.nha,1);
parfor i = 1:gg.nha
tmp = X2*spdiags(nEff_Ddiag(:,i),0,nEffAD.dnEff_size(1),nEffAD.dnEff_size(1))*X;
adjnEff(i) = tmp(:)'*adjA(:);

end


clear nEff_adi nEffAD nEff_Ddiag tmp;

%% Determine adjoint of Basal Slip
C_adi = struct('f', gg.S_h*C(:), 'dC',ones(gg.nha,1));   %Calculate main diagonal of D matrix (A = X2*D*X)
CAD = ism_dim_Ddiag_ADc(C_adi,rr.nEffn(:,j),aa,pp,gg,oo);
C_Ddiag = sparse(CAD.dC_location(:,1),CAD.dC_location(:,2), CAD.dC, CAD.dC_size(1), CAD.dC_size(2));


adjC = zeros(gg.nha,1);
parfor i = 1:gg.nha
tmp = X2*spdiags(C_Ddiag(:,i),0,CAD.dC_size(1),CAD.dC_size(1))*X;
adjC(i) = tmp(:)'*adjA(:);
end


if oo.hybrid
flag = 0;
C_adi = struct('f', gg.S_h*C(:), 'dC',ones(gg.nha,1));   %Calculate main diagonal of D matrix (A = X2*D*X)
CAD = ism_cslip_form_ADc(flag,F2,C_adi,aa,pp,gg,oo );
C_form = sparse(CAD.dC_location(:,1),CAD.dC_location(:,2), CAD.dC, CAD.dC_size(1), CAD.dC_size(2));

adjC = C_form'*adjC;
end

if oo.hybrid
adjC = adjC.*(1+(pp.c13*gg.S_h*Cb(:)).*F2).^-2; %Move from C effective to C basal
end

runC = runC + adjC;


clear C_adi CAD C_Ddiag tmp;


end

disp('Finished looping through picard iterations [Adjoint]')
rr.adjC = gg.S_h'*adjC;
rr.runC = gg.S_h'*runC;

end

