function [rr] = ism_adjAD_main(vv,rr ,aa,pp,gg,oo )
%% Shallow Stream Model 
% Inputs:
%   vv      struct containing initial solution variables
%   aa      prescribed fields, including inputs and boundary conditions
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs:

numPicIter = oo.pic_iter;                                  %Number of Picard Iterations
numAdjIter = oo.adj_iter;                                  %Number of Adjoint Iterations
numSaveIter = size(rr.uvn,2);                              %Number of saved picard iterations

runalpha = zeros(gg.nha,1);                                 %Variable accumalating the adjoint of alpha                   

%Initiate Constant matrices
X = [gg.du_x gg.dv_y; gg.du_x -gg.dv_y; gg.dhu_y gg.dhv_x; speye(gg.nua,gg.nua) sparse(gg.nua,gg.nva); sparse(gg.nva,gg.nua) speye(gg.nva,gg.nva)];
X2 = [gg.dh_x gg.dh_x gg.duh_y speye(gg.nua,gg.nua) sparse(gg.nva,gg.nua)'; gg.dh_y -gg.dh_y gg.dvh_x sparse(gg.nua,gg.nva)' speye(gg.nva,gg.nva)];

%Initial Adjoint of U [from cost function for first iteration]
adjU = rr.adjU; adjU_r = adjU;               

disp('Adjoint Calculation')


for adjIter = 1:numAdjIter;
disp(['Inverse Picard Iteration: ', num2str(adjIter)])

j = numSaveIter - adjIter;                  %Pass through saved iterations in reverse order, starting at second to last.

%% Initiate Variables
A_r = rr.An{j};                             %A matrix from current iteration (r indicates that it will be reduced during application of BC)
A_sp = (rr.An{j} ~= 0);                     %Sparsity pattern
b = zeros(gg.nua+gg.nva,1);          %Preallocate  b array

k = max(j-1,1);
if oo.hybrid, 
C = rr.Cbn(:,j)./(1 + (pp.c13*rr.Cbn(:,j)).*(rr.F2n(:,j))); 
Cp = rr.Cbn(:,k)./(1 + (pp.c13*rr.Cbn(:,k)).*(rr.F2n(:,k))); 
else
C = rr.Cbn(:,j);
Cp = rr.Cbn(:,k);
end;
clear k

%% Handle BC for A matrix and Velocity Array 

DEL = zeros(gg.nua+gg.nva,1);                   %Columns to delete
DEL2 = DEL;                                     %Rows to delete

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

clear A_r adjU_r adjU_r

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

adjA = reshape(-spdiags(b,0,gg.nua+gg.nva,gg.nua+gg.nva)*A_sp*spdiags(rr.uvn(:,j+1),0,gg.nua+gg.nva,gg.nua+gg.nva),(gg.nua+gg.nva)^2,1);
clear A_sp b tmp_a;


%% Determine Adjoints
%Blank variable for adjU, which is composed of components adjU_n and adjU_C
adjU = 0; 

%% Determine adjoint of Viscosity, then adjoint of uv
if oo.adjAD_uv
 
%Adjoint of Viscosity
nEff_adi = struct('f', rr.nEffn(:,j), 'dnEff',ones(gg.nha,1));   %Calculate main diagonal of D matrix (A = X2*D*X)
nEffAD = ism_dim_Ddiag_ADnEff(C,nEff_adi,aa,pp,gg,oo);
nEff_Ddiag = sparse(nEffAD.dnEff_location(:,1),nEffAD.dnEff_location(:,2), nEffAD.dnEff, nEffAD.dnEff_size(1), nEffAD.dnEff_size(2));


adjnEff = zeros(gg.nha,1);
matDim = nEffAD.dnEff_size(1);
matDim2 = (gg.nua+gg.nva)^2;
parfor i = 1:gg.nha
adjnEff(i) = reshape(X2*spdiags(nEff_Ddiag(:,i),0,matDim,matDim)*X,matDim2,1)'*adjA;
end

%Adjoint of Velocity from viscosity
U_adi = struct('f', rr.uvn(:,j), 'dU',ones(gg.nua+gg.nva,1));

if oo.hybrid
    UAD = ism_visc_diSADu(U_adi,rr.nEff_lyrsn{max(j-1,1)},Cp,aa,pp,gg,oo);
    U_visc = sparse(UAD.dU_location(:,1),UAD.dU_location(:,2), UAD.dU, UAD.dU_size(1), UAD.dU_size(2));
else
    UAD = ism_visc_ADu(U_adi,vv,aa,pp,gg,oo);
    U_visc = sparse(UAD.dU_location(:,1),UAD.dU_location(:,2), UAD.dU, UAD.dU_size(1), UAD.dU_size(2));
end

adjU_n = U_visc'*adjnEff; 
adjU = adjU + adjU_n;

%Adjoint of C from viscosity, alpha from C [Hybrid only]
if oo.hybrid
%Adjoint of C from viscosity
C_adi = struct('f', Cp, 'dC',ones(gg.nha,1));
CAD = ism_visc_diSADc(rr.uvn(:,j),rr.nEff_lyrsn{max(j-1,1)},C_adi,aa,pp,gg,oo);
U_C = sparse(CAD.dC_location(:,1),CAD.dC_location(:,2), CAD.dC, CAD.dC_size(1), CAD.dC_size(2));
adjC_n = U_C'*adjnEff; 

%Move from Ceff to C basal
flag = 1;
C_adi = struct('f', C, 'dC',ones(gg.nha,1));   
CAD = ism_cslip_form_ADc(flag,rr.F2n(:,j),C_adi,aa,pp,gg,oo );
C_form = sparse(CAD.dC_location(:,1),CAD.dC_location(:,2), CAD.dC, CAD.dC_size(1), CAD.dC_size(2));
adjC_n = C_form'*adjC_n;

%Determine adjoint of alpha from adjoint of C
alpha_adi = struct('f', rr.alpha, 'dalpha',ones(gg.nha,1));
AAD = ism_slidinglaw_ADa(alpha_adi,rr.uvn(:,j),rr.Cbn(:,max(j-1,1)),rr.F2n(:,j),vv,aa,pp,gg,oo);
C_alpha = sparse(AAD.dalpha_location(:,1),AAD.dalpha_location(:,2), AAD.dalpha, AAD.dalpha_size(1), AAD.dalpha_size(2));
adjalpha = C_alpha'*adjC_n;

runalpha = runalpha + adjalpha;
end

clear nEff_adi nEffAD nEff_Ddiag U_adi UAD U_visc adjnEff tmp matDim;
end

%% Determine adjoint of Basal Slip
if oo.adjAD_alpha
    
C_adi = struct('f', C, 'dC',ones(gg.nha,1));   %Calculate main diagonal of D matrix (A = X2*D*X)
CAD = ism_dim_Ddiag_ADc(C_adi,rr.nEffn(:,j),aa,pp,gg,oo);
C_Ddiag = sparse(CAD.dC_location(:,1),CAD.dC_location(:,2), CAD.dC, CAD.dC_size(1), CAD.dC_size(2));


adjC_A = zeros(gg.nha,1);
matDim = CAD.dC_size(1);
matDim2 = (gg.nua+gg.nva)^2;
parfor i = 1:gg.nha
adjC_A(i) = reshape(X2*spdiags(C_Ddiag(:,i),0,matDim,matDim)*X,matDim2,1)'*adjA;
end



%Move from C effective to C basal [Hybrid Only]
if oo.hybrid,                       
flag = 1;
C_adi = struct('f', C, 'dC',ones(gg.nha,1));   
CAD = ism_cslip_form_ADc(flag,rr.F2n(:,j),C_adi,aa,pp,gg,oo );
C_form = sparse(CAD.dC_location(:,1),CAD.dC_location(:,2), CAD.dC, CAD.dC_size(1), CAD.dC_size(2));

adjC_A = C_form'*adjC_A;
end


%Determine adjoint of alpha from adjoint of C
alpha_adi = struct('f', rr.alpha, 'dalpha',ones(gg.nha,1));
AAD = ism_slidinglaw_ADa(alpha_adi,rr.uvn(:,j),rr.Cbn(:,max(j-1,1)),rr.F2n(:,j),vv,aa,pp,gg,oo);
C_alpha = sparse(AAD.dalpha_location(:,1),AAD.dalpha_location(:,2), AAD.dalpha, AAD.dalpha_size(1), AAD.dalpha_size(2));
adjalpha = C_alpha'*adjC_A;

%Accumulate adjoints of alpha    
runalpha = runalpha + adjalpha;

%Determine adjoint of u from adjoint of C
U_adi = struct('f', rr.uvn(:,j), 'dU',ones(gg.nua+gg.nva,1));
UAD = ism_slidinglaw_ADu(rr.alpha,U_adi,rr.Cbn(:,max(j-1,1)),rr.F2n(:,j),vv,aa,pp,gg,oo);
C_U = sparse(UAD.dU_location(:,1),UAD.dU_location(:,2), UAD.dU, UAD.dU_size(1), UAD.dU_size(2));
adjU_C = C_U'*adjC_A;

adjU = adjU + adjU_C;
    

clear C_adi CAD C_Ddiag tmp;

end

%% Store Adjoints
rr.runalpha = runalpha;

adjU_r = adjU;
rr.adjU = adjU;

end

disp('Finished looping through picard iterations [Adjoint]')


end

