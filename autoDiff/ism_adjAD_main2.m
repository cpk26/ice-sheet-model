function [rr] = ism_adjAD_main(vv,rr ,aa,pp,gg,oo )
%% Shallow Stream Model 
% Inputs:
%   vv      struct containing initial solution variables
%   aa      prescribed fields, including inputs and boundary conditions
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs:

disp('Adjoint Backbone')
tic


numSaveIter = size(rr.uvn,2); 

fi = numSaveIter;               %Forward Iteration
ci = numSaveIter - 1;           %Current Iteration
pi = max(ci-1,1);               %Previous Iteration
    
    
%Initiate Constant matrices
X = [gg.du_x gg.dv_y; gg.du_x -gg.dv_y; gg.dhu_y gg.dhv_x; speye(gg.nua,gg.nua) sparse(gg.nua,gg.nva); sparse(gg.nva,gg.nua) speye(gg.nva,gg.nva)];
X2 = [gg.dh_x gg.dh_x gg.duh_y speye(gg.nua,gg.nua) sparse(gg.nva,gg.nua)'; gg.dh_y -gg.dh_y gg.dvh_x sparse(gg.nua,gg.nva)' speye(gg.nva,gg.nva)];

%Initial Adjoint of U [from cost function for first iteration]
adjU = rr.adjU; 
adjU_r = adjU;               

%% Initiate Variables
A_r = rr.An{ci};                             %A matrix from current iteration (r indicates that it will be reduced during application of BC)
A_sp = (rr.An{ci} ~= 0);                     %Sparsity pattern
b = zeros(gg.nua+gg.nva,1);          %Preallocate  b array


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

adjA = reshape(-spdiags(b,0,gg.nua+gg.nva,gg.nua+gg.nva)*A_sp*spdiags(rr.uvn(:,fi),0,gg.nua+gg.nva,gg.nua+gg.nva),(gg.nua+gg.nva)^2,1);
clear A_sp b tmp_a;


%% Determine Adjoints
%% Determine adjoint of Viscosity
if oo.adjAD_nEff
 
%Adjoint of Viscosity
nEff_adi = struct('f', rr.nEffn(:,ci), 'dnEff',ones(gg.nha,1));   %Calculate main diagonal of D matrix (A = X2*D*X)
nEffAD = ism_dim_Ddiag_ADnEff(rr.Cn(:,ci),nEff_adi,aa,pp,gg,oo);
nEff_Ddiag = sparse(nEffAD.dnEff_location(:,1),nEffAD.dnEff_location(:,2), nEffAD.dnEff, nEffAD.dnEff_size(1), nEffAD.dnEff_size(2));


adjnEff = zeros(gg.nha,1);
matDim = nEffAD.dnEff_size(1);
matDim2 = (gg.nua+gg.nva)^2;

parfor i = 1:gg.nha
adjnEff(i) = reshape(X2*spdiags(nEff_Ddiag(:,i),0,matDim,matDim)*X,matDim2,1)'*adjA;
end

rr.adjnEff = adjnEff;
end
clear nEff_adi nEffAD nEff_Ddiag U_adi UAD U_visc adjnEff tmp matDim;


%% Determine adjoint of Effective Drag
if oo.adjAD_C
    
C_adi = struct('f', rr.Cn(:,ci), 'dC',ones(gg.nha,1));   %Calculate main diagonal of D matrix (A = X2*D*X)
CAD = ism_dim_Ddiag_ADc(C_adi,rr.nEffn(:,ci),aa,pp,gg,oo);
C_Ddiag = sparse(CAD.dC_location(:,1),CAD.dC_location(:,2), CAD.dC, CAD.dC_size(1), CAD.dC_size(2));

adjC = zeros(gg.nha,1);
matDim = CAD.dC_size(1);
matDim2 = (gg.nua+gg.nva)^2;
parfor i = 1:gg.nha
adjC(i) = reshape(X2*spdiags(C_Ddiag(:,i),0,matDim,matDim)*X,matDim2,1)'*adjA;
end

rr.adjC = adjC;
clear C_adi CAD C_Ddiag tmp;

end


disp(['Finished:', num2str(toc)])

end


