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

if isfield(oo,'runC'), runC = vv.runC;               %Variable accumalating the adjoint of C
else runC = zeros(gg.nha,1); end;                    %initialize to vv.runC if provided


if strcmp(oo.pT, 'forward'); C = aa.C; end; %Problem Type
if strcmp(oo.pT, 'inverse'); C = vv.C; end;

X = [gg.du_x gg.dv_y; gg.du_x -gg.dv_y; gg.dhu_y gg.dhv_x; speye(gg.nua,gg.nua) sparse(gg.nua,gg.nva); sparse(gg.nva,gg.nua) speye(gg.nva,gg.nva)];
X2 = [gg.dh_x gg.dh_x gg.duh_y speye(gg.nua,gg.nua) sparse(gg.nva,gg.nua)'; gg.dh_y -gg.dh_y gg.dvh_x sparse(gg.nua,gg.nva)' speye(gg.nva,gg.nva)];


DEL = zeros(gg.nua+gg.nva,1);                   %Columns to delete


firstpass = 1;

disp('Adjoint Calculation')


for j = numIter:-1:1
disp(['Inverse Picard Iteration: ', num2str(j)])

%% Initiate Variables
A_r = rr.An{j};                             %A matrix from current iteration (r indicates that it will be reduced during application of BC)
Uf = rr.Un(:,j+1);                          %Velocity from iteration n+1 (forward iteration)
b = zeros(numel(Uf),1);                     %Preallocate  b array

if firstpass
    adjU = rr.adjU; adjU_r = adjU;               %Adjoint of U from cost function for first iteration
    firstpass = 0;    
    
else
    U_adi = struct('f', Uf, 'dU',ones(gg.nua+gg.nva,1));             %Otherwise calculate from adjoint of viscosity
    
    if oo.hybrid
    UAD = ism_visc_diAD(U_adi,rr.nEffn(:,j),C(:),aa,pp,gg,oo);    
    U_visc = sparse(UAD.dU_location(:,1),UAD.dU_location(:,2), UAD.dU, UAD.dU_size(1), UAD.dU_size(2));        
    else
    UAD = ism_visc_AD(U_adi,vv,aa,pp,gg,oo);    
    U_visc = sparse(UAD.dU_location(:,1),UAD.dU_location(:,2), UAD.dU, UAD.dU_size(1), UAD.dU_size(2));
    end
    
%     adjU = U_visc'*adjnEff;
%     imagesc(reshape(adjU(1:gg.nua),75,76)); colorbar; caxis([-1e-10,1e-10])
%     UAD = ism_visc_AD(U_adi,vv,aa,pp,gg,oo);    
%     U_visc = sparse(UAD.dU_location(:,1),UAD.dU_location(:,2), UAD.dU, UAD.dU_size(1), UAD.dU_size(2));
%     adjU = U_visc'*adjnEff;
%     figure
%     imagesc(reshape(adjU(1:gg.nua),75,76)); colorbar; caxis([-1e-10,1e-10])
%     figure;
%     imagesc(reshape(adjnEff,75,75));
%     
%     adjnEff = ones(75*75,1);
%     adjU_plain = U_visc'*adjnEff;
%     imagesc(reshape(adjU_plain(1:gg.nua),75,76)); colorbar; caxis([-1e-10,1e-10])
    
    adjU = U_visc'*adjnEff; adjU_r = adjU;
    clear U_adi UAD U_visc


end

%% Handle BC for A matrix and Velocity Array 

if any(gg.nperbc)

% tmp_a = [gg.S_u*(gg.nperbc_ugrid(:) < 0); gg.S_v*(gg.nperbc_vgrid(:) < 0)]; tmp_a = logical(tmp_a);
% tmp_b = [gg.S_u*(gg.nperbc_ugrid(:) > 0); gg.S_v*(gg.nperbc_vgrid(:) > 0)]; tmp_b = logical(tmp_b); 
% tmp_c = [gg.S_u*(gg.nperbc_ugrid(:) ~= 0); gg.S_v*(gg.nperbc_vgrid(:) ~= 0)]; tmp_c = logical(tmp_c);
% 
% adjU_r(tmp_c) = 0;
% adjU_r(tmp_a) = [];
% 
% A_r(:, tmp_b) = A_r(:, tmp_b) + A_r(:, tmp_a);
% DEL = DEL + [tmp_a];


tmp_a = [gg.S_u*(gg.nperbc_ugrid(:) ~= 0); gg.S_v*(gg.nperbc_vgrid(:) ~= 0)]; tmp_a = logical(tmp_a);
%adjU_r(tmp_a) = 0;
%A_r(:,tmp_a) = [];
A_r(tmp_a,:) = [];


clear tmp_a tmp_b tmp_c;
end

%% Intermediate step to determining the adjoint of A matrix
% DEL = logical(DEL);
% A_r(:,DEL) = [];
b_r = A_r'\adjU_r;

%% Apply BC to b array

if any(gg.nperbc)
% tmp_a = [gg.S_u*(gg.nperbc_ugrid(:) < 0); gg.S_v*(gg.nperbc_vgrid(:) < 0)]; tmp_a = logical(tmp_a);
% tmp_b = [gg.S_u*(gg.nperbc_ugrid(:) > 0); gg.S_v*(gg.nperbc_vgrid(:) > 0)]; tmp_b = logical(tmp_b); 
% 
% 
% b_r(tmp_a) = b_r(tmp_b);
%b = b_r;
tmp_a = [gg.S_u*(gg.nperbc_ugrid(:) ~= 0); gg.S_v*(gg.nperbc_vgrid(:) ~= 0)]; tmp_a = logical(tmp_a);
b(~tmp_a) = b_r;

% tmp_a = [gg.S_u*(gg.nperbc_ugrid(:) ~= 0); gg.S_v*(gg.nperbc_vgrid(:) ~= 0)]; tmp_a = logical(tmp_a);
% b(~tmp_a) = b_r;
clear tmp_a;
end   
  
%% Determine adjoint of A matrix
adjA = -b*Uf';  

clear b tmp_a;

%% Determine adjoint of Viscosity                                   
nEff_adi = struct('f', rr.nEffn(:,j), 'dnEff',ones(gg.nha,1));   %Calculate main diagonal of D matrix (A = X2*D*X)
nEffAD = ism_dim_Ddiag_ADnEff(C(:),nEff_adi,aa,pp,gg,oo);
nEff_Ddiag = sparse(nEffAD.dnEff_location(:,1),nEffAD.dnEff_location(:,2), nEffAD.dnEff, nEffAD.dnEff_size(1), nEffAD.dnEff_size(2));


adjnEff = zeros(gg.nha,1);
for i = 1:gg.nha
tmp = X2*spdiags(nEff_Ddiag(:,i),0,nEffAD.dnEff_size(1),nEffAD.dnEff_size(1))*X;
adjnEff(i) = tmp(:)'*adjA(:);

end


clear nEff_adi nEffAD nEff_Ddiag tmp;

%% Determine adjoint of Basal Slip
C_adi = struct('f', C(:), 'dC',ones(gg.nha,1));   %Calculate main diagonal of D matrix (A = X2*D*X)
CAD = ism_dim_Ddiag_ADc(C_adi,rr.nEffn(:,j),aa,pp,gg,oo);
C_Ddiag = sparse(CAD.dC_location(:,1),CAD.dC_location(:,2), CAD.dC, CAD.dC_size(1), CAD.dC_size(2));


adjC = zeros(gg.nha,1);
for i = 1:gg.nha
tmp = X2*spdiags(C_Ddiag(:,i),0,CAD.dC_size(1),CAD.dC_size(1))*X;
adjC(i) = tmp(:)'*adjA(:);
end

runC = runC + adjC;
%imagesc(reshape(runC,75,75))
clear C_adi CAD C_Ddiag tmp;


end

disp('Finished looping through picard iterations [Adjoint]')
rr.adjC = adjC;
rr.runC = runC;

end

