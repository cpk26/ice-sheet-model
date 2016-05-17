function [vv2] = ism_adjoint(vv,aa,pp,gg,oo )
%% Shallow Stream Model 
% Inputs:
%   vv      struct containing initial solution variables
%   aa      prescribed fields, including inputs and boundary conditions
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs:
%   vv2     struct containing new solution variables


%% Boundary Conditions    
[LHS, RHS] = ism_adjoint_fieldeq(vv,aa,pp,gg,oo);       %Field Equations
L = NaN(size(LHS,2),1);                               %Unmodified solution vector 

DEL = zeros(gg.nua+gg.nva,1);                   %Columns/rows to delete

if any(gg.nfxd(:))
DEL = DEL + [gg.S_u*gg.nfxd_ugrid(:); gg.S_v*gg.nfxd_vgrid(:)];  %lambda,mu = 0
end

if any(gg.nmgn(:))              %Ice Margin Nodes
tmp_a = gg.S_u*gg.nmgn_ugrid(:); tmp_a = logical(tmp_a);
tmp_b = [zeros(gg.nua,1); gg.S_v*gg.nmgn_vgrid(:)]; tmp_b = logical(tmp_b);

B = ones(gg.nua + gg.nva,1);
B(tmp_a) = 0; B(tmp_b) = 0;
RHS = RHS.*B;                   %Replace forcing at the edge

clear tmp_a tmp_b;
end

if any(gg.nperbc(:))    %Periodic BC nodes
    
tmp_a = [gg.S_u*(gg.nperbc_ugrid(:) < 0); gg.S_v*(gg.nperbc_vgrid(:) < 0)]; tmp_a = logical(tmp_a);
tmp_b = [gg.S_u*(gg.nperbc_ugrid(:) > 0); gg.S_v*(gg.nperbc_vgrid(:) > 0)]; tmp_b = logical(tmp_b); 

DEL = DEL + [tmp_a] + [tmp_b];

clear tmp_a tmp_b tmp_c tmp_d

end

DEL = logical(DEL);
    
LHS(:,DEL) = [];                                        %Apply boundary conditions
LHS(DEL,:) = [];   
RHS(DEL,:) = []; 

Lm = LHS\RHS;                                           %Solve

L(DEL) = 0;                                           %Return to L vector
L(isnan(L)) = Lm;


% if any(gg.nperbc)
% 
% tmp_a = [gg.S_u*(gg.nperbc_ugrid(:) < 0); gg.S_v*(gg.nperbc_vgrid(:) < 0)]; tmp_a = logical(tmp_a);
% tmp_b = [gg.S_u*(gg.nperbc_ugrid(:) > 0); gg.S_v*(gg.nperbc_vgrid(:) > 0)]; tmp_b = logical(tmp_b); 
% 
% L(tmp_a) = L(tmp_b);
% 
% clear tmp_a tmp_b;
% end


lambda = L(1:gg.nua);                             %lamba,mu fields
mu = L(gg.nua+1:end);

adjoint_norm = norm(RHS-LHS*Lm,oo.norm);

vv.lambda = lambda;
vv.mu = mu;
vv2=vv;

end

