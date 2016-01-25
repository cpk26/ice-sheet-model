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


%% Variables (Non-Dimensionalized)     
[LHS, RHS] = ism_adjoint_fieldeq(vv,aa,pp,gg,oo);      %Field Equations
lm = zeros(size(LHS,2),1);

bndry_ugrid = zeros(1,gg.nJ*(gg.nI+1));     %Boundary Conditions (lamdba = mu = 0)
bndry_vgrid = zeros(1,(gg.nJ+1)*gg.nI);

bndry_ugrid([1:gg.nJ,end-gg.nJ-1:end]) = 1;               %Lateral Boundaries (u-grid)
bndry_vgrid([1:gg.nJ+1:end,gg.nJ+1:gg.nJ+1:end]) = 1;   %Top/bottom Boundaries (v-grid)

bndryPts = logical([bndry_ugrid, bndry_vgrid]);        %Combine into one vector 
    
LHS(:,bndryPts) = [];   %Apply boundary conditions
lm_mod = LHS\RHS;       %Solve

lm(~bndryPts) = lm_mod; %Return to lm vector
lambda = gg.S_u'*lm(1:(gg.nI+1)*gg.nJ);    %lamba,mu fields
mu = gg.S_v'*lm((gg.nI+1)*gg.nJ+1:end);

vv.lambda = lambda;
vv.mu = mu;

vv2=vv;

end

