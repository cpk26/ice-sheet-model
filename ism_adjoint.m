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


%% Remap indices [from whole region, to masked area]

A = sum(gg.S_u); A2 = cumsum(A);            %U-grid
nbnd_uind = A2(gg.nbnd_uind);               %Boundary nodes 
nfxd_uind = A2(gg.nfxd_uind);               %Fixed nodes
nmgn_uind = A2(gg.nmgn_uind);               %Margin Nodes


vOff = sum(A);                              %offset to v values [number of u values]

A = sum(gg.S_v); A2 = cumsum(A);            %V-grid
nbnd_vind = A2(gg.nbnd_vind) + vOff;        %Boundary nodes 
nfxd_vind = A2(gg.nfxd_vind) + vOff;        %Fixed nodes 
nmgn_vind = A2(gg.nmgn_vind)+ vOff;         %Margin Nodes

%% Boundary Conditions    
[LHS, RHS] = ism_adjoint_fieldeq(vv,aa,pp,gg,oo);       %Field Equations
L = NaN(size(LHS,2),1);                               %Unmodified solution vector 

DEL = [];                   %Columns to delete
DEL2 = [];                  %Rows to delete

if any(gg.nfxd(:))
DEL = union(nfxd_uind,nfxd_vind);                       %lambda,mu = 0
end

if any(gg.nmgn(:))  %Ice Margin Nodes
B = ones(numel(RHS),1); 
B(nmgn_uind) = 1/(pp.x*gg.dx); B(nmgn_vind) = 1/(pp.x*gg.dy);

RHS = RHS.*B;           %Replace forcing at the edge
end
    
LHS(:,DEL) = [];                                        %Apply boundary conditions
LHS(DEL,:) = [];   
RHS(DEL,:) = []; 

Lm = LHS\RHS;                                           %Solve

L(DEL) = 0;                                           %Return to L vector
L(isnan(L)) = Lm;


lambda = L(1:vOff);                             %lamba,mu fields
mu = L(vOff+1:end);

adjoint_norm = norm(RHS-LHS*Lm,oo.norm);

vv.lambda = lambda;
vv.mu = mu;
vv2=vv;

end

