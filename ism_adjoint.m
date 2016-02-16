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
vOff = sum(A);                              %offset to v values [number of u values]

A = sum(gg.S_v); A2 = cumsum(A);            %V-grid
nbnd_vind = A2(gg.nbnd_vind) + vOff;        %Boundary nodes 

%% Boundary Conditions    
[LHS, RHS] = ism_adjoint_fieldeq(vv,aa,pp,gg,oo);       %Field Equations
L = NaN(size(LHS,2),1);                               %Unmodified solution vector 

DEL = union(nbnd_uind,nbnd_vind);                       %lambda,mu = 0
    
LHS(:,DEL) = [];                                        %Apply boundary conditions
LHS(DEL,:) = [];   
RHS(DEL,:) = []; 

Lm = LHS\RHS;                                           %Solve

L(DEL) = 0;                                           %Return to L vector
L(isnan(L)) = Lm;


lambda = L(1:vOff);                             %lamba,mu fields
mu = L(vOff+1:end);

vv.lambda = lambda;
vv.mu = mu;
vv.adjoint_norm = norm(RHS-LHS*Lm,oo.norm);
vv2=vv;

end

