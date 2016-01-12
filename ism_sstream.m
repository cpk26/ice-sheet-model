function [vv2] = ism_sstream(vv,aa,pp,gg,oo )
%% Shallow Stream Model 
% Inputs:
%   vv      struct containing initial solution variables
%   aa      prescribed fields, including inputs and boundary conditions
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs:
%   vv2     struct containing new solution variables
%   J       Jacobian matrix

numIter = 10;
iter_norm = zeros(numIter,1);
n = pp.n_Glen;                              %Ice Flow 


%% Variables (Non-Dimensionalized)
s = aa.s; s_diag = spdiags(s(:),0,gg.nIJ,gg.nIJ);      %Topography
h = aa.h; h_diag = spdiags(h(:),0,gg.nIJ,gg.nIJ);
Cslip = aa.C; Cslip_diag = spdiags(Cslip(:),0,gg.nIJ,gg.nIJ);

%Use gradient instead of gg.nddx/y 
%since periodic BC conditions do not apply
[Sx,Sy] = gradient(s, gg.dx, gg.dy);        
Sx = Sx(:); Sy = Sy(:); 

u = vv.u;                                       %Initial iterate velocity 
v = vv.v;


%% Picard Iterations
for j = 1:numIter

    [LHS, RHS] = ism_sstream_fieldeq(u,v,aa,pp,gg,oo);      %Field Equations
    U = zeros(size(LHS,2),1);                   %Unmodified velocity vector
    
    LHS(:, gg.u_BP1) = LHS(:, gg.u_BP1) + LHS(:, gg.u_BP2); % Apply Periodic BC
    LHS(:, gg.v_BP1) = LHS(:, gg.v_BP1) + LHS(:, gg.v_BP2); 
    LHS(:,[gg.u_BP2, gg.v_BP2]) = [];

    Um = LHS\RHS;               %Solve modified field equations
    
    U(gg.u_BP2) = NaN;          %Return to unmodified velocity vector             
    U(gg.v_BP2) = NaN;
    U(~isnan(U)) = Um;
    U(gg.u_BP2) = U(gg.u_BP1);
    U(gg.v_BP2) = U(gg.v_BP1);
    
    u = gg.S_u*U(1:(gg.nI+1)*gg.nJ);    %u,v velocity fields
    v = gg.S_v*U((gg.nI+1)*gg.nJ+1:end);
    
    iter_norm(j) = norm(RHS-LHS*Um,oo.norm); %iteration norm (using Um)
    
end

vv.u = u;
vv.v = v;
vv.iter_norm = iter_norm;

vv2=vv;

end


function [tbx, tby] = slidingLaw()

switch oo.sL
    case 'ismip'
        tbx = C.*u;
        tby = C.*v;
        
    case 'weertman'
        %To be implemented...
        %w = delta*(u.*Sx + v.*Sy);
        %Ub = sqrt(u.^2 + v.^2 + w.^2);
end


end


% for j = 1:numIter
% 
%     [LHS, RHS] = ism_sstream_fieldeq(u,v,aa,pp,gg,oo);
% 
%     E1 = gg.S_u_perBC;                                          %Add Periodic Boundary Conditions
%     E4 = gg.S_v_perBC;
%     E2 = zeros(size(E1,1), (gg.nJ+1)*gg.nI);
%     E3 = zeros(size(E4,1), gg.nJ*(gg.nI+1));
%     EE = [E1 E2; E3 E4];
%     
%     LL2 = [LHS;EE];                                              %Combine
%     
%     f1a = pp.c4*gg.c_hu*h_diag*Sx;                              %RHS SSA
%     f1b = pp.c4*gg.c_hv*h_diag*Sy;
%     f3 = zeros(size(EE,1),1);                                   %Periodic Boundary Conditions
%     
%     FF2 = [RHS ;f3];                                            %Combine
%     
%     tic
%     U = LL2\FF2;                                                %Solve
%     toc
%     u = gg.S_u*U(1:(gg.nI+1)*gg.nJ);
%     v = gg.S_v*U((gg.nI+1)*gg.nJ+1:end);
%     
%     iter_norm(j) = norm(FF2-LL2*U,oo.norm);
%     
% end
