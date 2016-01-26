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
sstream_norm = zeros(numIter,1);
n = pp.n_Glen;                              %Ice Flow 


%% Variables (Non-Dimensionalized)
s = aa.s; s_diag = spdiags(s(:),0,gg.nIJ,gg.nIJ);      %Topography
h = aa.h; h_diag = spdiags(h(:),0,gg.nIJ,gg.nIJ);
if strcmp(oo.pT, 'forward'); C = aa.C; end
if strcmp(oo.pT, 'inverse'); C = vv.C; end

%Use gradient instead of gg.nddx/y 
%since periodic BC conditions do not apply
[Sx,Sy] = gradient(s, gg.dx, gg.dy);        
Sx = Sx(:); Sy = Sy(:); 

u = vv.u;                                       %Initial iterate velocity 
v = vv.v;


%% Picard Iterations
for j = 1:numIter

    [LHS, RHS] = ism_sstream_fieldeq(u,v,C,aa,pp,gg,oo);      %Field Equations
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
    
    sstream_norm(j) = norm(RHS-LHS*Um,oo.norm); %iteration norm (using Um)
    
end

vv.u = u;
vv.v = v;
vv.sstream_norm = sstream_norm;

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

