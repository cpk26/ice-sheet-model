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
    U = Inf(size(LHS,2),1);                   %Unmodified velocity vector
    
    %% Remap indices [from whole region, to masked area]
    A = sum(gg.S_u); A2 = cumsum(A);  
    vOff = sum(A);                      %offset to v values [number of u values]
    nfxd_uind2 = A2(gg.nfxd_uind);
    
    A = sum(gg.S_v); A2 = cumsum(A); 
    nfxd_vind2 = A2(gg.nfxd_vind);
    
    
    
    
    %% Boundary Conditions
    
    %Need to figure out map from nfxd_u/v ind > new mapping
    %something like A = cumsum(S_u); indNew = A(indOld)??
    
    
    %Apply BC
    DEL = [];
    if any(gg.nfxd(:))
    RHS = RHS - LHS(:,nfxd_uind2)*aa.nfxd_uval;
    RHS = RHS - LHS(:,nfxd_vind2)*aa.nfxd_vval;    
    DEL = union(gg.nfxd_uind, gg.nfxd_vind);
    end
    
    if any(gg.nperbc(:))
    LHS(:, gg.nperbc_u1ind) = LHS(:, gg.nperbc_u1ind) + LHS(:, gg.nperbc_u2ind); % Apply Periodic BC
    LHS(:, gg.nperbc_v1ind) = LHS(:, gg.nperbc_v1ind) + LHS(:, gg.nperbc_v2ind); 
    DEL = union(DEL, [gg.nperbc_u2ind; gg.nperbc_v2ind]);
    end

    LHS(:,DEL) = [];
    
    %Solve 
    Um = LHS\RHS;               %Solve modified field equations
    
    %Return to original velocity vector
    U(DEL) = NaN;
    U(~isnan(U)) = Um;
    
    if any(gg.nfxd(:))
        U(nfxd_uind2) = aa.nfxd_uval;
        U(nfxd_vind2) = aa.nfxd_vval;
    end
    
    if any(gg.nperbc)
        U(gg.nperbc_u2ind) = U(gg.nperbc_u1ind);
        U(gg.nperbc_v2ind) = U(gg.nperbc_v1ind);
    end
    
    u = U(1:vOff);    %u,v velocity fields
    v = U(vOff+1:end);
    
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

