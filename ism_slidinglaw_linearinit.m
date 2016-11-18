function [alpha] = ism_slidinglaw_linearinit(B2,N,vv,aa,pp,gg,oo)
%% Field Equations for SSA velocities 
% Inputs:
%   v     	solution struct
%   aa      prescribed fields, including inputs and boundary conditions
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs: 
%   aa     


if ~isfield(oo,'slidinglaw'), oo.slidinglaw = 'linear'; end 


%% Linear
if strcmp(oo.slidinglaw, 'linear')                              
    alpha = B2;
    
else                                          
    
%% Non Linear

%Setup variables
N = max(N,0);
U = vv.U;


%% Basal velocities in the case of hybrid ice sheet model
if oo.hybrid
F2 = ism_falpha(2,vv.uv,vv.nEff_lyrs,vv,aa,pp,gg,oo );
tmpa = (1 + pp.c13*B2.*F2);
Ub = U./tmpa;
else
Ub = U;
end
  
%% Handle different sliding laws
% Test case for weertman/schoof, can delete
%   Ub = 100;
%   N = 1e7/pp.N;
%   B2 = 1e10;
%   

if strcmp(oo.slidinglaw, 'weertman')            %6a of Hewitt (2012)
    p = pp.p; q = pp.q;
    
    F = pp.c14 .*(N.^p .* Ub.^q);
    F = F .* (abs(Ub).^-1); 
    
    alpha = B2./F;

elseif strcmp(oo.slidinglaw, 'schoof')              %6b of Hewitt (2012)
    n = pp.n_Glen;
    F = pp.c15 .* N .* (Ub./ (pp.c16.*Ub + pp.c17.*N.^n)).^(1/n); 
    F = F .* (abs(Ub).^-1); 
    alpha = B2./F;
end

%% Reshape
alpha = reshape(gg.S_h'*alpha,gg.nJ,gg.nI);
N = reshape(gg.S_h'*N,gg.nJ,gg.nI);

% Fill regions where N (effective pressure) < eps
eps = 1e4/pp.phi;
mask = N<eps;
mask(gg.next) =0;

alpha(mask) = NaN;
alpha = inpaint_nans(alpha);

%% Reshape for output
alpha = gg.S_h*alpha(:);


end



end
