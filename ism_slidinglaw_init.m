function [aa] = ism_slidinglaw_init(B2, vv,aa,pp,gg,oo)
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
if isequal(oo.slidinglaw, 1)                              
    aa.B2 = B2;
    
else                                          
    
%% Non Linear

%Setup variables
N = max(gg.S_h*aa.N(:),0);
U = sqrt( (gg.c_uh*vv.u).^2 + (gg.c_vh*vv.v).^2 );

%% Basal velocities in the case of hybrid ice sheet model
if oo.hybrid
F2 = ism_falpha(2,vv.U,vv.nEff_lyrs,vv,aa,pp,gg,oo );
tmpa = (1 + pp.c13*B2(:).*F2);
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

if isequal(oo.slidinglaw, 2)            %6a of Hewitt (2012)
    p = pp.p; q = pp.q;
    
    F = pp.c14 .*(N.^p .* Ub.^q);
    F = F .* (abs(Ub).^-1); 
    
    mu = B2./F;

elseif isequal(oo.slidinglaw, 3)              %6b of Hewitt (2012)
    n = pp.n_Glen;
    F = pp.c15 .* N .* (Ub./ (pp.c16.*Ub + pp.c17.*N.^n)).^(1/n); 
    F = F .* (abs(Ub).^-1); 
    mu = B2./F;
end

%% Reshape
mu = reshape(gg.S_h'*mu,gg.nJ,gg.nI);


% Fill regions where N (effective pressure) < eps
eps = 1e4/pp.N;
mask = dd.N<eps;
mask(gg.next) =0;

mu(mask) = NaN;
mu = inpaint_nans(mu);
mu(gg.next)=0;

aa.mu = mu;

%% Save
aa.mu = mu;

end



end
