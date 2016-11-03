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
if strcmp(oo.slidinglaw, 'linear')                              
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

if strcmp(oo.slidinglaw, 'weertman')            %6a of Hewitt (2012)
    p = pp.p; q = pp.q;
    
    F = pp.c14 .*(N.^p .* Ub.^q);
    F = F .* (abs(Ub).^-1); 
    
    mu = B2./F;

elseif strcmp(oo.slidinglaw, 'schoof')              %6b of Hewitt (2012)
    n = pp.n_Glen;
    F = pp.c15 .* N .* (Ub./ (pp.c16.*Ub + pp.c17.*N.^n)).^(1/n); 
    F = F .* (abs(Ub).^-1); 
    mu = B2./F;
end

%% Reshape and save
mu = reshape(gg.S_h'*mu,gg.nJ,gg.nI);
aa.mu = mu;

end
