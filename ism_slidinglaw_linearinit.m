function [alpha] = ism_slidinglaw_linearinit(uv,B2,N,vv,aa,pp,gg,oo)
%% Field Equations for SSA velocities 
% Inputs:
%   v     	solution struct
%   aa      prescribed fields, including inputs and boundary conditions
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs: 
%   aa     


%% Linear
if isequal(oo.slidinglaw, 1)                              
    alpha = B2;
    
else                                          
    
%% Non Linear

%Setup variables
%Setup variables
u = uv(1:gg.nua);      
v = uv(gg.nua+1:end); 
U = ((gg.c_uh*u).^2 + (gg.c_vh*v).^2).^(1/2);
  
%% Basal velocities in the case of hybrid ice sheet model
evFac = 1;
if oo.hybrid
evFac = (1 + (pp.c13*B2).*vv.F2);
end

Ub = sqrt((U./evFac).^2 + pp.U_rp.^2);          %regularize, ensure > 0


if isequal(oo.slidinglaw, 2)            %6a of Hewitt (2012)
    p = pp.p; q = pp.q;
    
    F = pp.c14 .*(N.^p .* Ub.^q);
    F = F .* (Ub.^-1); 
    
    alpha = B2./F;

elseif isequal(oo.slidinglaw, 3)              %6b of Hewitt (2012)
    n = pp.n_Glen;
    F = pp.c15 .* N .* (Ub./ (pp.c16.*Ub + pp.c17.*N.^n)).^(1/n); 
    F = F .* (Ub.^-1); 
    alpha = B2./F;
end

%% Reshape
alpha = reshape(gg.S_h'*alpha,gg.nJ,gg.nI);
N = reshape(gg.S_h'*N,gg.nJ,gg.nI);

% Fill regions where N (effective pressure) is too low for an accurate
% estimate
eps = 1.5*pp.N_rp/pp.phi;
mask = N<eps;
mask(gg.next) =0;

alpha(mask) = NaN;
alpha = inpaint_nans(alpha);

%% Reshape for output
alpha = gg.S_h*alpha(:);


end



end
