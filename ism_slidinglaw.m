function [C] = ism_slidinglaw(vv,aa,pp,gg,oo)
%% Field Equations for SSA velocities 
% Inputs:
%   v     	solution struct
%   aa      prescribed fields, including inputs and boundary conditions
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs: 
%   C     


if ~isfield(oo,'slidinglaw'), oo.slidinglaw = 'linear'; end 


%% Linear

if strcmp(oo.slidinglaw, 'linear')                 
    if isfield(aa,'B2'), B2 = aa.B2; 
    elseif isfield(aa,'Cb'), B2 = aa.Cb;               %% Maintain backwards compatability
    else B2 = aa.C; end;
    
    C = B2;
else

%% Non Linear

%Setup variables
mu = gg.S_h*aa.mu(:);
N = max(gg.S_h*aa.N(:),0);
U = sqrt( (gg.c_uh*vv.u).^2 + (gg.c_vh*vv.v).^2 );

%% Basal velocities in the case of hybrid ice sheet model
if oo.hybrid

if strcmp(oo.pT, 'forward'),Cb = gg.S_h*aa.Cb;   %Current basal drag
else Cb = gg.S_h*vv.Cb; end  
   
F2 = ism_falpha(2,vv.U,vv.nEff_lyrs,vv,aa,pp,gg,oo ); 
tmpa = (1 + pp.c13*Cb(:).*F2);
Ub = U./tmpa;
else
Ub = U;
end

%% Handle different sliding laws
if strcmp(oo.slidinglaw, 'weertman')            %6a of Hewitt (2012)
    p = pp.p; q = pp.q;
    
    F = mu .* pp.c14 .*(N.^p .* Ub.^q);
    C = F .* (abs(Ub).^-1); 
        
elseif strcmp(oo.slidinglaw, 'schoof')              %6b of Hewitt (2012)
    n = pp.n_Glen;
    F = mu*pp.c15 .* N .* (Ub./ (pp.c16.*Ub + pp.c17.*N.^n)).^(1/n); 
    C = F .* (abs(Ub).^-1); 
    
end

%% Reshape 
C = reshape(gg.S_h'*C,gg.nJ,gg.nI);
end




end