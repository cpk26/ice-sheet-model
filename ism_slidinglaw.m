function [Cb2] = ism_slidinglaw(alpha,uv,Cb,F2,vv,aa,pp,gg,oo)
% Inputs:
%   U     velocities
%   Cb      basal drag
%   vv     	solution struct
%   aa      prescribed fields, including inputs and boundary conditions
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs: 
%   Cb2     updated basal drag    


%if ~isfield(oo,'slidinglaw'), oo.slidinglaw = 'linear'; end 


%% Linear
if strcmp(oo.slidinglaw,'linear')  
Cb2 = 1*alpha;    %1 is necessary for AD.

%% Non Linear
else

%Setup variables
u = uv(1:gg.nua);      
v = uv(gg.nua+1:end); 
U = ((gg.c_uh*u).^2 + (gg.c_vh*v).^2).^(1/2);


N = max(aa.N,0);
N = sqrt(N.^2 + pp.N_rp.^2);

%% Basal velocities in the case of hybrid ice sheet model
if oo.hybrid
tmpa = (1 + (pp.c13*Cb(:)).*F2);
Ub = U./tmpa;
else
tmpa = 1;
Ub = U./tmpa;
end

Ub = sqrt(Ub.^2 + pp.U_rp.^2);          %regularize, ensure > 0


% %% Handle different sliding laws
if strcmp(oo.slidinglaw,'weertman')               %6a of Hewitt (2012)
    p = pp.p; 
    q = pp.q;
    
    tmpa = alpha .*(N.^p .* Ub.^q);    
    tmpb = pp.c14 ./ Ub;
    Cb2 = tmpa.*tmpb;
    
elseif strcmp(oo.slidinglaw,'schoof')              %6b of Hewitt (2012)
    n = pp.n_Glen;
    tmpa = alpha .* N .* (Ub./(pp.c16.*Ub + pp.c17.*N.^n)).^(1/n); 
    tmpb = pp.c15 ./ Ub;
    Cb2 = tmpa.*tmpb;
        
end

end




end