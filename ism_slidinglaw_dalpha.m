function [dalpha] = ism_slidinglaw_dalpha(alpha,uv,Cb,F2,vv,aa,pp,gg,oo)
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
dalpha = ones(gg.nha,1);   


else

%% Non Linear

%Setup variables
u = uv(1:gg.nua);      
v = uv(gg.nua+1:end); 

N = max(aa.N,0);
U = ( (gg.c_uh*u).^2 + (gg.c_vh*v).^2 ).^(1/2);

%% Basal velocities in the case of hybrid ice sheet model
if oo.hybrid
tmpa = (1 + pp.c13*Cb(:).*F2);
Ub = U./tmpa;
else
Ub = U;
end


% %% Handle different sliding laws
if strcmp(oo.slidinglaw,'weertman')               %6a of Hewitt (2012)
    p = pp.p; 
    q = pp.q;
    
    F = pp.c14 .*(N.^p .* Ub.^q);
    dalpha = F .* (abs(Ub).^-1); 

        
elseif strcmp(oo.slidinglaw,'schoof')              %6b of Hewitt (2012)
    n = pp.n_Glen;
    F = pp.c15 .* N .* (Ub./ (pp.c16.*Ub + pp.c17.*N.^n)).^(1/n); 
    dalpha = F .* (abs(Ub).^-1); 
    
end

end




end