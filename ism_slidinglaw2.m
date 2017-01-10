function [Cb2] = ism_slidinglaw2(alpha,uv,Cb,F2,vv,aa,pp,gg,oo)
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
if isequal(oo.slidinglaw,1)  
Cb2 = 1.*alpha;    %1 is necessary for AD.

%% Non Linear
else
%Setup variables
u = uv(1:gg.nua);      
v = uv(gg.nua+1:end); 
U = ((gg.c_uh*u).^2 + (gg.c_vh*v).^2).^(1/2);

%% Basal velocities in the case of hybrid ice sheet model
evFac = 1;
if oo.hybrid
evFac = (1 + (pp.c13*Cb).*F2);
end

Ub = sqrt((U./evFac).^2 + pp.U_rp.^2);          %regularize


% %% Handle different sliding laws
if isequal(oo.slidinglaw,2)                 %6a of Hewitt (2012)
    p = pp.p; 
    q = pp.q;
    
    AA = alpha.*(aa.N.^p .* Ub.^q);    
    BB = pp.c14 ./ Ub;
    Cb2 = AA.*BB;
    
elseif isequal(oo.slidinglaw,3)                %6b of Hewitt (2012)
    n = pp.n_Glen;
    AA =  aa.N .* (Ub./(pp.c16.*Ub + pp.c17.*aa.N.^n)).^(1/n); 
    BB = (alpha.*pp.c15) ./ Ub;
    Cb2 = AA.*BB;
        
end


end



end