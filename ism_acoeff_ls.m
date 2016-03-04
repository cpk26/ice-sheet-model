
function [vv2] = ism_acoeff_ls(vv,aa, pp, gg, oo)
%% Line search to update inversion alpha coefficients 
% Based on code from C.T. Kelley, from SIAM's 'Solving Nonlinear Equations
% with Newton's Method'
% Inputs:
%   vv      struct containing initial solution variables
%   aa      prescribed fields, including inputs and boundary conditions
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs:
%   vv2     updated struct with new alpha coefficients

disp('Line Search Step')


%% Paramaters
vv2 = vv;               

maxarm = 15;            %Maximum number of steps (reductions)
alpha = 0.001;
iarm = 0;               %Counter of number of reductions
sigma = .25;            %Initial step size reduction
armflag = 0;            %Flag for exceeding max number of reductions

[mft0] = ism_inv_misfit(vv,aa,pp,gg, oo);  %Initial velocity misfit

%% calculate gradient, determine step size
[vv] = ism_adjoint(vv,aa,pp,gg,oo );
[vv] = ism_cslip_grad(vv, pp, gg, oo);


tau = mft0/norm(vv.agrad(:),2);  %Step coefficient.
tau = 1;
sdir = -vv.agrad / norm(vv.agrad(:),2);


%% Initial step
vv2.acoeff = vv.acoeff + tau*sdir;
vv2.C = ism_cslip_field(vv2, pp, gg, oo);       %Calculate new slipperiness field


[vv2] = ism_sia(aa.s,aa.h,vv2.C,vv2, pp,gg,oo); %SIA                                        
[vv2] = ism_sstream(vv2,aa,pp,gg,oo );          %SSA 

[mft] = ism_inv_misfit(vv2,aa,pp,gg, oo); %Current misfit

%% Iterate step size coefficient
while mft >= mft0 + alpha * tau * (vv.agrad(:)'*sdir(:));    %stopping conditions
    fprintf('Line search iteration: %i \n',iarm+1)
    %% Step Size
    tau = sigma*tau;

    %% Update vv2.coeff; keep the books on tau.
    vv2.acoeff = vv.acoeff + tau*sdir; %stepping from the original acoeff.
    %% Calculate misfit based on current step size
    vv2.C = ism_cslip_field(vv2, pp, gg, oo);      %Calculate new slipperiness field
    [vv2] = ism_sia(aa.s,aa.h,vv2.C,vv2,pp,gg,oo);  %Calculate corresponding velocities 
    [vv2] = ism_sstream(vv2,aa,pp,gg,oo );      %SSA 
    [mft] = ism_inv_misfit(vv2,aa,pp,gg, oo); %Current misfit
    
    %% Keep the books on the function norms.
    iarm = iarm+1;
    if iarm > maxarm
        disp('Line Search Failure, too many reductions ');
        vv2.armflag = 1;
        break
    end
    
    disp(iarm)
end
fprintf('Misfit: %i \n', mft)
disp('Completed')

end   

