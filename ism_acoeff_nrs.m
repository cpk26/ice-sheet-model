
function [vv2] = ism_acoeff_nrs(vv,aa, pp, gg, oo)
%% Newton-Raphson step to update inversion alpha coefficients 
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

disp('Newton-Raphson Step')


%% Paramaters
vv2 = vv;               

maxarm = 15;            %Maximum number of steps (reductions)
alpha = 0.01;
iarm = 0;               %Counter of number of reductions
sigma = .25;            %Initial step size reduction
armflag = 0;            %Flag for exceeding max number of reductions

[mft0] = ism_vel_misfit(vv.u,vv.v,aa,pp,gg, oo);  %Initial velocity misfit

tau = mft0/norm(vv.agrad(:),2);  %Step coefficient.
step = vv.agrad / norm(vv.agrad(:),2);

%% Initial step
vv2.acoeff = vv.acoeff - tau*step;
vv2.C = exp(idct2(vv2.acoeff));       %Calculate new slipperiness field


[vv2] = ism_sia(aa.s,aa.h,vv2.C,vv2, pp,gg,oo); %SIA                                        
[vv2] = ism_sstream(vv2,aa,pp,gg,oo );          %SSA 

[mft] = ism_vel_misfit(vv2.u,vv2.v,aa,pp,gg, oo); %Current misfit

%% Iterate step size coefficient
while mft >= mft0 - alpha * tau * norm(vv.agrad,2);    %stopping conditions
    fprintf('Line search iteration: %i \n',iarm+1)
    %% Step Size
    tau = sigma*tau;

    %% Update vv2.coeff; keep the books on tau.
    vv2.acoeff = vv.acoeff - tau*step; %stepping from the original acoeff.

    %% Calculate misfit based on current step size
    vv2.C = exp(idct2(vv2.acoeff));      %Calculate new slipperiness field
    [vv2] = ism_sia(aa.s,aa.h,vv2.C,vv2,pp,gg,oo);  %Calculate corresponding velocities 
    [vv2] = ism_sstream(vv2,aa,pp,gg,oo );      %SSA 
    [mft] = ism_vel_misfit(vv2.u,vv2.v,aa,pp,gg, oo); %Current misfit
    
    %% Keep the books on the function norms.
    iarm = iarm+1;
    if iarm > maxarm
        disp('Newton Step Failure, too many reductions ');
        vv2.armflag = 1;
        break
    end
    
    disp(iarm)
end
fprintf('Misfit: %i \n', mft)
disp('Completed')

end   

