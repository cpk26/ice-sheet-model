
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

if ~isfield(vv,'tau'), vv.tau = 1; end;        %Initial step coefficient

%% Paramaters
vv2 = vv;               

maxarm = 15;            %Maximum number of steps (reductions)
iarm = 0;               %Counter of number of reductions
sigma = .5;            %Initial step size reduction
armflag = 0;            %Flag for exceeding max number of reductions

[mft0] = ism_vel_misfit(vv.u,vv.v,aa,pp,gg, oo);  %Initial velocity misfit

tau = vv.tau;  %Step coefficient.
%step = mft0 * vv.agrad/norm(vv.agrad(:),2);
step = 0.1*mft0 * vv.agrad / norm(vv.agrad(:),2)^2;

%% Initial step
vv2.acoeff = vv.acoeff - tau*step;
vv2.C = exp(idct2(vv2.acoeff));       %Calculate new slipperiness field


[vv2] = ism_sia(aa.s,aa.h,vv2.C,vv2, pp,gg,oo); %SIA                                        
[vv2] = ism_sstream(vv2,aa,pp,gg,oo );          %SSA 

[mft] = ism_vel_misfit(vv2.u,vv2.v,aa,pp,gg, oo); %Current misfit

%% Iterate step size coefficient
%mft >= (1 - alpha*tau) * mft0

while mft >=  mft0;    %stopping conditions
    fprintf('Line search iteration: %i \n',iarm+1)
    %% Step Size
    tau = sigma*tau;

    %% Update vv2.coeff; keep the books on tau.
    %step = mft0 * vv.agrad/norm(vv.agrad(:),2);
    step = 0.1 * mft0 * vv.agrad / norm(vv.agrad(:),2)^2;
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
vv2.tau = tau; %Save tau as initial starting point for subsequent line search

end   

