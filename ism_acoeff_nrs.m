
function [vv2] = ism_acoeff_linesearch(vv,aa, pp, gg, oo)
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

maxarm = 25;            %Maximum number of steps (reductions)
iarm = 0;               %Counter of number of reductions
sigma = .5;            %Initial step size reduction
alpha = 1.d-4;          %Line search parameter
armflag = 0;            %Flag for exceeding max number of reductions
n_x = vv.n_x;           %Number of coefficients for basis in x/y dirs
n_y = vv.n_y;

[mft0] = ism_vel_misfit(vv.u,vv.v,aa,pp,gg, oo);  %Initial velocity misfit

tau = vv.tau; taum = 1; tauc = tau;     %Step coefficient. 'm': previous step; 'c': current step
step = mft0/norm(vv.agrad(:),1) * vv.agrad / norm(vv.agrad(:),1);


%% Initial step
vv2.acoeff = vv.acoeff - tau*step;
vv2.C = idct2(vv2.acoeff);       %Calculate new slipperiness field

[ii] = ism_sia(aa.s,aa.h,vv2.C, pp,gg,oo);  %Calculate corresponding velocities 
vv2.u = ii.u; vv2.v = ii.v;                 %SIA                                        
[vv2] = ism_sstream(vv2,aa,pp,gg,oo );      %SSA 

[mft] = ism_vel_misfit(vv2.u,vv2.v,aa,pp,gg, oo); %Current misfit
mft0_2 = mft0*mft0; mftc_2 = mft*mft; mftm_2 = mft*mft; %Misfit squared; 'm': previous step; 'c': current step


%% Iterate step size coefficient
%mft >= (1 - alpha*tau) * mft0

while mft >=  mft0;    %stopping conditions
    fprintf('Line search iteration: %i \n',iarm+1)
    %% Apply the three point parabolic model.
    tau = sigma*tau;
%     if iarm == 0
%         tau = sigma1*tau;
%     else
%         tau = parab3p(tauc, taum, mft0_2, mftc_2, mftm_2);
%     end

    %% Update vv2.coeff; keep the books on tau.
    step = mft0/norm(vv.agrad(:),1) * vv.agrad / norm(vv.agrad(:),1);
    vv2.acoeff = vv.acoeff - tau*step; %stepping from the original acoeff.
    taum = tauc;
    tauc = tau;

    %% Calculate misfit based on current step size
    vv2.C = idct2(vv2.acoeff).^2;      %Calculate new slipperiness field
    [ii] = ism_sia(aa.s,aa.h,vv2.C, pp,gg,oo);  %Calculate corresponding velocities 
    vv2.u = ii.u; vv2.v = ii.v;                 %SIA                                        
    [vv2] = ism_sstream(vv2,aa,pp,gg,oo );      %SSA 
    [mft] = ism_vel_misfit(vv2.u,vv2.v,aa,pp,gg, oo); %Current misfit
    
    %% Keep the books on the function norms.
    mftm_2 = mftc_2;
    mftc_2 = mft*mft;
    iarm = iarm+1;
    if iarm > maxarm
        disp(' Armijo failure, too many reductions ');
        armflag = 1;
        return;
    end
    
    disp(iarm)
end
fprintf('Misfit: %i \n', mft)
disp('Completed')
vv2.tau = tau; %Save tau as initial starting point for subsequent line search

end   

