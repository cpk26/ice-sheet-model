
function [cst,gradN] = ism_adjAD_optWrapper(acoeff,vv,aa, pp, gg, oo)
% Inputs:
%   acoeff  alpha coefficients determining basal cslip (C)
%   vv      struct containing initial solution variables
%   aa      prescribed fields, including inputs and boundary conditions
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs:
%   cst     cst function of misfit between observed and predicted velocities
%   grad    gradient of cost function w.r.t acoeff

%% Construct Basal Slip Field
vv.acoeff = reshape(acoeff,gg.nJ,gg.nI);   %Array=>matrix
vv.C = ism_cslip_field(vv, pp, gg, oo);    %Construct basal slipperiness

%% Initial velocites and viscosity

[vv] = ism_sia(aa.s,aa.h,vv.C,vv, pp,gg,oo);    %SIA

if oo.hybrid,                               %Setup Hybrid Approximation                           
Cb = vv.C;                                  %Basal Slipperiness
U = vv.U;
nEff = ism_visc(U,vv,aa,pp,gg,oo);          %Initial viscosity;

for j=[1:10]                                %Self consistent viscosity
F2 = ism_falpha(2,nEff,vv,aa,pp,gg,oo );    %Effective Basal Slipperiness
C = Cb(:)./(1 + Cb(:).*(gg.S_h'*F2));
nEff = ism_visc_di(U,nEff,gg.S_h*C(:),aa,pp,gg,oo); %Updated Viscosity
end   

vv.nEff = nEff;
end

%% DEISM

oo.savePicIter = 0;                                     %Depth Integrated Model
oo.pic_iter = 10; 
[vv, ~] = ism_deism(vv,aa,pp,gg,oo );           

if oo.hybrid
F1 = ism_falpha(1,vv.nEff,vv,aa,pp,gg,oo );          %Calculate F alpha factors [Hybrid]
F2 = ism_falpha(2,vv.nEff,vv,aa,pp,gg,oo );
else
F1 = ones(gg.nha,1);                                   %Placeholders for SSA
F2 = F1;
end

%% Cost

cst = ism_inv_cost(vv.U,gg.S_h*vv.C(:),F1,F2,vv,aa,pp,gg, oo);  %Current misfit

%% Return gradient if optimization routine requests it

if nargout > 1 % gradient required
    
    %% Depth Integrated Model
    %Second DEISM call, save parameters for adjoint method
    oo.savePicIter = 1;                                     %Save picard iterations
    oo.pic_iter = 1;                                        %Reduced number of picard iteration
    [vv, rr] = ism_deism(vv,aa,pp,gg,oo );   
    
    %% Adjoint method with automatic differentiation
    % Adjoint variable Uf*
    U_adi = struct('f', vv.U, 'dU',ones(gg.nua+gg.nva,1));
    UAD = ism_inv_cost_ADu(U_adi,gg.S_h*vv.C(:),F1,F2,vv,aa,pp,gg,oo);
    rr.adjU = UAD.dU;
    
    % For the hybrid approximation, determine adjoint state of Cf* 
    if oo.hybrid
    C_adi = struct('f', gg.S_h*vv.C(:), 'dC',ones(gg.nha,1));
    UAC = ism_inv_cost_ADc(vv.U,C_adi,F1,F2,vv,aa,pp,gg,oo);
    rr.adjC = UAC.dC;
    end
    
    %Determine adjoint variable C*
    rr = ism_adjAD_main(vv,rr,aa,pp,gg,oo );
    
    %Gradient of cst function w.r.t acoeff
    
    gradN = rr.runC.*(exp(acoeff(:)));
    gradN = gradN/max(abs(gradN));
    %gradN = grad./max(abs(grad));
    %gradN = grad./max(abs(grad));
    %imagesc(reshape(gradN,gg.nJ,gg.nI))
    
    %acoeff_orig = acoeff;
    %cst_orig = cst;
    %acoeff = acoeff_orig - .1*gradN;

end

end   

