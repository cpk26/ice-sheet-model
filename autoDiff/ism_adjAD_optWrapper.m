
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
if oo.hybrid, 
    vv.Cb = ism_cslip_field(vv, pp, gg, oo);    %Construct basal slipperiness
    C = vv.Cb;
else
    vv.C = ism_cslip_field(vv, pp, gg, oo); 
    C = vv.C;
end

%% Initial velocites and viscosity

[vv] = ism_sia(aa.s,aa.h,C,vv, pp,gg,oo);    %SIA

if oo.hybrid,                               %Setup Hybrid Approximation                           
Cb = vv.Cb;                                  %Basal Slipperiness
U = vv.U;
nEff = ism_visc(U,vv,aa,pp,gg,oo);          %Initial viscosity;
nEff_lyrs = repmat(nEff,1,oo.nl+1);

F2 = ism_falpha(2,U,nEff_lyrs,vv,aa,pp,gg,oo );
C = Cb(:)./(1 + (pp.c13*Cb(:)).*(gg.S_h'*F2));                   %Effective Basal Slipperiness

vv.C = C;
vv.F2 = F2;
vv.nEff = nEff;
vv.nEff_lyrs = nEff_lyrs;

else
U = vv.U;
nEff = ism_visc(U,vv,aa,pp,gg,oo);
vv.nEff = nEff;

end



%% DEISM

oo.savePicIter = 0;                                     %Depth Integrated Model
oo.pic_iter = 15; 
[vv, ~] = ism_deism(vv,aa,pp,gg,oo );           

if oo.hybrid
F1 = ism_falpha(1,vv.U,vv.nEff_lyrs,vv,aa,pp,gg,oo );          %Calculate F alpha factors [Hybrid]
F2 = ism_falpha(2,vv.U,vv.nEff_lyrs,vv,aa,pp,gg,oo );
cst = ism_inv_cost(vv.U,gg.S_h*vv.Cb(:),F1,F2,vv,aa,pp,gg, oo);  %Current misfit

else
cst = ism_inv_cost(vv.U,gg.S_h*vv.C(:),[],[],vv,aa,pp,gg, oo);  %Current misfit
end

%% Cost


%% Return gradient if optimization routine requests it

if nargout > 1 % gradient required
    
    %% Depth Integrated Model
    %Second DEISM call, save parameters for adjoint method
    oo.savePicIter = 1;                                     %Save picard iterations
    oo.pic_iter = 1;                                        %Reduced number of picard iteration
    [vv, rr] = ism_deism(vv,aa,pp,gg,oo ); 
    
    if oo.hybrid, 
    C = vv.Cb(:);                                                  %Use the correct Basal Drag
    F1 = ism_falpha(1,vv.U,vv.nEff_lyrs,vv,aa,pp,gg,oo );          %Update F1/F2 w/ new velocities/viscosities [Hybrid]
    F2 = ism_falpha(2,vv.U,vv.nEff_lyrs,vv,aa,pp,gg,oo );
    else C = vv.C(:); F1 = zeros(gg.nha,1); F2 = F1; end;
    
    %% Adjoint method with automatic differentiation
    % Adjoint variable Uf*
    U_adi = struct('f', vv.U, 'dU',ones(gg.nua+gg.nva,1));
    UAD = ism_inv_cost_ADu(U_adi,gg.S_h*C,F1,F2,vv,aa,pp,gg,oo);
    rr.adjU = UAD.dU;
    
    %Determine adjoint state of Cf* 
    C_adi = struct('f', gg.S_h*C, 'dC',ones(gg.nha,1));
    UAC = ism_inv_cost_ADc(vv.U,C_adi,F1,F2,vv,aa,pp,gg,oo);
    rr.adjC = UAC.dC;
    rr.runC = rr.adjC;
    
    %Determine adjoint variable C*
    rr = ism_adjAD_main(vv,rr,aa,pp,gg,oo );
    
    %Gradient of cst function w.r.t acoeff
    gradN = rr.runC.*(exp(acoeff(:)));
    
    %% For Manual Testing
    %gradN = gradN/max(abs(gradN));
    %imagesc(reshape(gradN,gg.nJ,gg.nI))
    %acoeff_orig = acoeff;
    %cst_orig = cst;
    %acoeff = acoeff_orig - .1*gradN;

end

end   

