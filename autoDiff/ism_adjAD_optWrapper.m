
function [cst,gradN] = ism_adjAD_optWrapper(acoeff,vv,aa, pp, gg, oo)
% Inputs:
%   vv      struct containing initial solution variables
%   aa      prescribed fields, including inputs and boundary conditions
%   pp      parameters
%   gg      grid and operators
%   oo      options
% Outputs:
%   vv2     updated struct with new alpha coefficients


vv.acoeff = reshape(acoeff,gg.nJ,gg.nI);   %Array=>matrix
vv.C = ism_cslip_field(vv, pp, gg, oo);    %Construct basal slipperiness

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

oo.savePicIter = 0;                                     %Depth Integrated Model
oo.pic_iter = 10; 
[vv, ~] = ism_deism(vv,aa,pp,gg,oo );           

F1 = ism_falpha(1,vv.nEff,vv,aa,pp,gg,oo );          %Calculate F alpha factors 
F2 = ism_falpha(2,vv.nEff,vv,aa,pp,gg,oo );

cst = ism_inv_cost(vv.U,gg.S_h*vv.C(:),F1,F2,vv,aa,pp,gg, oo);  %Current misfit

if nargout > 1 % gradient required
    
    %Depth Integrated Model
    %Second ism_deism call, such that our initial guess no longer is the
    %SIA
    oo.savePicIter = 1;                                     %Save picard iterations
    oo.pic_iter = 1;                                        %Reduced number of picard iteration
    [vv, rr] = ism_deism(vv,aa,pp,gg,oo );   
    
    %Call ism_AD_inv_Cst
    U_adi = struct('f', vv.U, 'dU',ones(gg.nua+gg.nva,1));
    UAD = ism_inv_cost_ADu(U_adi,gg.S_h*vv.C(:),F1,F2,vv,aa,pp,gg,oo);
    rr.adjU = UAD.dU;
    
    %imagesc(reshape(rr.adjU(1:gg.nua),75,76))
    if oo.hybrid
    C_adi = struct('f', gg.S_h*vv.C(:), 'dC',ones(gg.nha,1));
    UAC = ism_inv_cost_ADc(vv.U,C_adi,F1,F2,vv,aa,pp,gg,oo);
    rr.adjC = UAC.dC;
    end
    %imagesc(reshape(gg.S_h'*rr.adjC(:),gg.nJ,gg.nI))
    rr = ism_adjAD_main(vv,rr,aa,pp,gg,oo );
    
    grad = rr.runC.*(exp(acoeff(:)));
    gradN = grad./max(abs(grad));
    %gradN = grad./max(abs(grad));
    %imagesc(reshape(gradN,gg.nJ,gg.nI))
    
    %acoeff_orig = acoeff;
    %cst_orig = cst;
    %acoeff = acoeff_orig - .1*gradN;

end

end   

