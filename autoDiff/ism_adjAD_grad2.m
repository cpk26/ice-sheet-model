function [ gradN ] = ism_adjAD_grad2(acoeff,vv2,rr,aa,pp,gg,oo)


    %% Adjoint method with automatic differentiation
    rr.runalpha = [];
    
    numSaveIter = size(rr.uvn,2); 
    
    fi = numSaveIter;               %Forward Iteration
    ci = numSaveIter - 1;           %Current Iteration
    pi = max(ci-1,1);               %Previous Iteration
    
    %% Adjoint of final part of deism
      
    
    uv_adi = struct('f', rr.uvn(:,fi), 'duv',ones(gg.nua+gg.nva,1));
    C_adi = struct('f', rr.Cn(:,ci), 'dC',ones(gg.nha,1));
    alpha_adi = struct('f', vv2.alpha, 'dalpha',ones(gg.nha,1));
    ADOBJ = adi_deism_cost_AD(uv_adi,C_adi,alpha_adi,vv2,aa,pp,gg,oo);
    J_C = ADOBJ.dC;
    J_alpha = ADOBJ.dalpha;
    J_uv = ADOBJ.duv;
    
%     adjC_01 = zeros(gg.nha,1);
%     if oo.hybrid
%     C_adi = struct('f', rr.Cn(:,ci), 'dC',ones(gg.nha,1));
%     CAD = adi_deism_cost_ADc(rr.uvn(:,fi),C_adi,vv2.alpha,vv2,aa,pp,gg,oo);
%     adjC_01 = CAD.dC.^-1;
%     end
%     
%     alpha_adi = struct('f', vv2.alpha, 'dalpha',ones(gg.nha,1));
%     AAD = adi_deism_cost_ADa(rr.uvn(:,fi),rr.Cn(:,ci),alpha_adi,vv2,aa,pp,gg,oo);
%     adjalpha_01 = AAD.dalpha.^-1;
%     
%     
%     uv_adi = struct('f', rr.uvn(:,fi), 'duv',ones(gg.nua+gg.nva,1));
%     UAD = adi_deism_cost_ADu(uv_adi,rr.Cn(:,ci),vv2.alpha,vv2,aa,pp,gg,oo);
%     adjUV = UAD.duv.^-1;
%     rr.adjU = adjUV;  
    
    %% Determine adjoint variable nEff and C from uv
    oo.adjAD_C = 1;
    if oo.hybrid; oo.adjAD_nEff = 1;
    else oo.adjAD_nEff = 0; end
    
    rr.adjU = J_uv;
    [rr] = ism_adjAD_main2([],rr,aa,pp,gg,oo );
    
    adjnEff = rr.adjnEff;
    
    adjC_01 = J_C;
    adjC_02 = rr.adjC;
    adjC = adjC_01 + adjC_02;
    
    adjalpha_01 = J_alpha;
        
    alpha_adi = struct('f', vv2.alpha, 'dalpha',ones(gg.nha,1));
    AAD = adi_deism_C_ADa(alpha_adi,rr.uvn(:,ci),rr.Cn(:,pi),rr.F2n(:,ci),vv2,aa,pp,gg,oo);
    dalpha = sparse(AAD.dalpha_location(:,1),AAD.dalpha_location(:,2), AAD.dalpha, AAD.dalpha_size(1), AAD.dalpha_size(2));
    adjalpha_02 = dalpha'*adjC; 
    
    alpha_adi = struct('f', vv2.alpha, 'dalpha',ones(gg.nha,1));
    AAD = adi_deism_nEff_ADnEff(alpha_adi,rr.uvn(:,ci),rr.Cn(:,pi),rr.F2n(:,ci),vv2,aa,pp,gg,oo);
    dalpha = sparse(AAD.dalpha_location(:,1),AAD.dalpha_location(:,2), AAD.dalpha, AAD.dalpha_size(1), AAD.dalpha_size(2));
    adjalpha_03 = dalpha'*adjnEff;
    
    adjalpha = adjalpha_01 + adjalpha_02 + adjalpha_03;
    
    
    %% Gradient of cst function w.r.t acoeff
    gradN = (adjalpha.*exp(acoeff(:)));
    
end

