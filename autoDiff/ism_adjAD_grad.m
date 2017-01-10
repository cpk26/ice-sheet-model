function [ gradN ] = ism_adjAD_grad(acoeff,vv2,rr,aa,pp,gg,oo)


   
    C = vv2.Cb;  
    F1 = ones(gg.nha,1);
    F2 = ones(gg.nha,1);
    
    if oo.hybrid, 
    F1 = vv2.F1;       %Update F1/F2 w/ new velocities/viscosities [Hybrid]
    F2 = vv2.F2;
    end
    
    %% Adjoint method with automatic differentiation
    rr.runalpha = [];
    
    
    %% Adjoint of Cost Function
    
    %% Determine adjoint state of Cb_final* 

    C_adi = struct('f', C, 'dC',ones(gg.nha,1));
    UAC = ism_inv_cost_ADc(vv2.uv,C_adi,vv2.alpha,F1,F2,vv2,aa,pp,gg,oo);
    if isfield(UAC, 'dC'), adjC = UAC.dC;           %Hybrid vs SSA
    else adjC = zeros(gg.nha,1); end;
               
    
    %% Determine adjoint state of Alpha_final* 
    %Part A
    alpha_adi = struct('f', vv2.alpha, 'dalpha',ones(gg.nha,1));
    AAD = ism_slidinglaw_ADa(alpha_adi,vv2.uv,rr.Cbn(:,end-1),rr.F2n(:,end),vv2,aa,pp,gg,oo);
    C_alpha = sparse(AAD.dalpha_location(:,1),AAD.dalpha_location(:,2), AAD.dalpha, AAD.dalpha_size(1), AAD.dalpha_size(2));
    adjalpha1a = C_alpha'*adjC;
    
    %Part B
    alpha_adi = struct('f', vv2.alpha, 'dalpha',ones(gg.nha,1));
    AAD = ism_inv_cost_ADa(vv2.uv,C,alpha_adi,F1,F2,vv2,aa,pp,gg,oo);
    adjalpha1b = AAD.dalpha;
    
    %Combine
    adjalphaF =  adjalpha1a + adjalpha1b;
    
    
    %% Determine adjoint state of Uf* 
    %Part A
    if ~strcmp(oo.slidinglaw, 'linear')
    U_adi = struct('f', vv2.uv, 'dU',ones(gg.nua+gg.nva,1));
    UAD = ism_slidinglaw_ADu(vv2.alpha,U_adi,rr.Cbn(:,end-1),rr.F2n(:,end),vv2,aa,pp,gg,oo);
    C_U = sparse(UAD.dU_location(:,1),UAD.dU_location(:,2), UAD.dU, UAD.dU_size(1), UAD.dU_size(2));
    adjU_a = C_U'*adjC;
    else adjU_a=0; end
    %adjU_a = 0;
    
    %Part B
    U_adi = struct('f', vv2.uv, 'dU',ones(gg.nua+gg.nva,1));
    UAD = ism_inv_cost_ADu(U_adi,C,vv2.alpha,F1,F2,vv2,aa,pp,gg,oo);
    adjU_b = UAD.dU;
    
    %Combine
    adjUf = adjU_a + adjU_b;

    %% Construct AFPI (Goldberg, 2016) if option set, otherwise use adjUf
    w = adjUf;
    nIter = oo.adjAD_afpi_iter;
    
    if oo.adjAD_afpi                %AFPI
    
    %Set options
    oo.adjAD_alpha = 1;
    oo.adjAD_uv = 1;
    wNorm = zeros(1,nIter+1);
    wNorm(1) = norm(adjUf);
    
    for j=1:nIter
    rr.adjU = w;
    tic
    rr = ism_adjAD_main([],rr,aa,pp,gg,oo );    
    disp(['Elapsed Time (adjU): ', num2str(toc)])
    w = rr.adjU + adjUf;
    
    wNorm(j+1) = norm(w);
    end
    
    end 
    
    rr.adjU = w;  
    
    %% Determine adjoint variable Alpha_initial*
    oo.adjAD_alpha = 1;
    if oo.hybrid; oo.adjAD_uv = 1;
    else oo.adjAD_uv = 0; end
    
    tic
    rr = ism_adjAD_main([],rr,aa,pp,gg,oo );
    disp(['Elapsed Time (adjalpha): ', num2str(toc)])
    adjalphaI = rr.runalpha;    
    
    %% Gradient of cst function w.r.t acoeff
    adjalpha = adjalphaF + adjalphaI;
    gradN = (adjalpha.*exp(acoeff(:)));
    
end

