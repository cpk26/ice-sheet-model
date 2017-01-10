function [] = ism_adjAD_generate( vv,aa, pp, gg, oo )
%Automatically differentiate necessary functions using Adigator, specifying
%input sizes from current simulation

options = adigatorOptions('OVERWRITE',1);

%Cost Function w.r.t velocity, basal drag, and alpha
uv = adigatorCreateDerivInput([gg.nua+gg.nva,1], 'uv');
C = adigatorCreateDerivInput([gg.nha,1],'C');
alpha = adigatorCreateDerivInput([gg.nha,1],'alpha');
adigator('adi_deism_cost', {uv,C,alpha,vv,aa,pp,gg,oo },'adi_deism_cost_AD',options)
clear uv C alpha;

% %Cost Function w.r.t velocity
% oo.slidinglaw = 1;
% oo.inv_cst = 1;
% uv = adigatorCreateDerivInput([gg.nua+gg.nva,1], 'uv');
% C = adigatorCreateAuxInput([gg.nha,1]);
% alpha = adigatorCreateAuxInput([gg.nha,1]);
% adigator('adi_deism_cost', {uv,C,alpha,vv,aa,pp,gg,oo },'adi_deism_cost_ADu',options)
% clear uv C alpha;
% 
% %Cost Function w.r.t basal drag
% oo.slidinglaw = 1;
% oo.inv_cst = 1;
% C = adigatorCreateDerivInput([gg.nha,1],'C');
% uv = adigatorCreateAuxInput([gg.nua+gg.nva,1]);
% alpha = adigatorCreateAuxInput([gg.nha,1]);
% adigator('adi_deism_cost', {uv,C,alpha,vv,aa,pp,gg,oo },'adi_deism_cost_ADc',options)
% clear uv C alpha;
% 
% %Cost Function w.r.t alpha
% oo.slidinglaw = 1;
% oo.inv_cst = 1;
% alpha = adigatorCreateDerivInput([gg.nha,1],'alpha');
% uv = adigatorCreateAuxInput([gg.nua+gg.nva,1]);
% C = adigatorCreateAuxInput([gg.nha,1]);
% adigator('adi_deism_cost', {uv,C,alpha,vv,aa,pp,gg,oo },'adi_deism_cost_ADa',options)
% clear uv C alpha;



%Diagonal of A matrix w.r.t Basal slipperiness
C = adigatorCreateDerivInput([gg.nha,1], 'C');
nEff = adigatorCreateAuxInput([gg.nha,1]);
adigator('ism_dim_Ddiag', {C,nEff,aa,pp,gg,oo},'ism_dim_Ddiag_ADc',options)
clear C nEff;

%Diagonal of A matrix w.r.t depth integrated viscosity
nEff = adigatorCreateDerivInput([gg.nha,1], 'nEff');
C = adigatorCreateAuxInput([gg.nha,1]);
adigator('ism_dim_Ddiag', {C,nEff,aa,pp,gg,oo},'ism_dim_Ddiag_ADnEff',options)
clear C nEff;


%nEff Function w.r.t acoeff
alpha = adigatorCreateDerivInput([gg.nha,1],'alpha');
U = adigatorCreateAuxInput([gg.nua+gg.nva,1]);
C = adigatorCreateAuxInput([gg.nha,1]);
F2 = adigatorCreateAuxInput([gg.nha,1]);
adigator('adi_deism_nEff', {alpha,U,C,F2,vv,aa,pp,gg,oo },'adi_deism_nEff_ADnEff',options)
clear U C alpha;


%C Function w.r.t acoeff
alpha = adigatorCreateDerivInput([gg.nha,1],'alpha');
U = adigatorCreateAuxInput([gg.nua+gg.nva,1]);
C = adigatorCreateAuxInput([gg.nha,1]);
F2 = adigatorCreateAuxInput([gg.nha,1]);
adigator('adi_deism_C', {alpha,U,C,F2,vv,aa,pp,gg,oo },'adi_deism_C_ADa',options)
clear U C alpha;


end

