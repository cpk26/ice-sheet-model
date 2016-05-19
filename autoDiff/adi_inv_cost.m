

U = adigatorCreateDerivInput([11400,1], 'U');


adigator('ism_inv_cost', {U,vv,aa,pp,gg, oo},'myderiv')





%Testing Routine
% U = [u;v];
% z = gg.S_h*aa.s(:);
% z = gg.S_h*aa.b(:);
% C = gg.S_h*aa.C(:);
% ism_visc(z,U,nEff,C,aa,pp,gg,oo);