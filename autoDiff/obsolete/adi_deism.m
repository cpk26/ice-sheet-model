function [ J ] = adi_deism(acoeff,uv_init,SS,Cb_init,F2_init,vv,aa,pp,gg,oo )


alpha = ism_alpha_field(acoeff,vv, pp, gg, oo);
Cb = ism_slidinglaw2(alpha,uv_init,Cb_init,F2_init,vv,aa,pp,gg,oo);

evFac = 1;
if oo.hybrid
evFac = (1 + (pp.c13*Cb).*F2_init);
end

C = Cb./evFac;

[nEff] = adi_visc_di(uv_init,C,aa,pp,gg,oo); %Updated Viscosity



uv = adi_deism_placeholder(uv_init,SS,C,nEff,vv,aa,pp,gg,oo );
%[uv] = adi_deism_backbone(C,nEff,vv,aa,pp,gg,oo );


F1 = adi_falpha(1,uv,C,vv,aa,pp,gg,oo );
F2 = adi_falpha(2,uv,C,vv,aa,pp,gg,oo );
Cb = ism_slidinglaw3(alpha,uv,Cb,F2,vv,aa,pp,gg,oo);


[J] = ism_inv_cost2(uv,Cb,alpha,F1,F2,vv,aa,pp,gg, oo);    



end

