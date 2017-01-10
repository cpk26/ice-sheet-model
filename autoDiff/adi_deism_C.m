function [ C ] = adi_deism_C(alpha,uv_init,Cb_init,F2_init,vv,aa,pp,gg,oo )

Cb = ism_slidinglaw2(alpha,uv_init,Cb_init,F2_init,vv,aa,pp,gg,oo);

evFac = 1;
if oo.hybrid
evFac = (1 + (pp.c13*Cb).*F2_init);
end

C = Cb./evFac;


end

