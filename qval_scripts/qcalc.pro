function qcalc,lambda,sp,sigd=sigd
;based on Bouchy et al. 2001

if not keyword_set(sigd) then sigd = 0.

dv=deriv(lambda,sp)

w=lambda^2*dv^2/(sp+sigd^2)

qvalue=sqrt(total(w))/sqrt(total(sp))

return,qvalue
end
