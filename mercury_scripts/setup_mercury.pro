pro setup_mercury,file,mstar,nplanets=nplanets,alim=alim,mlim=mlim,run=run

restore,file,/v
mercpath='/astro/grads/sjg223/research/mercury6/'
examplepath=mercpath+'grid_ex/'
au = 1.496d11
bigG = 6.67259d-11 ;mks
if not keyword_set(nplanets) then nplanets=1
if nplanets gt 2 then begin
   print,'nplanets must be 1 or 2!'
   return
endif

if not keyword_set(alim) then alim=[min(as),max(as)]
if not keyword_set(mlim) then mlim=[min(ms),max(ms)]
;;find best fit
;junk = min(mmper, w)
;wa = w mod n_elements(mmper[*, 0])
;wm = w/n_elements(mmper[*, 0])
;best=pars[wa,wm,*]
;print,reform(best)

;mind=wm+10*(indgen(5)+5)
count=0
for i=0,n_elements(as)-1 do begin
   for j=0,n_elements(ms)-1 do begin
      count=count+1
      if as(i) ge alim[0] and as(i) le alim[1] and ms(j) ge mlim[0] and ms(j) le mlim[1] then begin
;for i=0,n_elements(mind)-1 do begin
         
         mfolder=mercpath+file+'_test'+strtrim(string(count),2)+'/';+'-long/'
         print,mfolder
         spawn,'mkdir '+mfolder
         
         mp=ms[j]*0.0009546d    ;solar masses
         ap=as[i]
;other params...
         ecc=es[i,j]
         q=ap*(1.-ecc)
         inc=90.0
         om=pars[i,j,3]
         lon_an=0.0
         T0=pars[i,j,1]
         if nplanets eq 2 then begin
            p2=pars[i,j,7]
            psec = p2*24.*3600.
            k2=pars[i,j,11]     ;m/s
            ecc2=pars[i,j,9]
            mp2=k2*sqrt(1.-ecc2^2)*psec^(1./3)*(mstar*1.99d30)^(2./3)/(2*!pi*bigG)^(1./3)/1.99d30
            ap2= psec^(2./3)*(bigG*mstar*1.99d30/(4*!pi^2))^(1./3)/au   
            q2=ap2*(1.-ecc2)
            inc2=inc
            om2=pars[i,j,10]
            lon_an2=lon_an
            T02=pars[i,j,8]
         endif
         
;write big.in
         get_lun,fid
         openw,fid,mfolder+'big.in'
         printf,fid,')O+_06 Big-body initial data  (WARNING: Do not delete this line!!)'
         printf,fid,' style (Cartesian, Asteroidal, Cometary) = Cometary'
         printf,fid,' epoch (in days) = 2453299.7'
         printf,fid,' PLANET1   m='+string(mp)+' r=1 d=1'
         printf,fid,string(q)+' '+string(ecc)+' '+string(inc)
         printf,fid,string(om)+' '+string(lon_an)+' '+string(T0)
         printf,fid,'0. 0. 0.'
         if nplanets eq 2 then begin
            printf,fid,' PLANET2   m='+string(mp2)+' r=1 d=1'
            printf,fid,string(q2)+' '+string(ecc2)+' '+string(inc2)
            printf,fid,string(om2)+' '+string(lon_an2)+' '+string(T02)  
            printf,fid,'0. 0. 0.'
         endif
         close,fid
         free_lun,fid
         
;write small.in
         spawn,'cp '+examplepath+'* '+mfolder
         

;write files.in
         get_lun,fid
         openw,fid,mfolder+'files.in'
         printf,fid,'big.in'
         printf,fid,'small.in'
         ;printf,fid,'param-long2.in'
         printf,fid,'param.in'
         printf,fid,'xv.out'
         printf,fid,'ce.out'
         printf,fid,'info.out'
         printf,fid,'big.dmp'
         printf,fid,'small.dmp'
         printf,fid,'param.dmp'
         printf,fid,'restart.dmp'
         close,fid
         free_lun,fid
      
         if keyword_set(run) then begin
            cd,mfolder
            spawn,'cd '+mfolder
            spawn,'g77 -o mercury6 mercury6_2.for'
            spawn,'./mercury6'
            cd,'/astro/grads/sjg223/research/rvorbit'
         endif
      endif
   endfor
endfor

end
