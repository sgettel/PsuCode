pro reformat_terraspec_dsst,dsstin,ord,wlt,sp0,wavl,rsp,orig=orig,ext=ext

;read original dsst for wavelength scale
yourpath='/pool/vonnegut0/mlb/rv/vel/sara'
;yourpath='/Volumes/RAID/rv/vel/sara'
dsstfile=yourpath+'/files/'+dsstin
restore,dsstfile,/v
g=where(dsst.ordt eq ord)

wavl=dblarr(700,n_elements(g))
rsp=dblarr(700,n_elements(g))
dst=dsst[g].dst

;were the telluric line strengths refit?
if keyword_set(orig) then $
  cdsstfile=yourpath+'/files/terraspec_output/'+dsstin+'_'+strtrim(string(ord),2)+'.dat.terra-restorels' else $
  cdsstfile=yourpath+'/files/terraspec_output/'+dsstin+'_'+strtrim(string(ord),2)+'.dat.terra-fitls'
;   cdsstfile=yourpath+'/files/terraspec_output/'+dsstin+'_'+strtrim(string(ord),2)+'.dat.terra0' else $
;      cdsstfile=yourpath+'/files/terraspec_output/'+dsstin+'_'+strtrim(string(ord),2)+'.dat.terra'

;if keyword_set(ext) then cdsstfile=cdsstfile+'-'+strtrim(string(ext),2)

;read telluric-cleaned dsst
readcol,cdsstfile,wvn,sp0,f='d,d'
wlt = (1./(wvn*1.d-8))
help,wlt,sp0
s=sort(wlt)
wlt=wlt(s)
sp0=sp0(s)

;linear approximation & cut into blocks
for i=0,n_elements(g)-1 do begin
   wavl[*,i]=dindgen(700)*dsst[g[i]].w1+replicate(dsst[g[i]].w0,700)
   ;rsp[*,i]=cspline(wlt,sp0,wavl[*,i])

   rsp[*,i]=interpol(sp0,wlt,wavl[*,i])
   
   

endfor

newdsst=dsst
newdsst(g).dst=rsp ;should I change the weights?

dsst=newdsst
if keyword_set(orig) then dsstfile2=dsstfile+'-restorels' else $
   dsstfile2=dsstfile+'-fitls'
if keyword_set(ext) then dsstfile2=dsstfile2+'-'+strtrim(string(ext),2)

if file_test(dsstfile2) eq 0 then begin
save,dsst,filename=dsstfile2
endif else print,dsstfile2+' exists'

;save,wavl,rsp,file='reformat_test.dat'


end
