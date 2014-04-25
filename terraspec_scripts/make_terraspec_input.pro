pro make_terraspec_input,obnm,ordob,twochip=twochip

yourpath='/pool/vonnegut0/mlb/rv/vel/sara'
;yourpath='/Volumes/RAID/rv/vel/sara'

if keyword_set(twochip) then begin
   obs_r = readfits('/pool/vonnegut0/mlb/rv/HET/reduced_sxw/ready/'+obnm+'.fits.gz', ext=1, /sil)
   obs_r=obs_r[*,0:18]
   obs_b = readfits('/pool/vonnegut0/mlb/rv/HET/reduced_sxw/ready/'+obnm+'.fits.gz', ext=2, /sil)
   obs_b=obs_b[*,0:32]
   cont_r = readfits('/pool/vonnegut0/mlb/rv/HET/reduced_sxw/ready/cont/'+obnm+'.cont.fits.gz', ext=1, /sil)
   cont_r=cont_r[*,0:18]
   cont_b = readfits('/pool/vonnegut0/mlb/rv/HET/reduced_sxw/ready/cont/'+obnm+'.cont.fits.gz', ext=2, /sil)
   cont_b=cont_b[*,0:32]
help,obs_r,obs_b
obs=transpose([transpose(obs_r),transpose(obs_b)])
cont=transpose([transpose(cont_r),transpose(cont_b)])

endif else begin ;assume blue chip only

   obs = readfits('/pool/vonnegut0/mlb/rv/HET/reduced_sxw/ready/'+obnm+'.fits.gz', ext=2, /sil)
   cont = readfits('/pool/vonnegut0/mlb/rv/HET/reduced_sxw/ready/cont/'+obnm+'.cont.fits.gz', ext=2, /sil)
endelse

ob = reform(obs[*,ordob])

tharpath = '/pool/vonnegut0/mlb/rv/vel/sara/files/iraf_wavecal/'

if keyword_set(twochip) then begin
   tharfile_r = tharpath+obnm+'_1.sp.wl'
   restore,tharfile_r,/v
   lambda_r=lambda[*,0:18] ;double check this...
   lambda_r=reverse(lambda_r,2) ;reverse IRAF orders to match REDUCE
   tharfile_b = tharpath+obnm+'_2.sp.wl'
   restore,tharfile_b,/v
   lambda_b=lambda[*,0:32]  ;make sure array has the right size
   lambda_b=reverse(lambda_b,2)
   lambda=transpose([transpose(lambda_r),transpose(lambda_b)])
endif else begin
   tharfile = tharpath+obnm+'_2.sp.wl'
   restore,tharfile,/v
   lambda=lambda[*,0:32]    
   lambda=reverse(lambda,2)
endelse
 
wli=lambda[*,ordob]

smcont=top(cont[*,ordob],100) ;smooth errors from mask
nob=ob/smcont ;normalize order
bad=where(nob gt 1.5 or nob le 0,nb)
if nb gt 0 then nob(bad)=median(nob(0:400))

cont2=top(nob,1.d5)          ;take out weird slope on blue end of order?
nob=nob/cont2

if keyword_set(twochip) then $
   tspecfilein=yourpath+'/files/terraspec_input/'+obnm+'_'+strtrim(string(ordob),2)+'.twochip.dat' else $ 
      tspecfilein=yourpath+'/files/terraspec_input/'+obnm+'_'+strtrim(string(ordob),2)+'.dat'

;print,file_test(tspecfilein)
if file_test(tspecfilein) eq 0 then begin ;make input for terraspec
   wvlnvac=wli
   airtovac,wvlnvac
   wn = (1./(wvlnvac*1.d-8))
   wn=reverse(wn)
   nob=reverse(nob)
   openw,fid,tspecfilein,/get_lun
   for i=0,n_elements(ob)-1 do begin
      printf,fid,wn(i),nob(i),format='(d12.6,1x,d12.6)'
   end
   free_lun,fid
endif else print,tspecfilein+' exists'

end
