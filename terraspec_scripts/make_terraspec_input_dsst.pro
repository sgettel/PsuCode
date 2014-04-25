pro make_terraspec_input_dsst,dsstin,ord,wmerge,smerge

yourpath='/Volumes/RAID/rv/vel/sara'
dsstfile=yourpath+'/files/'+dsstin
restore,dsstfile,/v
g=where(dsst.ordt eq ord)

wavl=dblarr(700,n_elements(g))
dst=dsst[g].dst

;linear approximation for each block...
for i=n_elements(g)-1,0,-1 do begin
   wavl[*,i]=dindgen(700)*dsst[g[i]].w1+replicate(dsst[g[i]].w0,700)
   if i eq n_elements(g)-1 then begin
      wmerge=wavl[*,i];fill in reddest block as is
      smerge=dst[*,i]    
   endif else begin
      x=where(wavl[*,i] lt wavl[0,i+1],count)
      wmerge=[wmerge,wavl[x,i]]
      smerge=[smerge,dst[x,i]]
   endelse
endfor

s=sort(wmerge)
wmerge=wmerge(s)
smerge=smerge(s)

outfile=yourpath+'/files/terraspec_input/'+dsstin+'_'+strtrim(string(ord),2)+'.dat'
print,outfile

if file_test(outfile) eq 0 then begin ;make input for terraspec
   wvlnvac=wmerge ;I2 wavescale in vacuum wavelengths!
   wn = (1./(wvlnvac*1.d-8))
   wn=reverse(wn)
   sp=reverse(smerge)
   openw,fid,outfile,/get_lun
   for i=0,n_elements(wn)-1 do begin
      printf,fid,wn(i),sp(i),format='(d12.6,1x,d12.6)'
   end
   free_lun,fid
endif else print,outfile+' exists'

end
