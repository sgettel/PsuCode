pro qvaluemodel,warr,sarr,nsp,wlist,qlist,type=type

;Calculate qvalues for model calibration and stellar spectra
;keep everything in vacuum wavelengths!

;this one does all the calculations and outputs plots

yourpath = '/Volumes/RAID/rv/vel/sara'
files = yourpath+'/files/'

case type of 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Iodine
  
   'iod': begin
   ;read FTS spectrum
      fts_atlas='ftshet70.dat'
      iodfile = strtrim(fts_atlas,2)
      rdfts,warr,sarr,5000,9000,dfn=iodfile
      help,warr,sarr
      
   ;convolve with rotation kernel
      lsf = lsf_rotate(0.05,1,velgrid=vel) ;1 Km/s, sampled every 50 m/s
      num_conv,sarr,lsf,nsp

   ;break into blocks...
      wmin=5000.
      wlist=fltarr(2000)
      qlist=fltarr(2000)
      for i=0,1999 do begin
         wmax=wmin+2.
         b=where(warr ge wmin and warr lt wmax,count)
         if i eq 0 then print,warr(b(0))/(warr(b(1))-warr(b(0)))
         if i eq 1999 then print,warr(b(0))/(warr(b(1))-warr(b(0)))
         if count gt 0 then begin
            ;qvalue=qcalc(warr(b),sarr(b))     
            qvalue=qcalc(warr(b),nsp(b)) 
            qlist(i)=qvalue
         endif else qlist(i)=0.
         wlist(i)=wmin
         wmin=wmin+2.
      endfor

      set_plot,'ps'
      device,file='qiod.ps'
      plot,wlist,qlist,xtitle='Wavelength (A)',ytitle='Q Value',title='Iodine',xrange=[5000,6000],charsize=1.5,xthick=3,ythick=3,charthick=3,thick=3
      device,/close
      set_plot,'x'
   end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Telluric

   'tel': begin
      ;read model spectrum
      readcol,'int1-1.0.terra',wn1,fl1
      readcol,'int2-1.0.terra',wn2,fl2
      readcol,'int3-1.0.terra',wn3,fl3
      readcol,'int4-1.0.terra',wn4,fl4
      readcol,'int5-1.0.terra',wn5,fl5
      readcol,'int6-1.0.terra',wn6,fl6
      
      ;trim off duplicates
      keep1=where(wn1 ge 11000. and wn1 lt 12500.)
      wn1=wn1(keep1)
      fl1=fl1(keep1)
      keep2=where(wn2 ge 12500. and wn2 lt 14000.)
      wn2=wn2(keep2)
      fl2=fl2(keep2)
      keep3=where(wn3 ge 14000. and wn3 lt 15500.)
      wn3=wn3(keep3)
      fl3=fl3(keep3)
      keep4=where(wn4 ge 15500. and wn4 lt 17000.)
      wn4=wn4(keep4)
      fl4=fl4(keep4)
      keep5=where(wn5 ge 17000. and wn5 lt 18500.)
      wn5=wn5(keep5)
      fl5=fl5(keep5)
      keep6=where(wn6 ge 18500. and wn6 le 20000.)
      wn6=wn6(keep6)
      fl6=fl6(keep6)

      wn=[wn1,wn2,wn3,wn4,wn5,wn6]
      fl=[fl1,fl2,fl3,fl4,fl5,fl6]

      wl=1/wn*1d8
      warr1=reverse(wl)
      sarr1=reverse(fl)

      readcol,'int1-0.5.terra',wn1,fl1
      readcol,'int2-0.5.terra',wn2,fl2
      readcol,'int3-0.5.terra',wn3,fl3
      readcol,'int4-0.5.terra',wn4,fl4
      readcol,'int5-0.5.terra',wn5,fl5
      readcol,'int6-0.5.terra',wn6,fl6

      ;trim off duplicates
      keep1=where(wn1 ge 11000. and wn1 lt 12500.)
      wn1=wn1(keep1)
      fl1=fl1(keep1)
      keep2=where(wn2 ge 12500. and wn2 lt 14000.)
      wn2=wn2(keep2)
      fl2=fl2(keep2)
      keep3=where(wn3 ge 14000. and wn3 lt 15500.)
      wn3=wn3(keep3)
      fl3=fl3(keep3)
      keep4=where(wn4 ge 15500. and wn4 lt 17000.)
      wn4=wn4(keep4)
      fl4=fl4(keep4)
      keep5=where(wn5 ge 17000. and wn5 lt 18500.)
      wn5=wn5(keep5)
      fl5=fl5(keep5)
      keep6=where(wn6 ge 18500. and wn6 le 20000.)
      wn6=wn6(keep6)
      fl6=fl6(keep6)

      wn=[wn1,wn2,wn3,wn4,wn5,wn6]
      fl=[fl1,fl2,fl3,fl4,fl5,fl6]

      wl=1/wn*1d8
      warr2=reverse(wl)
      sarr2=reverse(fl)
      
      ;convolve with rotation kernel
      lsf = lsf_rotate(0.05,1,velgrid=vel) ;1 Km/s, sampled every 50 m/s
      num_conv,sarr1,lsf,nsp1
      num_conv,sarr2,lsf,nsp2

      ;break into blocks...
      wmin=5000.
      wlist1=fltarr(2000)
      qlist1=fltarr(2000)
      wlist2=fltarr(2000)
      qlist2=fltarr(2000)
      for i=0,1999 do begin
         wmax=wmin+2.

         b=where(warr1 ge wmin and warr1 lt wmax,count)
         if count gt 0 then begin
            ;qvalue=qcalc(warr1(b),sarr1(b))      
            qvalue=qcalc(warr1(b),nsp1(b)) 
            qlist1(i)=qvalue
         endif else qlist1(i)=0.
         wlist1(i)=wmin

         b=where(warr2 ge wmin and warr2 lt wmax,count)
         if count gt 0 then begin
            ;qvalue=qcalc(warr2(b),sarr2(b))     
            qvalue=qcalc(warr2(b),nsp2(b)) 
            qlist2(i)=qvalue
         endif else qlist2(i)=0.
         wlist2(i)=wmin 

         wmin=wmin+2.
      endfor

      set_plot,'ps'
      device,file='qtel1.ps'
      plot,wlist1,qlist1,xtitle='Wavelength (A)',ytitle='Q Value',title='Tellurics 1.0',charsize=1.5,xthick=3,ythick=3,charthick=3,thick=3
      device,/close
      device,file='qtel2.ps'
      plot,wlist2,qlist2,xtitle='Wavelength (A)',ytitle='Q Value',title='Tellurics 0.5',charsize=1.5,xthick=3,ythick=3,charthick=3,thick=3
      device,/close
      set_plot,'x'

   end

;;;;;;;;;;;;;;;;;;;;;;;; Solar

   'sol': begin
      ;read Kurucz spectrum - vacuum, http://kurucz.harvard.edu/stars/sun/
      readcol,'kurucz_models/fsunallp.500000',nm,flux,f='d,d'

      ;change units for consistency with old Phoenix models
      warr=nm*10 ;Angstroms
      ;10^12 erg/s/cm^2/cm - changed wavelengths, integrated over 4pi sr
      sarr=flux*1e7/1e12*(4*!pi) 
      help,warr,sarr

      ;no convolution, this already has Vrot = 2 Km/s
      nsp=sarr

      ;break into blocks...
      wmin=5000.
      wlist=fltarr(2000)
      qlist=fltarr(2000)
      for i=0,1999 do begin
         wmax=wmin+2.
         b=where(warr ge wmin and warr lt wmax,count)
         if i eq 0 then print,warr(b(0))/(warr(b(1))-warr(b(0)))
         if i eq 1999 then print,warr(b(0))/(warr(b(1))-warr(b(0)))
         if count gt 0 then begin
            qvalue=qcalc(warr(b),nsp(b))     
            qlist(i)=qvalue
         endif else qlist(i)=0.
         wlist(i)=wmin
         wmin=wmin+2.
      endfor

      wmin=5000.
      wlist2=fltarr(40)
      qlist2=fltarr(40)
      for i=0,39 do begin
         wmax=wmin+100.
         ;print,wmin,wmax
         b=where(warr ge wmin and warr lt wmax,count)
         ;if i eq 0 then print,warr(b(0))/(warr(b(1))-warr(b(0)))
         if count gt 0 then begin
            ;qvalue=qcalc(warr(b),sarr(b))     
            qvalue=qcalc(warr(b),nsp(b)) 
            qlist2(i)=qvalue
         endif else qlist2(i)=0.
         wlist2(i)=wmin
         wmin=wmin+100.
      endfor

      set_plot,'ps'
      device,file='qsol.ps',/color
      plot,wlist,qlist,xtitle='Wavelength (A)',ytitle='Q Value',title='Solar',zrange=[5000,9000],charsize=1.5,xthick=3,ythick=3,charthick=3,thick=3
      oplot,wlist2,qlist2,color=200
      device,/close
      set_plot,'x'

   end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;GJ 277.1

   'md': begin

      ;using Ryan's reader
      stuff=read_btsettl('/Volumes/RAID/HZPF/stellar_models/BT_SETTL/M-0.0/lte038-5.0-0.0.BT-Settl.7',[5000,9000]) 
      warr=reform(stuff[0,*]) ;A
      sarr=reform(stuff[1,*])*1e8/1e12 ;interpolated, in 10^12 erg/s/cm^2/cm
      help,warr,sarr

      ;convolve with rotation kernel
      lsf = lsf_rotate(0.05,1,velgrid=vel) ;1 Km/s, sampled every 50 m/s
      num_conv,sarr,lsf,nsp

      ;break into blocks...
      wmin=5000.
      wlist=fltarr(2000)
      qlist=fltarr(2000)
      for i=0,1999 do begin
         wmax=wmin+2.
         ;print,wmin,wmax
         b=where(warr ge wmin and warr lt wmax,count)
         if i eq 0 then print,warr(b(0))/(warr(b(1))-warr(b(0)))
         if count gt 0 then begin
            ;qvalue=qcalc(warr(b),sarr(b))     
            qvalue=qcalc(warr(b),nsp(b)) 
            qlist(i)=qvalue
         endif else qlist(i)=0.
         wlist(i)=wmin
         wmin=wmin+2.
      endfor

      wmin=5000.
      wlist2=fltarr(200)
      qlist2=fltarr(200)
      for i=0,199 do begin
         wmax=wmin+100.
         ;print,wmin,wmax
         b=where(warr ge wmin and warr lt wmax,count)
         if i eq 0 then print,warr(b(0))/(warr(b(1))-warr(b(0)))
         if count gt 0 then begin
            ;qvalue=qcalc(warr(b),sarr(b))     
            qvalue=qcalc(warr(b),nsp(b)) 
            qlist2(i)=qvalue
         endif else qlist2(i)=0.
         wlist2(i)=wmin
         wmin=wmin+100.
      endfor

      set_plot,'ps'
      device,file='qmd.ps',/color
      plot,wlist,qlist,xtitle='Wavelength (A)',ytitle='Q Value',title='M0 Dwarf',xrange=[5000,9000],charsize=1.5,xthick=3,ythick=3,charthick=3,thick=3
      oplot,wlist2,qlist2,color=200
      device,/close
      set_plot,'x'

   end

endcase



end



