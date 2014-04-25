pro iraftowavenum,spname,order,starname

;needs to be run on mlb

;matrix[pix, order] - 2d array of observed flux values
;lambda[pix, order] - 2d array of observed wavelengths

imagename=strsplit(spname,'.sp.fits',/extract,/regex)
contname=imagename[0]+'.cont.fits'

; Apply ThAr Dispersion Correction 
orderm1=order-1   ;b/c idl indices start from zero!

	;Reads in information from header and removes non-dispcor information
	;file=dialog_pickfile()

	spec=mrdfits(spname,0,header)
        cont=mrdfits(contname,0,header2)

	indices=where(strmatch(header,'WAT2*'),num_indices)
	firstchar=11

	;Converts header strings into one line
	complete_string=''
	for index=0, num_indices-1 do begin
		stringa=header[indices[index]]
		stringb=strmid(stringa,firstchar,(strlen(stringa)-firstchar-1))
		complete_string=complete_string+stringb	
	endfor

	;Breaks header into individual strings for each order
	order_strings=strsplit(complete_string,'spec*',/extract)
	num_strings=n_elements(order_strings)-1

	;Break order strings into number strings, converts them to doubles,
	; and stores them in a matrix
	headernums=dblarr(21,num_strings)
	for index=3, num_strings do begin
		number_string=strsplit(order_strings[index],' ',/extract)
		num_entries=n_elements(number_string)
		for count=2, num_entries-1 do begin
			astring=number_string[count]
			if (count eq 2) then begin
			;Removes initial quotation mark
				alength=strlen(astring)
				astring=strmid(astring,1,alength-1)
			endif
			if (count eq (num_entries-1)) then begin
			;Removes final quotation mark
				alength=strlen(astring);
				astring=strmid(astring,0,alength-1)
				if (index eq num_strings) then begin
					alength=strlen(astring);
					astring=strmid(astring,0,alength-1)
				endif
			endif
			headernums[(count-2),(index-3)]=double(astring)
			;print,headernums[count-2,index-3]
		endfor
	endfor	

	;Create wavelength solutions from matrices
	;wavelength: headernums[3,*], number of pixels: headernums[5,*],
	;polynomial order: headernums[12,*],first pixel: headernums[13,*],
	;last pixel: headernums[14,*], Chebyshev coefficients: headernums[15-20,*]

	for index=0,num_strings-3 do begin
		;print, headernums[*,0]
		numpix=fix(headernums[5,index])
		pix=dindgen(numpix)+1d0
		if index eq 0 then lambda=replicate(0d0,numpix,(num_strings))
		px1=headernums[13,index]
		pxn=headernums[14,index]
		npoly=headernums[12,index]
		coeff=dblarr(npoly)
		coeff[0:5]=headernums[15:(15+npoly-1),index]

		;Construction of the wavelength given Chebyshev coefficients		
		for count=0, numpix-1 do begin
			x=dblarr(npoly)
			pom=dblarr(npoly)
			xn=(pix[count]-(pxn+px1)/2.0)/((pxn-px1)/2.0)
			x[0]=1d0
			x[1]=xn
			pom[0]=x[0]*coeff[0]
			pom[1]=x[1]*coeff[1]
			for i=2, npoly-1 do begin
				x[i]=2d0*xn*x[i-1]-x[i-2]
				pom[i]=x[i]*coeff[i]
			endfor 
			lambda[count,index]=total(pom)
		endfor
	endfor

        print,"Dispersion Correction applied"

; convert wavelengths to wavenumbers
lambdaord=lambda[*,orderm1]
airtovac,lambdaord

waveord=1./(lambdaord*1.d-8)
;lambdaord = 1.d8/waveord

specord=spec[*,orderm1]/cont[*,orderm1]

;make header for terraspec
m=where(strmatch(header,'MJD*'))
mjd=strmid(header(m),firstchar,19)
jd=double(mjd)+24500000.5d
baryvel,jd,2000,vh,vb

;;Old BC correction
;;RA and DEC values from header are not at all accurate
;ra=strmid(header(where(strmatch(header,'RA*'))),11,11) 
;dec=strmid(header(where(strmatch(header,'DEC*'))),11,11)

;ra=ten(ra)*15/!RADEG
;dec=ten(dec)*15/!RADEG

;bc0 = vb[0]*cos(dec)*cos(ra) + vb[1]*cos(dec)*sin(ra) + vb[2]*sin(dec) 
;print, bc0

jd=double(mjd)+24500000.5d
;better BC correction
barycor_wrapper,starname,jd,corr_save,hjd_save
bc=corr_save

; write ascii files!
fileout=imagename+'_o'+strtrim(string(order),2)+'.dat'
get_lun,fid
openw,fid,fileout

print,mjd,bc,0.0,'comment',format='(d12.6,1x,d12.6,1x,f8.4,1x,a8)'
printf,fid,mjd,bc,0.0,'comment',format='(d12.6,1x,d12.6,1x,f8.4,1x,a8)'

for i=0,4097 do begin
   printf,fid,waveord(i),' ',specord(i)
endfor

close,fid
free_lun,fid

end
