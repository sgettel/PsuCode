pro read_mercury_grid,phrase,acount=acount,mcount=mcount,stable=st

mpath='/astro/grads/sjg223/research/mercury6/'
dirs=file_search(mpath+phrase+'*')
;print,dirs

stable=intarr(acount,mcount)

for i=0, n_elements(dirs)-1 do begin
   result=dirs[i]+'/info.out'
   ;print,result

   num=strsplit(dirs[i],'test',/extract,/regex)
   num2=strsplit(num[1],'-',/extract)
   ;print,num2;,n_elements(num2)
   index=fix(num2[0])-1
   ;print,index
   if n_elements(num2) eq 1 then begin ;skip long integrations
      spawn,'\rm deleteme'
      spawn,'grep PLANET '+result+' > deleteme'
      close,4 & openr,4,'deleteme'
      
      if (eof(4) eq 0) then st=0 else st=1
                                ;print,st
   

      stable[index]=st
   endif
   
endfor



end
