pro read_terraspec_log,obnm,col1,col7,status

yourpath = '/pool/vonnegut0/mlb/rv/vel/sara'
;yourpath='/Volumes/RAID/rv/vel/sara'
logfile=yourpath+'/files/terraspec_params.log'
status=0

spawn,'\rm deleteme'
spawnst='grep -i -A 18 '+obnm+' '+strtrim(logfile,2)+ ' >  deleteme'
spawn,spawnst ;put that block of text in a temp file

spawn,'wc deleteme',result
   lines=(str_sep(strtrim(result(0),2),' '))(0)

if lines eq 0 then begin ;find other images from that night
   status=1
   spawn,'\rm deleteme'
   spawnst='grep -i -A 18 '+strmid(obnm,0,8)+' '+strtrim(logfile,2)+ ' >  deleteme'
   spawn,spawnst                
endif

close,4 & openr,4,'deleteme'
dum='?'
col1=-1 
col7=-1
col1list=-1
col7list=-1
while (eof(4) eq 0) do begin    ;begin reading log 
   readf,4,dum
   pos1=strpos(dum,'H2O')
   pos7=strpos(dum,'O2')
   if pos1 ne -1 then begin
      str1=strsplit(dum,/extract,count=c1)
      col1=float(str1(c1-1))
      print,col1
      if col1list(0) eq -1 then col1list=col1 else col1list=[col1list,col1]
   endif
   if pos7 ne -1 then begin
      str7=strsplit(dum,/extract,count=c7)
      col7=float(str7(c7-1))
      print,col7
      if col7list(0) eq -1 then col7list=col7 else col7list=[col7list,col7]
   endif
endwhile

col1=col1list[0] ;use the first entry - do something more clever here!
col7=col7list[0]

end
