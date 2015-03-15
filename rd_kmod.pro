pro rd_kmod,teff,logg,metal,model,header,tail,type=type
;+
;   Extracts a Kurucz model from his theoretical grid.
;   The 'old' grid was ftp'ed from CCP7@Armagh circa 1997.
;   The 'new' (k2odfnew and ak2odfnew) grids were downloaded from Kurucz's
;   web site in 2006.
;
;    IN: teff    - float - Effective temperature (K) (Only <= 10000 K)
;        logg    - float - log(g) log_10 of the gravity (cm s-2)
;        metal    - float - [Fe/H] = log N(Fe)/N(H) -  log N(Fe)/N(H)[Sun]
;
;    OUT: model    - fltarr(7,ntau) - RHOX,T,P,XNE,ABROSS,ACCRAD,VTURB
;         header     - strarr     - header
;         tail    - strarr     - tail
;
;    KEYWORDS: type  - by default, the k2odfnew grid is used ('type'
;                is internally set to 'odfnew') but this
;                keyword can be also set to 'old' or 'alpha'
;                to use the old models from CCP7, or the
;                ak2odfnew models ([alpha/Fe]=+0.4),respectively.
;
;   C. Allende Prieto, UT, May 1999
;               bug fixed to avoid rounoff errors, UT, April 2005
;               odfnew grids (type keyword), April 2006
;
;-

if n_params() lt 3 then begin
    print,'% RD_KMOD: use -- rd_kmod,teff,logg,metal,model,header,tail[,type=type]'
    return
endif

kpath='~/idl/idl_database/kurucz_models/'

availteff=findgen(27)*250+3500
availlogg=findgen(11)*.5+0.
availmetal=findgen(7)*0.5-2.5

if not keyword_set(type) then type='odfnew'
if type eq 'old' then availmetal=findgen(13)*0.5-5.0

v1=where(abs(availteff-teff) le .1)
v2=where(abs(availlogg-logg) le 0.001)
v3=where(abs(availmetal-metal) le 0.001)
if (v1(0) ne -1 and v2(0) ne -1 and  v3(0) ne -1) then begin
teff=availteff[v1[0]] & logg=availlogg[v2[0]] & metal=availmetal[v3[0]]
s1='a'
if (metal ge 0.) then begin
    s2='p'
endif else begin
     s2='m'
endelse
s3=string(metal,format='(2a)')
if (abs(s3) lt 1.) then begin
     s3=strcompress('0'+string(abs(metal*10),format='(i)'),/remove_all)
endif else begin
    s3=strcompress(string(abs(metal*10),format='(i)'),/remove_all)
endelse

case type of
   'old'  : s4='k2.dat'
   'alpha': s4='ak2odfnew.dat'
    else  : s4='k2odfnew.dat'
endcase

filename=kpath+s1+s2+s3+s4

text='texto'
teffstring=string(teff,format='(f7.0)')
loggstring=string(logg,format='(f8.5)')
header=''

;get_lun,u
u=13
openr,u,filename
readf,u,text
while (strpos(text,teffstring) eq -1 or strpos(text,loggstring) eq -1) do begin
    readf,u,text
endwhile
while (strpos(text,'READ') eq -1) do begin
    header=[header,text]
    readf,u,text
endwhile
header=[header,text]
header=transpose(header)
header=header(1:n_elements(header)-1)
po=strpos(text,'RHOX')-4

ntau=fix(strcompress(strmid(text,po,4),/remove_all))
if ((ntau eq 64 and type eq 'old') or (ntau eq 72)) then begin
    if type eq 'old' then model=dblarr(7,ntau) else model=dblarr(10,ntau)
endif else begin
    print,'% RD_KMOD: trouble! ntau and type do not match!'
    print,'% RD_KMOD: or ntau is neither 64 nor 72'
    stop
endelse

readf,u,model
tail1=''
tail2=''
readf,u,tail1
readf,u,tail2
close,u
tail=transpose([tail1,tail2])
free_lun,u
endif else begin
    print,'The requested values of ([Fe/H],logg,Teff) are NOT available'
endelse

end
