pro kmod,teff,logg,metal,outfile,ntau=ntau,type=type
;+
;    Extracts and if necessary interpolates (linearly) a kurucz model
;    from his grid.
;    The routine is intended for stars cooler than 10000 K.
;    The grid was ftp'ed from CCP7.
;
;    IN: teff    - float - Effective temperature (K)
;        logg    - float - log(g) log_10 of the gravity (cm s-2)
;        metal    - float - [Fe/H] = log N(Fe)/N(H) -  log N(Fe)/N(H)[Sun]
;
;    OUT: outfile     - string - name for the output file
;
;    KEYWORD: ntau    - returns the number of depth points in the output model
;
;             type  - by default, the k2odfnew grid is used ('type'
;                is internally set to 'odfnew') but this
;                keyword can be also set to 'old' or 'alpha'
;                to use the old models from CCP7, or the
;                ak2odfnew models ([alpha/Fe]=+0.4),respectively.
;
;
;    C. Allende Prieto, UT, May 1999
;               bug fixed, UT, Aug 1999
;               bug fixed to avoid rounoff errors, keyword ntau added
;                    UT, April 2005
;               bug fixed, assignment of the right tauscale to each
;                model (deltaT<1%), UT, March 2006
;               odfnew grids (type keyword), April 2006
;
;-
npar=n_params(0)
if (npar eq 0) then begin
    print,'kmod,teff,logg,[Fe/H],model[,ntau=ntau,type=type]'
    return
endif


availteff=findgen(27)*250+3500
availlogg=findgen(11)*.5+0.
availmetal=findgen(7)*0.5-2.5

if not keyword_set(type) then type='odfnew'
if type eq 'old' then availmetal=findgen(13)*0.5-5.0

if type eq 'old' then ntau=64 else ntau=72

v1=where(abs(availteff-teff) le .1)
v2=where(abs(availlogg-logg) le 0.001)
v3=where(abs(availmetal-metal) le 0.001)

if (teff le max(availteff) and teff ge min(availteff) and logg le max(availlogg) and logg ge min(availlogg) and metal ge min(availmetal) and metal le max(availmetal)) then begin
if (v1(0) ne -1 and v2(0) ne -1 and  v3(0) ne -1) then begin
    ; Direct extraction of the model
    teff=availteff[v1[0]] & logg=availlogg[v2[0]] & metal=availmetal[v3[0]]
    rd_kmod,teff,logg,metal,model,header,tail,type=type
    ntau=n_elements(model(0,*))

endif else begin


    ; Linear Interpolation
    teffimif=max(where(availteff le teff)) ; immediately inferior Teff
    loggimif=max(where(availlogg le logg)) ; immediately inferior logg
    metalimif=max(where(availmetal le metal)) ;immediately inferior [Fe/H]
    teffimsu=min(where(availteff ge teff)) ; immediately superior Teff
    loggimsu=min(where(availlogg ge logg)) ; immediately superior logg
    metalimsu=min(where(availmetal ge metal)) ;immediately superior [Fe/H]


    if type eq 'old' then ncols=7 else ncols=10

    grid=fltarr(2,2,2,ncols)
    tm1=availteff(teffimif)
    tp1=availteff(teffimsu)
    lm1=availlogg(loggimif)
    lp1=availlogg(loggimsu)
    mm1=availmetal(metalimif)
    mp1=availmetal(metalimsu)

    if (tp1 ne tm1) then begin
        mapteff=(teff-tm1)/(tp1-tm1)
    endif else begin
        mapteff=0.5
    endelse
    if (lp1 ne lm1) then begin
        maplogg=(logg-lm1)/(lp1-lm1)
    endif else begin
        maplogg=0.5
    endelse
    if (mp1 ne mm1) then begin
        mapmetal=(metal-mm1)/(mp1-mm1)
    endif else begin
        mapmetal=0.5
    endelse

    ; Reading the corresponding models

    for i=1,8 do begin

        if i eq 1 then rd_kmod,tm1,lm1,mm1,model,type=type,header,tail
        if i eq 2 then rd_kmod,tm1,lm1,mp1,model,type=type
        if i eq 3 then rd_kmod,tm1,lp1,mm1,model,type=type
        if i eq 4 then rd_kmod,tm1,lp1,mp1,model,type=type
        if i eq 5 then rd_kmod,tp1,lm1,mm1,model,type=type
        if i eq 6 then rd_kmod,tp1,lm1,mp1,model,type=type
        if i eq 7 then rd_kmod,tp1,lp1,mm1,model,type=type
        if i eq 8 then rd_kmod,tp1,lp1,mp1,model,type=type

        if (n_elements(model(0,*)) gt ntau) then begin
            m2=fltarr(ncols,ntau)
            m2(0,*)=interpol(model(0,*),ntau)
            for j=0,ncols-1 do begin
                m2(j,*)=interpol(model(j,*),model(0,*),$
                m2(0,*))
            endfor
        model=m2
        endif
        ;getting the tauross scale
        rhox=model(0,*)
        kappaross=model(4,*)
        tauross=fltarr(ntau)
        tauross(0)=rhox(0)*kappaross(0)
        for ii=1,ntau-1 do begin
            tauross(ii)=int_tabulated(rhox(0:ii),kappaross(0:ii))
        endfor

        case i of
            1: begin
                model1=model
                tauross1=tauross
            end
            2: begin
                model2=model
                tauross2=tauross
            end
            3: begin
                model3=model
                tauross3=tauross
            end
            4 :begin
                model4=model
                tauross4=tauross
            end
            5: begin
                model5=model
                tauross5=tauross
            end
            6: begin
                model6=model
                tauross6=tauross
            end
            7: begin
                model7=model
                tauross7=tauross
            end
            8: begin
                model8=model
                tauross8=tauross
            end
            else: message,'% KMOD: i should be 1--8!'
        endcase

    endfor

    model=fltarr(ncols,ntau) ; cleaning up for re-using the matrix

    ; defining the  mass (RHOX;gr cm-2) sampling
    tauross=tauross1 ; re-using the vector tauross
    bot_tauross=min([tauross1(ntau-1),tauross2(ntau-1),$
    tauross3(ntau-1),tauross4(ntau-1),$
    tauross5(ntau-1),tauross6(ntau-1),tauross7(ntau-1),tauross8(ntau-1)])
    top_tauross=max([tauross1(0),tauross2(0),tauross3(0),$
    tauross4(0),tauross5(0),tauross6(0),tauross7(0),tauross8(0)])
    tauross_new=interpol(tauross(where(tauross ge $
    top_tauross and tauross le  bot_tauross)),ntau)

    ; let's interpolate for every depth
    for i=1,ntau-1 do begin
        for j=0,ncols-1 do begin
        grid(0,0,0,j)=interpol(model1(j,1:ntau-1),tauross1(1:ntau-1),$
        tauross_new(i))
        grid(0,0,1,j)=interpol(model2(j,1:ntau-1),tauross2(1:ntau-1),$
        tauross_new(i))
        grid(0,1,0,j)=interpol(model3(j,1:ntau-1),tauross3(1:ntau-1),$
        tauross_new(i))
        grid(0,1,1,j)=interpol(model4(j,1:ntau-1),tauross4(1:ntau-1),$
        tauross_new(i))
        grid(1,0,0,j)=interpol(model5(j,1:ntau-1),tauross5(1:ntau-1),$
        tauross_new(i))
        grid(1,0,1,j)=interpol(model6(j,1:ntau-1),tauross6(1:ntau-1),$
        tauross_new(i))
        grid(1,1,0,j)=interpol(model7(j,1:ntau-1),tauross7(1:ntau-1),$
        tauross_new(i))
        grid(1,1,1,j)=interpol(model8(j,1:ntau-1),tauross8(1:ntau-1),$
        tauross_new(i))
        model(j,i)=$
        interpolate(grid(*,*,*,j),mapteff,maplogg,mapmetal)
        endfor
    endfor

        for j=0,ncols-1 do begin
           model(j,0)=model(j,1)*0.999
        endfor


    ; editing the header
    tmpstr=header(0)
    strput,tmpstr,string(teff,format='(f7.0)'),5
    strput,tmpstr,string(logg,format='(f8.5)'),21
    header(0)=tmpstr
    tmpstr1=header(1)
    tmpstr2=header(4)
    if (metal lt 0.0) then begin
        if type eq 'old' then begin
          strput,tmpstr1,string(abs(metal),format='("-",f3.1)'),18
        endif else begin
          strput,tmpstr1,string(abs(metal),format='("-",f3.1)'),8
        endelse
        strput,tmpstr2,string(10^metal,format='(f9.5)'),16
    endif else begin
        if type eq 'old' then begin
          strput,tmpstr1,string(abs(metal),format='("+",f3.1)'),18
        endif else begin
          strput,tmpstr1,string(abs(metal),format='("+",f3.1)'),8
        endelse
        strput,tmpstr2,string(10^metal,format='(f9.5)'),16
    endelse
    header(1)=tmpstr1
    header(4)=tmpstr2
    tmpstr=header(22)
    strput,tmpstr,string(ntau,format='(i2)'),11
    header(22)=tmpstr

endelse

; writing the outputfile
;get_lun,u
u=14
openw,u,outfile
for i=0,n_elements(header)-1 do begin
    printf,u,header(i)
endfor
if type eq 'old' then begin
    for i=0,ntau-1 do begin
        printf,u,model(*,i),$
        format='(E15.8,x,f8.1,5(x,E9.3))'
    endfor
endif else begin
    for i=0,ntau-1 do begin
        printf,u,model(*,i),$
        format='(E15.8,x,f8.1,8(x,E9.3))'
    endfor
endelse
for i=0,n_elements(tail)-1 do begin
    printf,u,tail(i)
endfor
close,u
free_lun,u
endif else begin
    print,'% KMOD:  The requested values of ([Fe/H],logg,Teff) fall outside
    print,'% KMOD:  the boundaries of the grid.'
    print,'% KMOD:  Temperatures higher that 10000 K can be reached, by modifying rd_kmod.'

endelse
end
