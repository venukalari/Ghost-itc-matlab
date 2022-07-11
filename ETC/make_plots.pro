function fullrange, x
return, max(x)-min(x)
end

function midrange, x
return, min(x) + 0.5*(max(x)-min(x))
end

pro run_ghost_etc,PRINT=print, DIR=dir
if not keyword_set(DIR) then dir = './'

;  sel = [0,0,0,1]
  sel = [0,1,0,0]
  sel=[1,1,1,0]
  
;############
; MAG LIMITS
;############

  mags = RANGE(9,18)
  waves = [3630,3750,4500,5500,9000,9500]
;  modes = ["SR",'SF','BS','BF','HR','HF','PRV']
  modes = ["SR",'SF','SVF','BS','BF','BVF','HR','HF','PRV']
  nmags = n_elements(mags)
  nwaves = n_elements(waves)
  nmodes = n_elements(modes)

  TargetSN = 30
  pre=''
  bright = 0 & pre=''
  

  openw, lun, dir+'/ghost_maglims_budgeted_sn'+STRTRIM(STRING(TargetSN),2)+pre+'.dat',/get_lun
  printf,lun,'Modes',waves,F='(a,x,'+STRING(nwaves,F='(i)')+'(i5,x))'
  if sel[0] eq 1 then begin
     print,'Computing GHOST magnitude limits'
     print

; BUDGETED MAGLIMITS
     print,'Budgeted throughput'
     outfile = dir+'/ghost_maglims_budgeted_sn'+STRTRIM(STRING(TargetSN),2)+pre+'.dat'
     openw, lun, outfile,/get_lun
     printf,lun,'Modes',waves,F='(a,x,'+STRING(nwaves,F='(i)')+'(i5,x))'
     for k=0,nmodes-1 do begin
        maglims = fltarr(nwaves,nmodes)
        for i =0,n_elements(waves)-1 do begin
           ghost_etc,wave=waves[i],action='F',mode=modes[k],flux=m,/ab,sn=TargetSN,texp=3600,/quiet,BRIGHT=bright,THROUGHPUT="B"
           maglims[i,k] = m
        endfor
        printf,lun,modes[k],maglims[*,k],F='(a,x,'+STRING(nwaves,F='(i)')+'(f5.2,x))'
     endfor
     free_lun, lun
     make_tables2_pdf,outfile,"``Budgeted'' Throughput (AB mag)"
     
; BUDGETED MAGLIMITS, upward port (no M3)
     print,'Budgeted throughput upward port'
     outfile = dir+'/ghost_maglims_budgeted_upward_sn'+STRTRIM(STRING(TargetSN),2)+pre+'.dat'
     openw, lun, outfile,/get_lun
     printf,lun,'Modes',waves,F='(a,x,'+STRING(nwaves,F='(i)')+'(i5,x))'
     for k=0,nmodes-1 do begin
        maglims = fltarr(nwaves,nmodes)
        for i =0,n_elements(waves)-1 do begin $
           ghost_etc,wave=waves[i],action='F',mode=modes[k],flux=m,/ab,sn=TargetSN,texp=3600,/quiet,BRIGHT=bright,THROUGHPUT="B", PORT=1
           maglims[i,k] = m
        endfor
        printf,lun,modes[k],maglims[*,k],F='(a,x,'+STRING(nwaves,F='(i)')+'(f5.2,x))'
     endfor
     free_lun, lun
     make_tables2_pdf,outfile,"``Budgeted'' Throughput on Upward Port (AB mag)"
     
; PREDICTED MAGLIMITS
;     print,'Predicted throughput'
;     openw, lun, dir+'/ghost_maglims_predicted_sn'+STRTRIM(STRING(TargetSN),2)+pre+'.dat',/get_lun
;     printf,lun,'Modes',waves,F='(a,x,'+STRING(nwaves,F='(i)')+'(i5,x))'
;     for k=0,nmodes-1 do begin
;        maglims = fltarr(nwaves,nmodes)
;        for i =0,n_elements(waves)-1 do begin
;           ghost_etc,wave=waves[i],action='F',mode=modes[k],flux=m,/ab,sn=TargetSN,texp=3600,/quiet,BRIGHT=bright,THROUGHPUT="P"
;           maglims[i,k] = m
;        endfor
;        printf,lun,modes[k],maglims[*,k],F='(a,x,'+STRING(nwaves,F='(i)')+'(f5.2,x))'
;     endfor
;     free_lun, lun

; REQUIREMENTS MAGLIMITS
     print,'Requirements throughput (as per REQ1003, 4110, etc.)'
     outfile = dir+'/ghost_maglims_requirements_sn'+STRTRIM(STRING(TargetSN),2)+pre+'.dat'
     openw, lun, outfile,/get_lun
     printf,lun,'Modes',waves,F='(a,x,'+STRING(nwaves,F='(i)')+'(i5,x))'
     for k=0,nmodes-1 do begin
        maglims = fltarr(nwaves,nmodes)
        for i =0,n_elements(waves)-1 do begin
           ghost_etc,wave=waves[i],action='F',mode=modes[k],flux=m,/ab,sn=TargetSN,texp=3600,/quiet,BRIGHT=bright,THROUGHPUT="R"
           maglims[i,k] = m
        endfor
        printf,lun,modes[k],maglims[*,k],F='(a,x,'+STRING(nwaves,F='(i)')+'(f5.2,x))'
     endfor
     free_lun, lun
     make_tables2_pdf,outfile,"``Requirements'' Throughput (AB mag)"

; SHORT FIBER MAGLIMITS
     print,'Requirements assuming SHORT FIBER'
     outfile = dir+'/ghost_maglims_shortfiber_sn'+STRTRIM(STRING(TargetSN),2)+pre+'.dat'
     openw, lun, outfile,/get_lun
     printf,lun,'Modes',waves,F='(a,x,'+STRING(nwaves,F='(i)')+'(i5,x))'
     for k=0,nmodes-1 do begin
        maglims = fltarr(nwaves,nmodes)
        for i =0,n_elements(waves)-1 do begin
           ghost_etc,wave=waves[i],action='F',mode=modes[k],flux=m,/ab,sn=TargetSN,texp=3600,/quiet,BRIGHT=bright,THROUGHPUT="S"
           maglims[i,k] = m
        endfor
        printf,lun,modes[k],maglims[*,k],F='(a,x,'+STRING(nwaves,F='(i)')+'(f5.2,x))'
     endfor
     free_lun, lun
     make_tables2_pdf,outfile,"``Short Fiber'' Throughput on Upward Port (AB mag)"
  endif
  
  if sel[1] eq 1 then begin
     
     if keyword_set(PRINT) then begin
        if not keyword_set(DIR) then dir = './'
        pfile = dir+'performance.ps'
        set_ps,/open,file=pfile,/landscape
        setplot_paper
        !p.charthick=6
        !x.thick=6 &   !y.thick=6
     endif

     ts = fltarr(n_elements(mags),n_elements(waves))
;     cols = bytscl(waves,top=220)+30
     cols = ['purple','cyan','blue','forest green','orange','red']
     loadct,39
     sn=50
     for i =0,n_elements(waves)-1 do begin
        for j=0,n_elements(mags)-1 do begin
           ghost_etc,wave=waves[i],action='T',mode='SR',flux=mags[j],/ab,sn=sn,texp=texp,through='R'
           ts[j,i] = texp
        endfor
        if i eq 0 then plot,mags,ts[*,i],/ylog,xtitle='AB Mag',ytitle='S/N='+STRING(sn,F='(i3)')+' Exposure Time (s)',/nodata,xrange=minmax(mags),yrange=[1,10.^5],/ysty,/xsty,title='GHOST Efficiency'
        if waves[i] eq 4500 then thick=10 else thick=4
;        oplot,mags,ts[*,i],col=cols[i],thick=thick
        oplot,mags,ts[*,i],col=cgcolor(cols[i]),thick=thick
        dx = !x.crange[1]-!x.crange[0]  &       dy = !y.crange[1]-!y.crange[0]
        x1 = !x.crange[0]+0.05*dx        &       y1 = 10^(!y.crange[0]+(0.9-i*0.05)*dy)
        XYOUTS,[x1],[y1],STRING(waves[i],F='(i5)')+'('+cgsymbol("angstrom")+')',col=cgcolor(cols[i]),charthick=thick<7
     endfor
     if keyword_set(PRINT) then set_ps,file=pfile,/landscape,/close else stop  
  endif

;===================================
;Comparison with other spectrographs
;===================================
  if sel[2] eq 1 then begin
     readcol,'~/Idl/GHOST_ETC/RefData/ffeige66.dat',wave,flux,/silent
     fstd = flux / 1.e16        ; get into erg/s/cm^2/A
;stop
;readcol,'~/Gemini/GHOST/ETC/feige66_002.ascii',wave,fstd,/silent
;readcol,'~/Gemini/GHOST/ETC/f66.oke',wave,fstd,/silent
;fstd = fstd * 10.^7/10.^6
;stop
     wavex = RANGE(max(wave),10000) ; extend the wavelength range
     fluxx = planck(wavex,34500)*(5.8d6)^2/(10.*3.086d16)^2  ; Using Teff from Petit+11 = 34500K, R~0.009Rsun, D=10pc (abs mag)
     wave = [wave,wavex]
     fstd = [fstd,fluxx]
     lambdas=[RANGE(3630,10000,30),4499,4500,4501]
;     lambdas=[RANGE(3630,10000,30)]
     lambdas = lambdas[sort(lambdas)]
     nl = n_elements(lambdas)
     sns = fltarr(nl)
     sns2 = sns
     sns3 = sns
     for i=0, nl-1 do begin
        ghost_etc,wave=lambdas[i],action='S',mode='SR',texp=3600,flux=interpol(fstd,wave,lambdas[i]),sn=sn,/quiet,seeing=0.4,through='R'
        sns[i] = sn
        ghost_etc,wave=lambdas[i],action='S',mode='SR',texp=3600,flux=interpol(fstd,wave,lambdas[i]),sn=sn,/quiet,seeing=0.4,through='B'
        sns2[i] = sn
;        ghost_etc,wave=lambdas[i],action='S',mode='SR',texp=3600,flux=interpol(fstd,wave,lambdas[i]),sn=sn,/quiet,seeing=0.4,through='P'
;        sns3[i] = sn
     endfor
     ; REference design wavelength point
     ghost_etc,wave=4500.,action='S',mode='SR',texp=3600,flux=interpol(fstd,wave,4500.),sn=sn,/quiet,seeing=0.4,through='R'
     snref = sn
;     if keyword_set(PRINT) then begin
;        if not keyword_set(DIR) then dir = './'
;        pfile = dir+'feige66_ghost.ps'
;        set_ps,/open,file=pfile,/landscape
;        setplot_paper
;        !p.charthick=6
;        !x.thick=6 &   !y.thick=6
;     endif
;     plot,lambdas/10.,sns,yrange=[0,1300],/ysty,xrange=[410,1020],/xsty,/nodata
;     oplot,lambdas/10.,sns,thick=10,col=cgcolor('red')
;     if keyword_set(PRINT) then set_ps,/close,file=pfile,/landscape else stop

;     readcol,'/Users/rmcdermid/Gemini/GHOST/ETC/GRACES/GRACEShigh.dat',ghl,ghf
     readcol,'/Users/rmcdermid/Idl/GHOST_ETC/RefData/GRACESlow.dat',gll,glf
;     readcol,'/Users/rmcdermid/Gemini/GHOST/ETC/GRACES/HIREShigh.dat',hhl,hhf
     readcol,'/Users/rmcdermid/Idl/GHOST_ETC/RefData/HIRESlow.dat',hll,hlf
     readcol,'/Users/rmcdermid/Idl/GHOST_ETC/RefData/UVESlow.dat',ull,ulf

     if keyword_set(PRINT) then begin
        if not keyword_set(DIR) then dir = './'
        pfile = dir+'ghost_comparison.ps'
        set_ps,/open,file=pfile,/landscape
        setplot_paper
        !p.charthick=6
        !x.thick=6 &   !y.thick=6
     endif

     plot,gll*10.,glf,xrange=[3000,10500],/nodata,/xsty,xtitle='Wavelength (Angstrom)',yrange=[1,2100],ytitle='S/N per resolution element in 3600s',$
          title='GHOST Performance Comparison',/ysty
     oplot,gll*10.,glf,col=cgcolor('blue'),psym=3
     oplot,hll*10.,hlf,col=cgcolor('green'),psym=3
     oplot,ull*10.,ulf,col=cgcolor('magenta'),psym=3
     oplot,lambdas,sns,thick=10,col=cgcolor('red')
     oplot,lambdas,sns2,thick=5,col=cgcolor('red'),linesty=2
;     oplot,lambdas,sns3,thick=5,col=cgcolor('red'),linesty=2
;     oplot, [4500.],[snref], thick=5, col=cgcolor('red'),psym=4

;dx = !x.crange[1]-!x.crange[0]  &       dy = !y.crange[1]-!y.crange[0]
;x1 = !x.crange[0]+0.6*dx       &       y1 = !y.crange[0]+0.9*dy
;oplot,[x1-0.05*dx,x1-0.01*dx],[y1,y1]+0.015*dy,col=cgcolor('red'),thick=10
;XYOUTS,[x1],[y1],"GHOST (R=50k)"
;y1 = y1-0.05*dy
;oplot,[x1-0.05*dx,x1-0.01*dx],[y1,y1]+0.015*dy,col=cgcolor('blue'),thick=5;,line=0
;XYOUTS,[x1],[y1],"GRACES (R=40k)"
;y1 = y1-0.05*dy
;oplot,[x1-0.05*dx,x1-0.01*dx],[y1,y1]+0.015*dy,col=cgcolor('green'),thick=5;,line=0
;XYOUTS,[x1],[y1],"HIRES (R=37.5k)"
;y1 = y1-0.05*dy
;oplot,[x1-0.05*dx,x1-0.01*dx],[y1,y1]+0.015*dy,col=cgcolor('magenta'),thick=5;,line=0
;XYOUTS,[x1],[y1],"UVES (R=31k)"

;oplot,[x1-0.1*dx,x1+0.35*dx,x1+0.35*dx,x1-0.1*dx,x1-0.1*dx],$
;      [y1-0.02*dy,y1-0.02*dy,y1+0.2*dy,y1+0.2*dy,y1-0.02*dy]

dx = !x.crange[1]-!x.crange[0]  &       dy = !y.crange[1]-!y.crange[0]
x1 = !x.crange[0]+0.35*dx       &       y1 = !y.crange[0]+0.92*dy
XYOUTS,[x1],[y1],"GHOST (R=50k)",col=cgcolor('red')
y1 = y1-0.05*dy
XYOUTS,[x1],[y1],"GRACES (R=40k)",col=cgcolor('blue')
y1 = y1-0.05*dy
XYOUTS,[x1],[y1],"HIRES (R=37.5k)",col=cgcolor('green')
y1 = y1-0.05*dy
XYOUTS,[x1],[y1],"UVES (R=31k)",col=cgcolor('magenta')

oplot,[x1-0.03*dx,x1+0.63*dx,x1+0.63*dx,x1-0.03*dx,x1-0.03*dx],$
      [y1-0.03*dy,y1-0.03*dy,y1+0.2*dy,y1+0.2*dy,y1-0.03*dy]

x1 = !x.crange[0]+0.735*dx       &       y1 = !y.crange[0]+0.92*dy
oplot,[x1-0.05*dx,x1-0.01*dx],[y1,y1]+0.015*dy,col=cgcolor('red'),thick=10
XYOUTS,[x1],[y1],"Requirements"
;y1 = y1-0.06*dy
;oplot,[x1-0.05*dx,x1-0.01*dx],[y1,y1]+0.015*dy,col=cgcolor('red'),thick=5,linesty=2
;XYOUTS,[x1],[y1],"Predicted"
y1 = y1-0.06*dy
oplot,[x1-0.05*dx,x1-0.01*dx],[y1,y1]+0.015*dy,col=cgcolor('red'),thick=5,linesty=2
XYOUTS,[x1],[y1],"Budgeted"
;y1 = y1-0.06*dy
;oplot,[x1-0.03*dx],[y1]+0.015*dy,col=cgcolor('red'),thick=5,psym=4
;XYOUTS,[x1],[y1],"Reference"

;y1 = y1-0.04*dy
;XYOUTS,[x1],[y1],"Contingency"
;oplot,[x1-0.05*dx,x1-0.01*dx],[y1,y1]+0.015*dy,col=cgcolor('blue'),thick=5;,line=0
;oplot,[x1-0.05*dx,x1-0.01*dx],[y1,y1]+0.015*dy,col=cgcolor('green'),thick=5;,line=0
;oplot,[x1-0.05*dx,x1-0.01*dx],[y1,y1]+0.015*dy,col=cgcolor('magenta'),thick=5;,line=0



     if keyword_set(PRINT) then set_ps,/close,file=pfile,/landscape else stop
  endif

  if sel[3] eq 1 then begin
     waves = [4500.,RANGE(3630.,10000.,200)]
     waves = waves[sort(waves)]
     nwaves = n_elements(waves)
     mags = fltarr(nwaves)
     for i=0,nwaves-1 do begin
        ghost_etc,wave=waves[i],action='F',mode='SF',texp=3600,flux=mag,sn=30,/quiet,through='R',/ab
        mags[i] = mag
        print_progress,i,nwaves
     endfor
     stop
  endif
  
end

pro plot_extinct, PRINT=print, DIR=dir
; GS: Patat et al. 2011, A&A, 527, AA91
; GN: Buton et al. 2013, AAp, 549, AA8
readcol,'~/Idl/GHOST_ETC/RefData/paranal_patat11.dat',lext,ext,/silent
readcol,'~/Idl/GHOST_ETC/RefData/MK_extinction_Buton.dat',mklext,mkext,F='(f,f)',/silent

if keyword_set(PRINT) then begin
   if not keyword_set(DIR) then dir = './'
   pfile = dir+'extinction.ps'
   set_ps,/open,file=pfile,/landscape
   setplot_paper
endif

thick=6
plot,lext,ext,xtitle='Wavelength ('+cgsymbol("angstrom")+')',ytitle='Extinction Coefficient, '+textoidl('k_\lambda')+' [mag/airmass]',thick=thick,$
     xrange=[3000,10500],/xsty, title='Atmospheric Extinction'
oplot,mklext,mkext,col=cgcolor('blue'),thick=thick

dx = !x.crange[1]-!x.crange[0]  &       dy = !y.crange[1]-!y.crange[0]
x1 = !x.crange[0]+0.3*dx       &       y1 = !y.crange[0]+0.8*dy

oplot,[x1-0.05*dx,x1-0.01*dx],[y1,y1],thick=thick
oplot,[x1-0.05*dx,x1-0.01*dx],[y1,y1]+0.07*dy,col=cgcolor('blue'),thick=thick
XYOUTS,[x1],[y1-0.01*dy],'Cerro Paranal (Patat et al., 2011)'
XYOUTS,[x1],[y1+0.06*dy],'Mauna Kea (Buton et al., 2013)'

if keyword_set(PRINT) then set_ps,/close,file=pfile,/landscape else stop
end

pro fit_reflect, PRINT=print, DIR=dir
; Fit reflectivity of Gemini mirrors using data from Vucina et al 2008
; Proc. of SPIE Vol. 7012 70122Q-1
; http://www.saao.ac.za/~dod/M5_Washing/SPIE7012_101_Gemini.pdf

if keyword_set(PRINT) then begin
   if not keyword_set(DIR) then dir = './'
   pfile = dir+'reflectivity.ps'
   set_ps,/open,file=pfile,/landscape
   setplot_paper
   thick=8
endif

readcol,'~/Idl/GHOST_ETC/RefData/reflectivity.dat',x,y,dy,/silent
readcol,'~/Idl/GHOST_ETC/RefData/Gemini_Reflectivity_2008-2015.dat',x2,y2,/silent,F='(f,x,x,x,x,f)'
readcol,'~/Idl/GHOST_ETC/RefData/GS_M1_Ag_Sample_10-8-2010.dat',x3,y3,/silent,F='(f,f)'
readcol,'~/Idl/GHOST_ETC/RefData/reflectivity_fresh.dat',lnmRef,Reflectivity,/silent
ploterror,x,y,dy,psym=4,/xsty,ysty=3,xrange=[3000,10500],xtitle='Wavelength ('+cgsymbol("angstrom")+')',ytitle='Reflectivity (%)',thick=thick,yrange=[70,100],title='Gemini Primary Mirror Reflectivity',/nodata
oplot,lnmRef*10d0,Reflectivity,col=cgcolor('red'),thick=thick
oplot,lnmRef*10d0,Reflectivity*0.95,col=cgcolor('red'),linesty=2,thick=thick
oploterror,x,y,dy,psym=4
oplot,x2*10.,y2,col=cgcolor('black')
oplot,x3*10.,y3,col=cgcolor('green')

horizontal,val=95,col=cgcolor('black')
vertical,val=4500,col=cgcolor('black')

readcol,'~/Idl/GHOST_ETC/RefData/boccas_reflectivity.dat',lnmRef,Reflectivity,/silent
oplot,lnmRef*10d0,smooth(Reflectivity,2),col=cgcolor('blue'),thick=thick

;res=robust_linefit(alog10(x),y,yfit)
;xx = RANGE(3000,30000,200)
;oplot,xx,res[1]*alog10(xx)+res[0],col=cgcolor('black'),thick=3

dx = !x.crange[1]-!x.crange[0]  &       dy = !y.crange[1]-!y.crange[0]
x1 = !x.crange[0]+0.25*dx       &       y1 = !y.crange[0]+0.35*dy

oplot,[x1,x1+0.05*dx],[y1,y1],col=cgcolor('blue'),thick=thick
XYOUTS,[x1+0.07*dx],[y1-0.01*dy],"Bare Silver (Boccas et al. 2006)"
y1 = y1-0.05*dy
oploterror,[x1+0.025*dx],[y1-0.02*dy],[0.02*dy],psym=4
XYOUTS,[x1+0.07*dx],[y1-0.01*dy],"Gemini South M1 after 4 years"
y1 = y1-0.05*dy
XYOUTS,[x1+0.07*dx],[y1-0.01*dy],"(Vucina et al., 2008)"
y1 = y1-0.05*dy
oplot,[x1,x1+0.05*dx],[y1,y1],col=cgcolor('red'),thick=thick
XYOUTS,[x1+0.07*dx],[y1-0.01*dy],'2010 Gemini South Witness Sample'
y1 = y1-0.05*dy
oplot,[x1,x1+0.05*dx],[y1,y1],col=cgcolor('red'),thick=thick,linesty=2
XYOUTS,[x1+0.07*dx],[y1-0.01*dy],'2010 Witness Sample * 0.95'
;oplot,[x1-0.05*dx,x1],[y1-0.04*dy]#[1,1],col=cgcolor('black'),thick=3
;XYOUTS,[x1+0.02*dx],[y1-0.05*dy],STRING(REVERSE(res),FORMAT='("y = ",f5.2,"log(x)",f+7.2)')

if keyword_set(PRINT) then set_ps,/close,file=pfile,/landscape else stop
end


pro plot_reflectivity, PRINT=print, DIR=dir
; Plot the various data from the Excel sheet provided by Gemini:
; "Gemini_Reflectivity_2008-2015.txt", Version 22 April 2005
;
; The spreadshet gives complete reflectivity curves for witness
; samples at GN and GS, plus 'IRIS' reflectometer values at 470nm,
; 530nm, 650nm, and 880nm taken at regular intervals on M1 and M2.

if keyword_set(PRINT) then begin
   if not keyword_set(DIR) then dir = './'
   pfile = dir+'reflectivity.ps'
   set_ps,/open,file=pfile,/landscape,xsize=26,ysize=10
   setplot_paper
   thick=4
   !p.charsize=1.3
endif

; Witness Samples
readcol,'~/Idl/GHOST_ETC/RefData/GN_M1_Ag_Sample_8-14-2013.dat',x,y,/silent,F='(f,x,x,x,x,f)'
readcol,'~/Idl/GHOST_ETC/RefData/GS_M1_Ag_Sample_10-8-2010.dat',x2,y2,/silent,F='(f,f)'

;==================================================================
; SCALING FOR M1 and M2
;
; These scalings were determined by eye by R.McDermid, 10/11/2015,
; fitting the lower bound of the two bluest measured points.
;
m1scale = 0.95
m2scale = 0.975
;
;==================================================================

; IRIS measurements
readcol,'~/Idl/GHOST_ETC/RefData/GN_M1_IRIS.dat',gnM1_47,gnM1_53, gnM1_65, gnM1_88, /silent,F='(x,f,f,f,f)', skipline=2
readcol,'~/Idl/GHOST_ETC/RefData/GN_M2_IRIS.dat',gnM2_47,gnM2_53, gnM2_65, gnM2_88, /silent,F='(x,f,f,f,f)', skipline=2
readcol,'~/Idl/GHOST_ETC/RefData/GS_M1_IRIS.dat',gsM1_47,gsM1_53, gsM1_65, gsM1_88, /silent,F='(x,f,f,f,f)', skipline=2
readcol,'~/Idl/GHOST_ETC/RefData/GS_M2_IRIS.dat',gsM2_47,gsM2_53, gsM2_65, gsM2_88, /silent,F='(x,f,f,f,f)', skipline=2

!p.color=cgcolor('black')
!p.background=cgcolor('white')
!p.multi=[0,2,1]
plot,x*10.,y,psym=4,/xsty,ysty=3,xrange=[3000,10500],xtitle='Wavelength ('+cgsymbol("angstrom")+')', $
          ytitle='Reflectivity (%)',thick=thick,yrange=[70,100],title='Gemini North Mirror Reflectivity',/nodata
oplot,x*10.,y,col=cgcolor('black'),thick=thick ; GN Sample

oploterror,4715.,midrange(gnM1_47*100.),0.5*fullrange(gnM1_47*100.),psym=4,col=cgcolor('blue'),errcol=cgcolor('blue'),thick=thick
oploterror,5315.,midrange(gnM1_53*100.),0.5*fullrange(gnM1_53*100.),psym=4,col=cgcolor('blue'),errcol=cgcolor('blue'),thick=thick
oploterror,6515.,midrange(gnM1_65*100.),0.5*fullrange(gnM1_65*100.),psym=4,col=cgcolor('blue'),errcol=cgcolor('blue'),thick=thick
oploterror,8815.,midrange(gnM1_88*100.),0.5*fullrange(gnM1_88*100.),psym=4,col=cgcolor('blue'),errcol=cgcolor('blue'),thick=thick
oploterror,4685.,midrange(gnM2_47*100.),0.5*fullrange(gnM2_47*100.),psym=5,col=cgcolor('red'),errcol=cgcolor('red'),thick=thick
oploterror,5285.,midrange(gnM2_53*100.),0.5*fullrange(gnM2_53*100.),psym=5,col=cgcolor('red'),errcol=cgcolor('red'),thick=thick
oploterror,6485.,midrange(gnM2_65*100.),0.5*fullrange(gnM2_65*100.),psym=5,col=cgcolor('red'),errcol=cgcolor('red'),thick=thick
oploterror,8785.,midrange(gnM2_88*100.),0.5*fullrange(gnM2_88*100.),psym=5,col=cgcolor('red'),errcol=cgcolor('red'),thick=thick

dx = !x.crange[1]-!x.crange[0]  &       dy = !y.crange[1]-!y.crange[0]
x1 = !x.crange[0]+0.15*dx       &       y1 = !y.crange[0]+0.35*dy
!p.charsize=1.1

oplot,[x1,x1+0.05*dx],[y1,y1],col=cgcolor('black'),thick=thick
XYOUTS,[x1+0.07*dx],[y1-0.01*dy],"Fresh Witness Sample"

y1 = y1-0.05*dy
oploterror,[x1+0.025*dx],[y1],[0.02*dy],psym=4,col=cgcolor('blue'),errcol=cgcolor('blue')
XYOUTS,[x1+0.07*dx],[y1-0.01*dy],"M1 Measured"

y1 = y1-0.05*dy
oploterror,[x1+0.025*dx],[y1],[0.02*dy],psym=5,col=cgcolor('red'),errcol=cgcolor('red')
XYOUTS,[x1+0.07*dx],[y1-0.01*dy],"M2 Measured"

!p.charsize=1.3
plot,x*10.,y,psym=4,/xsty,ysty=3,xrange=[3000,10500],xtitle='Wavelength ('+cgsymbol("angstrom")+')', $
          ytitle='Reflectivity (%)',thick=thick,yrange=[70,100],title='Gemini South Mirror Reflectivity',/nodata
oplot,x2*10.,y2,col=cgcolor('black'),thick=thick ; GS Sample
oplot,x2*10.,y2*m1scale,col=cgcolor('blue'),thick=thick ; Scaled to worst M1
oplot,x2*10.,y2*m2scale,col=cgcolor('red'),thick=thick  ; Scaled to worst M2

oploterror,4715.,midrange(gsM1_47),0.5*fullrange(gsM1_47),psym=4,col=cgcolor('blue'),errcol=cgcolor('blue'),thick=thick
oploterror,5315.,midrange(gsM1_53),0.5*fullrange(gsM1_53),psym=4,col=cgcolor('blue'),errcol=cgcolor('blue'),thick=thick
oploterror,6515.,midrange(gsM1_65),0.5*fullrange(gsM1_65),psym=4,col=cgcolor('blue'),errcol=cgcolor('blue'),thick=thick
oploterror,8815.,midrange(gsM1_88),0.5*fullrange(gsM1_88),psym=4,col=cgcolor('blue'),errcol=cgcolor('blue'),thick=thick
oploterror,4685.,midrange(gsM2_47),0.5*fullrange(gsM2_47),psym=5,col=cgcolor('red'),errcol=cgcolor('red'),thick=thick
oploterror,5285.,midrange(gsM2_53),0.5*fullrange(gsM2_53),psym=5,col=cgcolor('red'),errcol=cgcolor('red'),thick=thick
oploterror,6485.,midrange(gsM2_65),0.5*fullrange(gsM2_65),psym=5,col=cgcolor('red'),errcol=cgcolor('red'),thick=thick
oploterror,8785.,midrange(gsM2_88),0.5*fullrange(gsM2_88),psym=5,col=cgcolor('red'),errcol=cgcolor('red'),thick=thick

dx = !x.crange[1]-!x.crange[0]  &       dy = !y.crange[1]-!y.crange[0]
x1 = !x.crange[0]+0.15*dx       &       y1 = !y.crange[0]+0.35*dy
!p.charsize=1.1
oplot,[x1,x1+0.05*dx],[y1,y1],col=cgcolor('black'),thick=thick
XYOUTS,[x1+0.07*dx],[y1-0.01*dy],'Fresh Witness Sample'

y1 = y1-0.05*dy
oploterror,[x1+0.025*dx],[y1],[0.02*dy],psym=4,col=cgcolor('blue'),errcol=cgcolor('blue')
XYOUTS,[x1+0.07*dx],[y1-0.01*dy],"M1 Measured"

y1 = y1-0.05*dy
oploterror,[x1+0.025*dx],[y1],[0.02*dy],psym=5,col=cgcolor('red'),errcol=cgcolor('red')
XYOUTS,[x1+0.07*dx],[y1-0.01*dy],"M2 Measured"

y1 = y1-0.05*dy
oplot,[x1,x1+0.05*dx],[y1,y1],col=cgcolor('blue'),thick=thick
XYOUTS,[x1+0.07*dx],[y1-0.01*dy],'Assumed M1 & M3 Reflectivity'

y1 = y1-0.05*dy
oplot,[x1,x1+0.05*dx],[y1,y1],col=cgcolor('red'),thick=thick
XYOUTS,[x1+0.07*dx],[y1-0.01*dy],'Assumed M2 Reflectivity'

!p.multi=0
if keyword_set(PRINT) then set_ps,/close,file=pfile,/landscape else stop
end

pro plot_slitloss, PRINT=print, DIR=dir
; Plot the fraction offlux within the IFU as a function of seeing
; Based on the FracSpa subroutine
if keyword_set(PRINT) then begin
   if not keyword_set(DIR) then dir = './'
   pfile = dir+'slit_loss.ps'
   set_ps,/open,file=pfile,/landscape
   setplot_paper
endif

readcol,'~/Idl/GHOST_ETC/RefData/slit_loss.dat',s,f,/silent
plot,s,f,xtitle='Seeing FWHM',ytitle='Fraction of flux in IFU',psym=-4,symsize=2,xsty=3,title='Fraction of Flux Enclosed by IFU'

if keyword_set(PRINT) then set_ps,/close,file=pfile,/landscape else stop
end

pro plot_throughput, PRINT=print, DIR=dir

if keyword_set(PRINT) then begin
   if not keyword_set(DIR) then dir = './'
   pfile = dir+'throughput.ps'
   set_ps,/open,file=pfile,/landscape
   setplot_paper
   !p.charthick=6
   !x.thick=6 &   !y.thick=6
endif

readcol,'~/Idl/GHOST_ETC/RefData/ghost_throughput_master.dat',lnm,cup,cub,specp,specb,/silent
thick1=8
thick2=11
plot,lnm*10.,specp*100.,xtitle='Wavelength ('+cgsymbol("angstrom")+')',ytitle='Throughput (%)',thick=thick,xrange=[3000,10500],/xsty,title='Spectrograph Throughput',/nodata,yrange=[0,100]
oplot,lnm*10.,specp*100.,col=cgcolor('blue'),linesty=1, thick=thick1
oplot,lnm*10.,specb*100.,col=cgcolor('red'),linestyle=1, thick=thick1

oplot,lnm*10.,cup*100.,col=cgcolor('blue'),linesty=2, thick=thick1
oplot,lnm*10.,cub*100.,col=cgcolor('red'),linesty=2, thick=thick1

oplot,lnm*10.,cup*100.*specp,col=cgcolor('blue'),thick=thick2
oplot,lnm*10.,cub*100.*specb,col=cgcolor('red'),thick=thick2

dx = !x.crange[1]-!x.crange[0]  &       dy = !y.crange[1]-!y.crange[0]
x1 = !x.crange[0]+0.1*dx       &       y1 = !y.crange[0]+0.9*dy
oplot,[x1-0.05*dx,x1-0.01*dx],[y1,y1]+0.015*dy,col=cgcolor('black'),linesty=2,thick=thick1
XYOUTS,[x1],[y1],"Cass Unit"
x1 = !x.crange[0]+0.1*dx       &       y1 = !y.crange[0]+0.85*dy
oplot,[x1-0.05*dx,x1-0.01*dx],[y1,y1]+0.015*dy,col=cgcolor('black'),linesty=1,thick=thick1
XYOUTS,[x1],[y1],"Spectrograph+CCD"
x1 = !x.crange[0]+0.1*dx       &       y1 = !y.crange[0]+0.8*dy
oplot,[x1-0.05*dx,x1-0.01*dx],[y1,y1]+0.015*dy,col=cgcolor('black'),linesty=0,thick=thick2
XYOUTS,[x1],[y1],"Combined"

x1 = !x.crange[0]+0.6*dx       &       y1 = !y.crange[0]+0.9*dy
XYOUTS,[x1],[y1],"Predicted",col=cgcolor('blue')
x1 = !x.crange[0]+0.76*dx       &       y1 = !y.crange[0]+0.9*dy
XYOUTS,[x1],[y1],"/",col=cgcolor('black')
x1 = !x.crange[0]+0.79*dx       &       y1 = !y.crange[0]+0.9*dy
XYOUTS,[x1],[y1],"Budgeted",col=cgcolor('red')
x1 = !x.crange[0]+0.76*dx       &       y1 = !y.crange[0]+0.85*dy
XYOUTS,[x1],[y1],"Requirements",col=cgcolor('green'),align=0.5

oplot,[3630,3900],[8,8],col=cgcolor('green'),thick=7
oplot,[3900,9000],[13.5,13.5],col=cgcolor('green'),thick=7
oplot,[9000,10000],[2.1,2.1],col=cgcolor('green'),thick=7
oplot,[4500],[27],col=cgcolor('green'),thick=5,psym=4

if keyword_set(PRINT) then set_ps,/close,file=pfile,/landscape else stop

end

pro plot_throughput2, PRINT=print, DIR=dir

if keyword_set(PRINT) then begin
   if not keyword_set(DIR) then dir = './'
   pfile = dir+'throughput.ps'
   set_ps,/open,file=pfile,/landscape
   setplot_paper
   !p.charthick=6
   !x.thick=6 &   !y.thick=6
endif

readcol,'~/Idl/GHOST_ETC/RefData/GHOST_throughput_FDR.dat',lnm,totTPpred, totTPbudg,/silent

thick1=8
thick2=11
plot,lnm*10.,totTPbudg*100.,xtitle='Wavelength ('+cgsymbol("angstrom")+')',ytitle='Throughput (%)',thick=thick,xrange=[3000,10500],/xsty,title='Spectrograph Throughput',/nodata,yrange=[0,55]
;oplot,lnm*10.,totTPpred*100.,col=cgcolor('red'), thick=thick1
oplot,lnm*10.,totTPbudg*100.,col=cgcolor('blue'), thick=thick1

dx = !x.crange[1]-!x.crange[0]  &       dy = !y.crange[1]-!y.crange[0]
x1 = !x.crange[0]+0.1*dx       &       y1 = !y.crange[0]+0.9*dy
;XYOUTS,[x1],[y1],"Predicted",col=cgcolor('red')
;y1 = y1 - 0.07*dy
XYOUTS,[x1],[y1],"Budgeted",col=cgcolor('blue')
y1 = y1 - 0.07*dy
XYOUTS,[x1],[y1],"Requirements",col=cgcolor('green');,align=0.5

oplot,[3630,3750],[8,8],col=cgcolor('green'),thick=7
oplot,[3750,4500],[12.7,12.7],col=cgcolor('green'),thick=7
oplot,[4500,9000],[13.6,13.6],col=cgcolor('green'),thick=7
oplot,[9000,9500],[5.5,5.5],col=cgcolor('green'),thick=7
oplot,[9500,10000],[2.1,2.1],col=cgcolor('green'),thick=7
oplot,[4500],[27],col=cgcolor('green'),thick=5,psym=4

if keyword_set(PRINT) then set_ps,/close,file=pfile,/landscape else stop

end

pro sr_vs_faint

  mags = range(18,23)
  nmags = n_elements(mags)
  snsr = fltarr(nmags)
  snf = fltarr(nmags)
  for i=0,nmags-1 do begin
     ghost_etc,action='S',F=mags[i],/ab,sn=sn,throughput='R',texp=3600,mode='SR'
     snsr[i] = sn
     ghost_etc,action='S',F=mags[i],/ab,sn=sn,throughput='R',texp=3600,mode='SF'
     snf[i] = sn
  endfor

  plot,mags,snsr,psym=-4
  oplot,mags,snf,psym=-2,col=cgcolor('red')
  stop
end
  

pro make_tables_pdf, file, title
  if n_params() lt 2 then title="GHOST Throughput"
  dir = FILE_DIRNAME(file)
  filename = FILE_BASENAME(file,'.dat') 
  readcol,file,f='(a,f,f,f,f,f,f)',mode,mag363,mag375,mag450,mag900,mag950,mag1000,skipline=1

  modes = ['Standard Resolution','Standard Resolution Faint','Beam Switching','Beam Switching Faint',$
           'High Resolution','High Resolution Faint','Precision Radial Velocity']
  endline = ['\\','\\','\\','\\','\\','\\','\\\hline']
  colors = ['Orange!30','Orange!30','Orange!30','Orange!30','Cyan!30','Cyan!30','Cyan!30']
  
openw,lun,dir+'/'+filename+'.tex',/get_lun
printf,lun,'\documentclass{article}'
printf,lun,'\usepackage[usenames,dvipsnames]{xcolor}'
printf,lun,'\usepackage{tcolorbox}'
printf,lun,'\usepackage{tabularx}'
printf,lun,'\usepackage{array}'
printf,lun,'\usepackage{colortbl}'
printf,lun,'\tcbuselibrary{skins}'
printf,lun,'\usepackage{multirow}'
printf,lun,'\pagestyle{empty}'
printf,lun,'\usepackage{hhline}'
printf,lun,''
;printf,lun,'\tcbset{tab1/.style={enhanced,fonttitle=\bfseries,fontupper=\normalsize\sffamily,'
;printf,lun,'colback=yellow!10!white,colframe=red!50!black,colbacktitle=Yellow!40!white,'
;printf,lun,'coltitle=black,center title}}'
;printf,lun,''
printf,lun,'\begin{document}'
printf,lun,''
;printf,lun,'\begin{tcolorbox}[tab1,tabularx={l|X|X|X|X|X|X|},title='+title+',boxrule=0.5pt,width=14cm]'
printf,lun,'\begin{tabular}{|l|c|c|c|c|c|c|} \hline'
printf,lun,'\multicolumn{7}{|c|}{'+title+'} \\ \hline'
printf,lun,'\rowcolor{Yellow!30} &  \multicolumn{6}{c|}{{\bf Wavelength}} \\ \hhline{|>{\arrayrulecolor{Yellow}}->{\arrayrulecolor{black}}*{6}-}'
printf,lun,'\rowcolor{Yellow!30} \multirow{-2}{*}{{\bf Modes}} & 363nm   & 375nm     & {\bf 450nm}    & 550nm & 900nm     & 950nm       \\\hline'
FF = '(f5.2)'
for i=0,n_elements(modes)-1 do begin
   printf,lun,'\rowcolor{'+colors[i]+'}{\bf '+modes[i] +'}'+ ' & '+STRING(mag363[i],F=FF)+ ' & '+STRING(mag375[i],F=FF)+' & {\bf '+STRING(mag450[i],F=FF)+'} & '+STRING(mag900[i],F=FF)+' & '+STRING(mag950[i],F=FF)+' & '+STRING(mag1000[i],F=FF) + endline[i]
endfor
printf,lun,'\end{tabular}'
;printf,lun,''
;printf,lun,'\end{tcolorbox}'
printf,lun,''
printf,lun,'\end{document}'

free_lun, lun

; Run Latex
CD,dir
spawn,'latex '+dir+'/'+filename+'.tex'
spawn,'dvips '+dir+'/'+filename+'.dvi'
spawn,'ps2pdf '+dir+'/'+filename+'.ps'
spawn,'pdfcrop --margins 10 '+dir+'/'+filename+'.pdf'
file_move,dir+'/'+filename+'-crop.pdf',dir+'/'+filename+'.pdf',/over
file_delete,dir+'/'+filename+'.tex'
file_delete,dir+'/'+filename+'.dvi'
file_delete,dir+'/'+filename+'.ps'
file_delete,dir+'/'+filename+'.aux'
file_delete,dir+'/'+filename+'.log'
end

pro make_tables2_pdf, file, title
  if n_params() lt 2 then title="GHOST Throughput"
  dir = FILE_DIRNAME(file)
  filename = FILE_BASENAME(file,'.dat') 
  readcol,file,f='(a,f,f,f,f,f,f)',mode,mag363,mag375,mag450,mag900,mag950,mag1000,skipline=1

  modes = ['Two Target','Two Target Faint','Two Target V.Faint','Beam Switching','Beam Switching Faint','Beam Switching V.Faint',$
           'High Resolution','High Resolution Faint','Precision Radial Velocity']
;  endline = ['\\','\\','\\','\\','\\','\\','\\','\\','\\\hline']
  endline = ['\\\hline','\\\hline','\\\hline','\\\hline','\\\hline','\\\hline','\\\hline','\\\hline','\\\hline']
  colors = ['Orange!30','Orange!30','Orange!30','Orange!30','Orange!30','Orange!30','Cyan!30','Cyan!30','Cyan!30']
  
openw,lun,dir+'/'+filename+'.tex',/get_lun
printf,lun,'\documentclass{article}'
printf,lun,'\usepackage[usenames,dvipsnames]{xcolor}'
printf,lun,'\usepackage{tcolorbox}'
printf,lun,'\usepackage{tabularx}'
printf,lun,'\usepackage{array}'
printf,lun,'\usepackage{tabu}'
printf,lun,'\usepackage{colortbl}'
printf,lun,'\tcbuselibrary{skins}'
printf,lun,'\usepackage{multirow}'
printf,lun,'\pagestyle{empty}'
printf,lun,'\usepackage{hhline}'
printf,lun,'\let\xvline\vline'
printf,lun,'\begin{document}'
printf,lun,''
printf,lun,'\begin{tabular}{|l|c|c|c|c|c|c|} \hline'
;;printf,lun,'\begin{tabu}{|l|[2pt]c|c|c|c|c|c|} \hline'
printf,lun,'\rowcolor{Yellow!40} \multicolumn{7}{|c|}{'+title+'} \\ \hline'
printf,lun,'\rowcolor{Yellow!40} &  \multicolumn{6}{c|}{{\bf Wavelength}} \\ \hhline{|>{\arrayrulecolor{Yellow!40}}->{\arrayrulecolor{black}}*{6}-}'
printf,lun,'\rowcolor{Yellow!40} \multirow{-2}{*}{{\bf Modes}} & 363nm   & 375nm     & {\bf 450nm}    & 550nm & 900nm     & 950nm       \\\hline'
;;printf,lun,'\rowcolor{Yellow!40} \multirow{-2}{*}{{\bf Modes}} & 363nm   & 375nm     & {\bf 450nm}    & 550nm & 900nm     & 950nm       \\\tabucline[2pt]{-}'
FF = '(f5.2)'
for i=0,n_elements(modes)-1 do begin
   printf,lun,'\rowcolor{'+colors[i]+'}{\bf '+modes[i] +'}'+ ' & '+STRING(mag363[i],F=FF)+ ' & '+STRING(mag375[i],F=FF)+' & {\bf '+STRING(mag450[i],F=FF)+'} & '+STRING(mag900[i],F=FF)+' & '+STRING(mag950[i],F=FF)+' & '+STRING(mag1000[i],F=FF) + endline[i]
endfor
printf,lun,'\end{tabular}'
;;printf,lun,'\end{tabu}'
printf,lun,''
printf,lun,'\end{document}'

free_lun, lun

; Run Latex
CD,dir
spawn,'latex '+dir+'/'+filename+'.tex',latexout
spawn,'dvips '+dir+'/'+filename+'.dvi',dvipsout
spawn,'ps2pdf '+dir+'/'+filename+'.ps'
spawn,'pdfcrop --margins 10 '+dir+'/'+filename+'.pdf'
file_move,dir+'/'+filename+'-crop.pdf',dir+'/'+filename+'.pdf',/over
file_delete,dir+'/'+filename+'.tex'
file_delete,dir+'/'+filename+'.dvi'
file_delete,dir+'/'+filename+'.ps'
file_delete,dir+'/'+filename+'.aux'
file_delete,dir+'/'+filename+'.log'
end

pro make_etc_plots
cleanplot
dir = '~/Idl/GHOST_ETC/Outputs/'
plot_throughput2,/print, DIR=dir
run_ghost_etc,/print, DIR=dir
SPAWN,'convert '+dir+'ghost_comparison.pdf '+dir+'ghost_comparison.png'
plot_slitloss,/print,DIR=dir
plot_reflectivity, /print, DIR=dir
plot_extinct, /print,DIR=dir
toto=sky(4500,/print,dir=dir)
end
