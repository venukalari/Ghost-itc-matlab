PRO FracSpaGauss, Mode, Seeing, lmu, Dlens, FA
;
; DESCRIPTION
; Compute the fraction of the total flux enclosed in the GHOST IFU for given spectral
; mode, seeing FWHM, and wavelength (wavelength dependence not
; implemented, so this is currently a dummy argument)
;
; INPUTS
; Mode - Instrument mode:
;        SR - Standard Resolution
;        HR - High Resolution
; Seeing  - Seeing FWHM
; lmu  - Wavelength (in micron) of requested fractional flux (not used)
; Dlens - lens diameter, flat-to-flat in arcsec
; OUTPUTS
; FA   - Returned fraction of the flux enclosed within aperture, at given wavelength
;
; COMMENTS
;
; HISTORY:
; - Written by R.McDermid
;
;=======================================================================================

; First Create the PSF image. Assume 0.01" pixels. These
; parameters are OK for typica seeing values (0.2"-1.2" FWHM)
  scale = 0.01                  ; arcsec/pixel scale for PSF image
  npix = 500                    ; Size of PSF image
  x = scale * [[FINDGEN(npix)-npix/2] # REPLICATE(1,npix)]
  y = scale * [REPLICATE(1,npix) # [FINDGEN(npix)-npix/2]]
  r = SQRT(x^2+y^2)

  mof = PSF_GAUSSIAN(npix=npix,ndim=2,fwhm=Seeing/scale,/norm,centroid=npix/2)

; ####################
; Set up the IFU field
; ####################
  ; Hexagon parameters. a = half flat-to-flat = apothem, s = side length
  a = Dlens/2.0
  s = sqrt(3.)*a/2.

; Set up lens centres
  ; Standard Rsolution modes
  if Mode eq "SR" or Mode eq "SF" or Mode eq "SVF" or Mode eq "BS" or Mode eq "BF" or Mode eq "BVF" then begin
     xcen = [0.,0.,-2.*s,-2.*s,0.,2.*s,2.*s]
     ycen = [0.,2.*a, a, -1.*a, -2.*a, -1.*a, a]
  endif else if Mode eq "HR" or Mode eq "HF" or Mode eq "PRV" then begin
  ; High Resolution modes
     xcen = [0.,0.,-2.*s,-2.*s,0.,2.*s,2.*s,$ ; inner ring
             0., -2.*s, -4.*s, -4.*s, -4.*s, -2.*s, $
             0., 2.*s, 4.*s, 4.*s, 4.*s, 2.*s]
     ycen = [0.,2.*a, a, -1.*a, -2.*a, -1.*a, a, $
             4.*a, 3.*a, 2.*a, 0., -2.*a, -3.*a, -4.*a, $
             -3.*a, -2.*a, 0., 2.*a, 3.*a]
  endif

  ; Now cycle through each lens, finding the regions of the PSF it covers
  nhex = n_elements(xcen)
  width = Dlens
  for i=0,nhex-1 do begin
     xx = x + xcen[i]
     yy = y + ycen[i]
     w = where(yy lt width/2. and yy gt (-width/2.) and $
               yy lt (width-sqrt(3.)*xx) and yy gt (-width+sqrt(3.)*xx) and $
               yy lt (width+sqrt(3.)*xx) and yy gt (-width-sqrt(3.)*xx))

     print=0 ; Collect boundary coordinates for plotting?
     if print eq 1 then begin
        bpts = find_boundary(w,xsize=npix,ysize=npix)
        if i eq 0 then xb=scale*reform(bpts[0,*]) + min(x) else xb=[xb,reform(scale*bpts[0,*]) + min(x)]
        if i eq 0 then yb=scale*reform(bpts[1,*]) + min(y) else yb=[yb,reform(scale*bpts[1,*]) + min(y)]
        if i eq 0 then lens=REPLICATE(i,n_elements(bpts[0,*])) else lens = [lens,REPLICATE(i,n_elements(bpts[0,*]))]
     endif

     ; Accumulate the PSF pixels covered by the lenses
     if i eq 0 then wtot = w else wtot = [wtot,w]
  endfor

; Comput enclosed flux. The PSF function is already normalised, so it
; is just the total of the PSF pixels covered by the lenses
  FA = TOTAL(mof[wtot])

  if print eq 1 then begin      ; Show boundaries of lenslets on a plot
     loadct,3
     image_plot,alog10(mof/max(mof)),x,y,/asp,title=STRING(Seeing,F='(f4.2)')+'" seeing: EE='+STRING(FA,F='(f4.2)'),/nodata,/noerase
     for i=0,nhex-1 do begin
        ind = where(lens eq i)
        oplot,xb[ind],yb[ind],thick=3,col=cgcolor('green')
     endfor
     stop
  endif
  
end

;
;=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
;

PRO FracSpa, Mode, Seeing, lmu, Dlens, FA
;
; DESCRIPTION
; Compute the fraction of the total flux enclosed in the GHOST IFU for given spectral
; mode, seeing FWHM, and wavelength (wavelength dependence not
; implemented, so this is currently a dummy argument)
;
; INPUTS
; Mode - Instrument mode:
;        SR - Standard Resolution
;        HR - High Resolution
; Seeing  - Seeing FWHM
; lmu  - Wavelength (in micron) of requested fractional flux (not used)
; Dlens - lens diameter, flat-to-flat in arcsec
; OUTPUTS
; FA   - Returned fraction of the flux enclosed within aperture, at given wavelength
;
; COMMENTS
; Function is normalised by the formal integral, so fraction is not
; dependent on PSF image size.
;
; HISTORY:
; - Written by R.McDermid, using moffat.pro from M.Ireland. 20-Nov-2014
;
;=======================================================================================

; First Create the PSF image. Assume 0.01" pixels. These
; parameters are OK for typica seeing values (0.2"-1.2" FWHM)
  scale = 0.01                  ; arcsec/pixel scale for PSF image
  npix = 500                    ; Size of PSF image
  x = scale * [[FINDGEN(npix)-npix/2] # REPLICATE(1,npix)]
  y = scale * [REPLICATE(1,npix) # [FINDGEN(npix)-npix/2]]
  r = SQRT(x^2+y^2)

  mof = moffat(r/scale,Seeing/2./scale)

; ####################
; Set up the IFU field
; ####################
  ; Hexagon parameters. a = half flat-to-flat = apothem, s = side length
  a = Dlens/2.0
  s = sqrt(3.)*a/2.

; Set up lens centres
  ; Standard Rsolution modes
  if Mode eq "SR" or Mode eq "SF" or Mode eq "SVF" or Mode eq "BS" or Mode eq "BF" or Mode eq "BVF" then begin
     xcen = [0.,0.,-2.*s,-2.*s,0.,2.*s,2.*s]
     ycen = [0.,2.*a, a, -1.*a, -2.*a, -1.*a, a]
  endif else if Mode eq "HR" or Mode eq "HF" or Mode eq "PRV" then begin
  ; High Resolution modes
     xcen = [0.,0.,-2.*s,-2.*s,0.,2.*s,2.*s,$ ; inner ring
             0., -2.*s, -4.*s, -4.*s, -4.*s, -2.*s, $
             0., 2.*s, 4.*s, 4.*s, 4.*s, 2.*s]
     ycen = [0.,2.*a, a, -1.*a, -2.*a, -1.*a, a, $
             4.*a, 3.*a, 2.*a, 0., -2.*a, -3.*a, -4.*a, $
             -3.*a, -2.*a, 0., 2.*a, 3.*a]
  endif

  ; Now cycle through each lens, finding the regions of the PSF it covers
  nhex = n_elements(xcen)
  width = Dlens
  for i=0,nhex-1 do begin
     xx = x + xcen[i]
     yy = y + ycen[i]
     w = where(yy lt width/2. and yy gt (-width/2.) and $
               yy lt (width-sqrt(3.)*xx) and yy gt (-width+sqrt(3.)*xx) and $
               yy lt (width+sqrt(3.)*xx) and yy gt (-width-sqrt(3.)*xx))

     print=0 ; Collect boundary coordinates for plotting?
     if print eq 1 then begin
        bpts = find_boundary(w,xsize=npix,ysize=npix)
        if i eq 0 then xb=scale*reform(bpts[0,*]) + min(x) else xb=[xb,reform(scale*bpts[0,*]) + min(x)]
        if i eq 0 then yb=scale*reform(bpts[1,*]) + min(y) else yb=[yb,reform(scale*bpts[1,*]) + min(y)]
        if i eq 0 then lens=REPLICATE(i,n_elements(bpts[0,*])) else lens = [lens,REPLICATE(i,n_elements(bpts[0,*]))]
     endif

     ; Accumulate the PSF pixels covered by the lenses
     if i eq 0 then wtot = w else wtot = [wtot,w]
  endfor

; Comput enclosed flux. The PSF function is already normalised, so it
; is just the total of the PSF pixels covered by the lenses
  FA = TOTAL(mof[wtot])

  if print eq 1 then begin      ; Show boundaries of lenslets on a plot
     loadct,3
     image_plot,alog10(mof/max(mof)),x,y,/asp,title=STRING(Seeing,F='(f4.2)')+'" seeing: EE='+STRING(FA,F='(f4.2)'),/nodata,/noerase
     for i=0,nhex-1 do begin
        ind = where(lens eq i)
        oplot,xb[ind],yb[ind],thick=3,col=cgcolor('green')
     endfor
     stop
  endif
  
end


FUNCTION AB2Flux, AB, Lbda, SI=SI, Arcsec=Arcsec

; AB is the AB magnitude
; Lbda is in A
; Flux will be in erg/s/cm^2 if option SI is not set
; Flux will be in kg/m/s^3 if option SI is set
; Flux will be in arcsec-2 if option Sec is set
; Flux will be in A-1 if option A is set

c = 299792458.0 ; in m/s
l = Lbda*1.e-10 ; in m
asec = 1.0/206265.0 ; arcsec in radian


if ~Keyword_set(SI) then begin
    Flux = 10^(-0.4*(AB+48.60)-10)*c/l^2
endif else begin
    Flux = 10^(-0.4*(AB+48.60)-3)*c/l^2
endelse

if Keyword_set(Arcsec) then begin
    Flux = Flux/asec^2
endif

return, Flux

end

FUNCTION Flux2AB, Flux, Lbda, SI=SI, Arcsec=Arcsec

; Lbda is in A
; Flux is in erg/s/cm2 if option SI is not set
; Flux is in kg/m/sec3 if option SI is set
; Flux is in arcsec-2 if option Sec is set
; Flux is in A-1 if option A is set
; result is AB magnitude

c = 299792458.0 ; in m/s
l = Lbda*1.e-10 ; in m
asec = 1.0/206265.0 ; arcsec in radian

F = Flux
if Keyword_set(Arcsec) then begin
    F = Flux*asec^2
endif

if ~Keyword_set(SI) then begin
    AB = -73.60 -2.5*alog10(l^2*F/c)
endif else begin
    AB = -56.10 -2.5*alog10(l^2*F/c)
endelse

return, AB

end

function AB2Vega, mab, l
; Vega reference spectrum from CALSPEC (Bohlin 2007):
readcol,'~/Idl/GHOST_ETC/RefData/alpha_lyr_stis_004.ascii',lvega,fvega,/silent
mabVega = FLUX2AB(interpol(fvega,lvega,l),l)
return, mab-mabVega
end

;
;=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
;

PRO OnlyOne, a, b, c
;
; DESCRIPTION
; Short routine to check for self-consistency of input arrays. Only one of the
; three input arrays should have more than one element.
;
if (n_elements(a) gt 1 and (n_elements(b) gt 1 or n_elements(c) gt 1) or $
    n_elements(b) gt 1 and (n_elements(a) gt 1 or n_elements(c) gt 1) or $
    n_elements(c) gt 1 and (n_elements(a) gt 1 or n_elements(b) gt 1)) then begin
    print, 'ERROR: Only one array permitted among sn,flux,texp,wave'
    return
endif
end

function sky, lmu, PRINT=print, DIR=dir, BRIGHT=bright
;
; Return the sky flux from Hanuschik 2003.
; Optionally make a plot, via PRINT keyword.
; BRIGHT keyword scales the Hanuschik continuum flux according to
  
  readcol,'~/Idl/GHOST_ETC/RefData/hanuschik_cont_tab.dat',x1,x2,y,/silent
  res = poly_fit([x1,x2]/1.e4,[y,y],3,yfit=yfit)
  sky = res[3]*lmu^3 + res[2]*lmu^2 + res[1]*lmu + res[0]

  if keyword_set(BRIGHT) then begin
     wv = [3656.,4353.,5477.,6349.,8797.,10000.]
     dk = [22.5,22.8,21.8,20.9,19.7,19.3]
     brt = [19.3,19.1,18.5,18.2,17.3,16.9]
     delt = dk-brt
     scl = 10.^(0.4*delt)
     res2 = poly_fit(wv/1.e4,scl,3,yfit=yfit)
     brtcor = res2[3]*lmu^3 + res2[2]*lmu^2 + res2[1]*lmu + res2[0]
     sky = sky*brtcor
  endif
  
  if keyword_set(PRINT) then begin
     if not keyword_set(DIR) then dir = './'
     pfile = dir+'hanuschick.ps'
     set_ps,/open,file=pfile,/landscape
     setplot_paper

     plot,x1,y,/nodata,xrange=[3000,10500],/xsty,xtitle='Wavelength ' + $
          '('+cgsymbol("angstrom")+')', ytitle='Flux (10!U-16!N erg/s/cm!U2!N/'+cgsymbol("angstrom")+'/arcsec!U2!N)',$
          title='Sky Continuum Brightness',yrange=[0,100]
     for i=0,n_elements(x1)-1 do begin
        oplot,[x1[i],x2[i]],[y[i],y[i]],col=cgcolor('magenta'),thick=15
     endfor

     xplot = range(3000.,11000.)/1.e4
     yplot = res[3]*xplot^3 + res[2]*xplot^2 + res[1]*xplot + res[0]
     if keyword_set(BRIGHT) then yplot = res2[3]*xplot^3 + res2[2]*xplot^2 + res2[1]*xplot + res2[0]

     oplot,xplot*1.e4,yplot,col=cgcolor('black'),thick=5

     oplot,[3400,3800],[0.05,0.05],thick=5
     XYOUTS,[3900],[0.048],STRING(REVERSE(res),FORMAT='("y = ",f5.2,"x!U3!N",f+5.2,"x!U2!N",f+5.2,"x",f+5.2)')
     oplot,[3400,3800],[0.04,0.04],thick=15,col=cgcolor('magenta')
     XYOUTS,[3900],[0.038],"Mean sky flux (Hanuschik, 2003)"
     set_ps,/close,file=pfile,/landscape
     print,res
  endif

  return,sky*1.e-16

end

;
;=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
;

PRO Comp_SN, sn, flux, texp, ErrFrac, Nexp, Mode, nspec, nspat, npix, npix_org, SkyFact, FA, ScatLight
;
; DESCRIPTION
; Computates S/N for the given flux and integration time
;
; INPUTS
; sn - Returned S/N value
; flux - Flux of input source
; texp - Exposure time in seconds
; ErrFrac - Returned percentage contributions to total variance [obj,sky,RDN,DC]
; Nexp - Number of exposures, each with length texp
; Mode - GHOST observing mode
; nspec - Number of SPECTRAL pixels summed in computation
; nspat - Number of SPATIAL elements summed in computation
; npix - number of binned object pixels (for read noise)
; npix_org - number of unbinned object pixels (for dark current)
;
;========================================================================================

Common Share_ghost

if (Disp eq 1) then print,'------------------------------------------------------------------------'

;###############################################################
; Ponctual & Continuum Source
; Here we must use the fractional flux only in the spatial domain.
;###############################################################
;
; Working in flux or AB magnitudes?
;
if flag_AB eq 1 then FObj = AB2Flux(flux, lbda, /SI) $ ; convert AB in flux by A in SI unit
                else FObj = flux*1e7                   ; convert flux in SI unit
;
;
; COMPUTE Effective collecting power for Object (KO) and Sky (KS) entering the GHOST
; science object IFU. For the Object, this is the product of number of
; spectral pixels (nspec), spectral pixel size (DS), fraction
; of included flux (FA), throughput (Tghost), collecting area (GEM),
; extinction, and flux of given wavelength (l/hc).
; For the Sky, this is the product of nspec, DS, area of sky observed
; (DASky), Tghost, GEM and flux.
;
KO = double(nspec*DS*FA*Tghost*GEM*Extinct*l/hc)
KS = double(nspec*DS*DASky*Tghost*GEM*l/hc)
;
; COMPUTE counts from object and sky, applying collecting power to
; flux rate of source and integration time
;
ObjCounts = KO*FObj*texp
SkyCounts = KS*SkyFact*FSky*texp
;
; COMPUTE noise from detector
; NOTE: The sky is subtracted via peripheral sky fibers. The Poisson
; and read noise due to this is included by scaling the effective
; number of pixels (for read noise and dark current) and sky counts by
; the SkyFact factor, which is based on the number and area of sky
; fibers relative to the IFU.
;
RDN       = SkyFact*npix*RN^2
DKN       = SkyFact*npix_org*DC*texp
;
; COMPUTE S/N using standard formula and input object flux. This includes the noise
; contribution from Poisson statistics, sky component and CCD noise (read noise &
; dark current). ScatLight accounts crudely for scattered light. This
; scales the effective counts from the source only, and includes this
; additional flux in the Poisson noise calculation (so only adding
; noise, not signal).
;
sn = sqrt(Nexp)*ObjCounts/sqrt(ObjCounts*ScatLight + SkyCounts + RDN + DKN)
;
; Assemble values expresing the fractional errors
;
ErrFrac = [ObjCounts, SkyCounts, RDN, DKN]
ErrTot = TOTAL(ErrFrac)
ErrFrac = 100.0 * ErrFrac / ErrTot
;
; OUTPUT numbers if running in verbose mode.
;
if (Disp eq 1) then begin
   print,format='("Summed spaxels: ",i4, " fraction of object flux ",f4.2)', npix, FA
   if (nspec gt 1) then begin
      resolv = l/(nspec*DS)
      print,format='("Low R: ", f8.1," : ", f5.1," spectral pixels sum up")', resolv, nspec
   endif
   if flag_AB eq 1 then begin
      print,format='("Lbda: ",f4.2," um Object AB magnitude: ",f5.2)', lmu, flux
   endif else begin
      print,format='("Lbda: ",f4.2," um Object flux: ",e12.4," erg/s/cm2/A")',$
            lmu, flux
   endelse
   print_noise, ObjCounts, SkyCounts, RDN, DKN
endif

;####################
; Final output values
;####################

if (Disp eq 1) then print,format='("S/N: ",f8.3, " in ",i2, " exposure(s) of ", f8.1," sec")', $
                  sn, Nexp, texp

return

end

;
;=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
;

PRO Comp_T, sn, flux, texp, ErrFrac, Nexp, Mode, nspec, nspat, npix, npix_org, SkyFact, FA, ScatLight
;
; DESCRIPTION
; Compute integration time given flux and S/N.
; For more complete comments, see also Comp_SN above.

Common Share_ghost

if (Disp eq 1) then print,'------------------------------------------------------------------------'

if flag_AB eq 1 then FObj = AB2Flux(flux, lbda, /SI) $ ; convert AB in flux by A in SI unit
    else FObj = flux*1e7                               ; convert flux in SI unit
;
; COMPUTE Effective collecting power for Object (KO) and Sky (KS)
;
KO = double(nspec*DS*FA*Tghost*GEM*Extinct*l/hc)
KS = double(nspec*DS*DASky*Tghost*GEM*l/hc)
;
; Here the standard formula for computing S/N is rearranged for exposure time,
; resulting in a quadratic formula with the following coefficients:
;
a = Nexp*(KO*Fobj/sn)^2
b = ScatLight*KO*FObj + KS*SkyFact*FSky + SkyFact*npix_org*DC
c = SkyFact*npix*RN^2
;
; COMPUTE EXPOSURE TIME using quadratic equation
;
texp = (b + sqrt(b^2+4*a*c))/(2*a)
;
;
; Compute variance contributions
;
ObjCounts = KO*FObj*texp
SkyCounts = KS*SkyFact*FSky*texp
RDN = SkyFact*npix*RN^2
DKN = SkyFact*npix_org*DC*texp
ErrFrac = [ObjCounts, SkyCounts, RDN, DKN]
ErrTot = TOTAL(ErrFrac)
ErrFrac = 100.0 * ErrFrac / ErrTot
;
; OUTPUT numbers if running in verbose mode.
;
if (Disp eq 1) then begin
   if (nspec gt 1) then begin
      resolv = l/(nspec*DS)
      print,format='("Low R: ", f8.1," : ", f5.1," spectral pixels sum up")', resolv, nspec
   endif
   print,format='("Summed spaxels: ",i4, " fraction of object flux ",f4.2)', npix, FA
   if flag_AB eq 1 then begin
      print,format='("Lbda: ",f4.2," um Object AB magnitude: ",f5.2)', lmu, flux
   endif else begin
      print,format='("Lbda: ",f4.2," um Object flux: ",e12.4," erg/s/cm2/A")',$
            lmu, flux
   endelse
   print_noise, ObjCounts, SkyCounts, RDN, DKN              
endif


;####################
; Final output values
;####################

if (Disp eq 1) then print,format='("S/N: ",f8.3, " in ",i2, " exposure(s) of ", f8.1," sec")', $
                  sn, Nexp, texp

return

end

;
;=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
;

PRO Comp_F, sn, flux, texp, ErrFrac, Nexp, Mode, nspec, nspat, npix, npix_org, SkyFact, FA, ScatLight
;
; DESCRIPTION
; Compute flux given integration time and S/N
; For more complete comments, see also Comp_SN above.
;
Common Share_ghost

if (Disp eq 1) then print,'------------------------------------------------------------------------'
;
; COMPUTE Effective collecting power for Object (KO) and Sky (KS)
;
KO = double(nspec*DS*FA*Tghost*GEM*Extinct*l/hc)
KS = double(nspec*DS*DASky*Tghost*GEM*l/hc)
;
; As for Comp_T, solve the standard formula for S/N in terms of the object
; flux, giving a quadratic formula to solve.
;
a = Nexp*(KO*texp/sn)^2
b = KO*texp*ScatLight
c = KS*SkyFact*FSky*texp + SkyFact*npix_org*DC*texp + SkyFact*npix*RN^2
;
; COMPUTE the flux
;
Fobj = (b + sqrt(b^2+4*a*c))/(2*a)
;
; Give in requested system
;
if flag_AB eq 1 then flux = Flux2AB(Fobj, lbda, /SI) $ ; convert flux by A in SI unit in AB mag
    else flux = Fobj*1e-7                              ; convert flux in erg/s/cm-2/A
;
; Compute variance contributions
;
ObjCounts = KO*FObj*texp
SkyCounts = KS*SkyFact*FSky*texp
RDN = SkyFact*npix*RN^2
DKN = SkyFact*npix_org*DC*texp
ErrFrac = [ObjCounts, SkyCounts, RDN, DKN]
ErrTot = TOTAL(ErrFrac)
ErrFrac = 100.0 * ErrFrac / ErrTot
;
; OUTPUT numbers if running in verbose mode.
;
if (Disp eq 1) then begin
   if (nspec gt 1) then begin
      resolv = (l/(nspec*DS))
      print,format='("Low R: ", f10.1," : ", f5.1," spectral pixels sum up")', resolv, nspec
   endif
   print,format='("Summed spaxels: ",i4, " fraction of object flux ",f4.2)', npix, FA
   if flag_AB eq 1 then begin
      print,format='("Lbda: ",f4.2," um Object AB magnitude: ",f5.2)', lmu, flux
   endif else begin
      print,format='("Lbda: ",f4.2," um Object flux: ",e12.4," erg/s/cm2/A")',$
            lmu, flux
   endelse
   print_noise, ObjCounts, SkyCounts, RDN, DKN
endif


;####################
; Final output values
;####################

if (Disp eq 1) then print,format='("S/N: ",f8.3, " in ",i2, " exposure(s) of ", f8.1," sec")', $
                  sn, Nexp, texp

return

end

;
;=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
;

PRO print_noise, VObj, VSky, VRN, VDC
;
; DESCRIPTION
; Compute and print break down of noise variance contributions coming from the
; different components, all given as a fraction of the total.
;
; INPUTS
; VObj - Variance of object flux
; VSky - Variance of sky flux
; VRN - Variance of read noise
; VDC - Variance of dark current
;
;===================================================================================

Vtot = VObj + VSky + VRN + VDC

print, format='("Noise variance statistics - Obj:",f5.1,"% Sky:",f5.1,"% Detector:",f5.1,"% [RN:",f5.1, "% DC:",f5.1,"]")', $
    100*VObj/Vtot, 100*VSky/Vtot, 100*(VRN+VDC)/Vtot, 100*VRN/Vtot, 100*VDC/Vtot
print, format='("Noise variance statistics - Obj:",f7.1,"e- Sky:",f7.1,"e- Detector:",f7.1,"e- [RN:",f7.1, "e- DC:",f7.1,"e-]")', $
    VObj, VSky, (VRN+VDC), VRN, VDC

end


;
;=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
;


PRO GHOST_ETC, Action=Action, Sn=sn, Flux=flux, Texp=Texp, AB=AB, Wave=Wave, Nexp=Nexp, $
             Seeing=Seeing, Airmass=Airmass, Mode=Mode, ErrFrac=ErrFrac, Nspec=Nspec, Nspat=Nspat, $
             Quiet=Quiet, Help=Help, Throughput=Throughput, Site=Site, BRIGHT=bright, PORT=port
;
; TITLE
; ghost_etc
;
; DESCRIPTION
; Routine to compute, in a flexible way, the signal-to-noise ratio (S/N) of a GHOST
; observation of a point continuum source, given different seeing conditions, wavelengths, integration
; times and so on.
;
; INPUTS
; Action  - Determine S/N, exposure time or object flux (S/T/F default S)
; Sn      - Input/output S/N as appropriate
; Flux    - Input/output object flux in erg/s/cm2 (by A if continuum source) as appropriate
; Texp    - Integration time in sec of one exposure (default 1 hour)
; AB      - If set, then Flux is assumed to be in AB magnitudes
; Wave    - Reference wavelength (in A) at which to perform calculation (defaut I band 7900 A)
; Nexp    - Number of exposure (defaul 1)
; Seeing  - Seeing FWHM for a Moffat (beta=4) profile
; Airmass - Observation airmass (defaulted to 1)
; Mode    - GHOST mode (SR, SF, BS, BF, HR, HF, PRV)
; ErrFrac - Output percentage variance contributions from: object, sky, read noise and dark current
; Nspec   - Number of spectral pixel to sum up (defaulted to 1), only in continuum
; Nspat   - Number of spatial pixel to sum up (defaulted to 1), only in extended
; Quiet   - If set, no print and no plots
; Site    - Use extinction data for Gemini South ("GS") or North ("GN")
; Bright  - Make rough estimate of increased continuum background due to full moon
; Port    - Select port on the ISS. 1=upward looking, 0=side facing (default). If "S" throughput is selected, the short fiber forces port=1
;
; HISTORY
; V1.0 - Adapted from MUSE ETC, V1.7 (2007), written by R. McDermid and
;        R. Bacon.
; V2.0 - A few important changes:
;          - Corrected the slit length to 75 pixesl in SR/HR modes. 
;          - Changed resolution element to 3.5 pix SR, 0.6x3.5 pix HR
;          - Include option for different sky background
; V3.0 - Removed unused code for other source types
;      - Crude scattered light prescription is now included at top level
; V4.0 - Restructured the code for a cleaner distribution via GitHub.
;      - Added 'Very Faint' modes to incorporate 1x4 and 1x8 modes
;      - Updated Gemini reflectivity to Gemini spreadsheet data
; V5.0 - Updated assumed telescope aperture based on engineering data
;        from Gemin (Tom Hayward)
;      - Updated to use FDR throughput predictsio fom Ross, which
;        include some meaused elements (detector QE, FBPI fibre
;        measured throughput)
; V5.1 - Added 'Port' keyword and 'S' throughput to explore shorter fiber and
;        fixed upward port installation.
; V5.2 - Updated with new throughput data using measured Echelle performance.
;      - Removed 'Predicted' performance
;
         version = '5.2 - 11/10/16'
;
;========================================================================================

;
; Principle common block of variables
;
Common Share_ghost, lbda, l, lmu, Extinct, Tghost, FSky, flag_AB, hc, GEM, DS, RN, DC, DAObj, DASky, Disp


;                          ################################
;                          # SET DEFAULTS/KEYWORDS/INPUTS #
;                          ################################

; Version
if not Keyword_set(quiet) then print,'GHOST_ETC.pro version', version

; Parameter description
if Keyword_set(help) then begin
    print,'Usage GHOST_ETC, action=S|F|T, flux=flux, sn=sn, texp=texp, /ab, wave=7900, nexp=1,'
    print,'               Seeing=0.8, airmass=1, mode=WFM|NFM, nspec=1, nspat=1'
    return
endif

if not keyword_set(THROUGHPUT) then throughput = 'R' ; Assume PDR Requirements by default
if not Keyword_set(Action) then Action = 'S'   ; default action = Compute S/N
if not Keyword_set(Flux) then Flux = 1.e-18    ; default flux
if not Keyword_set(SN) then SN = 30.0          ; default S/N
if not Keyword_set(Texp) then Texp = 3600.0    ; default Integration time in sec
if not Keyword_set(Wave) then Wave = 4500.0    ; default wavelength in Angstrom
if not Keyword_set(Nexp) then Nexp = 1         ; default number of exposure
if not Keyword_set(Seeing) then Seeing = 0.8   ; default seeing FWHM (Moffat, beta=4.0)
if not Keyword_set(Airmass) then Airmass = 1.0 ; default airmass (zenith)
if not Keyword_set(nspec) then nspec = 1.      ; default number of spectrum pixel to sum up
if not Keyword_set(nspat) then nspat = 1.      ; default number of spatial pixel to sum up
if not Keyword_set(Mode) then Mode = 'SR'      ; default instrument mode
if Keyword_set(AB) then flag_AB = 1 else flag_AB = 0 ; check if AB magnitude flag is set
if not Keyword_set(SITE) then Site = 'GS'      ; GS = Cerro Paranal, GN = Mauna Kea
if not Keyword_set(PORT) then Port = 0         ; Deafult to side port, but can be over-ridden for non-short fiber throughput
if THROUGHPUT eq "S" then Port=1               ; Short fiber througput implies the upward port

; Ensure inputs are upper case
Mode   = strupcase(strtrim(Mode, 2))
Site   = strupcase(strtrim(Site, 2))

;####################################################################
; Check which mode we are in. This is chosen from:
;  S - Return the S/N for given input source and exposure time
;  T - Return the required exposure time given source and target S/N
;  F - Return the limiting flux for a demand S/N and expo. time
;####################################################################

Case Action of
    'S': begin
           OnlyOne, wave, flux, texp
           sn = fltarr(n_elements([wave,flux,texp])-2)
           errfrac = fltarr(4,n_elements(sn))
           if not Keyword_set(quiet) then print, 'Compute S/N'
         end
    'T': begin
           OnlyOne, wave, flux, sn
           texp = fltarr(n_elements([wave,flux,sn])-2)
           errfrac = fltarr(4,n_elements(texp))
           if not Keyword_set(quiet) then print, 'Compute Exposure time'
         end
    'F': begin
           OnlyOne, wave, texp, sn
           flux = fltarr(n_elements([wave,texp,sn])-2)
           errfrac = fltarr(4,n_elements(flux))
           if not Keyword_set(quiet) then print, 'Compute Flux'
         end
else: begin
        print, 'ERROR: Unknown Action (',Action,') should be S, T or F'
        return
      end
endcase


;#############################################
; Check if Instrument Mode parameter is valid
;#############################################

case Mode of
   'SR': if not Keyword_set(quiet) then print, 'GHOST Standard Resolution Mode'
   'SF': if not Keyword_set(quiet) then print, 'GHOST Standard Resolution Mode - Faint'
   'SVF': if not Keyword_set(quiet) then print, 'GHOST Standard Resolution Mode - Very Faint'
   'BS': if not Keyword_set(quiet) then print, 'GHOST Beam-Switch Mode'
   'BF': if not Keyword_set(quiet) then print, 'GHOST Beam-Switch Mode - Faint'
   'BVF': if not Keyword_set(quiet) then print, 'GHOST Beam-Switch Mode - Very Faint'
   'HR': if not Keyword_set(quiet) then print, 'GHOST High Resolution Mode'
   'HF': if not Keyword_set(quiet) then print, 'GHOST High Resolution Mode - Faint'
   'PRV': if not Keyword_set(quiet) then print, 'GHOST PRV mode'
else: begin
          print, 'Error in parameter Mode: ', Mode
          return
      end
endcase

;###################################
; Check if Site parameter is valid
;###################################

case Site of
   'GS': if not Keyword_set(quiet) then print, 'Assuming Gemini South extinction'
   'GN': if not Keyword_set(quiet) then print, 'Assuming Gemini North extinction'
   else: begin
             print, 'Error in parameter Site: ', Site
             return
         end
endcase

;###################################
; Principle constants and parameters
;###################################

; Telescope pupil area, in m^2
;GEM = 48.5425 ; GEMINI primary useful aperture in m^2
;GEM = 46.0 ; CoDR value. No reference
;GEM = 49.1 ; PDR value. Based on optical stops on the primary between 8 and 1.2m
;GEM = !PI*(4.0^2 - 1.0^2) ; = 47.1 CDR Based on 8m outer diameter and 2m diameter central blockage by baffles+secondary
GEM = 43.748    ; From Tom Hayward, email to RMcD and
                ; MC. This includes; the M2 mask, the Deployable Baffle which has a radius of 1.0 m when
                ; in the Visible position, and even the obscuration from the spider vanes.

hc = 6.626075510e-34*299792458.0 ; planck constant * light speed
pscale = 1.64                    ; Plate scale arcsec/mm

;=====================================================================
; Instrument parameters. Resolution is as per REQ.
;
; The following parameters are asumed, following analysis by Gordon Robertson, based on ZEMAX model by John Pazder:
; Resolution Element = 3.5 pix in SR, 2.2 pix in HR
; Pixels per lens is 4.12 in SR, 2.38 in HR
StandardRes = 50000.0 ; Standard resolution (SR)
HighRes     = 75000.0 ; High resolution (HR)
ReselSR     = 3.5     ; Resolution element in SR (in pixels)
ReselHR     = 2.2     ; Resolution element in HR (in pixels)
PixLensSR   = 4.12    ; Pixels per lens in SR
PixLensHR   = 2.38    ; Pixels per lens in HR
;=====================================================================

;=====================================================================
; GHOST MODES
; Modes are defined in the ConOps document. Main differences are
; detector binning, spatial sampling and resolution
case Mode of
   'SR': begin
      Res = StandardRes         ; Resolution R=lambda/Delta_lambda
      ResEl = ReselSR           ; Spectral resolution element in pixels
      nObj = 7.                 ; Number of Object Fibers
      nSky = 3.                 ; Number of Sky Fibers
      Dlens = 240.              ; Lens size in micron flat-to-flat
      nSlit = nobj * PixLensSR  ; Number of unbinned pixels along object slit
      xbin = 1.                 ; Detector binning in spectral direction
      ybin = 2.                 ; Detector binning in spatial direction
   end
   'SF': begin
      Res = StandardRes         ; Resolution R=lambda/Delta_lambda
      ResEl = ReselSR           ; Spectral resolution element in pixels
      nObj = 7.                 ; Number of Object Fibers
      nSky = 3.                 ; Number of Sky Fibers
      Dlens = 240.              ; Lens size in micron flat-to-flat
      nSlit = nobj * PixLensSR  ; Number of unbinned pixels along object slit
      xbin = 1.                 ; Detector binning in spectral direction
      ybin = 4.                 ; Detector binning in spatial direction
   end
   'SVF': begin
      Res = StandardRes         ; Resolution R=lambda/Delta_lambda
      ResEl = ReselSR           ; Spectral resolution element in pixels
      nObj = 7.                 ; Number of Object Fibers
      nSky = 3.                 ; Number of Sky Fibers
      Dlens = 240.               ; Lens size in micron flat-to-flat
      nSlit = nobj * PixLensSR  ; Number of unbinned pixels along object slit
      xbin = 2.                 ; Detector binning in spectral direction
      ybin = 4.                 ; Detector binning in spatial direction
   end
   'BS': begin
      Res = StandardRes         ; Resolution R=lambda/Delta_lambda
      ResEl = ReselSR           ; Spectral resolution element in pixels
      nObj = 7.                 ; Number of Object Fibers
      nSky = 7.                 ; Number of Sky Fibers
      Dlens = 240.               ; Lens size in micron flat-to-flat
      nSlit = nobj * PixLensSR  ; Number of unbinned pixels along object slit
      xbin = 1.                 ; Detector binning in spectral direction
      ybin = 2.                 ; Detector binning in spatial direction
   end
   'BF': begin
      Res = StandardRes         ; Resolution R=lambda/Delta_lambda
      ResEl = ReselSR           ; Spectral resolution element in pixels
      nObj = 7.                 ; Number of Object Fibers
      nSky = 7.                 ; Number of Sky Fibers
      Dlens = 240.              ; Lens size in micron flat-to-flat
      nSlit = nobj * PixLensSR  ; Number of unbinned pixels along object slit
      xbin = 1.                 ; Detector binning in spectral direction
      ybin = 8.                 ; Detector binning in spatial direction
   end
   'BVF': begin
      Res = StandardRes         ; Resolution R=lambda/Delta_lambda
      ResEl = ReselSR           ; Spectral resolution element in pixels
      nObj = 7.                 ; Number of Object Fibers
      nSky = 7.                 ; Number of Sky Fibers
      Dlens = 240.              ; Lens size in micron flat-to-flat
      nSlit = nobj * PixLensSR  ; Number of unbinned pixels along object slit
      xbin = 2.                 ; Detector binning in spectral direction
      ybin = 8.                 ; Detector binning in spatial direction
   end
   'HR': begin
      Res = HighRes             ; Resolution R=lambda/Delta_lambda
      ResEl = ReselHR           ; Spectral resolution element in pixels
      nObj = 19.                ; Number of Object Fibers
      nSky = 7.                 ; Number of Sky Fibers
      Dlens = 144.              ; Lens size in micron, flat-to-flat
      nSlit = nobj * PixLensHR  ; Number of unbinned pixels along object slit
      xbin = 1.                 ; Detector binning in spectral direction
      ybin = 2.                 ; Detector binning in spatial direction
   end
   'HF': begin
      Res = HighRes             ; Resolution R=lambda/Delta_lambda
      ResEl = ReselHR           ; Spectral resolution element in pixels
      nObj = 19.                ; Number of Object Fibers
      nSky = 7.                 ; Number of Sky Fibers
      Dlens = 144.              ; Lens size in micron, flat-to-flat
      nSlit = nobj * PixLensHR  ; Number of unbinned pixels along object slit
      xbin = 1.                 ; Detector binning in spectral direction
      ybin = 8.                 ; Detector binning in spatial direction
   end
   'PRV': begin
      Res = HighRes             ; Resolution R=lambda/Delta_lambda
      ResEl = ReselHR           ; Spectral resolution element in pixels
      nObj = 19.                ; Number of Object Fibers
      nSky = 7.                 ; Number of Sky Fibers
      Dlens = 144.              ; Lens size in micron, flat-to-flat
      nSlit = nobj * PixLensHR  ; Number of unbinned pixels along object slit
      xbin = 1.                 ; Detector binning in spectral direction
      ybin = 1.                 ; Detector binning in spatial direction
   end
endcase

;=====================================================================
; SLIT LOSSES
; Compute fraction of flux sampled by the IFU - equivalent to the slit losses
; This function only considers the geometry of the IFU, no other
; losses. This should be checked against Ross's detailed
; modelling (but still kept independent from IFU throughput)
;
;FracSpaGauss, Mode, Seeing, wave, pscale * Dlens/1.e3, FA  ; Assumes Gaussian seeing
FracSpa, Mode, Seeing, wave, pscale * Dlens/1.e3, FA ; Assumes Moffat profile

;=====================================================================
; SPATIAL SAMPLING
;
Dhex = (pscale*Dlens/1.e3)/206265.0 ; Diameter of hexagonal lens, flat to flat, in radians
DA = 2.*sqrt(3.)*(Dhex/2.)^2       ; Area of each lens in rad^2
DASky = nSky * DA                  ; Total area of sky fibres on sky (rad^2)
DAObj = nObj * DA                  ; Total area of object fibres on sky (rad^2)

;=====================================================================
; SPECTRAL SAMPLING
;
npix = (ResEl/xbin)*(nSlit/ybin) ; Total number of pixels summed for object, after binning is applied
npix_org = ResEl*nSlit           ; Original unbinned pixel number for dark current
SkyFact = 1.+nSky/nObj           ; Factor for additional pixels used for sky subtraction from simultaneous sky fibers
Nspec = ResEl                    ; Number of pixels in resolution element
DS = (Wave*1.e-10)/resel/res     ; Size of a spectral pixel in m

;print,npix,'Pixels per object IFU'
;print,npix*nSky/nObj,'Pixels per sky'
;=====================================================================
; DETECTOR NOISE PROPERTIES from requirements
;
RN = 4.0                        ; readout noise in e-
DC = 2.6/3600                      ; dark current in e-/sec

; Control amount of output
case 1 of
   Keyword_set(quiet) : Disp = 0
   n_elements(wave)*n_elements(flux)*n_elements(sn)*n_elements(texp) gt 30 : begin
      Disp = 0
      print, 'more than 30 elements computed ... print supressed'
   end
   else: Disp = 1
endcase

;                          ###################
;                          # BEGIN MAIN LOOP #
;                          ###################

n = 0 ; Initialize counter for output array

for i = 0, n_elements(Wave)-1 do begin

    lbda = wave[i] ; Assign easy variable name

    ; check wavelength is OK
    if ((lbda lt 3630) or (lbda gt 10000)) then begin
        print, 'Error Wavelength ', lbda, ' is outside GHOST limits (3630-10000)'
        return
    endif

    l = lbda*1.e-10 ; wavelength in m
    lmu = l*1.e6    ; wavelength in microns

    ;########################################################################
    ; ATMOSPHERIC EXTINCTION:
    ; Interpolate value from the extinction coefficient.
    ; GS: Patat et al. 2011, A&A, 527, AA91
    ; GN: Buton, C., Copin, Y., Aldering, G., et al. 2013, A&A, 549, AA8
    ;########################################################################

    CASE site OF
       'GS':    readcol,'~/Idl/GHOST_ETC/RefData/paranal_patat11.dat',lext,ext,/silent
       'GN':    readcol,'~/Idl/GHOST_ETC/RefData/MK_extinction_Buton.dat',lext,ext,F='(f,f)',/silent
    ENDCASE
    
    Extinct = interpol(10.^(-0.4*ext*Airmass),lext,lbda)

    ;########################################################################
    ; THROUGHPUT
    ; Read throughput data from Ross Zelhem. This has predicted and budgeted data for Cass
    ; Unit and Spectrograph. 10-Oct-2016 - includes measured echelle data
    ;########################################################################
    readcol,'~/Idl/GHOST_ETC/RefData/GHOST_throughput_Echelle.dat',lnm2,totTPbudg,/silent
    
    case THROUGHPUT of
       ; REQUIREMENTS for CDR
       'R': begin
          if not keyword_set(QUIET) then print,'Assuming throughput specified in Requirement 4110'
          CASE 1 OF
             (lbda ge 3630 and lbda lt 3750): Tghost = 0.08
             (lbda ge 3750 and lbda lt 4500): Tghost = 0.127
             (lbda eq 4500): Tghost = 0.27
             (lbda gt 4500 and lbda le 9000): Tghost = 0.136
             (lbda ge 9000 and lbda le 9500): Tghost = 0.055
             (lbda gt 9500 and lbda le 10000): Tghost = 0.021
             else: print,'wavelength out of range'
          ENDCASE
       end
       ; Total budgeted throughput
       'B': begin
          if not keyword_set(QUIET) then print,'Assuming BUDGETED throughput from Ross Zhelem, 10-10-2016'
          Tghost = interpol(totTPbudg,lnm2*10.,lbda)
       end
       'S': begin
          if not keyword_set(QUIET) then print,'Assuming SHORT FIBER throughput from Mick Edgar, 15-08-2016'
          readcol,'~/Idl/GHOST_ETC/RefData/GHOSTFullTransmittanceShortFiber.txt',lnm,totT,/silent
          Tghost = interpol(totT,lnm*10.,lbda)
       end
    endcase

    ;########################################################################
    ; GEMINI MIRROR REFLECITIVITY
    ;
    ; Various references were used:
    ;  CoDR - Flat reflecitivty of 95% for each of three mirrors
    ;  PDR - Boccas et al. 2006, Thin Solid Films, 502, 275
    ;  CDR/FDR - Reference data from Excel spreadsheet from Madeline Close, "Gemini_Reflectivity_2008-2015.txt", Version 22 April 2005
    ; Function for reflectivity of Gemini mirrors using data from Vucina et al 2008
    ; Proc. of SPIE Vol. 7012 70122Q-1
    ; http://www.saao.ac.za/~dod/M5_Washing/SPIE7012_101_Gemini.pdf
    ;########################################################################

    refdata = 'CDR' ; Switch for reflectivity assumption
;    if refdata eq 1 then M123 = (0.01*(14.09*alog10(lbda) + 36.94))^3 else M123 = 0.925^3
    CASE refdata OF
       'CoDR': M123 = 0.95^3
       'PDR': begin
          readcol,'~/Idl/GHOST_ETC/RefData/boccas_reflectivity.dat',lnmRef,Reflectivity,/silent ; Digitized from Boccas et al. 2006, Thin Solid Films, 502, 275. Used at PDR
          M123 = 0.955*(0.01*interpol(Reflectivity, lnmRef*10., lbda))^3
       end
       'CDR': begin
          readcol,'~/Idl/GHOST_ETC/RefData/GS_M1_Ag_Sample_10-8-2010.dat',lnmRef,Reflectivity,/silent,F='(f,f)' ; Data from Excel spreadsheet from Madeline Close, "Gemini_Reflectivity_2008-2015.txt", Version 22 April 2005
          M1 = 0.95*(0.01*interpol(Reflectivity, lnmRef*10., lbda))  ; Based on worst M1
          M2 = 0.975*(0.01*interpol(Reflectivity, lnmRef*10., lbda)) ; Based on worst M2
          if port eq 1 then M3 = 1.0 else M3 = M1                    ; For M3, assume worst M1. Skip M3 if upward looking port is set (port=1)
          M123 = M1*M2*M3                                            ; Combined surfaces
       end
    ENDCASE

    Tghost = Tghost * M123
    
    ;########################################################################
    ; SKY BRIGHTNESS
    ; Compute polynomial approximation of no OH Paranal Sky (in erg/s/A/cm2/arcsec2)
    ; The coefficients come from mathcad fit to data from
    ; Hanuschik R.W., 2003, A&A, 407, 1157
    ;########################################################################

;    bright = 1
    if not keyword_set(Quiet) then if keyword_set(BRIGHT) then print,'Assuming bright sky' else print,'Assuming dark sky conditions'
    FSky = sky(lmu, BRIGHT=bright);,/print)
    FSky = FSky*206265.0^2*1.e7 ; convert in SI
                                ;
    ;
    ;########################################################################
    ; SCATTERED LIGHT
    ; Include a crude scattered light component depending on the
    ; wavelength. This scales the effective counts from the source only,
    ; and includes this additional flux in the Poisson noise calculation
    ; (so only adding noise, not signal).
    ;########################################################################

    case 1 OF
       lbda le 3750: ScatLight = 1.05
       lbda gt 3750: ScatLight = 1.02
    endcase

    ;########################################################################
    ; COMPUTE requested parameter
    ;########################################################################

    case Action of
       'S': begin
               for j = 0, n_elements(flux)-1 do begin
                 for k = 0, n_elements(texp)-1 do begin
                    Comp_SN, snc, flux[j], texp[k], errstat, Nexp, Mode, $
                        nspec, nspat, npix, npix_org, SkyFact, FA, ScatLight
                    sn[n] = snc
                    errfrac[*,n] = errstat
                    n++
                 endfor
                endfor
            end
       'T': begin
               for j = 0, n_elements(flux)-1 do begin
                 for k = 0, n_elements(sn)-1 do begin
                    Comp_T, sn[k], flux[j], texpc, errstat, Nexp, Mode, $
                        nspec, nspat, npix, npix_org, SkyFact, FA, ScatLight
                    texp[n] = texpc
                    errfrac[*,n] = errstat
                    n++
                 endfor
                endfor
            end
        'F': begin
               for j = 0, n_elements(texp)-1 do begin
                 for k = 0, n_elements(sn)-1 do begin
                    Comp_F, sn[k], fluxc, texp[j], errstat, Nexp, Mode, $
                        nspec, nspat, npix, npix_org, SkyFact, FA, ScatLight
                    flux[n] = fluxc
                    errfrac[*,n] = errstat
                    n++
                 endfor
                endfor
            end
    endcase

endfor

;                             ##############################
;                             # PLOT OUTPUTS (if relevant) #
;                             ##############################

if n gt 1 and not Keyword_set(quiet) then begin
    case Action of
       'S': begin
               if (n_elements(wave) gt 1) then begin
                 plot,wave,sn,xtitle='Wavelength',ytitle='S/N',$
                  title='S/N at Flux: '+ string(flux) + ' & ' + 'Integ/exp. ' + $
                   string(texp) + 's [ ' + string (nexp) + ' exp]'
               endif
               if (n_elements(flux) gt 1) then begin
                 plot,flux,sn,xtitle='Flux',ytitle='S/N',$
                  title='S/N at '+ string(wave) + 'A & Integ/exp. ' + $
                  string(texp) + 's [ ' + string (nexp) + ' exp]'
               endif
               if (n_elements(texp) gt 1) then begin
                 plot,texp,sn,xtitle='Integ/exp. (sec) [' + string(nexp) + $
                    ' exp]',ytitle='S/N',$
                  title='S/N at '+ string(wave) + 'A & Flux ' + string(flux)
               endif
            end
       'F': begin
               if (n_elements(wave) gt 1) then begin
                 plot,wave,flux,xtitle='Wavelength',ytitle='Flux',$
                  title='Flux for S/N '+ string(sn) + ' & ' + 'Integ/exp. ' + $
                   string(texp) + 's [ ' + string (nexp) + ' exp]'
               endif
               if (n_elements(sn) gt 1) then begin
                 plot,sn,flux,xtitle='S/N',ytitle='Flux',$
                  title='Flux at '+ string(wave) + 'A & Integ/exp. ' + $
                  string(texp) + 's [ ' + string (nexp) + ' exp]'
               endif
               if (n_elements(texp) gt 1) then begin
                 plot,texp,flux,xtitle='Integ/exp. (sec) [' + string(nexp) + $
                    ' exp]',ytitle='Flux',$
                  title='Flux at '+ string(wave) + 'A & S/N ' + string(sn)
               endif
            end
       'T': begin
               if (n_elements(wave) gt 1) then begin
                 plot,wave,texp,xtitle='Wavelength',ytitle='Integ/exp. ' + $
                    's [' + string(nexp) + ']', $
                    title='Int. Time/exp for S/N '+ string(sn) + ' & Flux ' + string (flux)
               endif
               if (n_elements(sn) gt 1) then begin
                 plot,sn,texp,xtitle='S/N',ytitle='Integ/exp. ' + $
                    's [' + string(nexp) + ']', $
                  title='Int. Time/exp at '+ string(wave) + 'A & Flux ' + string(flux)
               endif
               if (n_elements(flux) gt 1) then begin
                 plot,texp,texp,xtitle='Flux ' + string(flux), ytitle='Integ/exp. ' + $
                    's [' + string(nexp) + ']', $
                  title='Int. Time/exp at '+ string(wave) + 'A & S/N ' + string(sn)
               endif
            end
    endcase
endif


end

;
;=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
;

