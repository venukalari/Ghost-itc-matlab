function [varargout] = comp_sn(varargin)
%
% File generated by IDL2Matlab 1.6 130501 %%

% %Initialization of parameters
  I2Mkwn=char('I2M_a1', 'I2M_a2', 'I2M_a3', 'I2M_a4', 'I2M_a5', 'I2M_a6', 'I2M_a7', 'I2M_a8', 'I2M_a9', 'I2M_b1', 'I2M_b2', 'I2M_b3', 'I2M_b4', 'I2M_pos');
  I2Mkwv={'sn', 'flux', 'texp', 'errfrac', 'nexp', 'mode', 'nspec', 'nspat', 'npix', 'npix_org', 'skyfact', 'fa', 'scatlight', 'I2M_pos'};
  sn=[]; flux=[]; texp=[]; errfrac=[]; nexp=[]; mode=[]; nspec=[]; nspat=[]; npix=[]; npix_org=[]; skyfact=[]; fa=[]; scatlight=[]; I2M_pos=[];
  I2M_lst={}; I2M_out=''; lv=length(varargin); if rem(lv,2) ~= 0, I2M_ok=0; else, I2M_ok=1;
  for I2M=1:2:lv; I2M_tmp=varargin{I2M}; if ~ischar(I2M_tmp); I2M_ok=0; break; end; I2Mx=strmatch(I2M_tmp,I2Mkwn); if length(I2Mx) ~=1; I2M_ok=0; break; end; eval([I2Mkwv{I2Mx} '=varargin{I2M+1};']); I2M_lst{(I2M+1)/2}=I2Mkwv{I2Mx}; end; end;
  if ~I2M_ok; for I2M=1:lv; eval([I2Mkwv{I2M} '=varargin{I2M};']); end; end;
  if ~isempty(I2M_pos); for I2M=1:length(I2M_pos); I2Ms=num2str(I2M); I2M_out=[I2M_out 'varargout{' I2Ms '}=' I2M_lst{I2M_pos(I2M)} '; ']; end; end;

% %End of parameters initialization

  % Creation of undeclared variables of functions parameters
  lbda=1; lmu=1; 


  %
  % description
  % computates s/n for the given flux and integration time
  %
  % inputs
  % sn - returned s/n value
  % flux - flux of input source
  % texp - exposure time in seconds
  % errfrac - returned percentage contributions to total variance [obj,sky,rdn,dc]
  % nexp - number of exposures, each with length texp
  % mode - ghost observing mode
  % nspec - number of spectral pixels summed in computation
  % nspat - number of spatial elements summed in computation
  % npix - number of binned object pixels (for read noise)
  % npix_org - number of unbinned object pixels (for dark current)
  %
  %========================================================================================
  global 

  if ((M2I_disp == 1))
    printt('------------------------------------------------------------------------');
    %###############################################################
    % ponctual & continuum source
    % here we must use the fractional flux only in the spatial domain.
    %###############################################################
    %
    % working in flux or ab magnitudes?
    %
  end%if
  if (flag_ab == 1)
    [fobj, lbda, flux] = ab2flux('I2M_a1', flux, 'I2M_a2', lbda, 'si', 1, 'I2M_pos', [2, 1]);  
  else
% convert AB in flux by A in SI unit

        fobj = flux .* 1e7;    % convert flux in si unit
    %
    %
    % compute effective collecting power for object (ko) and sky (ks) entering the ghost
    % science object ifu. for the object, this is the product of number of
    % spectral pixels (nspec), spectral pixel size (ds), fraction
    % of included flux (fa), throughput (tghost), collecting area (gem),
    % extinction, and flux of given wavelength (l/hc).
    % for the sky, this is the product of nspec, ds, area of sky observed
    % (dasky), tghost, gem and flux.
    %
  end%if
  ko = doubll(nspec .* ds .* fa .* tghost .* gem .* extinct .* l ./ hc);
  ks = doubll(nspec .* ds .* dasky .* tghost .* gem .* l ./ hc);
  %
  % compute counts from object and sky, applying collecting power to
  % flux rate of source and integration time
  %
  objcounts = ko .* fobj .* texp;
  skycounts = ks .* skyfact .* fsky .* texp;
  %
  % compute noise from detector
  % note: the sky is subtracted via peripheral sky fibers. the poisson
  % and read noise due to this is included by scaling the effective
  % number of pixels (for read noise and dark current) and sky counts by
  % the skyfact factor, which is based on the number and area of sky
  % fibers relative to the ifu.
  %
  rdn = skyfact .* npix .* rn.^2;
  dkn = skyfact .* npix_org .* dc .* texp;
  %
  % compute s/n using standard formula and input object flux. this includes the noise
  % contribution from poisson statistics, sky component and ccd noise (read noise &
  % dark current). scatlight accounts crudely for scattered light. this
  % scales the effective counts from the source only, and includes this
  % additional flux in the poisson noise calculation (so only adding
  % noise, not signal).
  %
  sn = sqrt(nexp) .* objcounts ./ sqrt(objcounts .* scatlight + skycounts + rdn + dkn);
  %
  % assemble values expresing the fractional errors
  %
  errfrac = d1_array(objcounts,skycounts,rdn,dkn);
  [errtot, errfrac] = total('I2M_a1', errfrac, 'I2M_pos', [1]);
  errfrac = 100.0 .* errfrac ./ errtot;
  %
  % output numbers if running in verbose mode.
  %
  if ((M2I_disp == 1))

    [fa, npix] = printt('format', '("Summed spaxels: ",i4, " fraction of object flux ",f4.2)', 'I2M_a1', npix, 'I2M_a2', fa, 'I2M_pos', [3, 2]);
    if ((nspec > 1))

      resolv = l ./ (nspec .* ds);
      [nspec, resolv] = printt('format', '("Low R: ", f8.1," : ", f5.1," spectral pixels sum up")', 'I2M_a1', resolv, 'I2M_a2', nspec, 'I2M_pos', [3, 2]);
    end%if

    if (flag_ab == 1)

      [flux, lmu] = printt('format', '("Lbda: ",f4.2," um Object AB magnitude: ",f5.2)', 'I2M_a1', lmu, 'I2M_a2', flux, 'I2M_pos', [3, 2]);
    
    else

      [flux, lmu] = printt('format', '("Lbda: ",f4.2," um Object flux: ",e12.4," erg/s/cm2/A")', 'I2M_a1', lmu, 'I2M_a2', flux, 'I2M_pos', [3, 2]);
    end%if

    [dkn, rdn, skycounts, objcounts] = print_noise('I2M_a1', objcounts, 'I2M_a2', skycounts, 'I2M_a3', rdn, 'I2M_a4', dkn, 'I2M_pos', [4, 3, 2, 1]);
  end%if

  %####################
  % final output values
  %####################
  if ((M2I_disp == 1))
    [texp, nexp, sn] = printt('format', '("S/N: ",f8.3, " in ",i2, " exposure(s) of ", f8.1," sec")', 'I2M_a1', sn, 'I2M_a2', nexp, 'I2M_a3', texp, 'I2M_pos', [4, 3, 2]);
  end%if
  if ~isempty(I2M_out),eval(I2M_out);end;return;

if ~isempty(I2M_out),eval(I2M_out);end;
 return;
% % end of function comp_sn
