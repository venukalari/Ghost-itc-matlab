function [SN,vars]=SNmain_specbin(lambda,AB,t_exp,ZD,resoln,seeing,SB,N_mirror,SR_sky,bin_spat,N_sub,plotQ)
%
%   Second routine for signal/noise calculation for GHOST.
%   ** This one does spectral binning x2 ** which is only available for SR
%
%   Input parameters:
%   -----------------
%   lambda      : wavelength(s) at which S/N is required. Vector or scalar. In *nm*
%   AB          : AB magnitude of object, assumed unresolved
%   t_exp       : exposure time, seconds. *Total* over N_sub exposures if N_sub > 1
%   ZD          : Zenith distance of object, degrees
%   resoln      : Resolution mode - 'SR' for standard 
%   seeing      : Seeing disc FWHM, arcseconds
%   SB          : Sky brightness class, character 'SB20', 'SB50', 'SB80', 'SBAny'
%   N_mirror    : Number of Gemini reflections: 2 for axial port, 3 for side port
%   SR_sky      : Number of sky microlenses: 3, 7 or 10
%   bin_spat    : Binning factor in spatial direction; 1 for no binning
%   N_sub       : Number of sub-exposures comprising the total t_exp. 1 for single exp
%   plotQ       : 1 for plot, 0 for no plot
%
%   Output parameters:
%   ------------------
%   SN          : Signal/noise per resolution element at each value of output lambda
%   vars        : array of variances, one value for each lambda
%                  row 1: object, row 2 sky, row 3 dark, row 4 readout
%
%                                   G. Robertson  7 August 2019. [GHOST 3 159]
%
%   Presets
%
    area_tel = 491000.0;  % Gemini primary mirror area, cm^2
    RON_red = 2.2;  % red CCD read noise e- rms (NRC rep Jun 2019; slow readout; av of 4)
    RON_blue = 2.25; % blue CCD read noise e- rms (NRC rep Jun 2019; slow readout; av of 4)
    dark_red = 0.825; % red dark noise, e- /pix /hr (NRC rep Jun 2019; av of 4)
    dark_blue = 1.175; % blue dark noise, e- /pix /hr (NRC rep Jun 2019; av of 4)   
    area_SR = 0.939; % SR IFU area, arcsec^2
    lambda_cross = 533.50; % wavelength of crossover from blue to red camera, nm
    Npix_SR = 18.9;  % Number of (unbinned) pixels in length of 1 star SR slit
%    
    lambda = lambda(:).';  % ensures row vector 
    [dim1,n_lambda] = size(lambda);
    assert(dim1 == 1,'lambda is not a scalar or vector!')   
    assert(min(lambda)>=360,'lambda value(s) below blue limit!')
    assert(max(lambda)<=1000,'lambda value(s) above red limit!')
    assert(N_mirror == 2 || N_mirror == 3,'Illegal number of telescope reflections!')
    assert(abs(round(N_sub) - N_sub)<10*eps && N_sub>0,'Number of subexposures must be integer!')
    assert(strcmp(resoln,'SR'),'Spectral binning is calculated for standard resolution only!')
    vars = zeros(4,n_lambda);
%    
    area_IFU = area_SR;
    slit_len = Npix_SR;
    assert(SR_sky == 3 || SR_sky == 7 || SR_sky == 10,'Illegal number of SR sky lenses!')
    sky_coeff = 1 + 7/SR_sky;  
%
%   Get number of signal photons from the object in each 'pixel' (at lambda values, which
%   won't in general be actual CCD pixel wavelengths, but will use pixel width near that
%   value, which is OK)
%   
    log_f_lambda = -AB/2.5 - 2*log10(lambda*10) - 2.408/2.5;  % get f_lambda
    f_lambda = 10.^log_f_lambda;  % in erg /s /cm^2 /Å
%
    rate0 = f_lambda.*lambda*1E-9/(6.6256E-27*2.99792E8);
               % photon rate incident on top of atmosphere, phot /cm^2 /Å /s
    rate1 = rate0*area_tel;  %     phot /Å /s  for whole telescope area
%    
    extinc_frac = Extinc_Paranal(lambda,ZD,0);  % Paranal extinction; fraction passed
    rate2 = rate1.*extinc_frac; % after atmos extinction; phot /Å /s
%    
    rate3 = rate2.*GS_reflectivity(lambda,N_mirror,0); % after telescope mirrors; phot /Å /s  
    rate4 = rate3*IFU_trans(resoln,seeing); % drop by IFU injection loss
    rate5 = rate4.*GHOST_cable_eta(lambda); % drop by Cass Unit and cable throughputs 
    rate6 = rate5.*sgr_throughput(lambda,0); % drop by sgr and detector throughputs    
    rate7 = rate6.*BlazeFunction(lambda);    % drop by individual order blaze functions
%    
    RD = 2*nmperpix(lambda)*10;  % Reciprocal dispersion, Å/pix for binned pixels
    rate8 = rate7.*RD*0.99;   % rate in phot /s /pixel; drop by 1% loss on slitview beamsplitter
%    
    object = rate8*t_exp;    % phot /binned pix over whole exposure time; summed over slit length
%
%   Assemble noise variances
%
    sky_flux = Sky_contin(lambda,SB,0);  % sky *continuum* flux, erg /s /Å /cm^2 /arcs^2
    sky_rate1 = sky_flux.*lambda*1E-9*area_tel*t_exp/(6.6256E-27*2.99792E8); % phot /Å /arcs^2
    sky_rate2 = sky_rate1.*GS_reflectivity(lambda,N_mirror,0);   % losses through the system
    sky_rate2a = sky_rate2.*extinc_frac; % Hanuschik sky brightness are for outside atmosphere(!)
    sky_rate3 = sky_rate2a.*GHOST_cable_eta(lambda).*sgr_throughput(lambda,0);
    sky_rate4 = sky_rate3.*BlazeFunction(lambda);  
    sky_rate5 = sky_rate4*area_IFU; % phot /Å after injection into IFU
    sky_rate6 = sky_rate5.*RD; % phot / binned pix in whole exposure time; summed over slit length
    sky_rate7 = sky_rate6*sky_coeff*0.99; % scale up to allow for sky subtraction noise; 1% BS loss
%    
    dark0 = (lambda < lambda_cross)*dark_blue + (lambda >= lambda_cross)*dark_red; %  e- /pix /hr  
    dark1 = dark0*t_exp/3600; %  e- /pix 
    dark2 = 2*dark1*slit_len;   %  e- / binned pix, summed over slit length
%
    ron_rms = (lambda < lambda_cross)*RON_blue + (lambda >= lambda_cross)*RON_red;
    ron_var = slit_len*N_sub*(ron_rms.^2)/bin_spat; % read noise variance for 1 (binned) 
                                    % spectral pixel, all spatial pixels, all subexposures
%
    var_total = object + sky_rate7 + dark2 + ron_var; % for 1 (binned) spectral pixel
    vars = [object; sky_rate7; dark2; ron_var];
%
    SN_pix = object./sqrt(var_total);     % S/N per pixel (spectral binned x2)
    B2_1 = B2lookup(lambda,resoln,1);     % sum(B^2) on the original unbinned wavelength grid
    B2 = 0.5*(1 + B2_1); % i.e. on av sum(B^2) will have half the excess over 1.0 cf unbinned  
%
    SN = SN_pix.*sqrt(B2); % S/N for resolution element
%       
    if plotQ
        subplot(2,1,1)
        plot(lambda,SN,'LineWidth',1.5)
        ylabel('Signal/noise')
        xlabel('Wavelength  /nm')
        xlim([363 1000])
        subplot(2,1,2)
        plot(lambda,object,'red','LineWidth',1.5)
        hold on
        plot(lambda,sky_rate7,'blue','LineWidth',1.5)
        plot(lambda,dark2,'black','LineWidth',1.5)
        plot(lambda,ron_var,'green','LineWidth',1.5)
        xlabel('Wavelength  /nm')
        ylabel('Variances')
        xlim([363 1000])
    end
return
end

