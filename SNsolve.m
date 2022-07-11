function timeSN = SNsolve(lambda,AB_low,AB_high,SN,ZD,resoln,seeing,N_mirror,SR_sky,bin_spat,N_sub,plotQ)
%
%   Calculate (and optional plot) time to given S/N vs AB magnitude, at a single wavelength
%   * No provision for spectral binning as yet *
%
%   Input parameters:
%   -----------------
%   lambda      : wavelength at which S/N is required. Single scalar only. In *nm*
%   AB_low      : Brightest AB magnitude of object, assumed unresolved
%   AB_high     : Faintest AB magnitude of object
%   SN          : desired target Signal/Noise
%   ZD          : Zenith distance of object, degrees
%   resoln      : Resolution mode - 'SR' for standard, 'HR' for high resolution
%   seeing      : Seeing disc FWHM, arcseconds
%   N_mirror    : Number of Gemini reflections: 2 for axial port, 3 for side port
%   SR_sky      : For SR only - number of sky microlenses: 3, 7 or 10
%   bin_spat    : Binning factor in spatial direction; 1 for no binning
%   N_sub       : Number of sub-exposures comprising the total t_exp. 1 for single exp
%   plotQ       : 1 for plot, 0 for no plot
%
%   Output parameters:
%   ------------------
%   timeSN      : Time to reach target signal/noise per resolution element at 
%                  each value of input AB mag
%
%                                   G. Robertson  1 August 2019. [GHOST 3 153]
%
%   Presets
%
    n_AB = 11;         % number of separate evaluations between AB limits
    timeSN = zeros(1,n_AB);
    SN_calc=zeros(1,n_AB);
    vars = zeros(4,n_AB);
    [dim1,dim2] = size(lambda);
    assert(dim1 == 1 && dim2 == 1,'lambda must be a single scalar for SN vs AB!')
    options = optimset('TolX',1e-4);
    
%
    for i = 1:n_AB
        AB = AB_low + (i-1)*(AB_high - AB_low)/(n_AB-1);
        log_t_init = -3.34 + 0.424*AB +2*log10(SN/50);  % initial guess (applic to order centres)
        t_init = 10.^log_t_init; 
%
%       Use Matlab facility fzero to solve for exposure time giving desired SN
% 
    myfun = @(lambda,AB,t_exp,ZD,resoln,seeing,N_mirror,SR_sky,bin_spat,N_sub,plotQ2) ...
    (SNmain(lambda,AB,t_exp,ZD,resoln,seeing,N_mirror,SR_sky,bin_spat,N_sub,plotQ2) - SN);
        plotQ2 = 0;
        fun = @(t_exp) myfun(lambda,AB,t_exp,ZD,resoln,seeing,N_mirror,SR_sky,bin_spat,N_sub,plotQ2);
        t = fzero(fun,[0.2*t_init 10*t_init],options);
%             
        [SN_calc(i), vars(:,i)]=SNmain(lambda,AB,t,ZD,resoln,seeing,N_mirror,SR_sky,bin_spat,N_sub,0);
        timeSN(i) = t;
    end
%
    if plotQ
        subplot(2,1,1)
        semilogy(linspace(AB_low,AB_high,n_AB),timeSN,'LineWidth',1.5)
        ylabel(['Exposure Time for S/N = ',num2str(SN)])
        xlabel(['AB magnitude at ',num2str(lambda),' nm'])
        subplot(2,1,2)
        semilogy(linspace(AB_low,AB_high,n_AB),vars(1,:),'red','LineWidth',1.5)
        hold on
        semilogy(linspace(AB_low,AB_high,n_AB),vars(2,:),'blue')
        semilogy(linspace(AB_low,AB_high,n_AB),vars(3,:),'black')
        semilogy(linspace(AB_low,AB_high,n_AB),vars(4,:),'green')
        xlabel(['AB magnitude at ',num2str(lambda),' nm'])
        ylabel('Variances')
    end
%        
return
end

