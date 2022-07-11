function SN = SNvsAB(lambda,AB_low,AB_high,t_exp,ZD,resoln,seeing,N_mirror,SR_sky,bin_spat,plotQ)
%
%   Calculate (and optional plot) SN vs AB magnitude, at a single wavelength
%   * No provision for spectral binning as yet *
%
%   Input parameters:
%   -----------------
%   lambda      : wavelength at which S/N is required. Single scalar only. In *nm*
%   AB_low      : Brightest AB magnitude of object, assumed unresolved
%   AB_high     : Faintest AB magnitude of object
%   t_exp       : exposure time, seconds
%   ZD          : Zenith distance of object, degrees
%   resoln      : Resolution mode - 'SR' for standard, 'HR' for high resolution
%   seeing      : Seeing disc FWHM, arcseconds
%   N_mirror    : Number of Gemini reflections: 2 for axial port, 3 for side port
%   SR_sky      : For SR only - number of sky microlenses: 3, 7 or 10
%   bin_spat    : Binning factor in spatial direction; 1 for no binning
%   plotQ       : 1 for plot, 0 for no plot
%
%   Output parameters:
%   ------------------
%   SN          : Signal/noise per resolution element at each value of input lambda
%
%                                   G. Robertson  29 July 2019. [GHOST 3 140]
%
%   Presets
%
    n_AB = 101;         % number of separate evaluations between AB limits
    SN = zeros(1,n_AB);
    vars = zeros(4,n_AB);
    [dim1,dim2] = size(lambda);
    assert(dim1 == 1 && dim2 == 1,'lambda must be a single scalar for SN vs AB!')
%
    for i = 1:n_AB
        AB = AB_low + (i-1)*(AB_high - AB_low)/(n_AB-1);
        [SN(i), vars(:,i)]=SNmain(lambda,AB,t_exp,ZD,resoln,seeing,N_mirror,SR_sky,bin_spat,0);
    end
%
    if plotQ
        subplot(2,1,1)
        semilogy(linspace(AB_low,AB_high,n_AB),SN,'LineWidth',1.5)
        ylabel('Signal/noise')
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

