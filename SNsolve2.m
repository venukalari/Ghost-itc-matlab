function AB_SN = SNsolve2(lambda,t_low,t_high,SN,ZD,resoln,seeing,SB,N_mirror,SR_sky,bin_spat,N_sub,plotQ)
%
%   Calculate (and optional plot) AB for given S/N at range of exp times, at a single wavelength
%   * No provision for spectral binning as yet *
%
%   Input parameters:
%   -----------------
%   lambda      : wavelength at which S/N is required. Single scalar only. In *nm*
%   t_low       : Shortest exposure time
%   t_high      : Longest exposure time
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
%   AB_SN      : AB magnitude for target signal/noise per resolution element at 
%                  each value of input exp time.
%
%                                   G. Robertson  7 August 2019. [GHOST 3 161]
%
%   Modified to call SNmain including sky brightness. 28 July 2020. [GHOST 5 4]
%
%   Presets
%
    if t_low == t_high
        n_time = 1;
    else
        n_time = 11;         % number of separate evaluations between exp time limits
    end
    AB_SN = zeros(1,n_time);
    SN_calc=zeros(1,n_time);
    vars = zeros(4,n_time);
    [dim1,dim2] = size(lambda);
    assert(dim1 == 1 && dim2 == 1,'lambda must be a single scalar for AB vs t_exp!')
    options = optimset('TolX',1e-4);
    
%
    for i = 1:n_time
        t_exp = t_low + (i-1)*(t_high - t_low)/(n_time-1 + 10*eps);
        AB_init = 2.36*(log10(t_exp) + 3.34 - 2*log10(SN/50)); % initial guess (applic to order centres)
%
%       Use Matlab facility fzero to solve for exposure time giving desired SN
% 
    myfun = @(lambda,AB,t_exp,ZD,resoln,seeing,SB,N_mirror,SR_sky,bin_spat,N_sub,plotQ2) ...
    (SNmain(lambda,AB,t_exp,ZD,resoln,seeing,SB,N_mirror,SR_sky,bin_spat,N_sub,plotQ2) - SN);
        plotQ2 = 0;
        fun = @(AB) myfun(lambda,AB,t_exp,ZD,resoln,seeing,SB,N_mirror,SR_sky,bin_spat,N_sub,plotQ2);
        AB = fzero(fun,[AB_init-4 AB_init+2],options);
%             
        [SN_calc(i), vars(:,i)]=SNmain(lambda,AB,t_exp,ZD,resoln,seeing,SB,N_mirror,SR_sky,bin_spat,N_sub,0);
        AB_SN(i) = AB;
        assert(abs(SN_calc(i) - SN)/SN < 0.01,'Did not reach target SN!')
    end
%
    if plotQ && n_time ~= 1
        subplot(2,1,1)
        plot(linspace(t_low,t_high,n_time),AB_SN,'LineWidth',1.5)
        grid on
        hold on
        ylabel('AB mag')
        xlabel(['Exposure time at ',num2str(lambda),' nm  /s'])
        axis_vec = axis;
        plot([3600 3600],[axis_vec(3) axis_vec(4)],'color',[0.75 0.75 0.75],'LineWidth',1.5)
        subplot(2,1,2)
        plot(linspace(t_low,t_high,n_time),vars(1,:),'red','LineWidth',1.5)
        hold on
        plot(linspace(t_low,t_high,n_time),vars(2,:),'blue','LineWidth',1.5)
        plot(linspace(t_low,t_high,n_time),vars(3,:),'black','LineWidth',1.5)
        plot(linspace(t_low,t_high,n_time),vars(4,:),'green','LineWidth',1.5)
        xlabel(['Exposure time at ',num2str(lambda),' nm  /s'])
        ylabel('Variances')
        axis_vec = axis;
        plot([3600 3600],[axis_vec(3) axis_vec(4)],'color',[0.75 0.75 0.75],'LineWidth',1.5)
    end
%        
return
end

