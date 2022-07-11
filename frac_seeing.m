function frac=frac_seeing(FWHM_min,FWHM_max,resoln)
%
%   Produce data and plot for IFU fraction transmitted vs seeing
%
%
%   Presets
%
    lambda = 500;
    flag = 9;
    ZD = 30; % not used
    N_mirror = 3; % not used
    n_values = 101; % number of seeing values plotted
%
    for i = 1:n_values
        seeing = FWHM_min + (i-1)*(FWHM_max - FWHM_min)/(n_values - 1);
        frac(i)=Throughput_main(lambda,flag,ZD,resoln,seeing,N_mirror,0);
    end
%
    plot(linspace(FWHM_min,FWHM_max,n_values),100*frac,'LineWidth',1.5)
    ylabel('Throughput  /%')
    xlabel('Seeing FWHM  /arcsec')
    axis([FWHM_min FWHM_max 0 100]);
%    
return
end

