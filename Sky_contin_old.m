function flux = Sky_contin(lambda, plotQ)
%
%   Sky continuum emission vs wavelength.
%   Currently based on data from Hanuschik A&A 407 1157 2003.
%
%   Input parameters:
%   -----------------
%   lambda_in   : single value or column vector of wavelength values, in *nm*
%   plotQ       : 1 for plot, 0 for no plot
%
%   Output parameters:
%   ------------------
%   flux  : vector, sky continuum flux, (erg /s /Å /cm^2 /arcs^2) at wavelengths of lambda
%
%                                                       JGR 28 June 2019
%
%   Presets
%
    Han_lam_l = [314 374 480 583 670 860];  % nm
    Han_lam_h =[376 486 577 679 856 1043];
    Han_flux = 1e-16*[0.17 0.14 0.09 0.10 0.08 0.07];
    lambda_offset = 650;
    lam_min = 315;
    lam_max = 1040;
% 
    lambda = lambda(:).';  % ensures row vector 
    [dim1,~] = size(lambda);
    assert(dim1 == 1,'lambda is not a scalar or vector!')
    assert(min(lambda)>=lam_min,'lambda value(s) below blue limit!')
    assert(max(lambda)<=lam_max,'lambda value(s) above red limit!')
%
%   Fitted parameters previously found for decaying exp
%
    betahat = 100*[-0.000056876550681   0.001533398242066   0.000693822733270 -3.765787850052728];
%    
%   Do 3rd order polynomial fit to the data points
%
%     xpoly = (Han_lam_l + Han_lam_h)/2 - lambda_offset;
%     ypoly = Han_flux * 1e16;
%     poly = polyfit(xpoly,ypoly,3);
%     flux = polyval(poly, lambda_in - lambda_offset)*1e-16;
%

%   
%   Try decaying exponential fit
%
%     beta_0 = [-0.005 0.2 0.05 -350];
%     [betahat] = nlinfit(xpoly,ypoly,@exp_fit,beta_0);
%
    flux2 = 1e-16*exp_fit(betahat,lambda - lambda_offset);
    if plotQ
        for i=1:6
            plot([Han_lam_l(i) Han_lam_h(i)],1e16*[Han_flux(i) Han_flux(i)],'red','LineWidth',2)
            hold on
        end
%        plot(lambda_in,flux,'green','LineWidth',2)
        plot(lambda,1e16*flux2,'blue','LineWidth',2)
    end    
    flux = flux2;
%    
    return   
end

