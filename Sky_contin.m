function flux = Sky_contin(lambda, SB, plotQ)
%
%   Sky continuum emission vs wavelength.
%   Currently based on data from Hanuschik A&A 407 1157 2003.
%
%   Input parameters:
%   -----------------
%   lambda_in   : single value or column vector of wavelength values, in *nm*
%   SB          : sky brightness class, character 'SB20', 'SB50', 'SB80', 'SBAny' only
%   plotQ       : 1 for plot, 0 for no plot
%
%   Output parameters:
%   ------------------
%   flux  : vector, sky continuum flux, (erg /s /Å /cm^2 /arcs^2) at wavelengths of lambda
%
%                                                       JGR 28 June 2019
%
%   Added option for sky brightness due to moon.        GHOST 4 164,  26 July 2020
%   Hanuschik data apply to SB20, BUT it turns out that they have to have extinction
%   applied (to be done elsewhere)
%
%   Presets
%
    Han_lam_l = [314 374 480 583 670 860];  % nm
    Han_lam_h =[376 486 577 679 856 1043];
    Han_flux = 1e-16*[0.17 0.14 0.09 0.10 0.08 0.07];
    SB_char = {'SB20', 'SB50', 'SB80', 'SBAny'};
    SB_val = [1 1.74 5.25 20.9]; % ratios wrt V = 21.3, 20.7, 19.5, 18 mag/sq arcsec
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
    beta = 100*[-0.000056876550681   0.001533398242066   0.000693822733270 -3.765787850052728];
%    
%   
    flux2 = 1e-16*(beta(2)*exp(beta(1)*(lambda - lambda_offset - beta(4))) + beta(3));
%  
    if plotQ
        for i=1:6
            plot([Han_lam_l(i) Han_lam_h(i)],1e16*[Han_flux(i) Han_flux(i)],'red','LineWidth',2)
            hold on
        end
        plot(lambda,1e16*flux2,'blue','LineWidth',2)
    end 
%
    k = find(strcmpi(SB_char,SB));
    if isempty(k)
        error('Unknown sky brightness designation!')
    end
    flux = flux2*SB_val(k);
%    
    if plotQ
        plot(lambda,1e16*flux,'green','LineWidth',2)
    end
%    
    return   
end

