function frac=Throughput_main(lambda,flag,ZD,resoln,seeing,N_mirror,plotQ)
%
%   Main routine for throughput calculation for GHOST.
%   Enable flexible selection of any combination of loss processes
%
%   Input parameters:
%   -----------------
%   lambda      : wavelength(s) at which S/N is required. Vector or scalar. In *nm*
%   flag:       : 1 - Cass unit, cable 
%                 2 - s/gr and CCD 
%                 3 - all of above
%                 4 - atmos extinction, IFU injection, telescope mirrors
%                 5 - all of above and requirements
%                 6 - extinction only
%                 7 - requirements only
%                 8 - reflection only
%                 9 - IFU only
%                10 - cable + Cass unit only
%                11 - sgr + blaze only
%   ZD          : Zenith distance of object, degrees
%   resoln      : Resolution mode - 'SR' for standard, 'HR' for high resolution
%   seeing      : Seeing disc FWHM, arcseconds
%   N_mirror    : Number of Gemini reflections: 2 for axial port, 3 for side port
%   plotQ       : 1 for plot, 0 for no plot
%
%   Output parameters:
%   ------------------
%   frac        : Throughput at each value of input lambda
%
%                                   G. Robertson  27 July 2019. [GHOST 3 147]
%
%   Presets
%     
    area_SR = 0.939; % SR IFU area, arcsec^2
    area_HR = 0.917; % HR IFU area, arcsec^2
    reqs = [363 375 375 450 450 900 900 950 950 1000; 8 8 12.7 12.7 13.6 13.6 5.5 5.5 2.1 2.1]; 
%    
    lambda = lambda(:).';  % ensures row vector 
    [dim1,n_lambda] = size(lambda);
    assert(dim1 == 1,'lambda is not a scalar or vector!')   
    assert(min(lambda)>=360,'lambda value(s) below blue limit!')
    assert(max(lambda)<=1000,'lambda value(s) above red limit!')
    assert(N_mirror == 2 || N_mirror == 3,'Illegal number of telescope reflections!')
%
%   Set up matrix for turning each loss process on or off. Rows of the matrix are the
%   'flag' values, columns are 1. Extinc 2. Reflec  3. IFU  4. cable  5. sgr  6. Blaze  7. Plot Reqs
%
    select = [0 0 0 1 0 0 0; ...
              0 0 0 0 1 1 0; ...
              0 0 0 1 1 1 0; ...
              1 1 1 0 0 0 0; ...
              1 1 1 1 1 1 1; ...
              1 0 0 0 0 0 0; ...
              0 0 0 0 0 0 1; ...
              0 1 0 0 0 0 0; ...    % ?? The ... missing on 141019. But still worked (!)
              0 0 1 0 0 0 0; ...
              0 0 0 1 0 0 0; ...
              0 0 0 0 1 1 0];
        
    rate1 = ones(1,n_lambda);  
%    
    extinc_frac = Extinc_Paranal(lambda,ZD,0)*select(flag,1) + rate1*~select(flag,1);   % fraction passed
%  
    reflec_frac = GS_reflectivity(lambda,N_mirror,0)*select(flag,2) + rate1*~select(flag,2); % after telescope mirrors;
%
    IFU_frac = IFU_trans(resoln,seeing)*select(flag,3) + rate1*~select(flag,3); % drop by IFU injection loss
%
    cable_frac = GHOST_cable_eta(lambda)*select(flag,4) + rate1*~select(flag,4); % drop by CU and cable throughputs 
%
    sgr_frac = sgr_throughput(lambda,0)*select(flag,5) + rate1*~select(flag,5); % drop by sgr and CCD throughputs 
%
    Blaze_frac = BlazeFunction(lambda)*select(flag,6) + rate1*~select(flag,6); % drop by order blaze functions
%    
    frac = extinc_frac.*reflec_frac.*IFU_frac.*cable_frac.*sgr_frac.*Blaze_frac;
%      
%       
    if plotQ 
        plot(lambda,100*frac,'LineWidth',1.5)
        ylabel('Throughput  /%')
        xlabel('Wavelength  /nm')
        axis([360 1000 0 100]);
        hold on
        if select(flag,7)
            plot(reqs(1,:),reqs(2,:),'green','LineWidth',1.5)
            plot([450],[27],'gd','MarkerSize', 5,'MarkerFaceColor','green')
        end 
    end
return
end

