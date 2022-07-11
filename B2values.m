function  B2values
%
%   Collect results re B^2 values from zmx_an2
%
%   GHOST 3 135  19 July 2019
%
%   Results from GHOST 3 18 and FWHM's from the zmx_an2 runs 
%
    SR_blue_B2 = [2.2205 2.3751 2.4524 2.2449 2.4404 2.5287 2.2177 2.4959 2.5602]';
    SR_blue_FWHM = [3.0354 3.2432 3.3448 3.0687 3.3407 3.4505 3.0481 3.399 3.5002]';
    SR_blue_order = [95 95 95 78 78 78 64 64 64]';
    SR_blue_lambda =[360.35 362.26 364.16 438.38 441.21 444.04 533.53 537.72 541.93]'; % eyeball est
%
    SR_red_B2 = [1.9907 2.4788 2.5355 2.0782 2.4463 2.5356 2.1819 2.411 2.5077]';
    SR_red_FWHM = [2.7424 3.3808 3.468 2.8521 3.3427 3.4632 2.9833 3.2905 3.4204]';
    SR_red_order = [36 36 36 45 45 45 66 66 66]';
    SR_red_lambda = [942.6771 955.9578 969.2385 756.2718 764.7692 773.2667 517.4839 521.4341 525.3843]';
%
%   HR results from GHOST 3 12 and 14; FWHMs and lambdas from last table of CST_BGT_007
%
    HR_blue_B2 = [1.4084 1.4918 1.5399 1.4197 1.5406 1.5898 1.4460 1.5650 1.6154]';
    HR_blue_FWHM = [1.848 1.933 1.991 1.868 2.007 2.058 1.918 2.018 2.118]';
    HR_blue_order = [95 95 95 78 78 78 64 64 64]';
    HR_blue_lambda =[360.347 362.25 364.16 438.38 441.21 444.03 533.52 537.72 541.92]';
%
    HR_red_B2 = [1.3325 1.5631 1.6019 1.3655 1.5367 1.5979 1.3915 1.5171 1.5739]';
    HR_red_FWHM = [1.772 2.029 2.102 1.814 2.000 2.091 1.842 1.968 2.046]';
    HR_red_order = [34 34 34 45 45 45 66 66 66]';
    HR_red_lambda = [997.29 1012.20 1027.10 756.26 764.76 773.25 517.47 521.43 525.38]';


    plot(HR_red_FWHM,HR_red_B2,'rd','MarkerFaceColor','red','MarkerSize',4) 
    hold on
    plot(HR_blue_FWHM,HR_blue_B2,'bd','MarkerFaceColor','blue','MarkerSize',4)   
%    plot(SR_blue_lambda,SR_blue_B2,'kd','MarkerFaceColor','black','MarkerSize',3)
%    plot(SR_blue_lambda,SR_blue_FWHM,'kd','MarkerFaceColor','black','MarkerSize',3)
    poly = polyfit([HR_red_FWHM' HR_blue_FWHM'],[HR_red_B2' HR_blue_B2'],1);
    X = 1.75:0.1:2.15;
    Y = polyval(poly,X);
    plot(X,Y,'black')
    poly

    return
    end

 