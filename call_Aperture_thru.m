function poly = call_Aperture_thru(Res,init,final,degree)
%
%   Quick routine to do multiple calls to Aperture_thru and plot fraction vs seeing.
%   Extend it to make the 'definitive' data set and get a polynomial fit to it, for
%   later use in online ETC.
%
%   Input parameters:
%   -----------------
%   Res     : character string: 'SR' for standard resolution, 'HR' for high
%   init    : lowest value of seeing FWHM (in arcsec) for calc and plot
%   final   : highest seeing value
%   degree  : degree of polynomial to fit
%
%   Output parameters:
%   ------------------
%   poly    : coeficients of best-fit polynomial
%
%                                       G. Robertson 20 June 2019
%
    n_step = 30;
    X = zeros(n_step,1);
    Y = zeros(n_step,1);
    Richard_x = 0.3:0.1:1.2; % points read off from Richard McDermid's plot in PDD
    Richard_y = [0.970,0.919,0.849,0.758,0.669,0.589,0.519,0.450,0.399,0.348];
%    
    for i = 1:n_step
        FWHM = init + (final - init)*(i - 1)/(n_step - 1);
        X(i) = FWHM;
        Y(i) = Aperture_thru(FWHM,Res,0);
    end
%    
    plot(X,Y)
    hold on
    plot(Richard_x,Richard_y,'green')
    axis([init,final,0,1])
%
%   Do polynomial fit to the data
%
    poly = polyfit(X,Y,degree);
    Y_fit = polyval(poly,X);
    plot(X,Y_fit,'red')
%    
    return

end

