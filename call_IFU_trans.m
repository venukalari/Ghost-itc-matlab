function  call_IFU_trans(Res,init,final)
%
%   Quick check of polynomial calc of IFU transmission. Provide a loop call of
%   IFU_trans
%
    n_step = 300;
    X = zeros(n_step,1);
    Y = zeros(n_step,1);
    Richard_x = 0.3:0.1:1.2; % points read off from Richard McDermid's plot in PDD
    Richard_y = [0.970,0.919,0.849,0.758,0.669,0.589,0.519,0.450,0.399,0.348];
%    
    for i = 1:n_step
        FWHM = init + (final - init)*(i - 1)/(n_step - 1);
        X(i) = FWHM;
        Y(i) = IFU_trans(Res,FWHM);
    end
%    
    plot(X,Y)
    hold on
    plot(Richard_x,Richard_y,'green')
    axis([init,final,0,1])    


end

