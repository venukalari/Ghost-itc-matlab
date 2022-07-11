function [frac_out] = Extinc_CTIO(lambda_in,ZD,plotQ)
%
%   Calculate attenuation due to atmospheric extinction, for Gemini South.
%   It currently has to use data from CTIO because that is the best we have. 
%   It is extrapolated beyond 871 nm, but extinction vs wavelength is fairly 
%   linear there.
%   The output is fractonal throughput, *not* magnitudes.
%
%   Input parameters:
%   -----------------
%   lambda_in   : single value or column vector of wavelength values, in *nm*
%   ZD          : zenith distance of the observation, in *degrees*
%   plotQ       : 1 for plot, 0 for no plot
%
%   Output parameters:
%   ------------------
%   frac_out    : Nx2 array -
%                   column 1 : wavelength in nm, same as lambda_in
%                   column 2 : fractional throughput of the atmosphere, range 0 - 1
%
%                                                       JGR 26 May 2019
    [N_in,dim2] = size(lambda_in);
    assert(dim2 == 1,'Input wavelength array is not a column vector!')
    assert(max(lambda_in) <= 1000,'Input wavelengths must be in nm and in range!')
    assert(min(lambda_in) >= 320,'Input wavelengths must be in nm and in range!')
    assert(ZD >= 0 && ZD < 90, 'ZD out of range!')
%
%   Extinction values for CTIO.
%   Data from Gutierrez-Moreno PASP 98 1208 1986 Table 1
%   These are magnitudes per air mass
%
    lambda = [320,325,330,335,340,345,350,357.1,363.6, 370.4,379,386.2,...
        403.6,416.7,425.5,446.4,456.6,467.5,478.5,500,513,526.3,542,555.6,...
        570,584,595,605.6,618,631,643.6,664,679,710,725,740,755,778,789,...
        799,809,818,828,837,870.8]';
    data_Jan = [1.063,.949,.870,.801,.746,.705,.669,.622,.581,.551,.502,...
        0.479,.414,.372,.352,.303,.286,.266,.254,.234,.226,.226,.220,.215,...
        .199,.194,.185,.189,.160,.173,.152,.146,.135,.102,.144,.085,.077,...
        .083,.071,.061,.073,.138,.124,.045,.076]';
    data_Apr = [.952,.844,.758,.697,.647,.613,.579,.535,.498,.466,.432,...
        .397,.343,.304,.283,.234,.219,.204,.189,.162,.157,.166,.138,.137,...
        .136,.110,.137,.115,.104,.122,.092,.076,.048,.070,.071,.034,...
        .039,.038,.043,.009,.025,.049,.013,.030,.007]';
    data_May = [.953,.819,.754,.695,.643,.600,.564,.524,.482,.449,.412,...
        .384,.322,.283,.263,.224,.210,.194,.185,.155,.144,.149,.139,.133,...
        .138,.133,.120,.121,.110,.114,.093,.074,.068,.062,.070,.046,.044,...
        .039,.064,.043,.029,.069,.052,.032,.061]';
    data_Jul = [1.045,.921,.835,.770,.718,.670,.633,.589,.550,.513,.470,...
        .443,.376,.338,.313,.267,.245,.226,.201,.185,.188,.184,.179,.179,...
        .179,.165,.164,.158,.140,.135,.119,.100,.104,.082,.099,.079,.076,...
        .061,.053,.067,.034,.048,.075,.032,.064]';
%
    data_Av = (data_Jan + data_Apr + data_May + data_Jul)/4;
%
%   Do a spline fit to 1 nm spaced points, over whole range. It goes wild at top end
%   but doesn't matter because not used there
%
    xgrid = (320:950)';  
    lambda_mf = lambda.^-4;
    xgrid_mf = xgrid.^-4;
    spline_data=spline(lambda_mf,data_Av,xgrid_mf);
%
    if plotQ
        plot(xgrid,spline_data)
        hold on
        plot(lambda,data_Av,'kd','MarkerFaceColor','black','MarkerSize',3)
    end
%
%   Do linear fit from 570 nm up to 950 but extend the line over whole range
% 
    xpolyin = lambda(27:end);  
    ypolyin = data_Av(27:end); 
    lin_poly = polyfit(xpolyin,ypolyin,1);
    lin_data = lin_poly(1)*xgrid + lin_poly(2);
%
    if plotQ
        plot(xgrid,lin_data)
    end
%
%   Now use a linear ramp of weights to change from spline to linear over range
%   570-750 nm
%
    weight = zeros(631,1);
    weight(251:431) = linspace(0,1,181)';
    weight(432:end) = 1;
    merge_data = (1 - weight).*spline_data + weight.*lin_data; 
%
%   Now use linear interpolation to find values at the desired lambda points
%
    extinc = interp1((320:950),merge_data,lambda_in,'linear');
%  
    if plotQ
        plot(lambda_in,extinc,'green','LineWidth',2)
    end
%
    frac_out = zeros(N_in,2);
    frac_out(:,1) = lambda_in;
    airmass = 1/cosd(ZD);
    frac_out(:,2) = 10.^(-0.4*airmass*extinc); 
%
    return
end

