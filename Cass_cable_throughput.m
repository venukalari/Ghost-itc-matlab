function frac_out = Cass_cable_throughput(lambda_in,plotQ)
%
%   Throughput of GHOST Cass unit, IFUs, and cable vs wavelength.
%   Based on data from GVR-5124.4 by J Bassett and Ross Zhelem, 20 May 2019
%   Wavelengths 350 - 1000 nm. 
%
%   Input parameters:
%   -----------------
%   lambda_in   : single value or column vector of wavelength values, in *nm*
%   plotQ       : 1 for plot, 0 for no plot
%
%   Output parameters:
%   ------------------
%   frac_out    : Nx2 array -
%                   column 1 : wavelength in nm, same as lambda_in
%                   column 2 : fractional throughput of the mirrors, range 0 - 1
%
%                                                       JGR 17 June 2019
%
%   Presets
%
    lam_min = 350;
    lam_max = 1000;
%    
    load('RefData','CassU_cable')
    in_data = CassU_cable;      % allow simple change of input file namne later ...
    [N_in,dim2] = size(lambda_in);
    assert(dim2 == 1,'Input wavelength array is not a column vector!')
    assert(max(lambda_in) <= lam_max,'Input wavelengths must be in nm and in range!')
    assert(min(lambda_in) >= lam_min,'Input wavelengths must be in nm and in range!')
%   
%   Do a spline fit to 1 nm spaced points, over whole range. 
%
    xgrid = (lam_min:lam_max)';  
    spline_data=spline(in_data(:,1),in_data(:,2),xgrid);
%
    if plotQ
        plot(xgrid,spline_data)
        hold on
        plot(in_data(:,1),in_data(:,2),'kd','MarkerFaceColor','black','MarkerSize',2)
    end
%
%   Now use linear interpolation to find values at the desired lambda points
%
    thruput = interp1((lam_min:lam_max),spline_data,lambda_in,'linear');
%
    frac_out = zeros(N_in,2);
    frac_out(:,1) = lambda_in;
    frac_out(:,2) = thruput; 
%
    if plotQ
        plot(lambda_in,frac_out(:,2),'green','LineWidth',2)
    end
%    
    return   
end

