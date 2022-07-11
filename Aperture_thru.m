function fraction_tx = Aperture_thru(FWHM,Res,plotQ)
%
%   GHOST IFU aperture throughput function. It calculates fractional throughput of
%   either the standard resolution array of 7 hex microlenses, or the high resolution
%   19-element array. It uses the Moffat function for modelling the seeing, and Monte
%   Carlo integration for the throughput of the aperture.
%
%   Included microlens array throughput factors. SR 0.93, HR 0.85.  9 Jan 2020
%   
%   Input parameters:
%   -----------------
%   FWHM    : FWHM of the seeing, arcsec (scalar)
%   Res     : character string: 'SR' for standard resolution, 'HR' for high
%   plotQ   : 1 for diagnostic plot, 0 for no plot (normal operation)
%
%                                       G. Robertson 18 June 2019
%
%   Presets
%
    n_step = 100;         % number of steps along each straight segment of hex boundaries
    n_Moffat = 1000;      % number of points in Moffat function LUT
    n_MC = 1e5;           % number of points in Monte Carlo integration over IFU area
    SR_flat_flat = 0.240; % SR microlens flat to flat dimension, mm
    HR_flat_flat = 0.144; % HR microlens flat to flat dimension, mm
    image_scale = 1.6397; % arcsec/mm on focal surface after the telecentricity lens
    beta = 4.0;           % recommended value of Moffat parameter for natural seeing
    SR_IFU_tx = 0.93;     % SR IFU lenslet throughput factor (Ross, GVR-5124.4)
    HR_IFU_tx = 0.85;     % Equivalent for HR
    FWHM = FWHM/image_scale; % convert seeing FWHM to mm on focal surface
%
    switch Res
        case 'SR'
%
%   Vertices and extrema on scale where length of one hex side = 1
%
    SR_vertices_x = [-0.5,0.5,1,2,2.5,2,2.5,2,1,0.5,-0.5,-1,-2,-2.5,-2,-2.5,-2,-1]';
    SR_vertices_y = sqrt(3)*[-1.5,-1.5,-1,-1,-0.5,0,0.5,1,1,1.5,1.5,1,1,0.5,0,-0.5,-1,-1]';
    SR_max_x = 1.05*2.5*SR_flat_flat/sqrt(3);  % already in mm on focal surface
    SR_max_y = 1.05*1.5*SR_flat_flat;          % ditto. Add a bit extra 
%    
    boundary_r_SR = zeros(18*n_step,1); % vertices get duplicate entries - is OK. This is SR.
    boundary_t_SR = zeros(18*n_step,1);
%
%   Create array giving SR boundary in polar coordinates. There are 18 sectors for SR
%
    for i_sector = 1:18
        if i_sector == 18
            i_end = 1;
        else
            i_end = i_sector + 1;
        end
        X1 = SR_vertices_x(i_sector);  % x coord of start of line segment
        Y1 = SR_vertices_y(i_sector);  % y coord of start of line segment
        X2 = SR_vertices_x(i_end);     % x coord of end of line segment
        Y2 = SR_vertices_y(i_end);     % y coord of end of line segment
%
        j = (1:n_step)';                        % steps along one segment of boundary
        X = X1 + (j-1)*(X2 - X1)/(n_step - 1);  % linear interp follows the boundary
        Y = Y1 + (j-1)*(Y2 - Y1)/(n_step - 1);
        r_lim = sqrt(X.^2 + Y.^2);
        theta = atan2d(Y,X); 
%
%       Load into the overall r,theta array
%
        boundary_r_SR(n_step*(i_sector-1)+1:n_step*i_sector) = r_lim;
        boundary_t_SR(n_step*(i_sector-1)+1:n_step*i_sector) = theta;
    end
%
%   Sort boundary array into order of increasing theta. Its input units correspond to one
%   hex flat having length = 1, and theta is in degrees, range -180 to 180.
%
    [boundary_t_sort,I] = sort(boundary_t_SR); 
    boundary_r_sort = boundary_r_SR(I)*SR_flat_flat/sqrt(3); % r sorted to align with its theta  
%                                                     and r converted to mm in focal surface
    if plotQ
        polar(boundary_t_sort*pi/180,boundary_r_sort)
        hold on
    end
%
%   Make an inline function for Moffat function vs r, so its use can be vectorised. Based
%   on formulas from Richard McDermid in GHOSD-06 PDD. 
%   Normalisation is to peak = 1, not integral = 1
%
    alpha = FWHM/(2*sqrt(2^(1/beta) - 1));
    PSF = @(r) (1 + (r/alpha).^2).^(-beta);      
%
%   Monte Carlo integration of Moffat seeing profile with the IFU's multiple hex boundary
%
    X_MC = random('Uniform',zeros(n_MC,1) -SR_max_x,SR_max_x); % prepare random points
    Y_MC = random('Uniform',zeros(n_MC,1) -SR_max_y,SR_max_y);
    in_MC = zeros(n_MC,1); % will be set to 1 if this MC point is inside IFU boundary
    integral_MC = zeros(n_MC,1); % set to PSF value if inside boundary
%
    r_MC = sqrt(X_MC.^2 + Y_MC.^2);         % convert MC points to polar coords
    theta_MC = atan2d(Y_MC,X_MC);
%
    theta_tol = max(abs(diff(boundary_t_sort)))/1.5; % Allow more search range than 1/2 max diff
%
%   It seems impossible to fully vectorise the comparison of MC points with boundary,
%   so have to use a for loop.
%
    for i = 1:n_MC
        index = find(abs(theta_MC(i) - boundary_t_sort) < theta_tol,1); 
        if r_MC(i) <= boundary_r_sort(index)
            in_MC(i) = 1;
            integral_MC(i) = PSF(r_MC(i));
        end
    end
%
    if plotQ
        plot (X_MC,Y_MC,'k.','MarkerFaceColor','black','MarkerSize',1)
    end
%   
%   Calculate the fraction of PSF power that was within the boundary. To do this,
%   have to allow for effective area element of each MC point, and integral of PSF.
%
    integral_tx = sum(integral_MC)*4*SR_max_x*SR_max_y/n_MC; 
    fraction_tx = SR_IFU_tx*integral_tx*(beta - 1)/(pi*alpha^2);
%     
%
        case 'HR'
%
%   Vertices and extrema on scale where length of one hex side = 1
%
    HR_vertices_x = [-0.5,0.5,1,2,2.5,3.5,4,3.5,4,3.5,4,3.5,2.5,2,1,0.5,-0.5,-1, ...
                    -2,-2.5,-3.5,-4,-3.5,-4,-3.5,-4,-3.5,-2.5,-2,-1]';
    HR_vertices_y = sqrt(3)*[-2.5,-2.5,-2,-2,-1.5,-1.5,-1,-0.5,0,0.5,1,1.5,1.5,2,2,2.5,...
                    2.5,2,2,1.5,1.5,1,0.5,0,-0.5,-1,-1.5,-1.5,-2,-2]';
    HR_max_x = 1.05*4*HR_flat_flat/sqrt(3);  % already in mm on focal surface
    HR_max_y = 1.05*2.5*HR_flat_flat;          % ditto. Add a bit extra          
%    
    boundary_r_HR = zeros(30*n_step,1); % vertices get duplicate entries - is OK. This is HR.
    boundary_t_HR = zeros(30*n_step,1); 
%
%   Create array giving HR boundary in polar coordinates. There are 30 sectors for HR
%
    for i_sector = 1:30
        if i_sector == 30
            i_end = 1;
        else
            i_end = i_sector + 1;
        end
        X1 = HR_vertices_x(i_sector);  % x coord of start of line segment
        Y1 = HR_vertices_y(i_sector);  % y coord of start of line segment
        X2 = HR_vertices_x(i_end);     % x coord of end of line segment
        Y2 = HR_vertices_y(i_end);     % y coord of end of line segment
%
        j = (1:n_step)';                        % steps along one segment of boundary
        X = X1 + (j-1)*(X2 - X1)/(n_step - 1);  % linear interp follows the boundary
        Y = Y1 + (j-1)*(Y2 - Y1)/(n_step - 1);
        r_lim = sqrt(X.^2 + Y.^2);
        theta = atan2d(Y,X); 
%
%       Load into the overall r,theta array
%
        boundary_r_HR(n_step*(i_sector-1)+1:n_step*i_sector) = r_lim;
        boundary_t_HR(n_step*(i_sector-1)+1:n_step*i_sector) = theta;
    end
%
%   Sort boundary array into order of increasing theta. Its input units correspond to one
%   hex flat having length = 1, and theta is in degrees, range -180 to 180.
%
    [boundary_t_sort,I] = sort(boundary_t_HR); 
    boundary_r_sort = boundary_r_HR(I)*HR_flat_flat/sqrt(3); % r sorted to align with its theta  
%                                                     and r converted to mm in focal surface
    if plotQ
        polar(boundary_t_sort*pi/180,boundary_r_sort)
        hold on
    end
%
%   Make an inline function for Moffat function vs r, so its use can be vectorised. Based
%   on formulas from Richard McDermid in GHOSD-06 PDD
%   Normalisation is to peak = 1, not integral = 1
%
    alpha = FWHM/(2*sqrt(2^(1/beta) - 1));
    PSF = @(r) (1 + (r/alpha).^2).^(-beta);     
%
%   Monte Carlo integration of Moffat seeing profile with the IFU's multiple hex boundary
%
    X_MC = random('Uniform',zeros(n_MC,1) -HR_max_x,HR_max_x); % prepare random points
    Y_MC = random('Uniform',zeros(n_MC,1) -HR_max_y,HR_max_y);
    in_MC = zeros(n_MC,1); % will be set to 1 if this MC point is inside IFU boundary
    integral_MC = zeros(n_MC,1); % set to PSF value if inside boundary
%
    r_MC = sqrt(X_MC.^2 + Y_MC.^2);         % convert MC points to polar coords
    theta_MC = atan2d(Y_MC,X_MC);
%
    theta_tol = max(abs(diff(boundary_t_sort)))/1.5; % Allow more search range than 1/2 max diff
%
%   It seems impossible to fully vectorise the comparison of MC points with boundary,
%   so have to use a for loop.
%
    for i = 1:n_MC
        index = find(abs(theta_MC(i) - boundary_t_sort) < theta_tol,1); 
        if r_MC(i) <= boundary_r_sort(index)
            in_MC(i) = 1;
            integral_MC(i) = PSF(r_MC(i));
        end
    end
%
    if plotQ
        plot (X_MC.*in_MC,Y_MC.*in_MC,'k.','MarkerFaceColor','black','MarkerSize',1)
    end
%   
%   Calculate the fraction of PSF power that was within the boundary. To do this,
%   have to allow for effective area element of each MC point, and integral of PSF.
%
    integral_tx = sum(integral_MC)*4*HR_max_x*HR_max_y/n_MC; 
    fraction_tx = HR_IFU_tx*integral_tx*(beta - 1)/(pi*alpha^2);
%        
        otherwise
            error('invalid resolution selection')
    end

%
    return
end

