function  SNplot(lambda1,lambda2,AB_low,AB_high,t_exp,ZD,resoln,seeing,N_mirror,SR_sky,bin_spat,N_sub,col)
                         
%
%   Make plot of GHOST S/N vs AB magnitude, with filled area for range from order
%   centre (maximum) sensitivity to order edge (minimum)
%
%                           G. Robertson    6 Aug 2019; GHOST 3 156
%
%   Input parameters:
%   -----------------
%   lambda1     : Wavelength of an order centre. Single scalar only. In *nm*
%   lambda2     : Wavelength of adjacent order edge. Single scalar only. In *nm*
%   AB_low      : Brightest AB magnitude of object
%   AB_high     : Faintest AB magnitude of object
%   t_exp       : exposure time, seconds. *Total* over N_sub exposures if N_sub > 1
%   ZD          : Zenith distance of object, degrees
%   resoln      : Resolution mode - 'SR' for standard, 'HR' for high resolution
%   seeing      : Seeing disc FWHM, arcseconds
%   N_mirror    : Number of Gemini reflections: 2 for axial port, 3 for side port
%   SR_sky      : For SR only - number of sky microlenses: 3, 7 or 10
%   bin_spat    : Binning factor in spatial direction; 1 for no binning
%   N_sub       : Number of sub-exposures comprising the total t_exp. 1 for single exp
%   col         : colour: 'blue', 'violet', 'orange', 'red'
%
%
%   Modified to create plots for different moon conditions; also streamline the
%   colour selections
%                                                    27 July 2020 GHOST 4 167
%
    colours = [1.0     0.0      0.0; ... %red        row 1
               0.0     1.0      0.0; ... % green         2
               0.0     0.0      1.0; ... % blue          3
               1.0     1.0      0.0; ... % yellow        4
               1.0     0.0      1.0; ... % magenta       5
               0.0     1.0      1.0; ... % cyan          6
               0.0     0.0      0.0; ... % black         7
               0.5     0.5      0.5; ... % medium grey   8
               0.67    0.0      1.0; ... % violet        9
               1.0     0.4      0.0; ... % orange       10
               0.5     0.0      0.0; ... % dark red     11
               0.0     0.5      0.0; ... % dark green   12
               0.8     0.5      0.8; ... % light violet 13
               0.5     0.5      0.9; ... % light blue   14
               1.0     0.6      0.35;... % light orange 15
               0.9     0.5      0.5];    % light red    16
%
    switch col
        case 'blue'
            col1 = colours(3,:);   % the strong one
            col2 = colours(14,:);  % the lighter one
        case 'violet'
            col1 = colours(9,:);
            col2 = colours(13,:);
        case 'orange'
            col1 = colours(10,:);
            col2 = colours(15,:);
        case 'red'
            col1 = colours(1,:);
            col2 = colours(16,:);
        otherwise
            error('Invalid colour specification!')
    end
%   
%   Presets and checks   
%
    n_AB = 101;         % number of separate evaluations between AB limits
    AB = zeros(1,n_AB);
    SN1 = zeros(1,n_AB);
    vars1 = zeros(4,n_AB);
    SN2 = zeros(1,n_AB);
    vars2 = zeros(4,n_AB);
    SB_char = {'SB20', 'SB50', 'SB80', 'SBAny'};
    [dim1,dim2] = size(lambda1);
    assert(dim1 == 1 && dim2 == 1,'lambda1 must be a single scalar for SNplot!')
    [dim1,dim2] = size(lambda2);
    assert(dim1 == 1 && dim2 == 1,'lambda2 must be a single scalar for SNplot!')
%
%   Overall loop - through the 4 sky brightness values
%
    for i_moon = 1:4
%
%   Loop through AB values
%
        for i = 1:n_AB
            AB(i) = AB_low + (i-1)*(AB_high - AB_low)/(n_AB-1);
            [SN1(i), vars1(:,i)]=SNmain(lambda1,AB(i),t_exp,ZD,resoln,seeing,...
                SB_char(i_moon),N_mirror,SR_sky,bin_spat,N_sub,0);
            [SN2(i), vars2(:,i)]=SNmain(lambda2,AB(i),t_exp,ZD,resoln,seeing,...
                SB_char(i_moon),N_mirror,SR_sky,bin_spat,N_sub,0);
        end
    %
        X = [AB AB(end:-1:1)];
        Y = [SN1 SN2(end:-1:1)];

        subplot(2,2,i_moon)
        semilogy(AB,SN1,'color',col1,'LineWidth',2)
        hold on
        semilogy(AB,SN2,'color',col1,'LineWidth',2)
        fill(X,Y,col2,'EdgeColor',col1)
        ylabel('Signal/noise')
        xlabel('AB magnitude')   
        grid on
        title(SB_char(i_moon))
    end
%
%
end

