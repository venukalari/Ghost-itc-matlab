function effic = GHOST_cable_eta(lambda)
%
%   Provide value(s) of throughput for GHOST Cass unit and science cable, based on
%   results from Ross Zhelem 21 May 2019. The throughput includes Cass unit
%   (telecentricity lens and ADC) and fibre cable + microlens arrays at both ends,
%   but not any aperture losses. This routine uses spline fits to interpolate.
%
%   Input parameters:
%   -----------------
%   
%   lambda_in  : scalar or vector of wavelengths for throughput calculation, in *nm*
%
%   Output parameters:
%   ------------------
%
%   effic      : row vector giving throughput value, range 0 - 1, same size as lambda
%
%   
    lambda = lambda(:).';  % ensures row vector 
    [dim1,~] = size(lambda);
    assert(dim1 == 1,'lambda is not a scalar or vector!')
    assert(min(lambda)>=350,'lambda value(s) below blue limit!')
    assert(max(lambda)<=1000,'lambda value(s) above red limit!')
%
    load('RefData','CassU_cable')  % the data from Ross's printout 210519
%    
    effic = spline(CassU_cable(:,1),CassU_cable(:,2),lambda);
    
%    
    return
end

