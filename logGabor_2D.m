function [filter] = logGabor_2D(im,nornt,ns,scalefac,sigma)


nscale          = ns;          % Number of wavelet scales.    
norient         = nornt;       % Number of filter orientations.
minWaveLength   = 4;           % Wavelength of smallest scale filter.    
r               = scalefac;    % Scaling factor between successive filters.    
sigmaOnf        = sigma;       % Ratio of the standard deviation of the
                               % Gaussian describing the log Gabor filter's
                               % transfer function in the frequency domain
                               % to the filter center frequency.    
dThetaOnSigma   = 1.5;   % Ratio of angular interval between filter orientations    
                         % and the standard deviation of the angular Gaussian
                         % function used to construct filters in the
                         % freq. plane.

filter = cell(nscale, norient);                               
thetaSigma = pi/norient/dThetaOnSigma;  % Calculate the standard deviation of the
                                        % angular Gaussian function used to
                                        % construct filters in the freq. plane.

[rows,cols] = size(im);

% Pre-compute some stuff to speed up filter construction
% Set up X and Y matrices with ranges normalised to +/- 0.5
% The following code adjusts things appropriately for odd and even values
% of rows and columns.
if mod(cols,2)
    xrange = [-(cols-1)/2:(cols-1)/2]/(cols-1);
else
    xrange = [-cols/2:(cols/2-1)]/cols; 
end

if mod(rows,2)
    yrange = [-(rows-1)/2:(rows-1)/2]/(rows-1);
else
    yrange = [-rows/2:(rows/2-1)]/rows; 
end

[x,y] = meshgrid(xrange, yrange);

radius = sqrt(x.^2 + y.^2);       % Matrix values contain *normalised* radius from centre.
theta = atan2(-y,x);              % Matrix values contain polar angle.
                                  % (note -ve y is used to give +ve
                                  % anti-clockwise angles)
                                  
radius = ifftshift(radius);       % Quadrant shift radius and theta so that filters
theta  = ifftshift(theta);        % are constructed with 0 frequency at the corners.
radius(1,1) = 1;                  % Get rid of the 0 radius value at the 0
                                  % frequency point (now at top-left corner)
                                  % so that taking the log of the radius will 
                                  % not cause trouble.

sintheta = sin(theta);
costheta = cos(theta);
clear x; clear y; clear theta;    % save a little memory

% Filters are constructed in terms of two components.
% 1) The radial component, which controls the frequency band that the filter
%    responds to
% 2) The angular component, which controls the orientation that the filter
%    responds to.
% The two components are multiplied together to construct the overall filter.

% Construct the radial filter components...

% First construct a low-pass filter that is as large as possible, yet falls
% away to zero at the boundaries.  All log Gabor filters are multiplied by
% this to ensure no extra frequencies at the 'corners' of the FFT are
% incorporated as this seems to upset the normalisation process when
% calculating phase congrunecy.
lp = lowpassfilter([rows,cols],.45,15);   % Radius .45, 'sharpness' 15

logGabor = cell(1,nscale);

% crr = 1;
for s = 1:nscale
    if size(r,2)> 1
%         crr = crr * r(s);
        wavelength = minWaveLength * r(s);
    else
        wavelength = minWaveLength * r^(s-1);        
    end
    fo = 1.0/wavelength;                  % Centre frequency of filter.
    logGabor{s} = exp((-(log(radius/fo)).^2) / (2 * log(sigmaOnf)^2));  
    logGabor{s} = logGabor{s}.*lp;        % Apply low-pass filter
    logGabor{s}(1,1) = 0;                 % Set the value at the 0 frequency point of the filter
                                          % back to zero (undo the radius fudge).
end

% Then construct the angular filter components...

spread = cell(1,norient);

for o = 1:norient
  angl = (o-1)*pi/norient;           % Filter angle.

  % For each point in the filter matrix calculate the angular distance from
  % the specified filter orientation.  To overcome the angular wrap-around
  % problem sine difference and cosine difference values are first computed
  % and then the atan2 function is used to determine angular distance.

  ds = sintheta * cos(angl) - costheta * sin(angl);    % Difference in sine.
  dc = costheta * cos(angl) + sintheta * sin(angl);    % Difference in cosine.
  dtheta = abs(atan2(ds,dc));                          % Absolute angular distance.
  spread{o} = exp((-dtheta.^2) / (2 * thetaSigma^2));  % Calculate the
                                                       % angular filter component.
end
% The main loop...

for o = 1:norient                    % For each orientation.
%   angl = (o-1)*pi/norient;           % Filter angle.
  for s = 1:nscale,                  % For each scale.
    filter{s,o} = logGabor{s} .* spread{o};   % Multiply radial and angular
  end
end