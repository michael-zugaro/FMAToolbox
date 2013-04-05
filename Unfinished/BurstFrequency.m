%BurstFrequency - Compute spike train auto-correlogram and autospectrum
%
%  USAGE
%
%    [ccg,t] = BurstFrequency(spikes,<options>)
%
%    spikes         spike spikes
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'binSize'     bin size in s (default = 0.01)
%     'duration'    duration in s of each xcorrelogram (default = 2)
%    =========================================================================
%
%  SEE
%
%    See also PlotBurstFrequencies.

% Copyright (C) 2012 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function stats = BurstFrequency(spikes,varargin)

% Check parameters
if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help BurstFrequency">BurstFrequency</a>'' for details).');
end
if ~isvector(spikes),
	error('Parameters ''spikes'' must be a vector (type ''help <a href="matlab:help BurstFrequency">BurstFrequency</a>'' for details).');
end
spikes = spikes(:);

% Compute autocorrelogram
[ccg,t] = CCG(spikes,1,'duration',10);

% Smooth and find maxima
smoothed = [t Smooth(ccg(:,1,1),2)];
maxima = t(IsExtremum(smoothed));

% Compute autospectrum
[s,f] = MTSpectrum(smoothed);