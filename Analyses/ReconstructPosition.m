function [stats,lambda,Px] = ReconstructPosition(positions,spikes,phases,varargin)

%ReconstructPosition - Bayesian reconstruction of positions from spike trains.
%
% Instantaneous positions are reconstructed using a Bayesian algorithm.
% Instantaneous population firing rates can be estimated either over fixed time
% windows, or over fractions of the theta cycle (or of any other brain rhythm).
% Similarly, positions will be reconstructed either over time or phase windows.
% The model is first trained on a subset of the data, then tested on the rest.
%
% USAGE
%
%    [stats,lambda,Px] = ReconstructPosition(positions,spikes,phases,<options>)
%
%    positions      linear or two-dimensional positions <a href="matlab:help samples">samples</a>, in [0..1]
%    spikes         list of (t,ID) couples (obtained via e.g. <a href="matlab:help GetSpikes">GetSpikes</a>,
%                   using numbered output) 
%    phases         optional unwrapped phase <a href="matlab:help samples">samples</a> of the LFP (see <a href="matlab:help Phase">Phase</a>)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'training'    time interval over which the model should be trained
%                   (see NOTE below for defaults)
%     'window'      length of the time or phase window (default = 0.020 s for
%                   time, and pi/3 for phases)
%     'type'        two letters (one for X and one for Y) indicating which
%                   coordinates are linear ('l') and which are circular ('c')
%                   - for 1D data, only one letter is used (default 'll')
%     'nBins'       firing curve or map resolution (default = [200 200])
%     'mode'        perform only training ('train'), only reconstruction
%                   ('test'), or both ('both', default)
%     'lambda'      to provide previously generated model ('test' mode)
%     'Px'          to provide previously generated model ('test' mode)
%    =========================================================================
%
%   OUTPUT
%
%     stats.positions     real position across time or phase windows
%     stats.spikes        cell firing vector across time or phase windows
%     stats.estimations   estimated position across time or phase windows
%     stats.errors        estimation error across time or phase windows
%     stats.average       average estimation error in each phase window
%     stats.windows       time windows (possibly computed from phases)
%     stats.phases        phase windows (empty for fixed time windows)
%
%     lambda              firing map for each unit
%     Px                  occupancy probability map
%
%   NOTE
%
%     Positions and spikes are interpreted differently depending on the mode:
%
%      - For 'train', all positions and spikes are used to train the model
%      - For 'test', all positions and spikes are used to test the model, and
%        positions are optional (e.g. for reconstruction during sleep)
%      - For 'both', the optional parameter 'training' can be used to indicate
%        the training interval (default = first half of the position data)
%

% Copyright (C) 2012-2015 by Michaël Zugaro, (C) 2012 by Karim El Kanbi (initial, non-vectorized implementation),
% (C) 2015 by Céline Drieu (separate training vs test), (C) 2015 by Ralitsa Todorova (log-exp fix)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Defaults
wt = 0.020; % default time window
wp = pi/3; % default phase window
window = [];
nBins = 200;
training = 0.5;
type = '';
nDimensions = 1;
mode = 'both';

% Optional parameter 'phases'
if nargin == 2,
	phases = [];
elseif nargin >= 3 && ischar(phases),
	varargin = {phases,varargin{:}};
	phases = [];
end

% Check number of parameters
if nargin < 2 || mod(length(varargin),2) ~= 0,
	builtin('error','Incorrect number of parameters (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).');
end

% Check parameter sizes
if ~isempty(positions) && ~isdmatrix(positions),
	builtin('error','Incorrect positions (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).');
end
if ~isdmatrix(spikes,'@2') && ~isdmatrix(spikes,'@3'),
	builtin('error','Incorrect spikes (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).');
end
if size(positions,2) >= 3,
	nDimensions = 2;
end
if ~isempty(phases) && ~isdmatrix(phases),
	builtin('error','Incorrect value for property ''phases'' (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).');
end

% Parse parameters
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		builtin('error',['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'training',
			training = varargin{i+1};
			if ~isdvector(training,'<'),
				builtin('error','Incorrect value for property ''training'' (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).');
			end
		case 'window',
			window = varargin{i+1};
			if ~isdscalar(window,'>0'),
				builtin('error','Incorrect value for property ''window'' (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).');
			end
		case 'show',
			show = varargin{i+1};
			if ~isastring(show,'on','off'),
				builtin('error','Incorrect value for property ''show'' (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).');
			end
		case 'nbins',
			nBins = varargin{i+1};
			if isiscalar(nBins),
				builtin('error','Incorrect value for property ''nBins'' (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).');
			end
		case 'type',
			type = lower(varargin{i+1});
			if (nDimensions == 1 && ~isastring(type,'cc','cl','lc','ll')) || (nDimensions == 2 && ~isastring(type,'ccl','cll','lcl','lll','ccc','clc','lcc','llc')),
				builtin('error','Incorrect value for property ''type'' (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).');
			end
		case 'mode',
			mode = lower(varargin{i+1});
			if ~isastring(mode,'both','train','test'),
				builtin('error','Incorrect value for property ''mode'' (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).');
			end
		case 'lambda',
			lambda = varargin{i+1};
			if ~isnumeric(lambda) || length(size(lambda)) ~= 3,
				builtin('error','Incorrect value for property ''lambda'' (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).');
			end
		case 'px',
			Px = varargin{i+1};
			if ~isdvector(Px),
				builtin('error','Incorrect value for property ''Px'' (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).');
			end
		otherwise,
			builtin('error',['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).']);
	end
end

% Defaults (window)
if isempty(window),
	if isempty(phases),
		window = wt;
	else
		window = wp;
		if ~isiscalar((2*pi)/window),
			builtin('error',['Incorrect phase window: not an integer fraction of 2pi (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).']);
		end
	end
end
% Defaults (training)
if isastring(mode,'train','both') && isdscalar(training),
	training = [-Inf positions(1,1)+training*(positions(end,1)-positions(1,1))];
end
% Defaults (type)
if isempty(type),
	if nDimensions == 2,
		type = 'lll';
	else
		type = 'll';
	end
end
% Defaults (nBins)
nBinsX = nBins(1);
if length(nBins) > 2,
	nBinsY = nBins(2);
else
	if nDimensions == 2,
		nBinsY = nBinsX;
	else
		nBinsY = 1;
	end
end
% Defaults (mode)
if isastring(mode,'train','both') && ( exist('lambda','var') || exist('Px','var') ),
	warning(['Inconsistent inputs, lambda and Px will be ignored in mode ''' mode ''' (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).']);
	clear('lambda');clear('Px');
end
% Convert from legacy format for backward compatibility with previous versions of the code (spikes)
if isdmatrix(spikes,'@3'),
	% List units, assign them an ID (number them from 1 to N), and associate these IDs with each spike
	% (IDs will be easier to manipulate than (group,cluster) pairs in subsequent computations)
	[units,~,i] = unique(spikes(:,2:end),'rows');
	nUnits = length(units);
	index = 1:nUnits;
	id = index(i)';
	spikes = [spikes(:,1) id];
	warning('Spikes were provided as Nx3 samples - this is now obsolete (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).');
	if ~strcmp(mode,'both'),
		builtin('error','Obsolete format can be used only when training and test are performed together (type ''help <a href="matlab:help ReconstructPosition">ReconstructPosition</a>'' for details).');
	end
else
	if strcmp(mode,'test'),
		nUnits = size(lambda,3);
	else
		nUnits = max(spikes(:,2));
	end
end

% TRAINING

if isastring(mode,'both','train'),

	% Split data (training vs test)
	if strcmp(mode,'train'),
		% Mode = 'train', use all data
		trainingPositions = positions;
		trainingSpikes = spikes;
	else
		% Mode = 'both', used info from 'training' parameter
		trainingPositions = Restrict(positions,training);
		trainingSpikes = Restrict(spikes,training);
	end
	
	% Compute average firing probability lambda for each unit (i.e. firing maps)
	for i = 1:nUnits,
		unit = trainingSpikes(:,2) == i;
		s = trainingSpikes(unit,1);
		map = Map(trainingPositions,s,'nbins',nBins,'smooth',5,'type',type);
		lambda(:,:,i) = map.z;
	end

	% Compute occupancy probability P(x) (i.e. normalized occupancy map)
	Px = map.time;
	Px = Px ./ sum(Px(:));
	
end

% TEST

if strcmp(mode,'train'),
	stats.estimations = [];
	stats.spikes = [];
   stats.errors = [];
   stats.average = [];
   stats.windows = [];
   stats.phases = [];
   return
end

% Split data (training vs test)
if strcmp(mode,'test'),
	% Mode = 'test', use all data
	testPositions = positions;
	testSpikes = spikes;
else
	% Mode = 'both', used info from 'training' parameter
	testPositions = positions(~InIntervals(positions,training),:);
	testSpikes = spikes(~InIntervals(spikes,training),:);
end

% Determine time windows (using unwrapped phases if necessary)
if ~isempty(phases),
	testPhases = phases(~InIntervals(phases,training),:);
	if ~isempty(testPositions),
		drop = testPhases(:,1) < testPositions(1,1);
		testPhases(drop,:) = [];
	end
	startPhase = ceil(testPhases(1,2)/(2*pi))*2*pi;
	stopPhase = floor(testPhases(end,2)/(2*pi))*2*pi;
	windows = (startPhase:window:stopPhase)';
	stats.phases = windows;
	windows = Interpolate(testPhases(:,[2 1]),windows);
	windows = [windows(1:end-1,2) windows(2:end,2)];
else
	stats.phases = [];
	if ~isempty(testPositions),
		windows = (testPositions(1,1):window:testPositions(end,1))';
	else
		windows = (testSpikes(1,1):window:testSpikes(end,1))';
	end
	windows = [windows(1:end-1) windows(2:end)];
end
nWindows = size(windows,1);

stats.estimations = nan(nBinsY,nBinsX,nWindows);
stats.spikes = zeros(nUnits,nWindows);
% Loop over data windows
for i = 1:nWindows,

	% Get spikes for this window
	s = Restrict(testSpikes,windows(i,:));

	if isempty(s),
		% No spikes: set uniform probability
		stats.estimations(:,:,i) = ones(nBinsY,nBinsX,1)/(nBinsX*nBinsY);
		continue;
	end

	% Population spike count vector
	stats.spikes(:,i) = Accumulate(s(:,2),1,nUnits);
	% To avoid 'for' loops, prepare for vector computation:
	% assign a spike count to each position and unit (3D array)
	n = reshape(repmat(stats.spikes(:,i),1,nBinsX*nBinsY)',nBinsY,nBinsX,nUnits);

	% For each cell i, compute P(ni|x) using a Poisson model. The direct formula is:
	% 	 Pnix = (dt*lambda).^n./factorial(n).*exp(-dt*lambda); 
	% However, large values of (dt*lambda).^n can create overflow erros, so instead we compute
	% the log and then take the exponential (fix by Ralitsa Todorova)
	dt = windows(i,2) - windows(i,1);
	Pnix = exp(n.*log(dt*lambda)-logfactorial(n)-dt*lambda);
	% Compute P(n|x) assuming independent probabilities across units (hmm...)
	% i.e. P(n|x) = product over i of P(ni|x)
	Pnx = prod(Pnix,3);

	% Compute P(n) = sum over x of P(n|x)*P(x)
	Pn = sum(sum(Pnx.*Px));

	% Compute P(x|n) = P(n|x)*P(x)/P(n)
	Pxn = Pnx .* Px / Pn;

	% Store result
	stats.estimations(:,:,i) = Pxn;

end
stats.estimations = squeeze(stats.estimations);
stats.windows = windows;

% Estimation error

stats.errors = [];
stats.average = [];
if nDimensions == 1 && ~isempty(testPositions),
	% Bin test positions and compute distance to center
	stats.positions = Interpolate(testPositions,windows(:,1));
	stats.positions(:,2) = Bin(stats.positions(:,2),[0 1],nBinsX);
	dx = (round(nBinsX/2)-stats.positions(:,2))';
	% Shift estimated position by the real distance to center
	stats.errors = CircularShift(stats.estimations(:,1:length(dx)),dx);
	% Average over one or more cycles
    if ~isempty(phases),
        k = 2*pi/window;
        n = floor(size(stats.errors,2)/k)*k;
        stats.average = reshape(stats.errors(:,1:n),nBins,k,[]);
        stats.average = nanmean(stats.average,3);
    end
else
	warning('Computation of estimation error not yet implemented for 2D environments');
end
end

function data = logfactorial(data);

% We compute log(n!) as the sum of logs, i.e. log(n!) = sum log(i) for i=1:n
% First determine the largest n in the array
m = max(data(:));
% Create a look-up vector of sum log(i) for each i up to the largest n
sums = [0 cumsum(log(1:m))];
% Look-up the value for each item in the array
data(:) = sums(data+1);
end
