function [ccv,t,tau,C] = CrossCovariance(groups1,groups2,varargin)

%CrossCovariance - Compute multiple cross-covariances.
%
%  USAGE
%
%    [ccv,t,tau,C] = CrossCovariance(groups1,groups2,<options>)
%
%    groups1      electrode groups located in PFC (e.g [1:8] or [1:4,9:12])
%    groups2      electrode groups located in HPC
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'binSize'     bin size in s (default = 0.01)
%     'duration'    duration in s of each xcorrelogram (default = 2)
%     'smooth'      smoothing size in bins (0 = no smoothing, default)
%     'intervals'   list of [start, stop] pairs defining time intervals used
%     'restrict'    'stim' or 'ripples' to restrict analysis around those 
%                   events. Can be used in combination with 'intervals'
%                   Default = 'none'
%     'mode'        define pairs to analyze ('HPC-PFC', 'HPC', 'PFC' or 'all')
%                   Default = 'HPC-PFC'
%     'alpha'       Value used to determine significantly correlated pairs
%                   Default = 0.05
%    =========================================================================
%
%  OUTPUT
%
%      ccv     value of cross-covariance
%      tau     peak lag time for a given cell pair
%      C       cross-covariance at the peak lag time
%
%  SEE
%
%    See also ShortTimeCCG, CCG

% Copyright (C) 2012 by Michaël Zugaro, Marie Goutierre
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
duration = 2;
binSize = 0.01;
smooth = 0;
intervals = [];
alpha = 0.05;
restrict = 'none';
groups =[];
mode = 'HPC-PFC';

% Check parameters
if nargin < 2,
  error('Incorrect number of parameters (type ''help <a href="matlab:help CrossCovariance">CrossCovariance</a>'' for details).');
end
if ~isdvector(groups1),
	error('Parameter ''groups1'' is not a vector (type ''help <a href="matlab:help CrossCovariance">CrossCovariance</a>'' for details).');
end
if ~isdvector(groups2),
	error('Parameter ''groups2'' is not a vector (type ''help <a href="matlab:help CrossCovariance">CrossCovariance</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help CrossCovariance">CrossCovariance</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'binsize',
			binSize = varargin{i+1};
			if ~isdscalar(binSize,'>0'),
				error('Incorrect value for property ''binSize'' (type ''help <a href="matlab:help CrossCovariance">CrossCovariance</a>'' for details).');
			end
		case 'duration',
			duration = varargin{i+1};
			if ~isdscalar(duration,'>0'),
				error('Incorrect value for property ''duration'' (type ''help <a href="matlab:help CrossCovariance">CrossCovariance</a>'' for details).');
			end
		case 'smooth',
			smooth = varargin{i+1};
			if ~isdscalar(smooth,'>=0'),
				error('Incorrect value for property ''smooth'' (type ''help <a href="matlab:help CrossCovariance">CrossCovariance</a>'' for details).');
			end
		case 'intervals',
			intervals = varargin{i+1};
			if ~isdmatrix(intervals) || size(intervals,2) ~= 2,
				error('Incorrect value for property ''intervals'' (type ''help <a href="matlab:help CrossCovariance">CrossCovariance</a>'' for details).');
			end
		case 'restrict',
			restrict = varargin{i+1};
			if ~isstring(restrict,'ripples','stim','none'),
				error('Incorrect value for property ''restrict'' (type ''help <a href="matlab:help CrossCovariance">CrossCovariance</a>'' for details).');
			end
		case 'mode',
			mode = varargin{i+1};
			if ~isstring(mode,'HPC-PFC','HPC','PFC','all'),
				error('Incorrect value for property ''mode'' (type ''help <a href="matlab:help CrossCovariance">CrossCovariance</a>'' for details).');
			end
		case 'alpha',
			alpha = varargin{i+1};
			if ~isdscalar(alpha,'>0'),
				error('Incorrect value for property ''alpha'' (type ''help <a href="matlab:help CrossCovariance">CrossCovariance</a>'' for details).');
			end
	end
end

% Get spikes times, ID numbers and groups numbers

if strcmp(mode,'HPC-PFC'),
	HPC = GetSpikeTimes([groups2' -1*ones(length(groups2),1)],'output','numbered');
	PFC = GetSpikeTimes([groups1' -1*ones(length(groups1),1)],'output','numbered');
	nPFC = max(PFC(:,2));
	ID = [PFC(:,2);HPC(:,2)+nPFC];
	groups = [2*ones(size(PFC,1),1);ones(size(HPC,1),1)];
	times = [PFC(:,1);HPC(:,1)];
	data = [times ID groups];
	data = sortrows(data);
	times = data(:,1);
	ID = data(:,2);
	groups = data(:,3);

elseif strcmp(mode,'HPC'),
	HPC = GetSpikeTimes([groups2' -1*ones(length(groups2,1))],'output','numbered');
	HPC = sortrows(HPC);
	ID = HPC(:,2);
	times = HPC(:,1);

elseif strcmp(mode,'PFC'),
	PFC = GetSpikeTimes([groups1' -1*ones(length(groups1,1))],'output','numbered');
	PFC = sortrows(PFC);
	ID = PFC(:,2);
	times = PFC(:,1);

elseif strcmp(mode,'all')
	spikes = GetSpikeTimes([1:16]' -1*ones(16,1),'output','numbered');
	spikes = sortrows(spikes);
	times = spikes(:,1);
	ID = spikes(:,2);
end

if strcmp(restrict,'ripples'),
	ripples = GetEvents('Ripple peak.*');
	if ~isempty(intervals),
		ripples = ripples(InIntervals(ripples,intervals),:);
	end
	ripplesWindow = [ripples-0.25 ripples+0.25];
	intervals = ripplesWindow;

elseif strcmp(restrict,'stim'),
	stim = GetEvents('0');
	if ~isempty(intervals),
		stim = stim(InIntervals(stim,intervals),:);
	end
	stimWindow = [stim-0.25 stim+0.25];
	intervals = stimWindow;
end


%Determine Quiet Periods
positions = GetPositions;
velocity = LinearVelocity(positions,5);
[periods,quiescence] = QuietPeriods(velocity,0.005,10,3);

% Restrict signal to quiet periods
ID = ID(InIntervals(times,periods));
if ~isempty(groups),
	groups = groups(InIntervals(times,periods));
end
times = times(InIntervals(times,periods));

% Optionnally, restrict signal to a certain period
if ~isempty(intervals),
	ID = ID(InIntervals(times,intervals));
	if ~isempty(groups),
		groups = groups(InIntervals(times,intervals));
	end
	times = times(InIntervals(times,intervals));
end


% Compute CCGs
[ccg,t] = CCG(times,ID,'groups',groups,'binSize',binSize,'duration',duration,'smooth',smooth);

if ~isempty(groups),
	Group1 = unique(ID(groups == 1));
	Group2 = unique(ID(groups == 2));

	nGroups = max(max(Group1),max(Group2));
else
	nGroups = max(ID);
end

% Determine mean Firing Rate for each neuron
FiringRate = zeros(nGroups,1);

if ~isempty(intervals)
	TotalTime = sum(intervals(:,2)-intervals(:,1));
else
	TotalTime = max(times)-min(times);
end

for i = 1:nGroups,
	FiringRate(i) = sum(ID == i)/TotalTime;
end


% Determine standardized cross-covariances
ccv = zeros(size(ccg));
tau = zeros(size(ccg,2),size(ccg,3));
C = zeros(size(ccg,2),size(ccg,3));

nPairs = size(ccg,2)*size(ccg,3);
disp(['Number of pairs: ' num2str(nPairs) ' pairs.']);

threshold = sqrt(2)*erfinv(1-(alpha/length(t)));

for i = 1:size(ccg,2),
	for j = 1:size(ccg,3),
		% Compute and normalize CCVs from CCGs
		if ~isempty(groups),
			FR = FiringRate(Group1(i))*FiringRate(Group2(j));
		else
			FR = FiringRate(i)*FiringRate(j);
		end

		EstimateCCV = ccg(:,i,j)/(binSize*TotalTime)-FR;
		ccv(:,i,j) = sqrt((binSize*TotalTime)/FR)*EstimateCCV;

		% Smooth it with a 3-bin boxcar
		data = ccv(:,i,j);
		top = flipud(data(1:size(ccg,1),:));
		bottom = flipud(data(end-size(ccg,1)+1:end,:));
		data = [top;data;bottom];
		tmp = filter([1 1 1],3,data);
		n = size(tmp,1);
		d = n - size(ccg,1);
		start = d/2+1;
		stop = start + size(ccg,1) - 1;
		ccv(:,i,j) = tmp(start:stop);

		% Keep only the significantly correlated pairs
		absCCV = abs(ccv(:,i,j));
		if ~any(absCCV>threshold),
			ccv(:,i,j) = zeros(size(ccg,1),1);
		end

		% Find the peak lag time and the peak value
		[maxValue,maxIndex] = max(ccv(:,i,j));
		tau(i,j) = t(maxIndex);
		C(i,j) = maxValue;
	end
end


Cov = C(:);
tau = tau(:);
tau = tau(Cov ~=0);

disp(['Number of significantly correlated pairs: ' num2str(length(tau)) ' pairs.']);


