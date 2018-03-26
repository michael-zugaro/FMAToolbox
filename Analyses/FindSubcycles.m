function [intervals,subcycle,cycle] = FindSubcycles(angles,nSubcycles)

%FindSubcycles - Find timestamps for subcycles of an oscillating signal.
%
%  USAGE
%
%    [intervals,subcycle,cycle] = FindSubcycles(angles,nSubcycles)
%
%    angles         phase <a href="matlab:help samples">samples</a>
%    nSubcycles     number of subcycles
%
%  EXAMPLE
%
%    % Split theta cycles in six (i.e. subcycles of size pi/3)
%    angles = Phase(FilterLFP(lfp,'passband','theta'));
%    [intervals,subcycle] = FindSubcycles(angles,6);
%
%  SEE
%
%    See also Phase.

% Copyright (C) 2016-2017 by Ralitsa Todorova, Micha??l Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Make sure angles are unwrapped
if range(angles) < 2 * pi + 1,
	angles = unwrap(angles);
end

% Window size
window = 2*pi/nSubcycles;

% Start at cycle number 0
angles(:,2) = angles(:,2) - floor(angles(1,2)/(2*pi))*2*pi;

% Find unwrapped phases corresponding to subcycle boundaries
startPhase = ceil(angles(1,2)/(2*pi))*2*pi;
stopPhase = floor(angles(end,2)/(2*pi))*2*pi;
intervals = (startPhase:window:stopPhase)';
intervals = [intervals(1:end-1) intervals(2:end)];

% Compute cycle and subcycle IDs
n = size(intervals,1);
cycle = ceil(nanmean(intervals,2)/(2*pi));
subcycle = CumSum(ones(n,1),[false;diff(cycle)~=0]);

% Convert boundaries from phases to time
% (ignore aberrant non-unique values = identical unwrapped phases for two timestamps)
[~,ok] = unique(angles(:,2));
angles = angles(ok,:);
intervals = interp1(angles(:,2),angles(:,1),intervals);
