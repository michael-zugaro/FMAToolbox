function PlotLinkage(linkage,varargin)

%PlotLinkage - Plot network clustering obtained e.g. by FCA analysis.
%
%  USAGE
%
%    PlotLinkage(linkage,<options>)
%
%  INPUT
%
%    linkage        network similarity matrix obtained using <a href="matlab:help FunctionalClustering">FunctionalClustering</a>
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'threshold'   minimum value for significant linkage (default = 1)
%    =========================================================================
%
%  SEE
%
%  See also FunctionalClustering.

% Copyright (C) 2015 by Ralitsa Todorova, Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
threshold = 1;

if nargin < 1,
	error('Incorrect number of parameters (type ''help <a href="matlab:help PlotLinkage">PlotLinkage</a>'' for details).');
end
if size(linkage,2) < 3 || ~isimatrix(linkage(:,1:2)) || ~isdvector(linkage(:,3)),
	error('Incorrect linkage (type ''help <a href="matlab:help PlotLinkage">PlotLinkage</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help PlotLinkage">PlotLinkage</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'threshold',
			threshold = varargin{i+1};
			if ~isdscalar(threshold),
				error('Incorrect value for property ''threshold'' (type ''help <a href="matlab:help PlotLinkage">PlotLinkage</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help PlotLinkage">PlotLinkage</a>'' for details).']);

	end
end

pairs = linkage(:,1:2);
linkage = linkage(:,3);

% We will renumber IDs from 1 to N because 'dendrogram' cannot handle missing IDs
% We first save the original IDs, then renumber them
original = pairs;
[id,~,i] = unique(pairs(:),'rows');
index = 1:length(id);
pairs(:) = index(i);

% Plot ranked dendrogram
ranks = 1:size(pairs,1);
[~,~,labels] = dendrogram([pairs ranks'],0);
% Update x labels to show the original IDs
tmp(pairs(:)) = original(:);
labels = tmp(labels);
set(gca,'xticklabel',labels,'ytick',[]);
% Adjust axes and plot threshold
yLim1 = ylim;
ylim([0 yLim1(2)]);
PlotIntervals([find(linkage<threshold,1,'first')-0.5 yLim1(2)],'direction','h');
xlabel('ID');
title('(grey area not significant)');

% Plot side panel
SideAxes(gca,'left',0.3);
plot(ranks,linkage,'k.:');
box off;
xlim(yLim1);
yLim2 = ylim;
ylim([yLim2(1) max([threshold+0.5 yLim2(2)])]);
PlotIntervals([find(linkage<threshold,1,'first')-0.5 yLim1(2)],'direction','v');
PlotHVLines(threshold,'h','r');
set(gca,'xticklabel',id);
xlabel('Pair #');
ylabel('linkage');
title('(grey area not significant)');
