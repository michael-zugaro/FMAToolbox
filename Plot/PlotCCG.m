function PlotCCG(t,ccg,varargin)

%PlotCCG - Plot auto/cross-correlograms of point processes.
%
%  USAGE
%
%    PlotCCG(t,ccg)
%
%    t              time bins obtained using <a href="matlab:help CCG">CCG</a>
%    ccg            data obtained using <a href="matlab:help CCG">CCG</a>
%
%  SEE
%
%    See also CCG, ShortTimeCCG.

% Copyright (C) 2010-2012 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if nargin < 2,
	error('Incorrect number of parameters (type ''help <a href="matlab:help PlotCCG">PlotCCG</a>'' for details).');
end

if ~isdvector(t),
	error('Incorrect times bins (type ''help <a href="matlab:help PlotCCG">PlotCCG</a>'' for details).');
end
if ~isdvector(ccg) && ndims(ccg) ~= 3,
	error('Incorrect cross-correlogram array (type ''help <a href="matlab:help PlotCCG">PlotCCG</a>'' for details).');
end
if length(t) ~= size(ccg,1),
	error('Inconsistent time bins and cross-correlogram data (type ''help <a href="matlab:help PlotCCG">PlotCCG</a>'' for details).');
end

n = size(ccg,2);
for i = 1:n,
	for j = i:n,
	%FigureIndex = g1 + nGroups*(nGroups-g2);
		subplot(n,n,i+n*(n-j));
		b = bar(t,ccg(:,i,j));
		xlim([min(t) max(t)]);
		if i == j,
			set(b,'EdgeColor','none','FaceColor','b');
		else
			set(b,'EdgeColor','k','FaceColor','k');
		end
		if i == 1,
			ylabel(int2str(j));
		end
	end
	if i == j,
		xlabel('Time (s)');
	else
		set(gca,'xtick',[]);
	end
	if j == n,
		title(int2str(i));
	end
end
