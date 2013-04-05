function str = mean2str(m,s,n,varargin)

%mean2str - Convert mean and SEM to string.
%
%  USAGE
%
%    s = mean2str(m,s,<options>)
%
%    m              mean
%    s              standard error of the mean (use <href="matlab:help sem">sem</a>)
%    n              optional number of observations
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'precision'   number of digits (default = 3)
%     'split'       split the output in two cells, e.g. for small figures
%                   (default = 'off')
%    =========================================================================

% Copyright (C) 2013 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
precision = 3;
split = 'off';

% Check parameters
if nargin < 2,
  error('Incorrect number of parameters (type ''help <a href="matlab:help mean2str">mean2str</a>'' for details).');
end
if ~isdscalar(m),
  error('Incorrect mean (type ''help <a href="matlab:help mean2str">mean2str</a>'' for details).');
end
if ~isdscalar(s),
  error('Incorrect SEM (type ''help <a href="matlab:help mean2str">mean2str</a>'' for details).');
end

% Optional number of observations
if nargin < 3,
	n = [];
elseif ischar(n),
	varargin = {n,varargin{:}};
	n = [];
elseif ~isiscalar(n,'>0'),
  error('Incorrect number of observations (type ''help <a href="matlab:help mean2str">mean2str</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
  if ~ischar(varargin{i}),
    error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help mean2str">mean2str</a>'' for details).']);
  end
  switch(lower(varargin{i})),
    case 'precision',
      precision = varargin{i+1};
      if ~isiscalar(precision,'>0'),
        error('Incorrect value for property ''precision'' (type ''help <a href="matlab:help mean2str">mean2str</a>'' for details).');
      end
    case 'split',
      split = lower(varargin{i+1});
      if ~isstring(split,'on','off'),
        error('Incorrect value for property ''split'' (type ''help <a href="matlab:help mean2str">mean2str</a>'' for details).');
      end
    otherwise,
      error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help mean2str">mean2str</a>'' for details).']);
  end
end

format = ['%.' int2str(precision) 'f'];
str = [sprintf(format,m) ' +- '  sprintf(format,s)];
if ~isempty(n),
	if strcmp(split,'on'),
		str = {str,['(N=' int2str(n) ')']};
	else
		str = [str ' (N=' int2str(n) ')'];
	end
end
