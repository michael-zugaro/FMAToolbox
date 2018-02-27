function m = Monochrome(color,n,varargin)

%Monochrome - Monochrome colormap (from white to a given color).
%
%  USAGE
%
%    m = Monochrome(color,n,<options>)
%
%    color          optional target color (default = black)
%    n              optional number of rows in output matrix (default = 100)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'hgamma'      gamma-like correction for hue (1 = no correction, default)
%    =========================================================================

% Copyright (C) 2009-2018 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
hgamma = 1;
type = 'linear';
n0 = 100;
color0 = [0 0 0];

% Optional parameters
if nargin == 0,
	% No parameter, use default color and n rows
	color = color0;
	n = n0;
elseif nargin == 1,
	if isdvector(color,'#3','>=0'),
		% One parameter = color, use default n rows
		n = n0;
	elseif isiscalar(color),
		% One parameter = n rows, use default color
		n = color;
		color = color0;
	else
		% One parameter = neither color nor n rows, error
		error('Incorrect parameter (type ''help <a href="matlab:help Monochrome">Monochrome</a>'' for details).');
	end
else
	if isdvector(color,'#3','>=0'),
		% First parameter is color
		if ischar(n),
			% Second is an option, use default n rows
			varargin = {n,varargin{:}};
			n = n0;
		elseif ~isiscalar(n),
			% Second is neither n rows nor an option, error
			error('Incorrect number of output rows (type ''help <a href="matlab:help Monochrome">Monochrome</a>'' for details).');
		end
	elseif isiscalar(color),
		% First parameter is n rows, use default color
		varargin = {n,varargin{:}};
		n = color;
		color = color0;
	else
		% First parameter is an option, use default color and n rows
		varargin = {color,n,varargin{:}};
		n = n0;
		color = color0;
	end
end

% Check number of parameters
if mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help Monochrome">Monochrome</a>'' for details).');
end


% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help Monochrome">Monochrome</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'hgamma',
			hgamma = varargin{i+1};
			if ~isdscalar(hgamma,'>=0'),
			error('Incorrect value for property ''hgamma'' (type ''help <a href="matlab:help Monochrome">Monochrome</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help Monochrome">Monochrome</a>'' for details).']);
	end
end

for i = 1:3,
	m(:,i) = linspace(1,color(i),n).^(1/hgamma);
end
