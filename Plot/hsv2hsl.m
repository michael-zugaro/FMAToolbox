function y = hsv2hsl(x)

%hsv2hsl - Convert hue-saturation-value colors to hue-saturation-luminance.
%
%  USAGE
%
%    y = hsv2hsl(x)
%
%    x              Nx3 RGB matrix (values in [0..1])

% Copyright (C) 2013 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Check number of parameters
if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help hsv2hsl">hsv2hsl</a>'' for details).');
end
if ~isdmatrix(x) || size(x,2) ~= 3 || any(x(:)<0) || any(x(:)>1),
  error('Incorrect RGB matrix (type ''help <a href="matlab:help hsv2hsl">hsv2hsl</a>'' for details).');
end

% Convert HSV to HSL
h = x(:,1);
s = x(:,2) .* x(:,3);
l = (2-x(:,2)).*x(:,3);
if l <= 1,
	s = s ./ l;
else
	s = s ./ (2-l);
end
l = l/2;
y = [h s l];
