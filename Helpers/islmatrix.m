%islmatrix - Test if parameter is a logical matrix (>= 2 columns).
%
%  USAGE
%
%    test = islmatrix(x)
%
%    x              parameter to test
%
%  EXAMPLES
%
%    % Test if x is a logical matrix
%    islmatrix(x)
%
%    % Special test: test if x is a 3-line logical matrix
%    islmatrix(x,'#3')
%
%    % Special test: test if x is a 2-column logical matrix
%    islmatrix(x,'@2')
%
%  NOTE
%
%    To be considered logical, the matrix should contain only values 0 and 1, but it
%    does not need to actually be of class 'logical' (class(x) could be e.g. 'double').
%
%  SEE ALSO
%
%    See also islscalar, islvector, isdmatrix, isdvector, isdscalar, isivector, isiscalar,
%    isastring.
%

% Copyright (C) 2010-2015 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function test = islmatrix(x,varargin)

% Check number of parameters
if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help islmatrix">islmatrix</a>'' for details).');
end

% Test: logical, two dimensions, two or more columns?
test = islogical(x) & length(size(x)) == 2 & size(x,2) >= 2;

% Optional tests
for i = 1:length(varargin),
	try
		if varargin{i}(1) == '#',
			if size(x,1) ~= str2num(varargin{i}(2:end)), test = false; return; end
		elseif varargin{i}(1) == '@',
			if size(x,2) ~= str2num(varargin{i}(2:end)), test = false; return; end
		end
	catch err
		error(['Incorrect test ''' varargin{i} ''' (type ''help <a href="matlab:help isdmatrix">islmatrix</a>'' for details).']);
	end
end


function test = islogical(x)

test = builtin('islogical',x) | all(x(:)==0|x(:)==1);
