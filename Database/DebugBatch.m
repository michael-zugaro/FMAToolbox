function DebugBatch(mfile,bfile,item)

%DebugBatch - Assign variables to help debug a batch job.
%
%  This function can be used to help debug a batch. It parses the batch file
%  and stores the results in a structure (similar to what StartBatch does
%  does it internally).
%
%  USAGE
%
%    batch = DebugBatch(mfile,bfile,item)
%
%    mfile          batch function (M-file name or function handle)
%    bfile          batch file listing the parameters for each iteration
%    item           item number in batch function
%
%  OUTPUT
%
%    Instantiates the variables of the batch function using the values in the
%    batch file.
%
%  SEE
%
%    See also StartBatch.
%

% Copyright (C) 2007-2013 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


% Check number of parameters
if nargin < 3,
	error(['Incorrect number of parameters (type ''help <a href="matlab:help DebugBatch">DebugBatch</a>'' for details).']);
end

% Batch function name
if isa(mfile,'function_handle'),
	mfileName = func2str(mfile);
else
	mfileName = mfile;
end

% Check batch file and function are valid
if ~isstring(bfile) || ~exist(bfile,'file'),
	error(['Batch file not found (type ''help <a href="matlab:help DebugBatch">DebugBatch</a>'' for details).']);
end
if isempty(which(mfileName)),
	error(['Batch function not found (type ''help <a href="matlab:help DebugBatch">DebugBatch</a>'' for details).']);
end

% Open batch file
f = fopen(bfile,'r');
if f == -1, error(['Could not open file ''' bfile '''.']); end

% Parse batch file
b = ParseBatch(bfile);

% Open batch function
if isa(mfile,'function_handle'),
	mfile = func2str(mfile);
end
f = fopen(which(mfile));

% Find first line containing the (uncommented) 'function' keyword
found = [];
while isempty(found),
	line = fgets(f);
	if line == -1, error('Could not find function definition in batch function.'); end
	line = regexprep(line,'%.*function.*','');
	found = regexp(line,'.*function[^(]*\(([^)]*)\).*');
end
fclose(f);

% Extract parameter names, and assign them in 'base' workspace
parameters = regexprep(line,'.*function[^(]*\(([^)]*)\).*','$1');
parameters = regexp(parameters,'[^,]*','match');
for i = 1:length(parameters),
	assignin('base',parameters{i},b.field{item,i});
end
