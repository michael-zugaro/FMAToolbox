function target = DBMerge(source,target,varargin)

%DBMerge - Merge databases.
%
%  USAGE
%
%    target = DBMerge(source,target,<options>)
%
%    source         database where the data is read
%    target         database where the data is added
%
%  NOTE
%
%    This function cannot handle overlapping databases.
%
%  SEE
%
%    See also DBCreate, DBDuplicate.
%

% Copyright (C) 2007-2016 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Make sure MyM is installed and functional
CheckMyM;

% Check parameters
if nargin < 2,
  error('Incorrect number of parameters (type ''help <a href="matlab:help DBMerge">DBMerge</a>'' for details).');
end
if ~isastring(source),
  error('Incorrect source database name (type ''help <a href="matlab:help DBMerge">DBMerge</a>'' for details).');
end
if ~isastring(target),
  error('Incorrect target database name (type ''help <a href="matlab:help DBMerge">DBMerge</a>'' for details).');
end

% Make sure both databases exist
current = DBUse;
try
	DBUse(source);
catch
  error(['Source database ''' source ''' not found (type ''help <a href="matlab:help DBDuplicate">DBDuplicate</a>'' for details).']);
end
if ~isempty(current), DBUse(current); end
try
	DBUse(target);
catch
  error(['Target database ''' target ''' not found (type ''help <a href="matlab:help DBDuplicate">DBDuplicate</a>'' for details).']);
end
if ~isempty(current), DBUse(current); end

% Confirm
disp(['Please make sure that the databases do not overlap, otherwise ''' target ''' will be damaged.']);
s = lower(input('Type ''merge'' to confirm: ','s'));
if ~strcmp(s,'merge'),
	disp('*** Cancelled ***');
	return
end

% Copy database contents
try
	h = mym(['insert into ' target '.' 'figures select * from ' source '.figures']);
	h = mym(['insert into ' target '.' 'variables select * from ' source '.variables']);
catch
   error('FMAToolbox:DBMerge:mergeDB',['Could not merge databases ''' source ''' and ''' target ''' (non-unique eid-name pairs?)']);
end

% Copy external storage if necessary
storage = DBExternalStoragePath;
sourceDirectory = [storage '/' source];
targetDirectory = [storage '/' target];
s = [sourceDirectory '/variables'];
t = [targetDirectory '/variables'];
if isdir(s),
	if ~isdir(t), mkdir(t); end
	if ~copyfile([s '/*'],t),
		error(['Could not copy external storage for variables.']);
	end
end
s = [sourceDirectory '/figures'];
t = [targetDirectory '/figures'];
if isdir(s),
	if ~isdir(t), mkdir(t); end
	if ~copyfile([s '/*'],t),
		error(['Could not copy external storage for figures.']);
	end
end
