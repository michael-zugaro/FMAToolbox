function new = DBDuplicate(old,new,varargin)

%DBDuplicate - Duplicate database.
%
%  USAGE
%
%    name = DBDuplicate(old,new)
%
%    old            database to duplicate
%    new            new database name (see NOTE below)
%
%  NOTE
%
%    Database names can include wildcards to indicate current date and time:
%
%      %y    year
%      %m    month
%      %d    day
%      %t    time
%
%  EXAMPLE
%
%    DBDuplicate('TestData_20130215','TestData_%y%m%d');
%
%  SEE
%
%    See also DBCreate, DBMerge.
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
  error('Incorrect number of parameters (type ''help <a href="matlab:help DBDuplicate">DBDuplicate</a>'' for details).');
end
if ~isastring(old),
  error('Incorrect old database name (type ''help <a href="matlab:help DBDuplicate">DBDuplicate</a>'' for details).');
end
if ~isastring(new),
  error('Incorrect new database name (type ''help <a href="matlab:help DBDuplicate">DBDuplicate</a>'' for details).');
end

% Insert date
new = InsertDate(new);

% Make sure old database exists
current = DBUse;
try
	DBUse(old);
catch
  error(['Database ''' old ''' not found (type ''help <a href="matlab:help DBDuplicate">DBDuplicate</a>'' for details).']);
end
% Create new database
DBCreate(new);
if ~isempty(current), DBUse(current); end

% Copy database
try
	h = mym(['insert into ' new '.' 'figures select * from ' old '.figures']);
	h = mym(['insert into ' new '.' 'variables select * from ' old '.variables']);
	% Copy external storage if necessary
	storage = DBExternalStoragePath;
	sourceDirectory = [storage '/' old];
	targetDirectory = [storage '/' new];
	if ~copyfile(sourceDirectory,targetDirectory),
	   error(['Could not copy external storage for database ''' new '''.']);
	end
catch
   error('FMAToolbox:DBDuplicate:copyDB',['Could not copy database ''' new '''.']);
end
