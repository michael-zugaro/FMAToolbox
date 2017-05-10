function storage = DBExternalStoragePath

%DBExternalStoragePath - Get path for database external storage space.
%
%  This is an 'internal' function used by FMAToolbox. You should not need
%  to use it, unless you are developping new functions for this toolbox.


% Copyright (C) 2016 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

global SETTINGS;
if ~isfield(SETTINGS,'dbExternalStoragePath'),
    error('External storage path for databases not set (type ''help <a href="matlab:help Database">Database</a>'' for details).');
end
storage = SETTINGS.dbExternalStoragePath;
if ~isdir(storage),
    error(['External storage path ''' storage ''' does not exist (type ''help <a href="matlab:help Database">Database</a>'' for details).']);
end
name = '.test';
if ~mkdir(storage,name),
    error(['External storage path ''' storage ''' is not writable (type ''help <a href="matlab:help Database">Database</a>'' for details).']);
end
rmdir([storage '/' name]);
