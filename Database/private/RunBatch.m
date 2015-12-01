%RunBatch - Run a batch job. This should *not* be called directly.
%
% Run a batch job. This function is called automatically by the
% batch timer upon expiration of the required delay.
%
%  USAGE
%
%    RunBatch(timer,event,b)
%
%    timer          Matlab timer object
%    event          Matlab timer event type
%    b              batch structure
%

% Copyright (C) 2007-2013 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function RunBatch(timer,event,b)

if b.hide,
	status = Hide('status');
	Hide('on');
	warning('off','FMAToolbox:Hide:FigureHidden');
end

progress = ['########## 0/' int2str(length(b.field)) ' = 0%% done [starting] ##########\n'];
fprintf(1,progress);
if b.log ~= -1,
	try
		fprintf(b.log,progress);
	catch
		fprintf(2,[' (Progress information could not be saved to log file)\n']);
	end
end

tic;
while true,
	% Get next item
	[b,item,line] = GetNextItem(b);
	if isempty(item), break; end
	% Start building command line
	i = 1;
	clear('args');
	while true,
		% Get next field
		[b,field] = GetNextField(b);
		if isempty(field), break; end
		args{i} = field;
		i = i + 1;
	end

	try
		% Call batch function
		outputs = cell(1,nargout(b.mfile));
		[outputs{:}] = feval(b.mfile,args{:});
		% Append results to 'UserData' field of timer object (this is a cell array)
		data = timer.UserData;
		if isempty(data),
			data = {outputs{:}};
		else
			data(end+1,:) = {outputs{:}};
		end
		timer.UserData = data;
	catch
		% Print error messages
		fprintf(2,['Batch: error processing item ' int2str(item) ' on line ' int2str(line) '\n']);
		e = lasterror;
		fprintf(2,[' ' e.message '\n']);
		for j = 1:length(e.stack),
			fprintf(2,[' Error in ==> ' e.stack(j).name ' at ' int2str(e.stack(j).line) '\n']);
			if strcmp(e.stack(j).name,b.mfile), break; end
		end
		% Log error messages
		if b.log ~= -1,
			try
				fprintf(b.log,['Batch: error processing item ' int2str(item) ' on line ' int2str(line) '\n']);
				fprintf(b.log,[' ' e.message '\n']);
				for j = 1:length(e.stack),
					fprintf(b.log,[' Error in ==> ' e.stack(j).name ' at ' int2str(e.stack(j).line) '\n']);
					if strcmp(e.stack(j).name,b.mfile), break; end
				end
			catch
				fprintf(2,[' (Error messages could not be saved to log file)\n']);
			end
		end
	end

	% Progress information
	t = toc;
	proportion = item / length(b.field);
	left = t*(1/proportion-1);
	if left > 24*3600,
		left = strrep(datestr(datenum(0,0,0,0,0,left),'dd-HH:MM:SS'),'-','d ');
	else
		left = datestr(datenum(0,0,0,0,0,left),'HH:MM:SS');
	end
	progress = sprintf(['########## ' int2str(item) '/' int2str(length(b.field)) ' = %.2f%%%% done [' left ' left] ##########\n'],100*proportion);
	fprintf(1,progress);
	if b.log ~= -1,
		try
			fprintf(b.log,progress);
		catch
			fprintf(2,[' (Progress information could not be saved to log file)\n']);
		end
	end
end

% Elapsed time
t = toc;
done = sprintf('Batch finished in %f s.\n',t);
fprintf(1,done);
if b.log ~= -1,
	try
		fprintf(b.log,done);
	catch
		fprintf(2,[' (Elapsed time could not be saved to log file)\n']);
	end
end


if b.hide,
	Hide(status);
	if strcmp(status,'off'), Hide('none'); end
	warning('on','FMAToolbox:Hide:FigureHidden');
end

