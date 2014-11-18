function ProcessBinary(inputName,outputName,nChannels,nOverlap,f,varargin)

%ProcessBinary - Process binary data file.
%
% Process binary data file, e.g. filter LFP file. This function loads the
% binary file segment by segment, calls a user-supplied function to process
% each data segment, and saves the result to a new file.
%
% If the function requires overlapping segments (e.g. to avoid edge effects),
% it will receive data in the form [o1;s;o2], where s is the segment to process,
% and o1 and o2 are portions of the previous and next segments (half the overlap
% size each). The function should return the processed segment, WITHOUT the
% overlapping portions.
%
%  USAGE
%
%    ProcessBinary(inputName,outputName,nChannels,nOverlap,f,p1,p2...)
%
%    inputName      binary input file
%    outputName     binary output file
%    nChannels      number of channels in the input file
%    nOverlap       nOverlap between successive reads (in # samples)
%    f              function handle
%    p1,p2...       optional additional parameters for f
%

% Copyright (C) 2004-2014 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if nargin < 4,
	error('Incorrect number of parameters (type ''help <a href="matlab:help ProcessBinary">ProcessBinary</a>'' for details).');
end
if rem(nOverlap,2) ~= 0,
	error('Overlap must be even (type ''help <a href="matlab:help ProcessBinary">ProcessBinary</a>'' for details).');
end
segmentLength = 2^16;
parameters = varargin;

% Open input and output files
inputFile = fopen(inputName,'r');
outputFile = fopen(outputName,'w');

% Process first segment
if nOverlap == 0,
	overlap = zeros(nChannels,0);
else
	% Flip first segment LR and prepend it to avoid edge effects
	overlap = fread(inputFile,[nChannels,nOverlap/2],'int16');
	overlap = fliplr(overlap);
	frewind(inputFile);
end
segment = fread(inputFile,[nChannels,segmentLength],'int16');
segment = [overlap,segment]';
processed = feval(f,segment,parameters{:});
fwrite(outputFile,processed,'int16');
overlap = segment(end-(nOverlap-1):end,:);

% Process subsequent segments
while ~feof(inputFile),
	segment = fread(inputFile,[nChannels,segmentLength],'int16');
	segment = [overlap;segment'];
	processed = feval(f,segment,parameters{:});
	fwrite(outputFile,processed,'int16');
	overlap = segment(end-(nOverlap-1):end,:);
end

% Process trailing unprocessed segment
if nOverlap ~= 0,
	% Flip last segment UD and append it to avoid edge effects
	segment = segment(end-(nOverlap-1):end,:);
	tail = flipud(segment(end-(nOverlap/2-1):end,:));
	segment = [segment;tail];
	processed = feval(f,segment,parameters{:});
	fwrite(outputFile,processed,'int16');
end

% Close input and output files
fclose(inputFile);
fclose(outputFile);

