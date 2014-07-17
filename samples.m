% Samples - List of timestamps and corresponding values
%
% In FMAToolbox, 'samples' refer to matrices in which the first column
% contains timestamps, and subsequent columns optionally contain values
% recorded at these timestamps.
%
% (Obviously, in functions implementing statistical tests, the word 'samples'
% is used in its usual sense).
%
%  EXAMPLES
%
%    The simplest case is the spike train of a single neuron. This can be des-
%    cribed as a list of timestamps:
%
%       t = [ 0.124
%             0.242
%             0.587
%             1.005
%             ...
%            15.239 ];
%
%    A single local field potential channel contains two columns, one for time
%    and one for voltages:
%
%       v = [ 0.0008  562
%             0.0016  451
%             0.0024  399
%             0.0032  410
%             0.0040  362
%             ...
%            15.0480 -109 ];
%
%    Multiple channels would simply require more columns:
%
%       v = [ 0.0008  562  520  680
%             0.0016  451  509  702
%             0.0024  399  384  716
%             0.0032  410  358  673
%             0.0040  362  359  644
%             ...
%            15.0480 -109 -150 -287 ];
