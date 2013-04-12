% Specialized plotting functions for FMAToolbox.
%
% Interactive commands
%
%   Browse                   - Interactively browse the data in one or all subplots in a figure.
%   UIAddLine                - Interactively add horizontal/vertical line to existing plot.
%   UIInPolygon              - Find points in interactively defined polygon zone.
%   UISelect                 - Interactively select polygon zone in an existing plot.
%
% General plotting functions
%
%   PlotSamples              - Plot samples (multiple time series).
%   PlotXY                   - Plot one column vs another in a matrix.
%   MultiPlotXY              - Plot two columns of each input matrix against each other.
%   PlotMean                 - Plot mean and confidence intervals.
%   PlotIntervals            - Plot vertical bars to show interval limits.
%   PlotHVLines              - Plot vertical (resp. horizontal) lines at listed x (resp. y).
%   PlotTicks                - Plot ticks.
%   PlotSlope                - Plot line with given slope.
%   PlotRepeat               - Plot repeated (tiled) copies of the data.
%   PlotCircularDistribution - Plot phase distribution and statistics.
%
% Specialized plotting functions
%
%   PlotCCG                  - Plot auto/cross-correlograms of point processes.
%   PlotColorCurves          - Plot a series of curves as an aggregated color map.
%   PlotColorMap             - Plot a color map.
%   PlotCSD                  - Plot current source density.
%   PlotDistribution2        - For two random variables X and Y, plot distributions and Y vs X.
%   PlotPhasePrecession      - Plot phase precession plots and maps.
%   PlotRippleStats          - Plot ripple descriptive stats.
%   PlotShortTimeCCG         - Plot time-varying auto/cross-correlograms of point processes.
%   PlotSpikeWaveforms       - Plot spike waveforms.
%   PlotSync                 - Plot successive occurrences of a multidimensional variable.
%
% Simple helper functions
%
%   AdjustAxes               - Adjust axes limits for all subplots
%   AdjustColorMap           - Adjust colormap for current figure, i.e. change gamma.
%   Bright                   - Bright colormap (similar to HSV or JET, but brighter).
%   clim                     - Get or set color scaling limits for current axes.
%   Hide                     - Hide (or show) existing or future figures, e.g. to speed up batch processing.
%   hsl2hsv                  - Convert hue-saturation-luminance colors to hue-saturation-value.
%   hsv2hsl                  - Convert hue-saturation-value colors to hue-saturation-luminance.
%
% Figure elements (insets, subplots, titles...)
%
%   Insets                   - Create insets in current axes.
%   SideAxes                 - Add side axes to existing axes.
%   SquareSubplot            - Layout subplots in a square arrangement.
%   Subpanel                 - Add an (invisible) panel container to a figure.
%   SplitTitle               - Split figure title over multiple lines
%
% Tables
%
%   TableFigure              - Create table figure.
%   HTML                     - Create HTML formatted string.
