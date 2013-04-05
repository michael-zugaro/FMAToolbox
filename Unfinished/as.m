function [a,b,c] = as(x,y)

[a,b,c] = AdaptiveSmooth(x,y,1/39.0625,3*1e6,logical(0),logical(0));
