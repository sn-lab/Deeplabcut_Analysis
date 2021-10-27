function d = get_dist(x1,y1,x2,y2)
% FUNCTION get_dist(p1,p2)
% calculate the distance between pairs of 2-D cartesian points 
%
% INPUTS:
% (x1,y1,x2,y2)

dx = x1-x2;
dy = y1-y2;
d = sqrt((dx.^2)+(dy.^2));