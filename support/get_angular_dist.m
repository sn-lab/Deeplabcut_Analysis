function angular_dist = get_angular_dist(angle1,angle2,units)
%FUNCTION angular_dist = get_angular_dist(angle1,angle2,units)
%
% finds the angular distance between two angles (or arrays of angles)
%
%INPUTS
%angle1, angle2: angles, or arrays of angles
%units: 'radians' or 'degrees'

assert(all(size(angle1)==size(angle2)),'angle arrays must be the same size')

switch units
    case 'radians'
        fullcircle = 2*pi;
        halfcircle = pi;
    case 'degrees'
        fullcircle = 360;
        halfcircle = 180;
    otherwise
        error('unit must be "radians" or "degrees"')
end


%convert all angles to 0-2pi or 0-360
angle1 = mod(angle1,fullcircle);
angle2 = mod(angle2,fullcircle);

%calculate distance between the 2 angles (or arrays of angles)
angular_dist = angle1-angle2;

%if the distance between 2 angles is >halfcircle, then it's closer in the other direction
angular_dist(angular_dist>halfcircle) = fullcircle - angular_dist(angular_dist>halfcircle);
