function [directionLabels] = velocityToDir(velocity)
% Converts the 3 column velocity matrix of X-Y-Z velocity to one of 8
% cardinal directions, where 1 = Right, 2 = Up-Right (proceed CCW)...8 =
% Down-Right. NOTE Y-VALUES ARE FLIPPED (down = +... so mult by -1)
%
% Input: velocity - 3 column vectors with X, Y, and Z velocity (Z not used)
%
% Output: direction - 1 column vector with direction 1-8, where 1 = Right,
% proceed CCW for 2-7, and 8 is Down-Right; 0 = no movement


% first, find max velocity and determine where is "moving" (> 0.05 max vel)
velMag = sqrt(sum(velocity(:,1:2).^2,2));
minVel = 0.05*max(velMag);
indicesToUse = velMag > minVel;
%unitVel = velocity(:,1:2)./velMag;

% now assign direction based on velocity
dirAngles = [0 45 90 135 180 225 270 315];
%dirAngles = [0 90 180 270];
numDirs = length(dirAngles);
dirAngles = [360+dirAngles 720+dirAngles 1080+dirAngles dirAngles -360+dirAngles -720+dirAngles -1080+dirAngles ]; % include more extended angles
directionLabels = zeros(length(velMag),1);
directions = 180/pi*(unwrap(angle(velocity(:,1)+1i*-velocity(:,2)))); % I don't like that I have to do it like this... also - sign on velocity(:,2) bc y-axis is flipped
[~,closestInds] = min(abs(dirAngles-directions),[],2);

% finally, for all above the minimum threshold velocity, relabel
directionLabels(indicesToUse) = mod(closestInds(indicesToUse)-1,numDirs)+1; % the weird modulus thing is to get the negative angles
end

