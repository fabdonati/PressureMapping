function [ W, P ] = volrotateP_r(V,P,angle)

r = P(1);
Q = [ P(2) P(3) ];
clear P

% Rotation around row-axis
for i = 1:size(V,1)
    I = squeeze(V(i,:,:));
    [ W(i,:,:), P] = imrotateP(I,Q,angle);
end

P = [ r P(1) P(2) ];