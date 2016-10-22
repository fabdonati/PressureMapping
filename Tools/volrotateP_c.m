function [ W, P ] = volrotateP_c( V, P, angle )

c = P(2);
Q = [ P(1) P(3) ];
clear P

% Rotation around column-axis
for i = 1:size(V,2)
    I = squeeze( V(:,i,:) );
    [ W(:,i,:), P ] = imrotateP( I, Q, angle );
end

P = [ P(1) c P(2) ];