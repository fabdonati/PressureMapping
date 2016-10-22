function [ W, P ] = volrotateP_h( V, P, angle )

h = P(3);
Q = [ P(1) P(2) ];
clear P

% Rotation around height-axis
for i = 1:size(V,3)
    I = squeeze( V(:,:,i) );
    [ W(:,:,i), P ] = imrotateP( I, Q, angle );
end

P = [ P(1) P(2) h ];