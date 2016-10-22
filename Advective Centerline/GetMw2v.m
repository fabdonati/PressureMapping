function Mw2v = GetMw2v(Mv2w)
% Function to compute the inverse transformation, world to voxel
% coordinates (Mw2v), from the direct transformation voxel to world
% coordinates (Mv2w)
%
% By Pablo Lamata. Oxford, 22/01/2013

% Theory in: http://en.wikipedia.org/wiki/Affine_transformation#Augmented_matrix


R = Mv2w(1:3,1:3);

IR = inv(R);

b = -IR * Mv2w(1:3,4);

Mw2v = Mv2w;
Mw2v(1:3,1:3) = IR;
Mw2v(1:3,4)   = b;