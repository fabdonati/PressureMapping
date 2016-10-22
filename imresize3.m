function J = imresize3( I, scale, method)

I = imresize( I, sqrt(scale), method );
I = permute( I, [3 2 1] );
I = imresize( I, sqrt(scale), method );
I = permute( I, [1 3 2] );
I = imresize( I, sqrt(scale), method );
J = permute( I, [2 3 1] );