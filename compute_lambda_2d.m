function lambda = compute_lambda_2d( V, N, B, dx, stencil )

Nx = size(V,2); Ny = size(V,3);
dS = dx(1) * dx(2); % Pixel area
lambda = 0;

% Looping over the image and computing viscous components ...
for j = 2:Ny-1
    if(max(B(:,j)) < 0.5) % checking for any contributions to the sum
        continue
    end
    for i = 2:Nx-1
        lambda = lambda + dot( mean2d(V,i,j,stencil), N ) * dS * B(i,j);
    end
end