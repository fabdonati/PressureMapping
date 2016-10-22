function lambda = compute_lambda_2d( V, N, B, dx )

dS = dx(1) * dx(2); % Pixel area
A = dS*sum(B(:)); % Cross section area

V = permute(V,[2 3 1]);
Vx = V(:,:,1);
Vy = V(:,:,2);
Vz = V(:,:,3);

V = mean( [ Vx(B>0) Vy(B>0) Vz(B>0)] )';

lambda = N*V*A;