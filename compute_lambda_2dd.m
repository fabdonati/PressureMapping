function lambda = compute_lambda_2dd( V, N, B, dx, stencil )

border = 1;

ro = find(max(B,[],2)==1);
ro = min(ro)-border : max(ro)+border;

co = find(max(B,[],1)==1);
co = min(co)-border : max(co)+border;

B = B(ro,co);
V = V(ro,co);

if(stencil == 'filtered')
    V = mean2dd(V);
end

Vp = sum( sum( V.*B ) );

dS = dx(1) * dx(2);
lambda = Vp*dS;