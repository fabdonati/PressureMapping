function [ J, P ] = imrotateP( I, P, angle )

[ rI, cI ] = size(I);

if not( ( rI-(floor(rI/2)*2) ) )
    I(end+1,:) = zeros(1,cI);
    rI = rI + 1;
end

if not( ( cI-(floor(cI/2)*2) ) )
    I(:,end+1) = zeros(rI,1);
    cI = cI + 1;
end

r = P(1);
c = P(2);

if c <= (cI-1)/2
    I = [  zeros( rI, cI-2*c + 1 )    I  ];
else
    I = [  I   zeros( rI, 2*c-1 - cI )  ];
end

cI = size(I,2);

if r <= (rI-1)/2
    I = [  zeros( rI-2*r + 1, cI );    I  ];
else
    I = [  I;   zeros( 2*r-1 - rI, cI )  ];
end

J = imrotate(I,angle);

[ rJ, cJ ] = size(J);

P = [ (rJ-1)/2+1 (cJ-1)/2+1 ];
