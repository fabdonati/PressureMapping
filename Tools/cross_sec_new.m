function [ points, r ] = cross_sec_new( points, P, v )

th = 0.15;

v = normr(v);
nump = size(points,1);

for i = 1 : nump
    d(i) = norm(points(i,:) - P);
end

for i = nump : -1 : 1
    
    if norm(points(i,:) - P) > min(d)*3
        points(i,:) = [];
    end
    
end

nump = size(points,1);

for i = nump : -1 : 1
    w = normr( points(i,:) - P );
    if abs(dot(v,w)) >= th
        points(i,:) = [];
    end
end

nump = size(points,1);

dist = zeros(nump,1);
for i = 1 : nump
    dist(i) = norm(points(i,:) - P);
end

r = mean(dist);

    