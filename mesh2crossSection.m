function [ r, O ] = mesh2crossSection( f, v, P, N, bDebug )

for i = 1 : size(f,1)
    edges((i-1)*3+1:(i-1)*3+3,:) = [ f(i,[1 2]); f(i,[1 3]); f(i,[2 3]) ];
end

edges = sort(edges,2);
edges = unique( sort(edges,2), 'rows');

j = 0;
intP = [];
for i = 1 : size(edges,1)
    A = v(edges(i,1),:);
    B = v(edges(i,2),:);
    if norm((A+B)/2 - P) <= 20
        if( dot(A-P,N)*dot(B-P,N) ) < 0;
            j = j+1;
            d = dot((P-B), N) / dot(A-B, N);
            intP(j,:) = d*(A-B) + B;
        end
    end
end
if bDebug
    hold on
    plot3( intP(:,1), intP(:,2), intP(:,3), 'b.', 'MarkerSize', 20 )
    drawnow
end
O = mean(intP);
r = mean( sqrt( sum( (intP - ones(size(intP,1),1)*O).^2,2) ) );