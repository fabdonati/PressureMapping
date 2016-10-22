function [ r, O ] = edges2crossSection( edges, v, P, N, bDebug )

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