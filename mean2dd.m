function V = mean2dd(V)

val = V;
for j = 2:size(V,3)-1
    for i = 2:size(V,2)-1
        % Final result (the average of the 2 estimates)
        val(:,i,j)  = 0.5 * V(:,i,j) +  0.125 * (   V(:, i+1, j  )  +  V(:, i-1, j  )  ...
                                                  + V(:, i  , j+1)  +  V(:, i  , j-1)  );
    end
end
V = val;