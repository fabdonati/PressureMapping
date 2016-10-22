function [ vol, ro, co, he ] = reduce_volume( vol, border )

if nargin < 2
    border = 10;
end

ro = find(max(max(vol,[],3),[],2)==1);
ro = max([min(ro)-border 1]) : min( [ max(ro)+border size(vol,1) ] );

co = find(max(max(vol,[],3),[],1)==1);
co = max([min(co)-border 1]) : min( [ max(co)+border size(vol,2) ] );

he = find(max(max(vol,[],2),[],1)==1);
he = max([min(he)-border 1]) : min( [ max(he)+border size(vol,3) ] );

vol = vol(ro,co,he);

end