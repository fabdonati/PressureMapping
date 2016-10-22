function [len, dl] = linelen(line)
nump = size(line,1);
len = 0;
for i = 1 : nump-1
    dl(i) = norm( line(i+1,:) - line(i,:) );
    len = len + dl(i);
end