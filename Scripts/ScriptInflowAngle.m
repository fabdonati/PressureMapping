
RootDir = 'F:\VictoriaPFO\';
cases = [8 10 11 13 14 16];
for i = 1:numel(cases)
    iCase = cases(i);
    angle(i) = GetInflow2SeptumAngle(RootDir,iCase);
end

fprintf('Results:\n')
for i = 1:numel(cases)
    iCase = cases(i);
    fprintf('Case %i: %1.2f\n',iCase,angle(i));
end


