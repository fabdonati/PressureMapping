function [CaseInfo] = ReadCaseFile(DirectoryEnsightOrCaseFile,RootName,bReport)

    if nargin<3
        bReport = 1;
    end
    bDebug = 0;
    
    if nargin==1
        filename = DirectoryEnsightOrCaseFile;
    else
        filename = [DirectoryEnsightOrCaseFile '/' RootName '.case'];
    end
    fid = fopen(filename,'r');
    CaseInfo = NaN;
    if fid==-1
        fprintf('ERROR! File could not be opened: %s\n',filename);
    else
        if(bDebug)
            fprintf('Reading CASE file %s ...\n',filename);   
        end
    
        % Get the root name from the GEO definition:
        string2find = 'model:';
        [bFound,tline,l] = findWordInFile(fid,string2find,filename);
        if(bFound)
            RootN = sscanf(tline(numel(string2find)+1:end),'%s',1);
        else
            fprintf('ERROR! Not possible to read the GEO definition in case file %s\n',filename);
            return;
        end
        % remove the .geo:
        RootName = RootN(1:end-4);
        % Get the number of steps:
        string2find = 'number of steps:';
        [bFound,tline,l] = findWordInFile(fid,string2find,filename);
        if bFound
            nSteps = sscanf(tline(numel(string2find)+1:end),'%i',1);
            if nSteps ==0
                % There must be at least one step for further calculations:
                nSteps = 1;
            end
        else
            nSteps = 1;
        end        
        % Get the first index of file:
        string2find = 'filename start number:';
        [bFound,tline,l] = findWordInFile(fid,string2find,filename);
        if(bFound)
            T0filename = sscanf(tline(numel(string2find)+1:end),'%f')';
        else
            T0filename = 0;
        end
        
        % Get the time values:
        string2find = 'time values:';
        [bFound,tline,l] = findWordInFile(fid,string2find,filename);
        if(bFound)
            TimeValues = sscanf(tline(numel(string2find)+1:end),'%f')';
            if numel(TimeValues)==0
                TimeValues = 0;
            end
        else
            TimeValues = 0;
        end
        
        while numel(TimeValues)<nSteps
            tline=deblank2(fgetl(fid));
            NewTimeValues = sscanf(tline,'%f')';
            TimeValues = [TimeValues NewTimeValues];
        end
        
        if numel(TimeValues)>1
            % Calculate deltaT:
            SC = 1e4; %Compare the deltaT up to this Scale 
            delta = (TimeValues(2:end)-TimeValues(1:end-1));
            deltaI = (round(SC*delta));
            deltaT = double(unique(deltaI))/SC;
            if numel(deltaT)>1
                if(bReport)
                    fprintf('WARNING! The increment of time during acquisitions (deltaT) is not constant!\n');
                end
                for iT=1:numel(deltaT)
                      nTimes(iT) = numel(find(deltaI/SC==deltaT(iT)));
                      if(bReport)
                            fprintf('   deltaT = %i in %i cases\n',deltaT(iT),nTimes(iT));
                      end
                end
                iTmax = find(nTimes == max(nTimes));
                temp = deltaT; clear deltaT;
                deltaT = temp(iTmax(1));
            end
        else
            deltaT = 0;
        end
        
        % Case files usually have time units in ms, and we want s:
        if deltaT>1
            scale = 0.001;
            TimeValues = TimeValues * scale;
            deltaT = deltaT * scale;
        end
        CaseInfo = struct;
        CaseInfo.RootName = RootName;
        CaseInfo.nSteps = nSteps;
        CaseInfo.T0filename = T0filename;
        CaseInfo.TimeValues=TimeValues;
        CaseInfo.deltaT = deltaT;
        fclose(fid);
    end