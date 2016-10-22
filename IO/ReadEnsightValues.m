function [values,ReadResult,Magnitude] = ReadEnsightValues(directory,fileName,Imageheader,options)
% Default options:
    dimensionOfValues=1;
    DiscardLines=0;
    bBinaryFormat = 0;
    if nargin>=4
        if isfield(options,'dimensionOfValues'), dimensionOfValues = options.dimensionOfValues; end
        if isfield(options,'DiscardLines'),      DiscardLines = options.DiscardLines; end
        if isfield(options,'bBinaryFormat'),     bBinaryFormat = options.bBinaryFormat; end
    end
    
% Internal parameters:    
    OrderConvention = 5;
    bDebug=1;
    ReadResult=0;
    
    ImSize = Imageheader.dim;

    filename = fullfile(directory, fileName);
    if(bBinaryFormat)
        fprintf('Reading binary file %s\n',filename);
        filecontents = loadPerNodeVariableFile(filename,prod(ImSize),dimensionOfValues);        
    else
        fid = fopen(filename,'r');
        if fid==-1
            fprintf('ERROR! File could not be opened: %s\n',filename);
            values = NaN;
            return;
        else
            fprintf('Reading file %s ...\n',filename);
        end
    
        % The first line contains the magnitude:
        Magnitude = fgetl(fid);
        % Make fid to point to the next line to the one starting in 'block':
        bFound = findWordInFile(fid,'block',filename); 
        % Get rid of an arbitrary number of lines afterwards:
        for i=1:DiscardLines
            tline = fgetl(fid);
        end
        nData = prod(ImSize);
        values = NaN * ones(ImSize(1),ImSize(2),ImSize(3),dimensionOfValues);
        %v2 = zeros(nData,dimensionOfValues);
    end
    % Reorder the data (if binary), or continue reading (if ASCII)
    
    for iDimension = 1:dimensionOfValues
        if(~bBinaryFormat)
            [data,count] = fscanf(fid,'%E',nData);
            if(bDebug), fprintf('    reading first=%f and last=%f\n',data(1),data(count)); end
            if count~=nData
                fprintf('ERROR while reading file %s! Number of data expected not found in dimension %i! (expected=%i, found=%i)\n',fileName,iDimension,nData,count);
                ReadResult=-1;
                fclose(fid);
                return;
            end
        else
            % Extract the data of this dimension:
            data = filecontents.data(:,iDimension);
        end
        %v2(:,iDimension) = data;
        for I=1:ImSize(1)
            for J=1:ImSize(2)
                for K=1:ImSize(3)
                    switch OrderConvention
                        case 1
                            iData = I + (J-1)*ImSize(1) + (K-1)*ImSize(2)*ImSize(1);
                        case 2
                            iData = I + (K-1)*ImSize(1) + (J-1)*ImSize(3)*ImSize(1);
                        case 3
                            iData = K + (J-1)*ImSize(3) + (I-1)*ImSize(3)*ImSize(2);
                        case 4
                            iData = J + (I-1)*ImSize(2) + (K-1)*ImSize(1)*ImSize(2);
                        case 5
                            iData = ImSize(1)-I+1 + (J-1)*ImSize(1) + (K-1)*ImSize(2)*ImSize(1);
                    end
                    values(I,J,K,iDimension) = data(iData);
                    if K==1
                       a=0; 
                    end
                end
            end
        end
    end
    if(~bBinaryFormat)
        fclose(fid);    
    end
end