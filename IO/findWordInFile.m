function [bWordFound,tline,l] = findWordInFile(fid,word,filename,bDebugNotFound)
    bDebug = 0;
    bWordFound = 0;
    EoF =0;
    tline=NaN;
    l=0;
    bFirstEnd = 1;
    if nargin<4
        bDebugNotFound = 0;
    end
    while(~bWordFound&~EoF);
        tline = deblank2(fgetl(fid)); l=l+1;
        if (bDebug), fprintf('  ... New line: %s\n',tline); end
        if tline==-1
            if bFirstEnd
                fclose(fid);
                fid = fopen(filename);
                bFirstEnd = 0;
            else
                EoF=1;
                if(bDebugNotFound)
                    fprintf('WARNING! the keyworkd "%s" was not found in the file %s\n',word,filename);
                end
            end
        end
        bWordFound = findWordInBeginningOfString(tline,word);
    end
end

function b  = findWordInBeginningOfString(string,word)
    b=0;
    nCharacters = numel(word);
    if numel(string)>=nCharacters
        if strcmp(string(1:nCharacters),word)
            b=1;
        end
    end
end