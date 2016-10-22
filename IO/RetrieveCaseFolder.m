function [RootDir,NameConventionIN,NameConventionOUT,CaseString,CaseStringPlotting,EnsightVersion] = RetrieveCaseFolder(CaseNumber,WhichCohort)
% Control version:
% - 20/12/2012: addition of default output (in order to use it in a general
% fashion).

global WhichDisk;
WhichDisk = 'E';

%--------------------------------------------------------------------------
% Default output:
RootDir = CaseNumber;
NameConventionIN  =1; % like "_**"
NameConventionOUT =1;     
if ischar(CaseNumber)
    [RootDir CaseString] = GetPath(CaseNumber);
    if isempty(CaseString)
        % in case the last character is the /:
        [foo CaseString] = GetPath(CaseNumber(1:end-1));
    end
else
    CaseString = num2str(CaseNumber);
end
CaseStringPlotting = CaseString;
EnsightVersion = 92;
%--------------------------------------------------------------------------
DataPackage = NaN;
if nargin<2
    WhichCohort = 'Alex';
end
    
if isempty(WhichCohort)
    WhichCohort = 'Alex';
end

if isnumeric(CaseNumber)
    switch WhichCohort
        case 'Alex'
            CaseString = num2str(CaseNumber,'%02i');
            switch CaseNumber
                case {1 2 79 80 82 83 84}
                    % These cases were processed with Ensight92
                    DataPackage = 1;
                    if CaseNumber==1, CaseString='BA'; end
                    if CaseNumber==2, CaseString='Di'; end  
                case {5 6 56 57 81 832}
                    % These cases were processed with Ensight91, Sebastian's
                    % computer by Tom
                    DataPackage = 2;
                    if CaseNumber==832, CaseString='83_2'; end  
                case {16 41 48 50 51 55 59 62 63 64 67 68 69 70 72 73 74 77 78 86 87}
                    DataPackage = 4;
                otherwise
                    DataPackage = 5;
            end
            switch CaseNumber
                case{5 6 56 57 81 832}
                    % Cases which required a different numbering convention:
                    EnsightVersion = 91;            
                otherwise
                    % do nothing
                    EnsightVersion = 92;
            end
        case 'Rachel'
            EnsightVersion = 92;                        
            CaseNumberORdirectory = sprintf('%i',CaseNumber);
            CaseString = CaseNumberORdirectory;
            DataPackage = 3;
        case 'BAV'
            EnsightVersion = 92;                        
            CaseNumberORdirectory = sprintf('OXBAV%03i',CaseNumber);
            CaseString = CaseNumberORdirectory;
            DataPackage = 6;
        otherwise
            fprintf('Data cohort not defined (path should be provided)\n')
    end
else
    if CaseString(1)=='R'
        % These are the cases from Rachel
        DataPackage = 3;
        CaseString = CaseString(2:end);
    else
        % removed 20/12/12
        %DataPackage = 1;
    end                
end

if(CaseNumber==0)
    DataPackage = 0;
end

    if isunix()
        Root = '/media/My Passport/data/';
    else
        Root = [WhichDisk ':\data'];
        Root2 = 'G:\PPEdata';
        Root3 = 'F:\PPE\BAV';
%         Root1 = 'H:\data';        
%         if isdir(Root1)
%             Root = Root1;
%         else
%             Root1 = 'F:\data'; 
%             if isdir(Root1)
%                 Root = Root1;
%             else
%                 Root2 = 'E:\data'; 
%                 if isdir(Root2)
%                     Root = Root2;
%                 else
%                     fprintf('No root dir of data found! Is passport connected?\n')
%                     return;
%                 end
%             end
%         end
    end
    switch DataPackage
        case 0
            RootDir = [Root '/PressureCases/First cases/MeinHerzDataProcessing/'];
            NameConventionIN  =1;
            NameConventionOUT  =2;
        case 1
            %RootDir = [Root 'data/PressureCases/AlexCases/CASEMS0'  CaseString '/'];
            RootDir = [Root '/OCMR1/CASEMS0'  CaseString '/'];
            NameConventionIN  =1;
            NameConventionOUT =2;
        case 2
            RootDir = [Root '/OCMR2/CASEMS0'  CaseString '/'];
            NameConventionIN  =1; % like "_**"
            NameConventionOUT =1;     
        case 3
            RootDir = [Root '/KCLcases/Rachel/CASE'  CaseString '/'];
            NameConventionIN  =1; % like "_**"
            NameConventionOUT =2;    
        case 4
            RootDir = [Root '/OCMRnew/CaseMS0'  CaseString '/'];
            NameConventionIN  =2; % like "****"
            NameConventionOUT =2;    
        case 5
            RootDir = [Root2 '/MS0'  CaseString '/'];
            NameConventionIN  =2; % like "****"
            NameConventionOUT =2;    
        case 6
            RootDir = fullfile(Root3,CaseString,'pressure_mapping');
            NameConventionIN  =2; % like "****"
            NameConventionOUT =2;    
    end

    CaseStringPlotting = CaseString;
    switch CaseString 
        case 'Di'
            CaseStringPlotting = 'AoD';
        case 'BA'
            CaseStringPlotting = 'BAV';
        case '79'
            CaseStringPlotting = 'MFS';
        case '82'
            CaseStringPlotting = 'HV1';
        case '83'
            CaseStringPlotting = 'HV2';
        case '06'
            CaseStringPlotting = 'HV3';
    end