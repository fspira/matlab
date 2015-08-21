%%% Script for synchronizing change in LC and acquisition by Zen.
%%% Author: Shalin Mehta, March 7, 2014.
%%%         Marine Biological Laboratory.

%% 0. Establish communication with zen and OpenPolScope 
% Before running this: Start Zen, start Zen Controller and run it.


zenPort=22500;
tcpipZen=tcpip('localhost',zenPort,'NetworkRole','client');
fopen(tcpipZen);

comPort=11;
%configFile='D:\Shared Data\Shalin\OpenPolScopeData\settings\20140121_40xOil1.4_488nm_ULMController.polset';
configFile='D:\Shared Data\Shalin\OpenPolScopeData\settings\Shalin20150223_633nm40xOil1.4NA_ULMController.polset';
%configFile='D:\Shared Data\Shalin\OpenPolScopeData\settings\20150421_100x1.4NA_488nm_ULMController.polset';

[ variLC,LCsettings ] = VariLCSetup( comPort, configFile);
% %% Establish communication with OpenPolScope
% status=setLC(tcpipOPS,'I45');
% %% Text communication with OpenPolScope
% commandfile='D:\Shared Data\Shalin\settings\LCCommands.txt';
% fileOPS=fopen(commandfile,'w+');
% fprintf(fileOPS,'%s','1');
% fclose(fileOPS);



%% 1 Setup acquisition parameters and run. 


%%%% User specified parameters. %%%%%%%%%%
PollDir='D:\Shared Data\Felix\20150429_PhalloidinAF488_mitotic\cell12T'; % Directory in which overview scan and the data is stored.
fileBase = 'cell1T';

TPoints=10000;


mkdir(PollDir)
  
    %T = timer('TimerFcn',@zenPoller,'Period',0.1); %execute the mycallback function every 0.1 seconds.
  
    %start(T);
  
    TimerH = 0.1; %for 0.7s
    %TimerH = 0.3; % for 1s
    initializeFolder = PollDir;
    counter = 0;
    
    field1 = 'timeStamp';
    field2 = 'fileName';
    field3 = 'LCOrientation';
    field4 = 'consitency';
    
    clear s
    
    
    for lauf = 1:TPoints
        %%% Save time stamp, file name, LC stage, consistency and a second
        %%% file with overall consistency 4 files per loop.
        % Initialize the current state of the folder, changes in folder
        % structure will be evaluated against the initialized state.
       
        %%%%
        %%%% I0
        %%%%
        
        lauf
        
        initializeFolder =dir(PollDir);
        setLC(variLC,'I0');
        disp('VariLC set I0')
        counter = counter +1; % used for debugging
        disp(counter)
        lastFileInFolder = initializeFolder(end).name;
        
        % Test whether the folder structure has changed while the LC was
        % changed
        [consistencyFlag] = consistencyCheck(PollDir,initializeFolder);
        if consistencyFlag == 1
            disp('Consistency check TRUE') % True = no change in folder structure
        end
        
        % flagOut is the change flag, if flag is 0 the folder structure is
        % unchanged, if set to 1 a new file was added to the folder.
        flagOut = 0;
        while flagOut == 0
            
            [flagOut, timeStamp,newFilename,initializeFolderOut]= zenPoller_oldNoNameCheck(PollDir,initializeFolder);
           flagOut
            pause(TimerH)
        end
        s(counter) = struct(field1, timeStamp, field2, newFilename, field3, 'I0',...
            field4, consistencyFlag);
        % Used for debugging
        % disp(newFilename)
        % disp(initializeFolderOut)
        
        %%%
        %%% I45
        %%%
        
        
        initializeFolder =dir(PollDir);
        setLC(variLC,'I45');
        disp('VariLC set I45')
        counter = counter +1;
        disp(counter)
        [consistencyFlag] = consistencyCheck(PollDir,initializeFolder);
        if consistencyFlag == 1
            disp('Consistency check TRUE')
        end
        
        lastFileInFolder = initializeFolder(end).name;
        
        flagOut = 0;
        while flagOut == 0
               [flagOut, timeStamp,newFilename,initializeFolderOut]= zenPoller_oldNoNameCheck(PollDir,initializeFolder);
           flagOut
            pause(TimerH)
        end
        s(counter) = struct(field1, timeStamp, field2, newFilename, field3, 'I45',...
            field4, consistencyFlag);
        % disp(newFilename)
        % disp(initializeFolder)
        
        
        %%%
        %%% I90
        %%%
        
        initializeFolder = dir(PollDir);
       setLC(variLC,'I90');
        disp('VariLC set I90')
        counter = counter +1;
        disp(counter)
        [consistencyFlag] = consistencyCheck(PollDir,initializeFolder);
        if consistencyFlag == 1
            disp('Consistency check TRUE')
        end
        
        lastFileInFolder = initializeFolder(end).name;
        
        flagOut = 0;
        while flagOut == 0
               [flagOut, timeStamp,newFilename,initializeFolderOut]= zenPoller_oldNoNameCheck(PollDir,initializeFolder);
            flagOut
            pause(TimerH)
        end
        s(counter) = struct(field1, timeStamp, field2, newFilename, field3, 'I90',...
            field4, consistencyFlag);
        %disp(newFilename)
        %disp(initializeFolderOut)
        
        
        %%%
        %%% I135
        %%%
        
        initializeFolder = dir(PollDir);
        setLC(variLC,'I135');
        disp('VariLC set I135')
        counter = counter +1;
        disp(counter)
        [consistencyFlag] = consistencyCheck(PollDir,initializeFolder);
        if consistencyFlag == 1
            disp('Consistency check TRUE')
        end
        
        lastFileInFolder = initializeFolder(end).name;
        
        flagOut = 0;
        while flagOut == 0
               [flagOut, timeStamp,newFilename,initializeFolderOut]= zenPoller_oldNoNameCheck(PollDir,initializeFolder);
            flagOut
            pause(TimerH)
        end
        s(counter) = struct(field1, timeStamp, field2, newFilename, field3, 'I135',...
            field4, consistencyFlag);
        %disp(newFilename)
        %disp(initializeFolderOut)
        
        initializeFolderOld = dir(PollDir);
        
        pause(28)
    end
%% Check consistency of the acquired imaged=s
    for checkrun = 1:length(s)
    
        consistencyTmp(checkrun) = s(checkrun).consitency;
    end
       
    if sum(consistencyTmp) < length(s)
        conistencyCheck = 'false'
    else
        conistencyCheck = 'true'
    end
    
    FramesParameter = length(s)/4;
    
    disp(['Consistency check = ', conistencyCheck]);
    disp([num2str(FramesParameter) ' frames were correctly acquired']);
    
    
%% 2 Process the data just acquired.
ChannelOfInterest=1;
confocalPaths={'D:\Shared Data\Felix\20150429_PhalloidinAF488_mitotic\cell12T'...
    };
%test = s

mkdir([char(confocalPaths) '\' 'analysis'])
% save( [char(confocalPaths) '\' 'analysis' '\' 'ConsistencyOutput.mat'],'s')
    %'D:\Shared Data\Felix\SiRActin\cell1T\',...
    %'D:\Shared Data\Felix\SiRActin\cell2T\',...
    %'D:\Shared Data\Felix\SiRActin\cell3T\',...
    % };
 % This is the time-stamp at the start of the experiment that Zen uses for
 % all exported files.
prefix={fileBase,...
   };
   % 'speedTest_2015_04_25__20_39_29__',...
   % 'speedTest_2015_04_25__20_42_40__',...
   % 'speedTest_2015_04_25__20_43_36__'};
Frames={1:17};
    %'D:\Shared Data\Felix\20150421_PhalloidinAF488_mitotic\cell2_ZStack\P001'...
    %'D:\Shared Data\Felix\20150421_PhalloidinAF488_mitotic\cell3_ZStack\P001'...
    %'D:\Shared Data\Felix\20150421_PhalloidinAF488_mitotic\cell4_ZStack\P001'...
    %'D:\Shared Data\Felix\20150422_PhalloidinAF488_mitotic\cell1_zStack\P001'...
    %'D:\Shared Data\Felix\20150422_PhalloidinAF488_mitotic\cell2_zStack\P001'...
    %'D:\Shared Data\Felix\20150422_PhalloidinAF488_mitotic\cell3_zStack\P001'...
 

 isopath='D:\Shared Data\Felix\20150428_PhalloidinAF488_mitotic\calib';
 normFactors=computeConfocalNormalizations(isopath,633,1.4,110);
%                                                   %wavelength, NA, pixsize
% wavelength and pixel size in nm.

normFactors=[1 1 1 1];
for idx=1:numel(confocalPaths)
    processConfocalPolData(confocalPaths{idx},'frames',Frames{idx},'acqMethod','polling','normFactors',normFactors,'Channel',ChannelOfInterest,'displayStatus',true,'prefix',prefix{idx});  
    exportConfocalPolData(confocalPaths{idx},'Zslices',1,'Frames',Frames{idx},'anisoCeiling',0.4,'avgCeiling',2000);
end