%%% Script for synchronizing change in LC and acquisition by Zen.
%%% Author: Shalin Mehta, March 7, 2014.
%%%         Marine Biological Laboratory.

%% 0. Establish communication with zen and OpenPolScope 
% Before running this: Start Zen, start Zen Controller and run it.

% 
% zenPort=22500;
% tcpipZen=tcpip('localhost',zenPort,'NetworkRole','client');
% fopen(tcpipZen);

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
PollDir='D:\Shared Data\Felix\20150502_SirActin_5Frame\cell8'; % Directory in which overview scan and the data is stored.
fileBase = 'cell1T';

TPoints=10000;


mkdir(PollDir)
  
    %T = timer('TimerFcn',@zenPoller,'Period',0.1); %execute the mycallback function every 0.1 seconds.
  
    %start(T);
  
    TimerH = 0.1; %for 0.7s
    %TimerH = 0.3; % for 1s
    initializeFolder = PollDir;
    counter_I0 = 0;
    counter_I45 = 0;
    counter_I90 = 0;
    counter_I135 = 0;
    counter_I0_Bleach = 0;
    counter = 0;
    field1 = 'timeStamp';
    field2 = 'fileName';
    field3 = 'LCOrientation';
    field4 = 'consistency';
    
    clear s
    
        setLC(variLC,'I0');
        disp('VariLC set I0')
       
        for lauf = 1:TPoints
            counter = counter+1;
            initializeFolder =dir(PollDir);
             consistencyFlag = 1;
            %%%%%%%%%
            %%%% I0
            %%%%%%%%%
            
            % Check for changes in folder structure during LC-change
            
%             [consistencyFlag] = consistencyCheck(newFilename,expectedFileBase);
%             if consistencyFlag == 1
%                 disp('Consistency check TRUE')
%             else
%                 disp('Consistency check FALSE')
%             end
            
            % flagOut is the change flag, if flag is 0 the folder structure is
            % unchanged, if set to 1 a new file was added to the folder.
            
            flagOut = 0;
            expectedFileBase = 'I0_';
            counter_I0 = counter_I0 +1;
            while flagOut == 0
                pause(TimerH)
                [flagOut, timeStamp,newFilename,initializeFolderOut]= zenPoller_5Frame(PollDir,initializeFolder,expectedFileBase,counter_I0);
            end
            disp(['Frame: ', num2str(counter_I0)])
            s(counter) = struct(field1, timeStamp, field2, newFilename, field3, 'I0',...
                field4, consistencyFlag);
            
            %%%%%%%%
            %%% I45
            %%%%%%%%
            
            counter = counter+1;
            initializeFolder =dir(PollDir);
            setLC(variLC,'I45');
            disp('VariLC set I45')
            
            
%             [consistencyFlag] = consistencyCheck(newFilename,expectedFileBase);
%             if consistencyFlag == 1
%                 disp('Consistency check TRUE')
%             else
%                 disp('Consistency check FALSE')
%             end
            
            flagOut = 0;
            expectedFileBase = 'I45_';
            counter_I45 = counter_I45 +1;
            while flagOut == 0
                pause(TimerH)
                [flagOut, timeStamp,newFilename,initializeFolderOut]= zenPoller_5Frame(PollDir,initializeFolder,expectedFileBase,counter_I45);
                
            end
             disp(['I45 frame: ', num2str(counter_I45)])
            s(counter) = struct(field1, timeStamp, field2, newFilename, field3, 'I45',...
                field4, consistencyFlag);
            
            
            %%%%%%%%
            %%% I90
            %%%%%%%%
            counter = counter+1;
            initializeFolder = dir(PollDir);
            setLC(variLC,'I90');
            disp('VariLC set I90')
            
            
%             [consistencyFlag] = consistencyCheck(newFilename,expectedFileBase);
%             if consistencyFlag == 1
%                 disp('Consistency check TRUE')
%             else
%                 disp('Consistency check FALSE')
%             end
            
            flagOut = 0;
            expectedFileBase = 'I90_';
            counter_I90 = counter_I90 +1;
            while flagOut == 0
                pause(TimerH)
                [flagOut, timeStamp,newFilename,initializeFolderOut]= zenPoller_5Frame(PollDir,initializeFolder,expectedFileBase,counter_I90);
                
            end
            
             disp(['I90 frame: ', num2str(counter_I90)])
            s(counter) = struct(field1, timeStamp, field2, newFilename, field3, 'I90',...
                field4, consistencyFlag);
            
            
            
            %%%%%%%%
            %%% I135
            %%%%%%%%
            counter = counter+1;
            initializeFolder = dir(PollDir);
            setLC(variLC,'I135');
            disp('VariLC set I135')
           
%             [consistencyFlag] = consistencyCheck(newFilename,expectedFileBase);
%             
%             if consistencyFlag == 1
%                 disp('Consistency check TRUE')
%             else
%                 disp('Consistency check FALSE')
%             end
            
            flagOut = 0;
            expectedFileBase = 'I135_';
            counter_I135 = counter_I135 +1;
            
            while flagOut == 0
                pause(TimerH)
                [flagOut, timeStamp,newFilename,initializeFolderOut]= zenPoller_5Frame(PollDir,initializeFolder,expectedFileBase,counter_I135);
                
                
            end
             disp(['I135 frame: ', num2str(counter_I135)])
            s(counter) = struct(field1, timeStamp, field2, newFilename, field3, 'I135',...
                field4, consistencyFlag);
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Bleach correction LC0
            %%%%%%%%%%%%%%%%%%%%%%%%%
            counter = counter+1;
            initializeFolder = dir(PollDir);
            setLC(variLC,'I0');
            disp('VariLC set I0')
            
            
%             [consistencyFlag] = consistencyCheck(newFilename,expectedFileBase);
%             if consistencyFlag == 1
%                 disp('Consistency check TRUE')
%             else
%                 disp('Consistency check FALSE')
%             end
            
            flagOut = 0;
            expectedFileBase = 'I0Bleach_';
            counter_I0_Bleach = counter_I0_Bleach +1;
            while flagOut == 0
                pause(TimerH)
                [flagOut, timeStamp,newFilename,initializeFolderOut]= zenPoller_5Frame(PollDir,initializeFolder,expectedFileBase,counter_I0_Bleach);
            end
            disp(['I0 bleach frame: ', num2str(counter_I0_Bleach)])
            s(counter) = struct(field1, timeStamp, field2, newFilename, field3, 'I0b',...
                field4, consistencyFlag);
            initializeFolderOld = dir(PollDir);

            
            pause(25)
        end
        %% Check consistency of the acquired imaged=s
    for checkrun = 1:length(s)
   
        consistencyTmp(checkrun) = s(checkrun).consistency;
    end
       
    if sum(consistencyTmp) < length(s)
        conistencyCheck = 'FALSE'
    else
        conistencyCheck = 'TRUE'
    end
    
    FramesParameter = length(s)/5;
    
    disp(['Consistency check = ', conistencyCheck]);
    disp([num2str(FramesParameter) ' frames were correctly acquired']);
    
    
%% 2 Process the data just acquired.
ChannelOfInterest=1;
confocalPaths={'D:\Shared Data\Felix\20150502_SirActin_5Frame\cell2'...
    };
%test = s

mkdir([char(confocalPaths) '\' 'analysis'])
 save( [char(confocalPaths) '\' 'analysis' '\' 'ConsistencyOutput.mat'],'s')
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
    processConfocalPolData5(confocalPaths{idx},'frames',Frames{idx},'acqMethod','polling','normFactors',normFactors,'Channel',ChannelOfInterest,'displayStatus',true,'prefix',prefix{idx});  
    exportConfocalPolData(confocalPaths{idx},'Zslices',1,'Frames',Frames{idx},'anisoCeiling',0.4,'avgCeiling',2000);
end