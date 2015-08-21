%%% Script for synchronizing change in LC and acquisition by Zen.
%%% Author: Shalin Mehta, March 7, 2014.
%%%         Marine Biological Laboratory.

%% 0. Establish communication with zen and OpenPolScope 
% Before running this, 
% 1. Start Zen, start Zen Controller and run it.
% 2. Start OpenPolScope (Micro-Manager), select LC settings file, enable
% TCP/IP listner.

zenPort=22500;
tcpipZen=tcpip('localhost',zenPort,'NetworkRole','client');
fopen(tcpipZen);

comPort=11;
%configFile='D:\Shared Data\Shalin\OpenPolScopeData\settings\20140121_40xOil1.4_488nm_ULMController.polset';
configFile='D:\Shared Data\Shalin\OpenPolScopeData\settings\Shalin20150223_633nm40xOil1.4NA_ULMController.polset';
%configFile='D:\Shared Data\Shalin\OpenPolScopeData\settings\20150421_100x1.4NA_488nm_ULMController.polset';

[ variLCObj,LCsettings ] = VariLCSetup( comPort, configFile);
% %% Establish communication with OpenPolScope
% status=setLC(tcpipOPS,'I45');
% %% Text communication with OpenPolScope
% commandfile='D:\Shared Data\Shalin\settings\LCCommands.txt';
% fileOPS=fopen(commandfile,'w+');
% fprintf(fileOPS,'%s','1');
% fclose(fileOPS);



%% 1 Setup acquisition parameters and run. 


%%%% User specified parameters. %%%%%%%%%%
datadir='D:\Shared Data\Felix\20150502_SirActin_5Frame\calib256px'; % Directory in which overview scan and the data is stored.
configLabel='100x_633nm_Gasp_12perc_800_4av_zoom3_2_256px_calib';
nFrames=10;
DelayBetweenT=0;
acqTime=4; % Acquisition time in seconds.
%%%%%%%%%%


if(exist(datadir,'dir'))
    overwrite=questdlg('The directory exists. Do you want to overwrite?');
else
    overwrite='Yes';
end

switch(overwrite)
    case('Yes')
    mkdir(datadir);
    setupstatus=setZenForPolStack(tcpipZen);
    case 'No'
        setupstatus=false;
        msgbox('Check the name of the directory to which the data should be written.','Acuisition aborted.');
end 


% If the setupstatus is true/1, proceed with acquisition.
if(setupstatus)
   acquireZenPolStackOverTime(tcpipZen,variLCObj,datadir,acqTime, 'Config',configLabel,'nFrames',nFrames,'DelayBetweenT',DelayBetweenT);
end


%% 3. Run this if you kill the acquisition.
waitforZen(tcpipZen,2);
%% 4. Clean up or not at the end.
fclose(tcpipZen);
fclose(tcpipOPS);

%% 5 Process the data just acquired.
ChannelOfInterest=1;
%Zslices={1:21,1:19,1:18,1:14,1:44,1:48,1:52};
Zslices= {1:7};
confocalPaths={'D:\Shared Data\Felix\20150424_PhalloidinAF488_mitotic\FS41\cell2_zStack\P001'...
     };
    %'D:\Shared Data\Felix\20150421_PhalloidinAF488_mitotic\cell2_ZStack\P001'...
    %'D:\Shared Data\Felix\20150421_PhalloidinAF488_mitotic\cell3_ZStack\P001'...
    %'D:\Shared Data\Felix\20150421_PhalloidinAF488_mitotic\cell4_ZStack\P001'...
    %'D:\Shared Data\Felix\20150422_PhalloidinAF488_mitotic\cell1_zStack\P001'...
    %'D:\Shared Data\Felix\20150422_PhalloidinAF488_mitotic\cell2_zStack\P001'...
    %'D:\Shared Data\Felix\20150422_PhalloidinAF488_mitotic\cell3_zStack\P001'...
 

Frames=1;
isopath='D:\Shared Data\Felix\20150424_PhalloidinAF488_mitotic\calib\';
normFactors=computeConfocalNormalizations(isopath,488,1.4,80);


for idx=1:numel(confocalPaths)
    processConfocalPolData(confocalPaths{idx},'normFactors',normFactors,'ZSlices',Zslices{idx},'Frames',Frames,'Channel',ChannelOfInterest,...
        'displayStatus',true,'suffix',['_' configNames{1}]);  
    exportConfocalPolData(confocalPaths{idx},'Zslices',Zslices{idx},'Frames',Frames,'anisoCeiling',0.5,'avgCeiling',160,'suffix',['_' configNames{1}]);
end