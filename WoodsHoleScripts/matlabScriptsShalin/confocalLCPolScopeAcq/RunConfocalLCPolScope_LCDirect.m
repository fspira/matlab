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
%configFile='D:\Shared Data\Shalin\OpenPolScopeData\settings\Shalin20150223_633nm40xOil1.4NA_ULMController.polset';
configFile='D:\Shared Data\Shalin\OpenPolScopeData\settings\20150421_100x1.4NA_488nm_ULMController.polset';

[ variLCObj,LCsettings ] = VariLCSetup( comPort, configFile);
% %% Establish communication with OpenPolScope
% status=setLC(tcpipOPS,'I45');
% %% Text communication with OpenPolScope
% commandfile='D:\Shared Data\Shalin\settings\LCCommands.txt';
% fileOPS=fopen(commandfile,'w+');
% fprintf(fileOPS,'%s','1');
% fclose(fileOPS);



%% 1 Setup acquisition parameters. 


%%%% User specified parameters. %%%%%%%%%%
datadir='D:\Shared Data\Felix\20150501_Alexa488Phalloidin\cell8'; % Directory in which overview scan and the data is stored.
acquireZStack=false; % Set false if want to image current slice.
acquireTiles=false;
DelayBetweenZ=0;% Delay in seconds.
DelayBetweenXY=5;% Delay in seconds.
%configNames={'488nm','561nm'};
%PolConfig=[true,false];

configNames={'100x_488nm_12percLL'};
PolConfig=true;
%%%%%%%%%%


% Read current X,Y,Z from Zen.
[xCurrent,yCurrent,zCurrent]=readZenXYZ(tcpipZen);

% (If doing tile scan)Read tile positions from overview scan.
% This feature does not work with multi-Pos dialog. Because the multiple
% position images store the absolute stage position in metadata as seen in the LCD display. Whereas,
% zencontroller uses relative positions of the stage. To make this work,
% one should set zen and LCD to zero at the same location.

if(acquireTiles)
    [tileScanFile,pathtile]=uigetfile([datadir '\..\*.lsm'],'Select LSM file with overview tile scan.');
    r=bfGetReader([pathtile tileScanFile]); 
    nTiles=r.getSeriesCount();
    xPositions=NaN(1,nTiles);
    yPositions=NaN(1,nTiles);
    zStart=NaN(1,nTiles);
    
    for idS=1:nTiles
        xposkey=['X position for position #' int2str(idS)];
        yposkey=['Y position for position #' int2str(idS)];
        zposkey=['Z position for position #' int2str(idS)];
        
        xPositions(idS)=r.getMetadataValue(xposkey)+xCurrent;
        yPositions(idS)=r.getMetadataValue(yposkey)+yCurrent; % For some reason, the yPosition read from metadata needs to be flipped to match the value in the interface. Double-check.
        zStart(idS)=r.getMetadataValue(zposkey)+zCurrent;
    end
    
    % The metadata in the file is the offset from center. It is important
    % to add the current positions.
%     xPositions=xPositions+xCurrent;
%     yPositions=yPositions+yCurrent;
else
    nTiles=1;
    xPositions=xCurrent;
    yPositions=yCurrent;
    zStart=zCurrent;
end

% (If doing z stack) read z positions from zen.
if(acquireZStack)
    dz=readZendz(tcpipZen);
    nz=readZenNz(tcpipZen);
    zPositions=bsxfun(@plus,zStart(:),(0:1:nz-1)*dz); % Arrange zStart at each position along first dimension, and individual depths along second dimension. 
else
    zPositions=repmat(zCurrent,[length(xPositions) 1]);
end


%%% Make sure these readings match with readings seen on the zen interface
%%% AND on the LCD display.
xPositions
yPositions
zPositions

%% 2. Acquire.
% In Zen, switch off Z and Tile before running this section.
% Check if directory exists and we want to overwrite the data.




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
   dlmwrite([datadir 'zpositions.txt'],zPositions,'delimiter','\t');
   dlmwrite([datadir 'xPositions.txt'],[1:nTiles; xPositions],'delimiter','\t');
   dlmwrite([datadir 'yPositions.txt'],[1:nTiles; yPositions],'delimiter','\t');    
   acquireZenPolStack(tcpipZen,variLCObj,xPositions,yPositions,zPositions,datadir,DelayBetweenXY,DelayBetweenZ,'Configs',configNames,'PolConfig',PolConfig);
end


%%% Bring stage back to original position.
zstatus=setZenXYZ( tcpipZen,xCurrent,yCurrent,zCurrent);  


%% 3. Run this if you kill the acquisition.
waitforZen(tcpipZen,2);
%% 4. Clean up or not at the end.
fclose(tcpipZen);
fclose(tcpipOPS);

%% 5 Process the data just acquired.
ChannelOfInterest=1;
%Zslices={1:21,1:19,1:18,1:14,1:44,1:48,1:52};
Zslices= {1};
confocalPaths={'D:\Shared Data\Felix\20150501_Alexa488Phalloidin\cell8\P001'...
     };
    %'D:\Shared Data\Felix\20150421_PhalloidinAF488_mitotic\cell2_ZStack\P001'...
    %'D:\Shared Data\Felix\20150421_PhalloidinAF488_mitotic\cell3_ZStack\P001'...
    %'D:\Shared Data\Felix\20150421_PhalloidinAF488_mitotic\cell4_ZStack\P001'...
    %'D:\Shared Data\Felix\20150422_PhalloidinAF488_mitotic\cell1_zStack\P001'...
    %'D:\Shared Data\Felix\20150422_PhalloidinAF488_mitotic\cell2_zStack\P001'...
    %'D:\Shared Data\Felix\20150422_PhalloidinAF488_mitotic\cell3_zStack\P001'...
 

    mkdir([char(confocalPaths{1}) '\analysis'])
    
Frames=1;
isopath='D:\Shared Data\Felix\20150501_Alexa488PhalloidinLL\calib';
normFactors=computeConfocalNormalizations(isopath,488,1.4,80);


for idx=1:numel(confocalPaths)
    processConfocalPolData(confocalPaths{idx},'normFactors',normFactors,'ZSlices',Zslices{idx},'Frames',Frames,'Channel',ChannelOfInterest,...
        'displayStatus',true,'suffix',['_' configNames{1}]);  
    exportConfocalPolData(confocalPaths{idx},'Zslices',Zslices{idx},'Frames',Frames,'anisoCeiling',0.5,'avgCeiling',160,'suffix',['_' configNames{1}]);
end