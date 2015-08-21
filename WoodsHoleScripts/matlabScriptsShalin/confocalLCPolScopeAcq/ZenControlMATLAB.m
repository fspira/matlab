%%% Script for contorlling acquisition by Zen 2011/2012.
%%% Author: Shalin Mehta, March 7, 2014.
%%%         Marine Biological Laboratory.

%% 0. Establish communication with zen and OpenPolScope 
% Before running this, 
% 1. Start Zen, start Zen Controller and run it.


zenPort=22500;
tcpipZen=tcpip('localhost',zenPort,'NetworkRole','client');
fopen(tcpipZen);




%% 1 Setup acquisition parameters.

%%%% User specified parameters. %%%%%%%%%%
datadir='D:\Shared Data\Shalin\20140606Re-Calibration40x561nm\';
acquireZStack=true; % Set false if want to image current slice.
acquireTiles=true; % Set false if want to image just the current view.
DelayBetweenZ=0;% Delay in seconds.
DelayBetweenXY=5;% Delay in seconds.

%%%%%%%%%%


% Read current X,Y,Z from Zen.
[xCurrent,yCurrent,zCurrent]=readZenXYZ(tcpipZen);

% (If doing tile scan)Read tile positions from overview scan.
if(acquireTiles)
    tileScanFile=uigetfile([datadir '*.lsm'],'Select LSM file with overview tile scan.');
    r=bfGetReader([datadir tileScanFile]); 
    nTiles=r.getSeriesCount();
    xPositions=NaN(1,nTiles);
    yPositions=NaN(1,nTiles);
    for idS=1:nTiles
        xposkey=['X position for position #' int2str(idS)];
        yposkey=['Y position for position #' int2str(idS)];
        
        xPositions(idS)=r.getMetadataValue(xposkey);
        yPositions(idS)=r.getMetadataValue(yposkey);
    end
    
    % The metadata in the file is the offset from center. It is important
    % to add the current positions.
    xPositions=xPositions+xCurrent;
    yPositions=yPositions+yCurrent;
else
    nTiles=1;
    xPositions=xCurrent;
    yPositions=yCurrent;
end

% (If doing z stack) read z positions from zen.
if(acquireZStack)
    dz=readZendz(tcpipZen);
    [zStart,zEnd]=readZenZrange(tcpipZen);
    zPositions=zStart:dz:zEnd;
else
    zPositions=zCurrent;
end


%% 2. Acquire.

% Check if directory exists and we want to overwrite the data.
if(exist(datadir,'dir'))
    overwrite=questdlg('The directory exists. Do you want to overwrite?');
else
    overwrite='Yes';
end

switch(overwrite)
    case('Yes')
    mkdir(datadir);
    setupstatus=setZenForPolStack(tcpipZen); %Essentially turn off time and Z. Replace by suitable 'setup function'.

    case 'No'
        setupstatus=false;
        msgbox('Check the name of the directory to which the data should be written.','Acuisition aborted.');
end 

% If the setupstatus is true/1, proceed with acquisition.
if(setupstatus)
   acquireZenPolStack(tcpipZen,tcpipOPS,xPositions,yPositions,zPositions,datadir,DelayBetweenXY,DelayBetweenZ);
   dlmwrite([datadir 'zpositions.txt'],[1:length(zPositions); zPositions],'delimiter','\t');
   dlmwrite([datadir 'XTiles.txt'],[1:nTiles; xPositions],'delimiter','\t');
   dlmwrite([datadir 'YTiles.txt'],[1:nTiles; xPositions],'delimiter','\t');

end


%%% Bring stage back to original position.
zstatus=setZenXYZ( tcpipZen,xCurrent,yCurrent,zCurrent);  

%% 3 Clean up
fclose(tcpipZen);
