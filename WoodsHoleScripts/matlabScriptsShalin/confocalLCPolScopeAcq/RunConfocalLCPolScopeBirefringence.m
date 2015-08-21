%%% Script for synchronizing change in LC and acquisition by Zen.
%%% Author: Shalin Mehta, October 21, 2014.
%%%         Marine Biological Laboratory.

%% 0. Establish communication with zen and OpenPolScope 
% Before running this, 
% 1. Start Zen, start Zen Controller and run it.
% 2. Start OpenPolScope (Micro-Manager), select LC settings file, enable
% TCP/IP listner.
clear tcpipZen tcpipOPS;
zenPort=22500;
tcpipZen=tcpip('localhost',zenPort,'NetworkRole','client');
fopen(tcpipZen);

opsPort=56789;
tcpipOPS=tcpip('localhost',opsPort,'NetworkRole','client');
fopen(tcpipOPS);



%% 1 Setup acquisition parameters. Run once for each position.


%%%% User specified parameters. %%%%%%%%%%
datadir='D:\Shared Data\Shalin\zebrafishConfocalBirefringence\25xOilSomite\'; % Directory in which overview scan and the data is stored.
acquireZStack=false; % Set false if want to image current slice.
acquireTiles=false; % Set false if want to image just the current view.
DelayBetweenZ=0;% Delay in seconds.
DelayBetweenXY=5;% Delay in seconds.

%%%%%%%%%%


% Read current X,Y,Z from Zen.
[xCurrent,yCurrent,zCurrent]=readZenXYZ(tcpipZen);

% (If doing tile scan)Read tile positions from overview scan.
if(acquireTiles)
    [tileScanFile, pathname]=uigetfile([datadir '*.lsm'],'Select LSM file with overview tile scan.');
    r=bfGetReader([pathname '/' tileScanFile]); 
    nTiles=r.getSeriesCount();
    xPositions=NaN(1,nTiles);
    yPositions=NaN(1,nTiles);
    zStart=NaN(1,nTiles);
    
    for idS=1:nTiles % Make sure that tile-data has some overlap. Only then metadata records the positions of individual tiles.
        xposkey=['X position for position #' int2str(idS)];
        yposkey=['Y position for position #' int2str(idS)];
        zposkey=['Z position for position #' int2str(idS)];
        
        xPositions(idS)=r.getMetadataValue(xposkey);
        yPositions(idS)=-r.getMetadataValue(yposkey); % For some reason, the yPosition read from metadata needs to be flipped to match the value in the interface. Double-check.
        zStart(idS)=zCurrent;
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
    zPositions=zCurrent;
end

xPositions
yPositions
zPositions

configNames={'actinBiref'};
PolConfig=[true];
 setupstatus=setZenForPolStack(tcpipZen);
if(exist(datadir,'dir'))
    errordlg('Output directory exists, please change the name');
else 
    mkdir(datadir);
end

%% 2. Acquire.



% If the setupstatus is true/1, proceed with acquisition.
if(setupstatus)
   dlmwrite([datadir 'zpositions.txt'],[1:length(zPositions); zPositions],'delimiter','\t');
   dlmwrite([datadir 'xPositions.txt'],[1:nTiles; xPositions],'delimiter','\t');
   dlmwrite([datadir 'yPositions.txt'],[1:nTiles; yPositions],'delimiter','\t');    
   acquireZenPolStack(tcpipZen,tcpipOPS,xPositions,yPositions,zPositions,datadir,DelayBetweenXY,DelayBetweenZ,'Configs',configNames,'PolConfig',PolConfig);
end


%%% Bring stage back to original position.
%zstatus=setZenXYZ( tcpipZen,xCurrent,yCurrent,zCurrent);  

%% 3. Run this if you cancel the acquisition.
waitforZen(tcpipZen,2);
waitforZen(tcpipOPS,2);

%% 4. Run the analysis: set parameters

calibPath='D:\Shared Data\Shalin\cellulose\20141014_Bara_Calibration\P001\';
analysisPath='D:\Shared Data\Shalin\cellulose\20141014_Bara_Pollen tube 0.005%-4\P001\';
nSlices=45;
colorCeiling=0.3;
%% 5. Calibration.

ItoSMatrix=[0.5 0.5 0.5 0.5; 1 0 -1 0; 0 1 0 -1];
NA=1.4; %0
Wavelength=0.580; PixSize=0.104; % In microns

[I0Iso,I45Iso,I90Iso,I135Iso,~,blackiso]=readconfocalPolData(calibPath);


sigma=0.1*Wavelength/NA;
sigmaPix=sigma/PixSize;
FiltGauss=fspecial('gaussian',round(7*sigmaPix),sigmaPix);

I0IsoFilt=imfilter(I0Iso,FiltGauss,'replicate','same');
I45IsoFilt=imfilter(I45Iso,FiltGauss,'replicate','same');
I90IsoFilt=imfilter(I90Iso,FiltGauss,'replicate','same');
I135IsoFilt=imfilter(I135Iso,FiltGauss,'replicate','same');

eq45=I0IsoFilt./I45IsoFilt;
eq90=I0IsoFilt./I90IsoFilt;
eq135=I0IsoFilt./I135IsoFilt;

normFactors={ones(size(eq45)), eq45, eq90, eq135};


[OrientationBeforeNormalization, AnisotropyBeforeNormalization, IntensityBeforeNormalization]=...
    ComputeFluorAnisotropy(I0Iso,I45Iso,I90Iso,I135Iso,...
    'anisotropy','BlackLevel',0,'normFactors',[1 1 1 1]);

[OrientationAfterNormalization, AnisotropyAfterNormalization, IntensityAfterNormalization]=...
    ComputeFluorAnisotropy(I0Iso,I45Iso,I90Iso,I135Iso,...
    'anisotropy','BlackLevel',0,'normFactors',normFactors);

hwideIso=togglefig('isotropic slide',1); colormap gray;
set(hwideIso,'Position',[100 100 2500 1000],'defaultaxesfontsize',15);

ha=imagecat(OrientationBeforeNormalization,OrientationAfterNormalization, OrientationBeforeNormalization,OrientationAfterNormalization, AnisotropyBeforeNormalization, AnisotropyAfterNormalization, IntensityBeforeNormalization,IntensityAfterNormalization,...
    'equal','colorbar','off');

[countsIso,levelsIso]=hist((180/pi)*OrientationBeforeNormalization(:),0:0.5:180);
axes(ha(3)); stem(levelsIso,countsIso);  title('Histogram of orientaiton before normalization');
axis tight; xlim([0 180]); xlabel('Orientation');
 
[countsIsoEq,levelsIsoEq]=hist((180/pi)*OrientationAfterNormalization(:),0:0.5:180);
axes(ha(4)); stem(levelsIsoEq,countsIsoEq);  title('Histogram of orientaiton after normalization');
axis tight; xlim([0 180]); xlabel('Orientation');


%% 6. Analyze data.
[I0, I45, I90, I135, I0bleach, black]=readconfocalPolData(analysisPath,'Zslices',1:nSlices,'registerPol',false,'bin',1,'channel',1);
        
        
[orient, aniso, avg]=ComputeFluorAnisotropy(I0,I45,I90,I135,...
'anisotropy','BlackLevel',0,'normFactors',normFactors);

aniso(aniso<0)=0;
aniso(aniso>1)=1;

polstack=cat(3,aniso,orient,avg,I0,I45,I90,I135,I0bleach);
writePolStack(polstack,[analysisPath '_Z.tif'],'PolStack','anisoCeiling',1,'bitDepth',8);
writePolStack(polstack,[analysisPath '_Z.tif'],'OrientationMap','colorCeiling',colorCeiling,'bitDepth',8,'invertIntensity',false,'legend',true);