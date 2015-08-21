%% Characterize response of the light-path of the confocal using LC and the linear polarizer.

%% Step-1: Connect with LC and Zen.
zenPort=22500;
tcpipZen=tcpip('localhost',zenPort,'NetworkRole','client');
fopen(tcpipZen);

comPort=11;
palleteFile='D:\Shared Data\Shalin\OpenPolScopeData\settings\20140121_40xOil1.4_488nm_ULMController.polset';
[ variLCObj,LCsettings ] = VariLCSetup( comPort, palleteFile);
datadir=;
%% Step-2: Set linear polarizer's orientation, iterate over LC retardances, and record images. 

% Data format.
% Intensity.I(orient)=[]%Dim 1: LC-A, Dim 2: LC-B.
thisOrient=0;
fieldName=['I' num2str(thisOrient,'%03u')];
LCRange=0.2:0.1:1.2;
I=zeros(length(LCRange),length(LCRange));
for LCA=1:numel(LCRange);
    for LCB=1:numel(LCRange);
        LCcommand=['L ' num2str(LCRange(LCA)) num2str(LCRange(LCB))];
        setLC(variLCObj,LCcommand);
        acquireZenImg(tcpipZen,datadir,filename);
        % Read the file.
        I=mean(data);
    end
end
Intensity.(fieldName)=I;
%% Step-3: Compare intensity maps at different orientations of linear polarizer.