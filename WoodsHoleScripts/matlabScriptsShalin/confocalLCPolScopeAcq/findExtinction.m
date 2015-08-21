%% Script to quickly find extinction setting.
datadir='D:\Shared Data\Shalin\zencontrollerTest\';

%% Setup communication.
zenPort=22500;
tcpipZen=tcpip('localhost',zenPort,'NetworkRole','client');
fopen(tcpipZen);

opsPort=56789;
tcpipOPS=tcpip('localhost',opsPort,'NetworkRole','client');
fopen(tcpipOPS);


%% Run first iteration: step size 0.1

AStart=0.15;  BStart=0.15;
AEnd=1.1; BEnd=1.1;
LCA=AStart:0.1:AEnd;
LCB=BStart:0.1:BEnd;
I=zeros(length(LCA),length(LCB));
setupstatus=setZenForPolStack(tcpipZen);

for idA=1:numel(LCA)
    for idB=1:numel(LCB)
       setLC(tcpipOPS,['L' num2str(LCA) num2str(LCB)]);
       acquireZenImg(tcpipZen,datadir,'Ext',0);
       ImgCurr=bfopen([datadir 'Ext' '_Z' num2str(0) '.lsm']);
       ImgCurr=ImgCurr{1}{1};
        I(idA,idB)= idxA+idxB;%mean(ImgCurr(:)); 
    end
end
%% Run second iteration: step size 0.01

%% Run third iteration: step size 0.001

