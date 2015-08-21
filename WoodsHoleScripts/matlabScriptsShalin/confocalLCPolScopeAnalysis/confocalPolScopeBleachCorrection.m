%%% Test bleach correction
confocalPath='/media/sanguine/backup/shalin/confocalPolScope/siRactin/20150501_SirActin_5Frame/SeqNames3/';

r0=bfGetReader([confocalPath 'I0_t1.lsm']);
r45=bfGetReader([confocalPath 'I45_t1.lsm']);
r90=bfGetReader([confocalPath 'I90_t1.lsm']);
r135=bfGetReader([confocalPath 'I135_t1.lsm']);
r0b=bfGetReader([confocalPath 'I0bleach_t1.lsm']);

I0=bfGetPlane(r0,2);
I45=bfGetPlane(r45,2);
I135=bfGetPlane(r135,2);
I90=bfGetPlane(r90,2);
I0bleach=bfGetPlane(r0b,2);

[orient, aniso, avg]=ComputeFluorAnisotropy(I0,I45,I90,I135,'anisotropy');
[orientB, anisoB, avgB]=ComputeFluorAnisotropy(I0,I45,I90,I135,'anisotropy','I0bleach',I0bleach);

fig=togglefig('BleachTest');
colormap gray;
imagecat(orient, aniso, avg, orientB, anisoB, avgB,'link','equal','colorbar',fig); 

%% Bleach artificially to test the effect and correction.
bleachRatio=1.05; %I0/I0bleach=1.05.
BleachConstant=(1/4)*log(bleachRatio);
I135b=I135*exp(-BleachConstant);
I90b=I90*exp(-2*BleachConstant);
I45b=I45*exp(-3*BleachConstant);
I0b=I0*exp(-4*BleachConstant);

[orient, aniso, avg]=ComputeFluorAnisotropy(I0,I45,I90,I135,'anisotropy');
[orientUnCorr, anisoUnCorr, avgUnCorr]=ComputeFluorAnisotropy(I0,I45b,I90b,I135b,'anisotropy');
[orientCorr, anisoCorr, avgCorr]=ComputeFluorAnisotropy(I0,I45b,I90b,I135b,'anisotropy','I0bleach',I0b);

fig=togglefig('BleachTest');
colormap gray;
imagecat(orientUnCorr-orient, anisoUnCorr-aniso, avgUnCorr-avg, orientCorr-orient, anisoCorr-aniso, avgCorr-avg, 'link','equal','colorbar',fig); 

%% Call for computing corrected anisotropy, orientation, average from average of a ROI.

% Compute calibration parameters before starting.
isopath='D:\Shared Data\Felix\20150428_PhalloidinAF488_mitotic\calib';
normFactors=computeConfocalNormalizations(isopath,633,1.4,110);
 
[orientCorr, anisoCorr, avgCorr]=ComputeFluorAnisotropy(I0,I45b,I90b,I135b,'anisotropy','I0bleach',I0b,'bleachROI',1,'normFactors',normFactors);