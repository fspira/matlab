%% Generate hyperstack out of individual slices
datadir='D:\Shared Data\Shalin\20140306George_MantleMusclePhalloidinAF633\FOV2\';
zpositions=dlmread([datadir 'zpositions.txt']);
nz=zpositions(1,end);

