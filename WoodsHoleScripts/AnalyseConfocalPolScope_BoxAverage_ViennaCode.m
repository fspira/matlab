%%%%% anisotropy Linescan
%%
clear all

distvec = linspace(0,1,256)';
red = zeros(256,3);
red(:,1) = distvec;

green = zeros(256,3);
green(:,2) = distvec;

blue = zeros(256,3);
blue(:,3) = distvec;
addpath('/Users/spira/Documents/Matlab_scripte/Image_Processing_utils')

addpath('/Users/spira/Documents/Matlab_scripte/tiffIO')
addpath('/Users/spira/Documents/Matlab_scripte/')
addpath('/Users/spira/Desktop/Desktop/LifeactCherry_GlGPIEgfp/131204')
addpath('/Users/spira/Documents/MATLAB_scripte/ImageProcessing/Utilities')

addpath('/Users/spira/Desktop/programme/calculateAnisotropieRatio')
addpath('/Users/spira/Desktop/programme/calculateanisotropyRatio')
addpath('/Users/spira/Desktop/programme/centerOfMass')
addpath('/Users/spira/Desktop/programme/curvature')
addpath('/Users/spira/Desktop/programme/determineAngle')
addpath('/Users/spira/Desktop/programme/staging')
addpath('/Users/spira/Desktop/programme/tools')
addpath('/Users/spira/Desktop/programme/WoodsHoleScripts')

javaaddpath '/Applications/Fiji.app/jars/ij-1.49h.jar'
addpath('/Applications/Fiji.app/scripts')
addpath('/Applications/Fiji.app/plugins/')
addpath('/Users/spira/Desktop/programme/WoodsHoleScripts/matlabScriptsShalin')
%addpath('/Users/spira/Desktop/2015_WoodsHole/matlabScriptsShalin')

%addpath('/Users/spira/Desktop/programme/miji/mij.jar')
%addpath('/Users/spira/Desktop/programme/miji/ij.jar')

javaaddpath '/Applications/MATLAB_R2012a_Student.app/java/mij.jar'

javaaddpath '/Applications/MATLAB_R2012a_Student.app/java/ij.jar'


 %MIJ.start('/Applications/Fiji.app')


curdir = pwd;

origDir = '/Users/spira/Desktop/5FramePoleScope/20150503_SirActin_5Frame_dish2/cell2'

fileList= dir([curdir '/' '*I0_*'])

isopath='/Users/spira/Desktop/5FramePoleScope/calib256';
%normFactors=computeConfocalNormalizations(isopath,633,1.4,110);


clear avgFile




for lauf = 1 : length(fileList)
 
   tmpFile = tiffread30(fileList(lauf).name);
   tmpFileFluor = tmpFile.data{1};
    
    avgFile(:,:,lauf) =   tmpFileFluor;
   
end

tiffwrite_mat(avgFile, [origDir '/' 'analysis' '/' 'I0Channel.tif'])

cd(origDir)

fileListOrig = dir([origDir '/'  '*I*'])



  
[m n p] =size(avgFile)

AnalysisEnd =p
firstFrame =1
lastFrame = 33;
anaOnset =10;


saveFileName = 'cell1Analysis';

%%% load images

lengthDirectory = length(fileListOrig)






fileListOrig0 = dir([origDir '/'  '*I0_*'])
fileListOrig45 = dir([origDir '/'  '*I45_*'])
fileListOrig90 =  dir([origDir '/'  '*I90_*'])
fileListOrig135 = dir([origDir '/'  '*I135_*'])
fileListOrigB = dir([origDir '/'  '*I0Bleach_*'])
fileListOrigH2B = dir([origDir '/' '*H2B_*'])



          
  
p = length(fileListOrig0)

AnalysisEnd =p
%firstFrame = 22;
%lastFrame = 42;
%anaOnset =5;


%%% load images

lengthDirectory = length(fileListOrig0)

            clear listDirectory sortOut
            clear sortedList sortArray sortIdx

            for subrun = 1:length(fileListOrig0)-1
                laufName = fileListOrig0(subrun+1).name
                listDirectory{subrun} =  laufName
                suffixStrSort = findstr( laufName,'.lsm');
                prefixStrSort = findstr( laufName,'_');


                sortArray(subrun) = str2num(laufName(prefixStrSort+1:suffixStrSort-1));
            end

            [sortArray0 listDirectory0] = doSortDirectory(fileListOrig0)
            [sortOut0 sortIdx0] = sort(sortArray0)
           
            [sortArray45 listDirectory45] = doSortDirectory(fileListOrig45)
            [sortOut45 sortIdx45] = sort(sortArray45)
           
            [sortArray90 listDirectory90] = doSortDirectory(fileListOrig90)
            [sortOut90 sortIdx90] = sort(sortArray90)
           
            [sortArray135 listDirectory135] = doSortDirectory(fileListOrig135)
            [sortOut135 sortIdx135] = sort(sortArray135)
           
            [sortArrayB listDirectoryB] = doSortDirectory(fileListOrigB)
            [sortOutB sortIdxB] = sort(sortArrayB)
           
            [sortArrayH2B listDirectoryH2B] = doSortDirectory(fileListOrigH2B)
            [sortOutH2B sortIdxH2B] = sort(sortArrayH2B)
           

            sortedList0{1} = fileListOrig0(1).name
            sortedList45{1} = fileListOrig45(1).name
            sortedList90{1} = fileListOrig90(1).name
            sortedList135{1} = fileListOrig135(1).name
            sortedListB{1} = fileListOrigB(1).name
            sortedListH2B{1} = fileListOrigH2B(1).name

            for subrun = 1:length(sortIdx0)
                sortedList0{subrun+1} = listDirectory0{sortIdx0(subrun)};
                sortedList45{subrun+1} = listDirectory45{sortIdx45(subrun)};
                sortedList90{subrun+1} = listDirectory90{sortIdx90(subrun)};
                sortedList135{subrun+1} = listDirectory135{sortIdx135(subrun)};
                sortedListB{subrun+1} = listDirectoryB{sortIdxB(subrun)};
                sortedListH2B{subrun+1} = listDirectoryH2B{sortIdxH2B(subrun)};
            end

            
            Frame = length(fileListOrig0);
            p = length(fileListOrig0);
            
            
            
  
for lauf = 1:p
    I0Name(:,:,:,lauf) = tiffread30(char(sortedList0((lauf))));
    I45Name(:,:,:,lauf) = tiffread30(char(sortedList45((lauf))));
    I90Name(:,:,:,lauf) = tiffread30(char(sortedList90((lauf))));
    I135Name(:,:,:,lauf) = tiffread30(char(sortedList135((lauf))));
    IBName(:,:,:,lauf) = tiffread30(char(sortedListB((lauf))));
    IH2BName(:,:,:,lauf) = tiffread30(char(sortedListH2B((lauf))));

end

for lauf = 1:p
 
    I0Tmp = I0Name(:,:,:,lauf);
    I0Tmp = cat(3,I0Tmp.data);
    I0File(:,:,lauf) =  I0Tmp{1};
    
    I45Tmp = I45Name(:,:,:,lauf);
    I45Tmp = cat(3,I45Tmp.data);
    I45File(:,:,lauf) = I45Tmp{1};
    
    I90Tmp = I90Name(:,:,:,lauf);
    I90Tmp = cat(3,I90Tmp.data);
    I90File(:,:,lauf) = I90Tmp{1};
    
    I135Tmp = I135Name(:,:,lauf);
    I135Tmp = cat(3,I135Tmp.data);
    I135File(:,:,lauf) = I135Tmp{1};
    
     IBTmp = IBName(:,:,lauf);
    IBTmp = cat(3,IBTmp.data);
    IBFile(:,:,lauf) = IBTmp{1};
    
    
    IH2BTmp = IH2BName(:,:,lauf);
    IH2BTmp = cat(3,IH2BTmp.data);
    IH2BFile(:,:,lauf) = IH2BTmp{1};
 

end


I0Segment = doSegmentImage(I0File);


I45Segment = doSegmentImage(I45File);

I90Segment = doSegmentImage(I90File);
I135Segment = doSegmentImage(I135File);
IBSegment = doSegmentImage(IBFile);


for lauf = 1:p

    averageImage(:,:,lauf)= (I0File(:,:,lauf) + I45File(:,:,lauf)+ I90File(:,:,lauf)...
         + I135File(:,:,lauf)+ IBFile(:,:,lauf))/5;

end







for lauf = 1:p
   
    
     greenStackRGB(:,:,:,lauf) = ind2rgb(normalizedImage(averageImage(:,:,lauf)),green);
     redStackRGB(:,:,:,lauf) = ind2rgb(normalizedImage(IH2BFile(:,:,lauf)),red);
      blueStackRGB(:,:,:,lauf) = ind2rgb(normalizedImage(IH2BFile(:,:,lauf)),blue);
     imgMerge(:,:,:,lauf) = greenStackRGB(:,:,:,lauf) + redStackRGB(:,:,:,lauf)+blueStackRGB(:,:,:,lauf);
end


tiffwrite_mat(averageImage, [origDir '/' 'analysis' '/' 'averageImage.tif'])

tiffwrite_mat(I0File, [origDir '/' 'analysis' '/' 'I0File.tif'])
tiffwrite_mat(I45File, [origDir '/' 'analysis' '/' 'I45File.tif'])
tiffwrite_mat(I90File, [origDir '/' 'analysis' '/' 'I90File.tif'])
tiffwrite_mat(I135File, [origDir '/' 'analysis' '/' 'I135File.tif'])
tiffwrite_mat(IBFile, [origDir '/' 'analysis' '/' 'IBleachFile.tif'])
tiffwrite_mat(IH2BFile, [origDir '/' 'analysis' '/' 'H2BFile.tif'])


tiffwrite_RBG(imgMerge, [origDir '/' 'analysis' '/' 'imgMerge.tif'])


%imshow(averageImage(:,:,1),[]);



voxelX = getfield(I0Name(:,:,:,lauf),'lsm','VoxelSizeX');
voxelX_mum = voxelX*1000000;

timeInterval = 30;

anaTime = anaOnset*timeInterval;
timeMax = (p*timeInterval+p) - anaTime;
anaTime = 0;
 
 timeVec = 1:timeInterval:p*timeInterval;
            
timeVec = timeVec-timeVec(anaOnset);


%voxelX_mum
%voxelX_H2B

%%%% Select first channel from file


%%%% Crop the file

    clear I0Tmp
    I0Tmp = I0File(:,:,firstFrame:lastFrame);
    clear I0File
    I0File = I0Tmp;
    clear I0Tmp
    
    clear I45Tmp
    I45Tmp = I45File(:,:,firstFrame:lastFrame);
    clear I45File
    I45File = I45Tmp;
    clear I45Tmp
    
    clear I90Tmp
    I90Tmp = I90File(:,:,firstFrame:lastFrame);
    clear I90File
    I90File = I90Tmp;
    clear I90Tmp
    
    clear I135Tmp
    I135Tmp = I135File(:,:,firstFrame:lastFrame);
    clear I135File
    I135File = I135Tmp;
    clear I135Tmp
    
 

%%

%%
%% Background correction I135File

[I0File I45File I90File I135File I135Fileb] =  bkgCorrection5FrameConfocalPolscope(I0File,  I45File,I90File, I135File, IBFile);
%%%% average files

          for lauf = 1:p
       
       

                [BW_Box BWI0 BWI45 BWI90 BWI135 BWI135B x y...
                     I0Avg I45Avg I90Avg I135Avg I0BleachAvg] = doDrawAnalysisBox5Orientations(I0File(:,:,lauf),I0File(:,:,lauf), I45File(:,:,lauf),...
                    I90File(:,:,lauf), I135File(:,:,lauf), I0File(:,:,lauf) ,voxelX_mum)



                xFurrow1{lauf} = x;
                yFurrow1{lauf} = y;

       


               %%%% Parse data

               I0_Furrow(lauf) =(I0Avg(1) + I0Avg(2))/2;
               I0_Pole(lauf) = (I0Avg(3) + I0Avg(4))/2;

               I45_Furrow(lauf) = (I45Avg(1) + I45Avg(2))/2;
               I45_Pole(lauf) = (I45Avg(3) + I45Avg(4))/2;

               I90_Furrow(lauf) = (I90Avg(1) + I90Avg(2))/2;
               I90_Pole(lauf) = (I90Avg(3) + I45Avg(4))/2;

               I135_Furrow(lauf) = (I135Avg(1) + I135Avg(2))/2;
               I135_Pole(lauf) = (I135Avg(3) + I135Avg(4))/2;

               I0Bleach_Furrow(lauf) = (I0BleachAvg(1) + I0BleachAvg(2))/2;
               I0Bleach_Pole(lauf) = (I0BleachAvg(3) + I0BleachAvg(4))/2;

           
              I0_Store{lauf} = I0Avg;
              I45_Store{lauf} = I45Avg;
              I90_Store{lauf} = I90Avg;
              I135_Store{lauf} = I135Avg;
              I0Bleach_Store{lauf} = I0BleachAvg;
              BW_BoxStore{lauf} =    BW_Box;

              BWI0_Store{lauf} = BWI0;
              BWI45_Store{lauf} = BWI45;
              BWI90_Store{lauf} = BWI90;
              BWI135_Store{lauf} = BWI135;
              BWI135B_Store{lauf} = BWI135B;
              x_Store{lauf} = x;
              y_Store{lauf} = y;
       

          end

          for lauf = 1:p
              x_StoreTmp = cell2mat(x_Store{lauf});
              y_StoreTmp = cell2mat(y_Store{lauf});
              
               cleavageFurrowDiameter(lauf) = (pdist2([x_StoreTmp(1),y_StoreTmp(1)],[x_StoreTmp(2),y_StoreTmp(2)])) * voxelX_mum
          end
          
          plot(cleavageFurrowDiameter)

      
          

plot(I0_Furrow)
            



arg.Zslices=1; 
arg.Frames=1;
arg.Channel=1; % Select the channel to process.
%arg.Position=1;
arg.registerPol=false; % Do pol-channels need registration?
arg.bin=1;
arg.normFactors=[1 1 1 1];
arg.suffix=''; %Suffix after the slice name.
arg.prefix=''; %Prefix before the time-stamp for data acquired by polling when files are written by Zen.
arg.BlackLevel=0;
arg.bitDepth=NaN; %Bit-depth has no effect on how data is read. It only affects the format in which computed results are written.
arg.displayStatus=true;
arg.acqMethod='zen-controller';
%arg=parsepropval(arg,varargin{:});
clear anisoStore
clear avgStore
for lauf = 1:AnalysisEnd
    
  %  [orient, aniso, avg]=...
  %      ComputeFluorAnisotropy(I0_Furrow(lauf),I45_Furrow(lauf),I90_Furrow(lauf),I135_Furrow(lauf),...
  %      'anisotropy','BlackLevel',arg.BlackLevel,'normFactors',arg.normFactors);
    [orientCorr, anisoCorr, avgCorr]=ComputeFluorAnisotropy(I0_Furrow(lauf),I45_Furrow(lauf),I90_Furrow(lauf),I135_Furrow(lauf),...
        'anisotropy','I0bleach',I0_Furrow(lauf),'bleachROI',1,'normFactors',normFactors);

    %orientStore{lauf} = orient;
    anisoStore(lauf) = anisoCorr;
    avgStore(lauf) = avgCorr;
    
end

for lauf = 1:p
    ItoSMatrix=[0.5 0.5 0.5 0.5; 1 0 -1 0; 0 1 0 -1];
    S=ItoSMatrix*[I0_Furrow(lauf) I45_Furrow(lauf) I90_Furrow(lauf) I135_Furrow(lauf)]';

    aniso=sqrt(S(3,:).^2 + S(2,:).^2)./S(1,:);
         aniso=(1+aniso)./(1-aniso); %Then Ratio
         aniso(aniso<1)=NaN;
         
         ratioStore(lauf) = aniso
         
end


for lauf = 1:p
    ItoSMatrix=[0.5 0.5 0.5 0.5; 1 0 -1 0; 0 1 0 -1];
    S=ItoSMatrix*[I0_Pole(lauf) I45_Pole(lauf) I90_Pole(lauf) I135_Pole(lauf)]';

    aniso=sqrt(S(3,:).^2 + S(2,:).^2)./S(1,:);
         aniso=(1+aniso)./(1-aniso); %Then Ratio
         aniso(aniso<1)=NaN;
         
         ratioStorePole(lauf) = aniso
         
end



for lauf = 1:p
  S=ItoSMatrix*[I0_Furrow(lauf) I45_Furrow(lauf) I90_Furrow(lauf) I135_Furrow(lauf)]';
     aniso=sqrt(S(3,:).^2 + S(2,:).^2)./S(1,:);

   
         
         anisoStore(lauf) = aniso
         
end



h= figure
plot(ratioStore,'r')
hold on
plot(ratioStorePole,'b')

%axis([timeVec(1) 1000 0 25])  
%yL = get(gca,'YLim');
%line([ingressionTime ingressionTime],yL,'Color','r');
 
legend('Ratio cleavage furrow','Ratio pole')

 
xlabel ('Time [s]','FontSize', 16);
ylabel('Ratio [A.U.]','FontSize', 16);
title('Polarization ratio polscope','FontSize', 16);

print(h,'-dpdf', ['_polscopeRatio.pdf']);%tifCurvetifFilename);



for lauf = 1:AnalysisEnd
    
[orientCorr, anisoCorr, avgCorr]=ComputeFluorAnisotropy(I0,I45b,I90b,I135b,'anisotropy','I0bleach',I0b,'bleachROI',1,'normFactors',normFactors);
   % [orient, aniso, avg]=...
    %    ComputeFluorAnisotropy(I0_Pole(lauf),I45_Pole(lauf),I90_Pole(lauf),I135_Pole(lauf),...
    %    'ratio','BlackLevel',arg.BlackLevel,'normFactors',arg.normFactors);
    
    %orientStore{lauf} = orient;
    anisoStorePole(lauf) = anisoCorr;
    avgStorePole(lauf) = avgCorr;
    
end

ratioFurrowPole = ratioStore./ratioStorePole

figure(1)
figure(2)
plot(ratioFurrowPole)
hold on

plot(ratioStore)

ratioStore = anisoStore
%
close all
plot(anisoStorePole)
plot(avgStorePole)
plot(avgStore)


plot(avgAverage)

[AX,H1,H2] = plotyy(1:AnalysisEnd,anisoStore,1:AnalysisEnd,avgStore);

[AX,H1,H2] = plotyy(1:AnalysisEnd,anisoStore,1:AnalysisEnd,cleavageFurrowDiameter);

         h = figure(1) 
        plot((0-(newDistance-1):newDistance)*voxelX_mum,linescanSelector{lauf}')
       
     
      
        hold on
          set(gca,'FontSize',16,'FontName', 'Arial')
      plot(x_axOrig_store{lauf},y_fittedOrig_store{lauf},'r')
       axis([-15 +15 -1 1]) 
  title('Anisotropy SiR-actin green/red Exp. 2433','FontSize', 16,'FontName', 'Arial')
   %title('Gaussian fit to SiR width at the cleavage furrow. 2466, SiRActin','FontSize', 16)
    xlabel ('Distance [µm]','FontSize', 16,'FontName', 'Arial');
    ylabel('Ratio green/red  [A.U.]' ,'FontSize', 16,'FontName', 'Arial');
    legend([num2str(FWHMOrigStore(lauf)) 'µm' ' at ' num2str(timeVec(lauf)) 's' ])
        pause(0.2)
       
       
      % print(h,'-dpdf', [curdir '\' 'greenFit' '\' tifFilename, 'frame_', num2str(lauf) ,'_GaussFit_5percAVG.pdf']);%tifCurvetifFilename);
   
      print(h,'-dpdf', [curdir '\' ,'redFit' '\' , 'frame_', num2str(lauf) ,'_GaussFit_FWHM_anisotropy.pdf']);%tifCurvetifFilename);
   close all
    
        
   
   cd(curdir)
   
       save([saveFileName,'.mat'])
   
            dummyVariable = (1:p)';
            dummyVariable(:) = 0;
           

            saveVariables = {};
            saveVariables{1} = timeVec';
            saveVariables{2} = cleavageFurrowDiameter';
            saveVariables{3} =  avgStore';
            
            saveVariables{4} =   avgStorePole';
             saveVariables{5} =  anisoStore';
            saveVariables{6} =   anisoStore';
            
           
            
           
            
            clear csvData;
            csvData=saveVariables;
    
            header= ['time, cleavageFurrowDiameter, intensityCleavageFurrow, intensityPole',...
                ',ansioCleavageFurrow, anisoPole']
            
            outid = fopen([saveFileName,'Analysis.csv'], 'w+');

            fprintf(outid, '%s', header);
            fclose(outid);
            dlmwrite ([saveFileName,'Analysis.csv'],csvData,'roffset',1,'-append')
            
            
   
   
   
   
   
   
   
   
   [anisoPath,orientPath,avgPath]=writeComputedChannels(aniso,orient,avg,directory,idZ,idF,arg.Channel,arg.suffix,bitDepth);
        filenames=cat(1,filenames,{anisoPath,orientPath,avgPath}');

        
        
        
        
        
        
    function  [anisoPath,orientPath,avgPath]=writeComputedChannels(aniso,orient,avg,directory,Slice,Frame,Channel,suffix,bitDepth)
Zstr=num2str(Slice,'%03u'); Tstr=num2str(Frame,'%04u'); Chstr=num2str(Channel,'%u');
anisoPath=[directory '/analysis/I1-aniso' '_Z' Zstr '_T'  Tstr '_Ch' Chstr suffix '.tif'];
orientPath=[directory '/analysis/I2-orient' '_Z' Zstr '_T'  Tstr '_Ch' Chstr suffix '.tif'];
avgPath=[directory '/analysis/I3-avg' '_Z' Zstr '_T'  Tstr '_Ch' Chstr  suffix '.tif'];

    switch(bitDepth)
        case 'uint8'
            aniso=(2^8-1)*aniso;
            orient=(180/pi)*orient;
            imwrite(uint8(aniso),anisoPath);
            imwrite(uint8(orient),orientPath);
            imwrite(uint8(avg),avgPath);
            
        case 'uint16'
            aniso=(2^16-1)*aniso;
            orient=100*(180/pi)*orient;
            imwrite(uint16(aniso),anisoPath);
            imwrite(uint16(orient),orientPath);
            imwrite(uint16(avg),avgPath);
            
        case 'single'
            saveastiff(anisoPath,single(aniso));
            saveastiff(orientPath*(180/pi),single(orient));
            saveastiff(avgPath,single(avg));
    end
       end