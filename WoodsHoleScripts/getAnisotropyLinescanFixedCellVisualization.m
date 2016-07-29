%%%%% getanisotropyLinescanforFixedCells

%%%%%% Generate nice anisotropy colormap for 5 frame poslcope data - load
%%%%%% workspace of the respective cell load workspace first

curdir = pwd;
clear newLinescan
%fileIdx = 19
lauf =1
   
    cc = regionprops(t3Store,'all');
    B = bwboundaries(t3Store,8);
    boundary = B{1};
    
    bStore{lauf} = B{1};
    
    MIJ.createImage(imgIntensity(:,:,lauf));
    MIJ.run('setLine8');
    MIJ.setRoi( [ boundary(:,2)'; boundary(:,1)'], ij.gui.Roi.POLYLINE);
  
    k = waitforbuttonpress
    clear logData
    
    MIJ.run('getSplineCoordinates')
    logData = MIJ.getLog;
    MIJ.run('closeLogWindow')
    logDataTmp = char(logData);
    logDataTmp1 = strread(logDataTmp, '%s');
    
    xLogTmp = logDataTmp1(2:4:end);
    yLogTmp = logDataTmp1(3:4:end);
    
    
    xLog = zeros(length(xLogTmp),1,'double');
    yLog = zeros(length(yLogTmp),1,'double');
    
    for subrun = 1 : length(xLogTmp)
        
        xLog(subrun) = round(str2num(xLogTmp{subrun}));
        yLog(subrun) = round(str2num(yLogTmp{subrun}));
        
    end
    
    newLinescan{lauf} = [xLog yLog]
    
    
    
    
    MIJ.run('closeAllWindows');


            AnalysisEnd = 1;
                
for lauf = 1:AnalysisEnd %length(redStack)
        
    
    
    
     
    %%%% Load XY coordinates of the two flanking regions into the tmp
    %%%% variable which is later used for analysis. This variable will be
    %%%% updated every loop
         
         flank1Tmp = newLinescan{lauf};
        
         
         %%%% Load the red image
         MIJ.createImage(imgAniso(:,:,lauf));
         MIJ.run('setLine8');%%%% increase linescan width
         MIJ.setRoi( [ flank1Tmp(:,1)'; flank1Tmp(:,2)'], ij.gui.Roi.POLYLINE);%%%% Draw line into the image
         
         
         %%%%% Update InterpolatedLine
          clear logData
          MIJ.run('getSplineCoordinates')%%%% Get interpolated coordinates of the line
          logData = MIJ.getLog;%%%%% Reading MIJI log file to get coordinates
          
          MIJ.run('closeLogWindow')%%%% Close log file
         
          MIJ.run('getLinescanRed');%%%% Get intensities under the line
         yRed1{lauf} = MIJ.getColumn('y');%%%% Intensity values
         xRed1{lauf} = MIJ.getColumn('x');%%%% Running variable for each pixel can be discarded
        logDataTmp = char(logData);%%%%% convert values to characters - important to search for a string in the next step
        logDataTmp1 = strread(logDataTmp, '%s'); %%%% Sort log data using regular expression
        
         xLogTmp = logDataTmp1(2:4:end); %%%% Save X coordinates of the interpolated linescan
         yLogTmp = logDataTmp1(3:4:end);%%% Save Y coordinates of the interpolatd linescan
         
             
        xLog = zeros(length(xLogTmp),1,'double'); %%%% Pre allocate X and Y variables
        yLog = zeros(length(yLogTmp),1,'double');
         
         for subrun = 1 : length(xLogTmp)%%%% Save XY coordinates of the interpolated linescan
             
             xLog(subrun) = round(str2num(xLogTmp{subrun}));
             yLog(subrun) = round(str2num(yLogTmp{subrun}));
             
         end
         
         flank1StoreNew{lauf} = [xLog yLog]%%% Save coordinates into store array
       
        
             MIJ.run('closeAllWindows');
          flank1Tmp = flank1StoreNew{lauf};
            % flank1Tmp = flank1StoreNew{lauf};
         
             %%%%%% Get lineprofiles from the green channel
         MIJ.createImage(imgAniso(:,:,lauf));
         MIJ.run('setLine8');
         MIJ.setRoi( [ flank1Tmp(:,1)'; flank1Tmp(:,2)'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yAniso{lauf} = MIJ.getColumn('y');
         xAniso{lauf} = MIJ.getColumn('x');
         MIJ.run('closeAllWindows');
         
         
         
        
         
         
end



slidingWindow =7




cmap = jet;
 cmap(1,3) = 0;
 
for lauf = 1:p




file1Tmp = imgAniso(:,:,lauf);
file1Tmp(:,:) = 0;
file1Tmp = double(file1Tmp);
file2Tmp = file1Tmp;

cFurrowContour1 =  flank1Tmp(:,1);
rFurrowContour1 =  flank1Tmp(:,2);


anisoOutPlot1 = yAniso{lauf};

for i = 1:length(rFurrowContour1)
   
    file1Tmp(rFurrowContour1(i),cFurrowContour1(i)) = anisoOutPlot1(i);
    
end


%file1Norm = normalizedImage(file1Tmp);
file1Norm = file1Tmp;
imshow(file1Norm,[])


 
 
 [res1,res_tmp1]=splineDilate(file1Tmp,file1Norm,[cFurrowContour1 rFurrowContour1], [2,4], 3);
 

 
 
resNorm = normalizedImage(res1);
resStoreNorm(:,:,lauf) = resNorm;
resStore(:,:,lauf) = res1;
 
imshow(res1,[])
colormap(cmap)



pause(0.5)

end

%imwrite(resNorm,cmap,'normColormap')

tiffwrite_mat(resStoreNorm, 'Norm_linewidth8_slidingWindow7')


tiffwrite_mat(resStore, [curdir,'/','Orig_linewidth8_slidingWindow7'])

[m n p] = size(imgAniso)
imgIntensity

resNorm


[mm nn pp] = size(imgMid)

imgSave(:,:,1) = imgIntensity;
imgSave(:,:,2) = imgAniso;
imgSave(:,:,3) = resNorm;
imgSave(:,:,4) = imgMid(:,:,1);

 save([I0File,'_Linescan.mat'])
 
 
 tiffwrite_mat(imgSave,[fileIdx,'_mergedFiles.tif'])

clear all