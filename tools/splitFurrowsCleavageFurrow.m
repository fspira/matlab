function [xCoord,yCoord,subRegionIntersect,selectedRoi] = splitFurrows(L2,x,y,boundaryRoi, orientationSet,sW,imgOrig,cc)    

%L2 = im2bw(L2);


%boundaryRoi = boundary;
%cc.Centroid(1,1),cc.Centroid(1,2),x,y,

%img =L2;
%cc = s1;
%sW = 400;
%img = redImg(:,:,lauf);
%x=xFurrow1{lauf}
%y=yFurrow1{lauf}
%cNew1Ref = cNew1Ref{lauf}
%rNew1Ref = rNew1Ref{lauf}
%sW = 180

%cNewStore1 = cNew1Ref;
%rNewStore1 = rNew1Ref;

%cNew1Ref = cNewStore1;
%rNew1Ref = rNewStore1;

        % lauf = 1
        
        k=1;
       [px py] = doSelectEllipsePoint(cc,orientationSet ,k) 
         px = sqrt(px^2);
           py = sqrt(py^2);
       px = [px cc.Centroid(1,1)]
        
          py = [py cc.Centroid(1,2)]
          
          L2= uint8(L2).*128;
          L1 = L2;
          L1(:,:)=0;
         tmp1 = doConnectTwoPoints(px,py,L1);
         tmp1 = uint8(tmp1).*128;
         L3 = tmp1+L2;
         
         %%%%%% Only show pixels with intensity of 255 grey values
         
         L4 = L3 - 254;
         L3 = L4 .*255;
         
         
        
   
         MIJ.createImage(L3(:,:,1));
         MIJ.run('setLine8');
         MIJ.setRoi( [ boundaryRoi(:,2)'; boundaryRoi(:,1)'], ij.gui.Roi.POLYLINE);
         
        
         
          MIJ.run('getSplineCoordinates')
         logData = MIJ.getLog;
       
         logDataTmp = char(logData);
         logDataTmp1 = strread(logDataTmp, '%s');
        
         xLogTmp = logDataTmp1(2:4:end);
         yLogTmp = logDataTmp1(3:4:end);
        
        
         
        xLog = zeros(length(xLogTmp),1,'double');
        yLog = zeros(length(yLogTmp),1,'double');
         
         for lauf = 1 : length(xLogTmp)
             
             xLog(lauf) = str2num(xLogTmp{lauf});
             yLog(lauf) = str2num(yLogTmp{lauf});
             
         end
         MIJ.run('getLinescanRed');
         
         
         yRedAxis1 = MIJ.getColumn('y');
         xRedAxis1 = MIJ.getColumn('x');
         
         
         
          MIJ.run('closeResultsWindow')
          MIJ.closeAllWindows
    
          
          regionIntersect = find(yRedAxis1==max(yRedAxis1))
         
   
   
   
   boundaryRoi = [round(xLog) round(yLog)];
         
   
      
      if regionIntersect < (sW/2)+1
          catVar = (sW/2)-regionIntersect;
         
          selectedRoiTmp1 = boundaryRoi(length(boundaryRoi)-catVar:length(boundaryRoi),:);
          selectedRoiTmp2 = boundaryRoi(1:regionIntersect+(sW/2),:);
          
      
          
         selectedRoi = cat(1, selectedRoiTmp1, selectedRoiTmp2);
         %%%%% Test whether the selected window fits into the length of the
         %%%%% ROI. If window size is longer than the ROI, Values from the
         %%%%% beginning have to be concatenated at the end. 
      elseif length(boundaryRoi) < regionIntersect+(sW/2)
          catVar = length(boundaryRoi)-regionIntersect;
          catVar = (sW/2)-catVar-1;
          
          
          selectedRoiTmp1 = boundaryRoi(regionIntersect-(sW/2):length(boundaryRoi),:);
          selectedRoiTmp2 = boundaryRoi(1:catVar,:);
         
          selectedRoi = cat(1, selectedRoiTmp1, selectedRoiTmp2);
      else
           selectedRoi = boundaryRoi(regionIntersect-(sW/2):regionIntersect+(sW/2),:);
      end
      
      %%%%%%% This part identifies the intersecting element from the
      %%%%%%% subroi. The resulting value will be the center pixel of the
      %%%%%%% ingressing furrow and also the pole
     % MIJ.closeAllWindows
    %    MIJ.createImage(L3(:,:));
    %     MIJ.run('setLine8');
    %  MIJ.setRoi( [ boundaryRoi(300:350,1)';  boundaryRoi(300:350,2)'], ij.gui.Roi.POLYLINE);
    %  MIJ.setRoi( [ selectedRoi(:,1)'; selectedRoi(:,2)'], ij.gui.Roi.POLYLINE);
    %  MIJ.run('getLinescanRed');
    %  yRedAxis1 = MIJ.getColumn('y');
    %  xRedAxis1 = MIJ.getColumn('x');
       %  MIJ.run('closeAllWindows');
    %   MIJ.run('closeResultsWindow')
    %     MIJ.closeAllWindows
    %     subRegionIntersect = find(yRedAxis1==max(yRedAxis1))
    subRegionIntersect = 150;
         xCoord = selectedRoi(subRegionIntersect,2);
         yCoord = selectedRoi(subRegionIntersect,1);
         
         
         %%% MeasureLinescan from average Image S1 and S2 channel overlay
        % MIJ.createImage(imgOrig(:,:));
        % MIJ.run('setLine8');
     % MIJ.setRoi( [ boundaryRoi(:,1)';  boundaryRoi(:,2)'], ij.gui.Roi.POLYLINE);
     % MIJ.setRoi( [ selectedRoi(:,1)'; selectedRoi(:,2)'], ij.gui.Roi.POLYLINE);
     % MIJ.run('getLinescanRed');
     % yRedAxis1 = MIJ.getColumn('y');
     % MIJ.run('closeResultsWindow')
     % MIJ.closeAllWindows
         
         
         
         
      
end
  