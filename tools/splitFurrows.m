function [xCoord,yCoord,subRegionIntersect,selectedRoi, yRedAxis1] = splitFurrows(img,x,y,cNew1Ref,rNew1Ref, sW,imgOrig)    
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

%cNew1Ref = cNewStore;
%rNew1Ref = rNewStore;

        % lauf = 1
            boundaryRoi =[cNew1Ref rNew1Ref];
              

%%%%%%% The program has fit an ellipse to the ROI and this part of the
%%%%%%% program detects the intersection from the center of the ellipse in
%%%%%%% an angle to the ROI. This is done to chop the original ROI into
%%%%%%% subrois containing each flank of the ingressing furrow and the two
%%%%%%% poles. The last part of this function determines the mid point of
%%%%%%% each subroi.
%img= redImg(:,:,1);

 
    img(:,:)=0;


    img = uint16(img);  

    img(round(y),round(x)) = 64000;
   
   
  %  MIJ.closeAllWindows
   
         MIJ.createImage(img(:,:));
         MIJ.run('setLine8');
         MIJ.setRoi( [  cNew1Ref'; rNew1Ref'], ij.gui.Roi.POLYLINE);
         
         %%%%% reads the spline coordinates and opens a Log window, the log
         %%%%% window values are then converted into a variable
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
      %   MIJ.run('closeAllWindows');
         
       
         
       %  interpolatedSpline = [round(Log(:,2)) round(Log(:,3))];
         regionIntersect = find(yRedAxis1==max(yRedAxis1))
         
         
       
        
         
         %%%%%% Select a certain window for later analysis. In case the
         %%%%%% intersecting pixel is close to the end or beginning of the
         %%%%%% line, the ROI has to be concatenated. This part tests
         %%%%%% whether the center is close to start or end and if so, will
         %%%%%% concatenate the roi.
     
         boundaryRoi = [xLog  yLog];
      
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
        MIJ.createImage(img(:,:));
         MIJ.run('setLine8');
     % MIJ.setRoi( [ boundaryRoi(:,1)';  boundaryRoi(:,2)'], ij.gui.Roi.POLYLINE);
      MIJ.setRoi( [ selectedRoi(:,1)'; selectedRoi(:,2)'], ij.gui.Roi.POLYLINE);
      MIJ.run('getLinescanRed');
      yRedAxis1 = MIJ.getColumn('y');
      xRedAxis1 = MIJ.getColumn('x');
       %  MIJ.run('closeAllWindows');
       MIJ.run('closeResultsWindow')
         MIJ.closeAllWindows
         subRegionIntersect = find(yRedAxis1==max(yRedAxis1))
         xCoord = selectedRoi(subRegionIntersect,2);
         yCoord = selectedRoi(subRegionIntersect,1);
         
         
         %%% MeasureLinescan from average Image S1 and S2 channel overlay
         MIJ.createImage(imgOrig(:,:));
         MIJ.run('setLine8');
     % MIJ.setRoi( [ boundaryRoi(:,1)';  boundaryRoi(:,2)'], ij.gui.Roi.POLYLINE);
      MIJ.setRoi( [ selectedRoi(:,1)'; selectedRoi(:,2)'], ij.gui.Roi.POLYLINE);
      MIJ.run('getLinescanRed');
      yRedAxis1 = MIJ.getColumn('y');
      MIJ.run('closeResultsWindow')
      MIJ.closeAllWindows
         
         
         
         
      
end
  