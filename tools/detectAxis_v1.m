function [xCoord,yCoord,regionIntersect,selectedRoi] = detectAxis(img,cc,x,y,boundary, sW)    
%cc.Centroid(1,1),cc.Centroid(1,2),x,y,

%img =L2;
%cc = s1;
%sW = 400;


%%%%%%% The program has fit an ellipse to the ROI and this part of the
%%%%%%% program detects the intersection from the center of the ellipse in
%%%%%%% an angle to the ROI. This is done to chop the original ROI into
%%%%%%% subrois containing each flank of the ingressing furrow and the two
%%%%%%% poles. The last part of this function determines the mid point of
%%%%%%% each subroi.

close all

    img = uint16(img);  

    imshow(img(:,:))
     hold on
     plot([cc.Centroid(1,1),x], [cc.Centroid(1,2),y])
     hLine = imline(gca,[cc.Centroid(1,1),x], [cc.Centroid(1,2),y]); 
     binaryImage1 = hLine.createMask();
     normImg1 = normalizedImage(binaryImage1);
     
    
    
  
    
   
    
        
         MIJ.createImage(normImg1(:,:));
         MIJ.run('setLine8');
         MIJ.setRoi( [ boundary(:,2)'; boundary(:,1)'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yRedAxis1 = MIJ.getColumn('y');
         xRedAxis1 = MIJ.getColumn('x');
         MIJ.run('closeAllWindows');
         
         regionIntersect = find(yRedAxis1==max(yRedAxis1))
         
         
         %%%%%% Select a certain window for later analysis. In case the
         %%%%%% intersecting pixel is close to the end or beginning of the
         %%%%%% line, the ROI has to be concatenated. This part tests
         %%%%%% whether the center is close to start or end and if so, will
         %%%%%% concatenate the roi.
     
      
      if regionIntersect < (sW/2)+1
          catVar = (sW/2)-regionIntersect;
         
          selectedRoiTmp1 = boundary(length(boundary)-catVar:length(boundary),:)
          selectedRoiTmp2 = boundary(1:regionIntersect+(sW/2),:)
          
      
          
         selectedRoi = cat(1, selectedRoiTmp1, selectedRoiTmp2)
         
      elseif length(boundary) < regionIntersect+(sW/2)
          catVar = length(boundary)-regionIntersect;
          catVar = (sW/2)-catVar-1;
          
          
          selectedRoiTmp1 = boundary(regionIntersect-(sW/2):length(boundary),:)
          selectedRoiTmp2 = boundary(1:catVar,:)
          
          selectedRoi = cat(1, selectedRoiTmp1, selectedRoiTmp2)
      else
           selectedRoi = boundary(regionIntersect-(sW/2):regionIntersect+(sW/2),:)
      end
      
      %%%%%%% This part identifies the intersecting element from the
      %%%%%%% subroi. The resulting value will be the center pixel of the
      %%%%%%% ingressing furrow and also the pole
      
        MIJ.createImage(normImg1(:,:));
         MIJ.run('setLine8');
      MIJ.setRoi( [ selectedRoi(:,2)'; selectedRoi(:,1)'], ij.gui.Roi.POLYLINE);
      MIJ.run('getLinescanRed');
      yRedAxis1 = MIJ.getColumn('y');
      xRedAxis1 = MIJ.getColumn('x');
         MIJ.run('closeAllWindows');
         
         subRegionIntersect = find(yRedAxis1==max(yRedAxis1))
         xCoord = selectedRoi(subRegionIntersect,2);
         yCoord = selectedRoi(subRegionIntersect,1);
      
end
  