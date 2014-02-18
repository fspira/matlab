function [furrowCenter1,furrowCenter2, curveCenter1, curveCenter2]= doDetectCenterOfRoi(redStackMid ,splineFitOutline,Furrow1,Furrow2,Curve1,Curve2)

    lauf =1
   % redStackMid = redStack(:,:,lauf);
    
    cNew1Ref = splineFitOutline(:,1);
    rNew1Ref = splineFitOutline(:,2);
    
    
    cMark1 = Furrow1(1);
    rMark1 = Furrow1(2);
    
    cMark2 =  Furrow2(1);
    rMark2 =  Furrow2(2);
     
  
    
    cPoleMark1 =   Curve1(1);
    rPoleMark1 =   Curve1(2);
    
    
    cPoleMark2 =   Curve2(1);
    rPoleMark2 =   Curve2(2);
    
    
         redStackMid(:,:) = 0;
         redStackMid(rMark1,cMark1) = 65000;
         MIJ.createImage(redStackMid(:,:));
         MIJ.run('setLine12');
         MIJ.setRoi( [ rNew1Ref'; cNew1Ref'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yRedMid1{lauf} = MIJ.getColumn('y');
         xRedMid1{lauf} = MIJ.getColumn('x');
         MIJ.run('closeAllWindows');
          
        
          
         redMax1(lauf) = find(yRedMid1{lauf}==max(yRedMid1{lauf}))
         
         redStackMid(:,:) = 0;
         redStackMid(rMark2,cMark2) = 65000;
         MIJ.createImage(redStackMid(:,:));
         MIJ.run('setLine12');
         MIJ.setRoi( [ rNew1Ref'; cNew1Ref'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yRedMid2{lauf} = MIJ.getColumn('y');
         xRedMid2{lauf} = MIJ.getColumn('x');
         MIJ.run('closeAllWindows');
         
          
         redMax2(lauf) = find(yRedMid2{lauf}==max(yRedMid2{lauf}))
         redStackMid(:,:) = 0;
         
           
         %%%%%%% Max at some non furrow point at the spline
         
         redStackMid(:,:) = 0;
         redStackMid(rPoleMark1,cPoleMark1) = 65000;
         MIJ.createImage(redStackMid(:,:));
         MIJ.run('setLine12');
         MIJ.setRoi( [ rNew1Ref'; cNew1Ref'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yPoleMid1{lauf} = MIJ.getColumn('y');
         xPoleMid1{lauf} = MIJ.getColumn('x');
         MIJ.run('closeAllWindows');
          
         redPoleMax1(lauf) = find(yPoleMid1{lauf}==max(yPoleMid1{lauf}))
         
         redStackMid(:,:) = 0;
         redStackMid(rPoleMark2,cPoleMark2) = 65000;
         MIJ.createImage(redStackMid(:,:));
         MIJ.run('setLine12');
         MIJ.setRoi( [ rNew1Ref'; cNew1Ref'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yPoleMid2{lauf} = MIJ.getColumn('y');
         xPoleMid2{lauf} = MIJ.getColumn('x');
         MIJ.run('closeAllWindows');
         
         redPoleMax2(lauf) = find(yPoleMid2{lauf}==max(yPoleMid2{lauf}))
 
         furrowCenter1 = redMax1
         furrowCenter2 = redMax2
         
         curveCenter1 =  redPoleMax1
         curveCenter2 = redPoleMax2;
         

