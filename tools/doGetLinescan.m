function [ S1Linescan, S2Linescan] = doGetLinescan(S1Norm, S2Norm,splineFitOutline)

lauf = 1

         cNew1Ref = splineFitOutline(:,1);
         rNew1Ref = splineFitOutline(:,2);

         MIJ.createImage(S1Norm(:,:,lauf));
         MIJ.run('setLine8');
         MIJ.setRoi( [ rNew1Ref'; cNew1Ref'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yRed1{lauf} = MIJ.getColumn('y');
         xRed1{lauf} = MIJ.getColumn('x');
         MIJ.run('closeAllWindows');
         
         MIJ.createImage(S2Norm(:,:,lauf));
         MIJ.run('setLine8');
         MIJ.setRoi( [ rNew1Ref'; cNew1Ref'], ij.gui.Roi.POLYLINE);
         MIJ.run('getLinescanRed');
         yRed2{lauf} = MIJ.getColumn('y');
         xRed2{lauf} = MIJ.getColumn('x');
         MIJ.run('closeAllWindows');
         
         S1Linescan = [yRed1,xRed1];
         S2Linescan = [yRed2,xRed2];
         