function[Furrow1 Furrow2 Curve1 Curve2]= doPickFurrowSecondCurve(imgMerge)
       lauf = 1;
                MIJ.createImage(imgMerge(:,:,lauf));
                MIJ.run('8-bit');
              
                %MIJ.run('Stack to RGB');
               
                MIJ.setRoi( [50;50], ij.gui.Roi.POINT); 
                k = waitforbuttonpress 
                coords =    MIJ.getRoi(1);
                cFurrow1(lauf) = coords(1);
                rFurrow1(lauf) = coords(2);
                 
                MIJ.setRoi( [50;50], ij.gui.Roi.POINT); 
                k = waitforbuttonpress 
                coords =    MIJ.getRoi(1);
                cFurrow2(lauf) = coords(1);
                rFurrow2(lauf) = coords(2);
                
                MIJ.setRoi( [50;50], ij.gui.Roi.POINT); 
                k = waitforbuttonpress 
                coords =    MIJ.getRoi(1);
                cCurve1(lauf) = coords(1);
                rCurve1(lauf) = coords(2);
                
                MIJ.setRoi( [50;50], ij.gui.Roi.POINT); 
                k = waitforbuttonpress 
                coords =    MIJ.getRoi(1);
                cCurve2(lauf) = coords(1);
                rCurve2(lauf) = coords(2);
                
                 MIJ.run('closeAllWindows');

                
                Furrow1 = [cFurrow1 rFurrow1];
                Furrow2 = [cFurrow2 rFurrow2];
                
                Curve1 =  [cCurve1 rCurve1];
                Curve2 =  [cCurve2 rCurve2];