
function [cCoordsNew,rCoordsNew] = setMarker(cCoords, rCoords)     



              
               % MIJ.selectWindow('Import from Matlab');
                MIJ.setRoi( [cCoords;rCoords], ij.gui.Roi.POINT); 
                k = waitforbuttonpress 
                coords =    MIJ.getRoi(1);
                cCoordsNew = coords(1);
                rCoordsNew = coords(2);
              
end

            

