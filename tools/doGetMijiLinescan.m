function [cNew1Ref rNew1Ref] = doGetMijiLinescan(imgMerge)

[mm nn pp] = size(imgMerge)
pp=1;
for lauf = 1:pp
    
    
cSort1{lauf}  = 50 ;
rSort1{lauf} = 50;


    MIJ.createImage(imgMerge(:,:,:,lauf));
    MIJ.run('8-bit');
    MIJ.run('Stack to RGB');
   % MIJ.selectWindow('Import from Matlab');
    MIJ.setRoi( [cSort1{lauf}';rSort1{lauf}'], ij.gui.Roi.POLYLINE);
    k = waitforbuttonpress 
    MIJ.run('saveMacro');
    sROI= ReadImageJROI('Roi.zip');
    tmp = sROI{:,:}.mnCoordinates;
    cNew1Ref{lauf} = tmp(:,1);
    rNew1Ref{lauf} = tmp(:,2);
    
    MIJ.run('closeAllWindows');
    
    

end