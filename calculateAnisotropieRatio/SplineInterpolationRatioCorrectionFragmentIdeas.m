 
    for subsubrun = 1:9
        spline = cscvn(I1_coords(:,subsubrun:10:end));
        fnplt(spline,'or',2.5); hold on

        t = 1:9:max(spline.breaks);
        cv = fnval(spline, t);
        cdv = fnval(fnder(spline), t);
        quiver(cv(1,:),cv(2,:), cdv(1,:),cdv(2,:));

        cvSub{subsubrun} = cv;
        cvdSub{subsubrun} = cdv;

    end
        concA = cvSub{1}
        concB = cvSub{2}
        concC = cvSub{3}
        concD = cvSub{4}
        concE = cvSub{5}
        concF = cvSub{6}
        concG = cvSub{7}
        concH = cvSub{8}
        concI = cvSub{9}
        
        
        counter =1;
        
    for lauf = 1:length(concA)
        
       
       
        cvSConcTmp(:,counter:counter+8) = [concA(:,lauf), concB(:,lauf), concC(:,lauf),concD(:,lauf),concE(:,lauf),concF(:,lauf),concG(:,lauf),concH(:,lauf),concI(:,lauf)]
        
        counter = counter + 9;
        
    end
    
    
        concA =  cvdSub{1}
        concB =  cvdSub{2}
        concC =  cvdSub{3}
        concD =  cvdSub{4}
        concE =  cvdSub{5}
        concF =  cvdSub{6}
        concG =  cvdSub{7}
        concH =  cvdSub{8}
        concI =  cvdSub{9}
        
        
        counter =1;
        
    for lauf = 1:length(concA)
        
       
       
        cvdConcTmp(:,counter:counter+8) = [concA(:,lauf), concB(:,lauf), concC(:,lauf),concD(:,lauf),concE(:,lauf),concF(:,lauf),concG(:,lauf),concH(:,lauf),concI(:,lauf)]
        
        counter = counter + 9;
        
    end
    
      
    
    %%%%%%% Identifiy points within the spline that are present in the
    %%%%%%% original image
    
    flank1Shift = flank1'
    
    firstIdx = find(flank1Shift(1,lauf)  == round(cvSConcTmp(1,:))) 
    
    secondIdx = find(flank1Shift(2,lauf) == round(cvSConcTmp(2,:)))
    
   
    idxCounter = 1;
    for subrun = 1:length(firstIdx)
        
        newIdx = find(firstIdx(subrun) == secondIdx)
        
        if isempty(newIdx) ==1
        else
            newIdxStore(idxCounter) =  find(firstIdx(subrun) == secondIdx)
            idxCounter = idxCounter +1;
        end
    
    end
    