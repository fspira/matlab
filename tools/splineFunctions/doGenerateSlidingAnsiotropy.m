function [anisoOut] = doGenerateSlidingAnsiotropy(I0,I45,I90,I135,slidingWindow,ItoSMatrix)
%lauf = 1
% slidingWindow =1;
% y0_tmp = yI0_2{1}
% y45_tmp = yI45_2
% y90_tmp = yI90_2
% y135_tmp = yI135_2


y0_tmp = I0;
y45_tmp = I45;
y90_tmp = I90;
y135_tmp = I135;


%%%%% Format the strings to account for boundary pixels - in this case I
%%%%% use padding

y0_tmp_padding = paddingWindow(y0_tmp',slidingWindow);
y45_tmp_padding = paddingWindow(y45_tmp',slidingWindow);
y90_tmp_padding = paddingWindow(y90_tmp',slidingWindow);
y135_tmp_padding = paddingWindow(y135_tmp',slidingWindow);

increment = 1+((slidingWindow-1)/2);
rangeSelector = (slidingWindow-1)/2;

%slidingWindow = 5

for i=1:length(y0_tmp_padding)-(slidingWindow-1)
    
    
    
    yI0_Interpolated(i) = mean(y0_tmp_padding(increment-rangeSelector:increment+rangeSelector));
    yI45_Interpolated(i) = mean(y45_tmp_padding(increment-rangeSelector:increment+rangeSelector));
    yI90_Interpolated(i) = mean(y90_tmp_padding(increment-rangeSelector:increment+rangeSelector));
    yI135_Interpolated(i) = mean(y135_tmp_padding(increment-rangeSelector:increment+rangeSelector));
    
    increment = increment+1;
    
end




for i = 1:length(yI135_Interpolated)
    S=ItoSMatrix*[yI0_Interpolated(i) yI45_Interpolated(i) yI90_Interpolated(i) yI135_Interpolated(i)]';
    aniso=sqrt(S(3,:).^2 + S(2,:).^2)./S(1,:);
    
    
    
    anisoOut(i) = aniso;
    
end

end
% lauf =1
% 
% flank2Tmp = flank2StoreNew{lauf};
% file1Tmp = I0File(:,:,lauf);
% fileTmp(:,:) = 0;
% fileTmp = double(fileTmp);
% file1Tmp = double(file1Tmp);
% cFurrowContour2 =  flank2Tmp(:,1);
% rFurrowContour2 =  flank2Tmp(:,2);
% 
% 
% for lauf = 1:length(anisoStore)
%    
%     file1Tmp(cFurrowContour2(lauf),rFurrowContour2(lauf)) = anisoStore(lauf);
%     
% end
% 
% file1Norm = normalizedImage(file1Tmp);
% 
% imshow(file1Tmp,[])
% 
% cmap = jet;
%  cmap(1,3) = 0;
% 
% imshow(splineDilate(fileTmp,file1Norm,flank2Tmp, [1,3], 3), cmap)
% 
% 
