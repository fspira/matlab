function t3Store = doSegmentImageCircleMarker(S1)
test  = S1;%redStack(:,:,1);
%test  = I_crop;
figure(1)
imshow(test(:,:,1),[]);
%marker = testNorm;
testNorm = test;
marker = test;
marker(:,:) = 0




fh = figure(1);
 title('Mark the Background')
  [xx,yy] = ginput(3)
close(fh);
                
for lauf = 1:3

                marker(round(xx(lauf)),round(yy(lauf))) = 255;
                
end



%h = fspecial('gaussian', 5, 5) 
redStackGauss = test%imfilter(S1(:,:,:),h);
imshow(redStackGauss(:,:,1),[])

[m n p]= size(testNorm);

for lauf = 1:p
             testNorm = normalizedImage(redStackGauss(:,:,lauf));

                mp = imimposemin(testNorm,marker);
                L2 = watershed(mp);

                %mpReg = imregionalmin(mp);
                t3 = false(size(testNorm));
                t3(L2 == 0 ) = 0;
                t3Store(:,:,lauf) = t3;
               % imshow(mpReg(:,:,1),[]);
                
end
                
              %tiffwrite_mat(t3Store, [folderName,'_SegmentGauss.tif']);
    
close all

