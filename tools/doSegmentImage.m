function t3Store = doSegmentImage(S1)
test  = S1;%redStack(:,:,1);
figure(1)
imshow(test(:,:,:),[]);
%marker = testNorm;
testNorm = normalizedImage(test);
marker = false(size(testNorm));


for lauf = 1:2

fh = figure(1);
 title('Mark the Background')
frapBkg = roipoly(testNorm(:,:,1));
close(fh);
                

                marker(frapBkg) = 255;
                
end

testNorm = normalizedImage(S1(:,:,1));


h = fspecial('gaussian', 10, 10) 
redStackGauss = imfilter(S1(:,:,:),h);
imshow(redStackGauss(:,:,1),[])

[m n p]= size(S1);

for lauf = 1:p
             testNorm = normalizedImage(redStackGauss(:,:,lauf));

                mp = imimposemin(testNorm,marker);
                L2 = watershed(mp);

                %mpReg = imregionalmin(mp);
                t3 = false(size(testNorm));
                t3(L2 == 0 ) = 255;
                t3Store(:,:,lauf) = t3;
               % imshow(mpReg(:,:,1),[]);
                
end
                
              %tiffwrite_mat(t3Store, [folderName,'_SegmentGauss.tif']);
    
close all

