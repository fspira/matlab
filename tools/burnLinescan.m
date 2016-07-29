function imgOut = burnLinescan(img,burnLinesacan)

for funcRun = 1:length(burnLinesacan)
    
    img(burnLinesacan(funcRun,2),burnLinesacan(funcRun,1)) = 255;
end

imgOut = img;
