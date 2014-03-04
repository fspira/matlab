function [filenameCenter, filenamePole,  chromatinDistance_center,aniso_vertical_center,  ratio_vertical_center, predictedTime_vertical_center,fileName_vertical_center,...
     aniso_horizontal_center,ratio_horizontal_center,predictedTime_horizontal_center, fileName_horizontal_center...
      chromatinDistance_pole, aniso_vertical_pole,ratio_vertical_pole, predictedTime_vertical_pole, fileName_vertical_pole...
     aniso_horizontal_pole,ratio_horizontal_pole, predictedTime_horizontal_pole, fileName_horizontal_pole] = doLoadFiles(refDir,loadPath)

%refDir = refdir_3
%loadPath = loadPath_3
filenameCenter = {};
filenamePole = {};
centerData = {};
poleData = {};

for r = 3:length({refDir.name})
    
  
         fileList = getfield(refDir,{r},'name')
            
                
                
          
            
    
                if findstr(fileList, '.DS_Store')
                    
                      else

                        suffixCenter = 0;
                        suffixPole = 0;


                        if suffixCenter <= findstr(fileList,'Center.')
                            suffixCenter = findstr(fileList,'Center') 
                        else
                             suffixCenter = 5;
                        end

                        if suffixPole <= findstr(fileList,'Pole.') 
                            suffixPole = findstr(fileList,'Pole') 
                        else
                            suffixPole = 5;
                        end

                        if fileList(suffixCenter:suffixCenter+5) == ('Center') 


                                filenameCenter{length(filenameCenter)+1} = fileList;

                        elseif fileList(suffixPole:suffixPole+3) == ('Pole')

                            filenamePole{length(filenamePole)+1} = fileList;

                        end
                end
            
end
            
       
         counterHorizontal = 1;
         counterVertical = 1;

      for lauf = 1:length(filenameCenter)
            
             load([loadPath ,'/', filenameCenter{lauf}])
            
            
              f = fieldnames(s);
            

           

                if s.(f{1}).orientation.orientation > 45 

                        chromatinDistance_center(counterVertical)= s(1).(f{1}).chromatin_distance.chromatin_distance


                        aniso_vertical_center(counterVertical) = s(1).(f{1}).averageRatioD_Vertical.averageRatioD_Vertical  




                        ratio_vertical_center(counterVertical) = s(1).(f{1}).averageRatio_Vertical.averageRatio_Vertical 





                        predictedTime_vertical_center(counterVertical) = s.(f{1}).predicted_time.predicted_time
                        fileName_vertical_center{counterVertical} = f{1};

                       
                        counterVertical = counterVertical+1

                else

                      aniso_horizontal_center(counterHorizontal) = s(1).(f{1}).averageRatioD_Horizontal.averageRatioD_Horizontal

                      ratio_horizontal_center(counterHorizontal) = s(1).(f{1}).averageRatio_Horizontal.averageRatio_Horizontal


                        predictedTime_horizontal_center(counterHorizontal) = s.(f{1}).predicted_time.predicted_time

                       fileName_horizontal_center{counterHorizontal} = f{1};
                       
                        counterHorizontal = counterHorizontal+1


                end
       
   clear s

            
            
      end
            
                 
        
      
      counterHorizontal = 1;
      counterVertical = 1;    
      for subRun = 1:length(filenamePole)
            
           PoleData{subRun} = load([loadPath ,'/', filenamePole{subRun}])
           
          load([loadPath ,'/', filenamePole{subRun}])
           

      f = fieldnames(s);
            
                if s.(f{1}).orientation.orientation > 45 


                        chromatinDistance_pole(counterVertical)= s(1).(f{1}).chromatin_distance.chromatin_distance
                        aniso_vertical_pole(counterVertical) = s(1).(f{1}).averageRatioD_Vertical.averageRatioD_Vertical  
                        ratio_vertical_pole(counterVertical) = s(1).(f{1}).averageRatio_Vertical.averageRatio_Vertical

                        predictedTime_vertical_pole(counterVertical) = s.(f{1}).predicted_time.predicted_time
                        fileName_vertical_pole{counterVertical} = f{1};

                       
                        counterVertical = counterVertical+1

                else


                      aniso_horizontal_pole(counterHorizontal) = s(1).(f{1}).averageRatioD_Horizontal.averageRatioD_Horizontal

                      ratio_horizontal_pole(counterHorizontal) = s(1).(f{1}).averageRatio_Horizontal.averageRatio_Horizontal

                        predictedTime_horizontal_pole(counterHorizontal) = s.(f{1}).predicted_time.predicted_time

                       fileName_horizontal_pole{counterHorizontal} = f{1};
                       
                        counterHorizontal = counterHorizontal+1


                end
            

        clear s

            
      end
            
            
            
            
end

