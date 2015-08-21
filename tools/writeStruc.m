

for jj = 1:4
c = s(rr).(headers{jj}); 
str = ''; 
if sz(jj,1)<ii 
str = repmat(',',1,sz(jj,2)); 
%% Added code 
if sz(jj,2) ==0 %MOdified 
str = repmat(',',1,1); 
end 
%%Finish added code 
else 
if isnumeric(c) 
for kk = 1:sz(jj,2)
end
end
end
end
jj = 1:4
c = s(1).(headers{jj}); 
d = size(c); 
if isnumeric(c) 
str = [num2str(c(ii,kk)),',']; 
elseif d(2)==1 
str = ['"',c(ii),'",']; 
else 
str = ['"',c{ii,kk},'",']; 
end

 saveVariables = {};

            saveVariables{1} = timeVec(1:p)';%time{1};
            saveVariables{2} = orientation_Store(1:p)';
            
clear saveVariable_1 saveVariable_2 saveVariable_3 saveVariable_4
  
            
            names = fieldnames(s)
            
            header = ['timeStamp','filename','LCOrientation','consistency'];
          
            
            
            [ii pp]= size(s)
            
            for jj = 1:pp
               
                saveVariable_1(jj) = getfield(s(jj), 'timeSTamp')
                saveVariable_2{jj} = getfield(s(jj), 'fileName')
                saveVariable_3(jj) = getfield(s(jj), 'LCOrientation')
                saveVariable_4(jj) = getfield(s(jj), 'consistencyCheck')
            end

            
            saveVariables{1} =  saveVariable_1
            saveVariables{2} =  saveVariable_2
            saveVariables{3} =  saveVariable_3
            saveVariables{4} =  saveVariable_4
            
            
            
            clear csvData;
            csvData=saveVariables;
    
           
             
            outid = fopen(['Analysis.csv'], 'w+');

            fprintf(outid, '%s', header);
            fclose(outid);
            dlmwrite (['Analysis.csv'],csvData,'roffset',1,'-append')
            
            
         A(1,:) = {'Year','Month','Day','Hour','Weekday','kW'}   

    