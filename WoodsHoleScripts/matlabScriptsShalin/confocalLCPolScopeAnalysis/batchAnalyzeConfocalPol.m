ChannelOfInterest=1;
Zslices={1};

confocalPaths={'D:\Shared Data\Felix\20150422_PhalloidinAF488_mitotic\cell6\P001'...
     };
   % 'D:\Shared Data\Felix\20150422_PhalloidinAF488_mitotic\cell6\P001'...
     
    %'D:\Shared Data\Felix\20150422_PhalloidinAF488_mitotic\cell1_zStack\P001'...{1:44
    %'D:\Shared Data\Felix\20150422_PhalloidinAF488_mitotic\cell2_zStack\P001'...1:48
    %'D:\Shared Data\Felix\20150422_PhalloidinAF488_mitotic\cell3_zStack\P001'...1:52
   
    %D:\Shared Data\Felix\20150421_PhalloidinAF488_mitotic\cell1_ZStack\P001'...{1:21
    %'D:\Shared Data\Felix\20150421_PhalloidinAF488_mitotic\cell2_ZStack\P001'...1:19
    %'D:\Shared Data\Felix\20150421_PhalloidinAF488_mitotic\cell3_ZStack\P001'...1:18
    
   
        
     %%% 8 perc
    %'D:\Shared
    %Data\Felix\20150423_PhalloidinAF488_mitotic\cell3_zStack\P001',...{1:10};
    %'D:\Shared Data\Felix\20150423_PhalloidinAF488_mitotic\cell4_zStack\P001'...{1:9}
    %'D:\Shared Data\Felix\20150423_PhalloidinAF488_mitotic\cell8_zStack\P001'... {1:20}
    
configNames={'100x_488nm_6perc'};
Frames=1;
isopath='D:\Shared Data\Felix\20150422_PhalloidinAF488_mitotic\Calib_20150421_lessIntensity';
normFactors=computeConfocalNormalizations(isopath,488,1.4,80);


for idx=1:numel(confocalPaths)
    processConfocalPolData(confocalPaths{idx},'normFactors',normFactors,'ZSlices',Zslices{idx},'Frames',Frames,'Channel',ChannelOfInterest,...
        'displayStatus',true,'suffix',['_' configNames{1}]);  
    exportConfocalPolData(confocalPaths{idx},'Zslices',Zslices{idx},'Frames',Frames,'anisoCeiling',0.5,'avgCeiling',160,'suffix',['_' configNames{1}]);
end

%%%%% 20150423_PhalloidinAF488_mitotic cell6 
%%
ChannelOfInterest=1;
%Zslices={1:21,1:19,1:18,1:14,1:44,1:48,1:52};
Zslices= {4:16};
confocalPaths={'D:\Shared Data\Felix\20150423_PhalloidinAF488_mitotic\cell9_zStack\P001'...
      };
  %  'D:\Shared Data\Felix\20150423_PhalloidinAF488_mitotic\cell6_zStack\P001'...{1:12}
 % 'D:\Shared Data\Felix\20150423_PhalloidinAF488_mitotic\cell7_zStack\P001'... {1:15}
   

configNames={'100x_488nm_8perc'};
Frames=1;
isopath='D:\Shared Data\Felix\20150422_PhalloidinAF488_mitotic\Calib_20150421_lessIntensity';
normFactors=computeConfocalNormalizations(isopath,488,1.4,80);


for idx=1:numel(confocalPaths)
    processConfocalPolData(confocalPaths{idx},'normFactors',normFactors,'ZSlices',Zslices{idx},'Frames',Frames,'Channel',ChannelOfInterest,...
        'displayStatus',true,'suffix',['_' configNames{1}]);  
    exportConfocalPolData(confocalPaths{idx},'Zslices',Zslices{idx},'Frames',Frames,'anisoCeiling',0.5,'avgCeiling',160,'suffix',['_' configNames{1}]);
   
end


%% folder 20150424


ChannelOfInterest=1;
%Zslices={1:21,1:19,1:18,1:14,1:44,1:48,1:52};
Zslices= {1:13,1:53};
confocalPaths={'D:\Shared Data\Felix\20150424_PhalloidinAF488_mitotic\cell2_zStack\P001'...
    'D:\Shared Data\Felix\20150424_PhalloidinAF488_mitotic\cell3_zStack\P001'...
      };
  %  'D:\Shared Data\Felix\20150423_PhalloidinAF488_mitotic\cell6_zStack\P001'...{1:12}
 % 'D:\Shared Data\Felix\20150423_PhalloidinAF488_mitotic\cell7_zStack\P001'... {1:15}
   

configNames={'100x_488nm_8perc'};
Frames=1;
isopath='D:\Shared Data\Felix\20150424_PhalloidinAF488_mitotic\calib';
normFactors=computeConfocalNormalizations(isopath,488,1.4,80);


for idx=1:numel(confocalPaths)
    processConfocalPolData(confocalPaths{idx},'normFactors',normFactors,'ZSlices',Zslices{idx},'Frames',Frames,'Channel',ChannelOfInterest,...
        'displayStatus',true,'suffix',['_' configNames{1}]);  
    exportConfocalPolData(confocalPaths{idx},'Zslices',Zslices{idx},'Frames',Frames,'anisoCeiling',0.5,'avgCeiling',160,'suffix',['_' configNames{1}]);
   
end