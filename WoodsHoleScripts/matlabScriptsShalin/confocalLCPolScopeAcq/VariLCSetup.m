function [ variLCObj,LCsettings ] = VariLCSetup( COMno, palleteFile)


% Find a VISA-Serial object at specified COM port.
InstrStr=['ASRL' int2str(COMno) '::INSTR'];
variLCObj= instrfind('Type', 'visa-serial', 'RsrcName', InstrStr, 'Tag', '');

% Create the VISA-Serial object if it does not exist
% otherwise use the object that was found.
if isempty(variLCObj)
    variLCObj = visa('NI', 'ASRL11::INSTR');
else
    fclose(variLCObj);
end

% Connect to instrument object, obj1.
fopen(variLCObj);
% 
% Configure instrument object, obj1.
set(variLCObj, 'Timeout', 1.0);
set(variLCObj, 'Terminator', {'CR','CR'});

% When setting the LC - read echo.
% When querying the LC - read echo and then response.

% set 'brief mode' This omits the command character from the response.
fprintf(variLCObj, 'B 1'); % command version
echo= fscanf(variLCObj, '%s\n'); % read echo response
% response=fscanf(variLCObj, '%s\n');
% display([echo, response]);

% Get version and display.
fprintf(variLCObj, 'V?'); % command version
echo = fscanf(variLCObj, '%s\n'); % read echo response
response=fscanf(variLCObj, '%s\n');
display(['Communicating with controller version:' response]);

% fprintf(variLCObj, 'L?'); % command version
% echo = fscanf(variLCObj, '%s\n'); % read echo response
% response=fscanf(variLCObj, '%s\n');
% display([echo, response]);

LCsettings=loadjson(palleteFile);

% Set wavelegth.
fprintf(variLCObj, ['W ' int2str(LCsettings.wavelength)]); % command version
echo = fscanf(variLCObj, '%s\n'); % read echo response
% response=fscanf(variLCObj, '%s\n');
% display([echo, response]);


% Set all five palletes.
DefinePallete(0,LCsettings.lcPalElsD(1,:));
DefinePallete(1,LCsettings.lcPalElsD(2,:));
DefinePallete(2,LCsettings.lcPalElsD(3,:));
DefinePallete(3,LCsettings.lcPalElsD(4,:));

% Read the calibration and display.
fprintf(variLCObj, 'W?'); % command version
echo = fscanf(variLCObj, '%s\n'); % read echo response
response=fscanf(variLCObj, '%s\n');
display(['Set wavelength:' response]);
activatePallete(0);
activatePallete(1);
activatePallete(2);
activatePallete(3);
activatePallete(0);

% VariLC Commands Basic - Refer to blue folder next to Amit's desk (3rd shelf) for expanded list
%
% V ? - query voltage
% B 1 - sets to short mode reply (no command term in response - easier
% parsing in some cases). B1 removes the preceding command letter.
% L 0.25 0.5
% L ? - query retardance
% D 1 - defines Palette 1 based on current retardance values
% P 1 - sets LC to Pallette 1
% C 1 - clears Palette 1
% E 3 - exercise for 3 cycles
% W 550 - sets to 550nm wavelength calibration
% W ? - query wavelength
% R 1 - resets LC error light


    function DefinePallete(pno,values)
        
        % Clear any existing error.
        fprintf(variLCObj,'R 1');
        echo = fscanf(variLCObj, '%s\n');
        
        % First current retardance to the specified values.
        fprintf(variLCObj,['L ' num2str(values(1)) ' ' num2str(values(2)) ]);
        echo = fscanf(variLCObj, '%s\n'); % read echo 
        
        % Check if these settings have caused an error.
        fprintf(variLCObj,'R?');
        echo = fscanf(variLCObj, '%s\n');
        response= fscanf(variLCObj, '%s\n');
        if(str2double(response))
            error(['Setting LCs to [' num2str(values(1)) ' ' num2str(values(2)) '] causes an error#' response '.']);
        else %If not define the pallete.
            fprintf(variLCObj, ['D ' int2str(pno)]); 
            echo = fscanf(variLCObj, '%s\n'); % read echo 
        end  
    end

    function LC=activatePallete(pno)
         fprintf(variLCObj, ['P ' int2str(pno)]); 
         echo = fscanf(variLCObj, '%s\n'); % read echo 
         %response= fscanf(variLCObj, '%s\n'); 
         %display([echo, response]);
         pause(0.05); %Pause for 50ms.
         fprintf(variLCObj, ['L?']); 
         echo = fscanf(variLCObj, '%s\n'); % read echo 
         response= fscanf(variLCObj, '%s\n'); 
         LC=cell2mat(textscan(response,'%f %f'));
         display(['Pallete#' int2str(pno)  ' ['  num2str(LC(1)) ' ' num2str(LC(2)) '].' ]);    
    end

end

