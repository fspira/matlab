% Find a VISA-Serial object.
obj1 = instrfind('Type', 'visa-serial', 'RsrcName', 'ASRL11::INSTR', 'Tag', '');

% Create the VISA-Serial object if it does not exist
% otherwise use the object that was found.
if isempty(obj1)
    obj1 = visa('NI', 'ASRL11::INSTR');
else
    fclose(obj1);
    obj1 = obj1(1);
end

% Connect to instrument object, obj1.
fopen(obj1);
% 
% Configure instrument object, obj1.
set(obj1, 'Timeout', 1.0);
set(obj1, 'Terminator', {'CR','CR'});

% Communicating with instrument object, obj1.
% When reading - read echo and then response.
fprintf(obj1, 'V?'); % command version
data1 = fscanf(obj1, '%s\n'); % read echo response
data1 = fscanf(obj1, '%s\n'); % read answer response


display(['VariLC response: ', data1]);

% Flush the data in the input buffer.
flushinput(obj1);


% When writing do the single read.
fprintf(obj1, 'B0'); % command version
data1 = fscanf(obj1, '%s\n'); % read echo response
%data1 = fscanf(obj1, '%s\n'); % read answer response

fprintf(obj1, 'L?'); % command version
data1 = fscanf(obj1, '%s\n'); % read echo response
data1 = fscanf(obj1, '%s\n'); % read answer response

display(['VariLC response: ', data1]);

% Disconnect from instrument object, obj1.
%fclose(obj1);

% Clear the visa object
%clear obj1

% VariLC Commands Basic - Refer to blue folder next to my desk (3rd shelf) for expanded list
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

