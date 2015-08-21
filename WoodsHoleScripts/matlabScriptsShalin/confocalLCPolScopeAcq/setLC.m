function [ status ] = setLC( variLC,state)

switch(class(variLC))
    case 'tcpip' %If communicating with micro-manager over tcpip.
        %char(10) is LF character expected as delimiter by OPS.
        switch(state)
            case 'I0'
                fwrite(variLC,['0' char(10)]);
            case 'I135'
                fwrite(variLC,['1' char(10)]);
            case 'I90'
                fwrite(variLC,['2' char(10)]);
            case 'I45'
                fwrite(variLC,['3' char(10)]);
            case 'I0bleach'
                fwrite(variLC,['4' char(10)]);
            otherwise
                fwrite(variLC,[state char(10)]);
        end

        pause(0.1); % Wait 200ms for LC to settle, tcpip buffer to be ready.
        msg=fread(variLC,variLC.BytesAvailable);
        msg=char(msg)';

        if(isempty(strfind(msg,'SUCCESS')))
            status=0;
        else 
            status=1;
        end
        
    case 'visa' %If communicating directly with the LC over serial port.
        switch(state)
            case 'I0'
                fprintf(variLC,'P 0');
                echo= fscanf(variLC, '%s\n'); % read echo 
            case 'I135'
                fprintf(variLC,'P 1');
                echo= fscanf(variLC, '%s\n'); % read echo 
            case 'I90'
                fprintf(variLC,'P 2');
                 echo= fscanf(variLC, '%s\n'); % read echo 
            case 'I45'
                fprintf(variLC,'P 3');
                 echo= fscanf(variLC, '%s\n'); % read echo 
            case 'I0bleach'
                fprintf(variLC,'P 0');
                 echo= fscanf(variLC, '%s\n'); % read echo 
            otherwise
                fprintf(variLC,state);
                 echo= fscanf(variLC, '%s\n'); % read echo 
        end
        
        %display([ state ' ' echo]);
        % Let the LC settle.
        pause(0.1); %100ms.
        
        
end

end

