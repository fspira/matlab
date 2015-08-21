function [ status ] = getLC( variLC,command)

switch(class(variLC))
    case 'tcpip' %If communicating with micro-manager over tcpip.
        %char(10) is LF character expected as delimiter by OPS.
      
        fwrite(variLC,[command char(10)]);
      
        pause(0.2); % Wait 200ms for LC to settle, tcpip buffer to be ready.
        msg=fread(variLC,variLC.BytesAvailable);
        msg=char(msg)';

        if(isempty(strfind(msg,'SUCCESS')))
            status=0;
        else 
            status=1;
        end
        
    case 'visa' %If communicating directly with the LC over serial port.

        fprintf(variLC,command);
        echo= fscanf(variLC, '%s\n'); % read echo 
        status=fscanf(variLC,'%s\n'); 
end

end

