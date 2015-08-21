function [ status ] = waitforZen( tcpipZen,timeout )
% Wait for response from zen until SUCCESS, REJECT or TIMEOUT.
% SUCCESS: status=1, REJECT: status=0, TIMEOUT: status=-1.

BinMesg=[];
a=tic;
while(1) % Keep reading.
    
    while(~tcpipZen.BytesAvailable) % Wait until bytes available to read, but only up to time out.
        t=toc(a);
        if(t>timeout)
            status=-1;
            return;
        end
    end
    
    BinMesg=[BinMesg; fread(tcpipZen,tcpipZen.BytesAvailable)];  
    Mesg=char(BinMesg)';
    if(~isempty(strfind(Mesg,'SUCCESS')))
        status=1;
        return;
    end
    
    if(~isempty(strfind(Mesg,'REJECTED')))
        status=0;
        return;
    end
    
end

end

