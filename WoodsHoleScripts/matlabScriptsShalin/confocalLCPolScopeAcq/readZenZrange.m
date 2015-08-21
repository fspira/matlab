function [ zStart,zEnd ] = readZenZrange( tcpipZen )
msgdelim=';;';
successdelim=':';
timeout=1;
fwrite(tcpipZen,['-get_zstack_range' msgdelim]);
pause(timeout); %100ms delay for microscope for response to arrive. Do not use waitforZen, because it discards the message.
zMesg=fread(tcpipZen,tcpipZen.BytesAvailable);
zMesg=char(zMesg)';
zMesgBreak=textscan(zMesg,'%s','delimiter',successdelim); 
zRange=textscan(zMesgBreak{1}{2},'%f');
zRange=zRange{:};
zStart=zRange(1); zEnd=zRange(2); 

end

