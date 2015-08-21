function [ dz ] = readZendz( tcpipZen )
msgdelim=';;';
successdelim=':';
timeout=1;

fwrite(tcpipZen,['-get_zstack_dz' msgdelim]);

pause(timeout);

dzMesg=fread(tcpipZen,tcpipZen.BytesAvailable);
dzMesg=char(dzMesg)';
dzMesgBreak=textscan(dzMesg,'%s','delimiter',successdelim); 
dz=textscan(dzMesgBreak{1}{2},'%f');
dz=dz{:};
end

