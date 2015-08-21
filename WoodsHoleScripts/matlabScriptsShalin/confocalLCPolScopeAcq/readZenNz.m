function [ nz ] = readZenNz( tcpipZen )
msgdelim=';;';
successdelim=':';
timeout=1;

fwrite(tcpipZen,['-get_zstack_nz' msgdelim]);

pause(timeout);

Mesg=fread(tcpipZen,tcpipZen.BytesAvailable);
Mesg=char(Mesg)';
MesgBreak=textscan(Mesg,'%s','delimiter',successdelim); 
nz=textscan(MesgBreak{1}{2},'%f');
nz=nz{:};
end

