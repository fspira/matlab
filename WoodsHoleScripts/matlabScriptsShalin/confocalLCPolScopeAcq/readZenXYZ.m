function [ x,y,z ] = readZenXYZ( tcpipZen )
msgdelim=';;';
successdelim=':';
timeout=1;
fwrite(tcpipZen,['-get_stage_xyz' msgdelim]);

pause(timeout);
xyzMesg=fread(tcpipZen,tcpipZen.BytesAvailable);
xyzMesg=char(xyzMesg)';
xyzMesgBreak=textscan(xyzMesg,'%s','delimiter',successdelim); 
xyz=textscan(xyzMesgBreak{1}{2},'%f');
xyz=xyz{:};
x=xyz(1); y=xyz(2); z=xyz(3);

end

