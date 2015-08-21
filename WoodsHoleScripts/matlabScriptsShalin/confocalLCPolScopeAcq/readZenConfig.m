function [ config ] = readZenConfig( tcpipZen )
    msgdelim=';;';
    successdelim=':';
    timeout=1;

    fwrite(tcpipZen,['-get_config_name' msgdelim]);

    pause(timeout);

    % Read the mesg and covert in characters.
    Mesg=fread(tcpipZen,tcpipZen.BytesAvailable);
    Mesg=char(Mesg)';

    % Separate the message based on success delimiter
    MesgBreak=textscan(Mesg,'%s','delimiter',successdelim); 

    % Break further based on delimiter for fields ';'
    MesgBreak2=textscan(MesgBreak{1}{2},'%s','delimiter',';');

    % Extract the config name.
    config= MesgBreak2{1}{1};
end

