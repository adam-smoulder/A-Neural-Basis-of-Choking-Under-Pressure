function [dateTimeString] = grabDateTimeString()
%%% Very quick function just to get a simple datetime string to be used for
%%% appending to filenames
%%%
%%% Adam Smoulder, 9/10/18 (edit 5/5/22)
%%%

% thedate = date;
% time = clock;
% 
% % ensure each of the clock values has 2 digits
% timeString = {num2str(time(4)), num2str(time(5)), num2str(round(time(6)))}; % hour/minute/second
% for i = 1:length(timeString)
%     if length(timeString{i})<2
%         timeString{i} = ['0' timeString{i}];
%     end
% end; clear i
% 
% dateTimeString = [thedate timeString{:}];


% Easier way, though had to google forever to find it...
dateTimeString = datestr(now,'yyyymmdd_HHMMSS');

end

