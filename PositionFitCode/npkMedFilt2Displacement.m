function [ DSCx_filtered, DSCy_filtered ] = npkMedFilt2Displacement( DSCx, DSCy, ampthresh, anglethresh )
% Nathanael Kazmierczak, 02042020
%
% Because the displacement field is so dependent on x and y components
% being synchronized, this function does a median filter simultaneously in
% both x and y directions.
%
% Basic xy filtering works for the angle but greatly screws up the
% amplitude. I'm not completely sure why that is, but the half-hexagon
% might be implicated.

DSCamp = (DSCx.^2 + DSCy.^2).^0.5;
DSCangle = atan2(DSCy,DSCx);

ampfilt = medfilt2(DSCamp,[3 3],'symmetric');
anglefilt = medfilt2(DSCangle,[3 3],'symmetric');
sub_amp = abs(DSCamp - ampfilt) > ampthresh;
sub_angle = abs(DSCangle - anglefilt) > anglethresh;  % This needs to be more clever

% For the angle, because 0 and pi are equivalent, we need to ensure the
% wraparound is okay
DSCangle_phased = DSCangle;
DSCangle_phased(DSCangle_phased > pi/2) = DSCangle_phased(DSCangle_phased > pi/2) - pi;
anglephasedfilt = medfilt2(DSCangle_phased,[3 3],'symmetric');
% % % undo the phasing [no]
% anglephasedfilt(anglephasedfilt < 0) = anglephasedfilt(anglephasedfilt < 0) + pi;
sub_angle_phased = abs(DSCangle - anglephasedfilt) > anglethresh;
sub_angle_comp = sub_angle & sub_angle_phased;  
% logic here is that it has to be far away from the median under both the
% wraparound and the non wraparound, or else that means one interpretation
% of the wraparound made it work.

allsub = sub_amp | sub_angle_comp;

DSCamp_filtered = DSCamp;
DSCangle_filtered = DSCangle;
DSCamp_filtered(allsub) = ampfilt(allsub);
DSCangle_filtered(allsub) = anglefilt(allsub);
 
% DSCx_filtered = xfilt;
% DSCy_filtered = yfilt;

DSCx_filtered = DSCamp_filtered.*cos(DSCangle_filtered);
DSCy_filtered = DSCamp_filtered.*sin(DSCangle_filtered);
% DSCx_filtered(allsub) = xfilt(allsub);
% DSCy_filtered(allsub) = yfilt(allsub);



% % % % % xfilt = medfilt2(DSCx,[3 3],'symmetric');
% % % % % yfilt = medfilt2(DSCy,[3 3],'symmetric');
% % % % % sub_x = abs(DSCx - xfilt) > thresh;
% % % % % sub_y = abs(DSCy - yfilt) > thresh;
% % % % % allsub = sub_x | sub_y;
% % % % % 
% % % % % DSCamp = (DSCx.^2 + DSCy.^2).^0.5;
% % % % % DSCangle = atan2(DSCy,DSCx);
% % % % % 
% % % % % DSCamp_filtered = DSCamp;
% % % % % DSCangle_filtered = DSCangle;
% % % % % DSCamp_filtered(allsub) = 
% % % % %  
% % % % % % DSCx_filtered = xfilt;
% % % % % % DSCy_filtered = yfilt;
% % % % % 
% % % % % DSCx_filtered = DSCx;
% % % % % DSCy_filtered = DSCy;
% % % % % % DSCx_filtered(allsub) = xfilt(allsub);
% % % % % % DSCy_filtered(allsub) = yfilt(allsub);

end

