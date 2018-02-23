function [basaltet]=eulerori(eulert)

rotLG=eulerLG(eulert(1),eulert(2),eulert(3));

lz=[0 0 1]';
ly=[0 1 0]';
lx=[1 0 0]';

vbasal=[0 0 1]';

basaltet(1)=acosd(dot(rotLG*vbasal,lx)); 
if (basaltet(1)>90)
    basaltet(1)=180-basaltet(1);
end % basal misori with x

basaltet(2)=acosd(dot(rotLG*vbasal,ly)); 
if (basaltet(2)>90); 
    basaltet(2)=180-basaltet(2); 
end; %basal misori with y

basaltet(3)=acosd(dot(rotLG*vbasal,lz)); 
if (basaltet(3)>90); 
    basaltet(3)=180-basaltet(3); 
end; % basal misor with 