%% This subbroutine determines common grains across layers. 
%% compared datam1 & datam2
function [idpair, idp1f]=commongrainacrosslayers(datam1,datam2, Layer1, Layer2, TolDist, TolTotalMis, TolBasMis)
idp1=0;
for i=1:size(datam1,1)
    cmc=datam1(i,4:6);
    euler1=datam1(i,7:9);
    for j=1:size(datam2,1)
        dist1=((cmc(1)-datam2(j,4))^2+...
               (cmc(2)-datam2(j,5))^2+...
               (cmc(3)-datam2(j,6))^2)^0.5;
        
        if dist1 < TolDist
            euler2=datam2(j,7:9);
            misorab=totalmisor(euler1,euler2);
            if (misorab(1)<TolTotalMis & misorab(2) < TolBasMis)
                idp1=idp1+1;
                
                idpair(idp1,1)=Layer1; % 1st layer 
                idpair(idp1,2)=Layer2; % 2nd layer
                
                idpair(idp1,3)=i; 
                idpair(idp1,4)=j;
                idpair(idp1,5)=dist1;
                idpair(idp1,6:7)=misorab(1:2);  
                
            end            
        end          
    end    
end
idp1f=idp1;
if (idp1f==0)
    idpair=NaN;
end