function rotate=eulerLG(Phi1,Phi,Phi2)

%%%%%%%%%Modfied fpor Bruker!
%% euler angles should be from sample to crystal but at the end because of applying a transpose
%% transforming the local coordinate to the global
    phi1=Phi1*pi/180.d0;
	phi= Phi*pi/180.d0;
	phi2=Phi2*pi/180.d0;
	rotate(1,1)=cos(phi1)*cos(phi2)-sin(phi1)*sin(phi2)*cos(phi);
	rotate(1,2)=sin(phi1)*cos(phi2)+cos(phi1)*sin(phi2)*cos(phi);
	rotate(1,3)=sin(phi2)*sin(phi);
	rotate(2,1)=-cos(phi1)*sin(phi2)-sin(phi1)*cos(phi2)*cos(phi);
	rotate(2,2)=-sin(phi1)*sin(phi2)+cos(phi1)*cos(phi2)*cos(phi);
	rotate(2,3)=cos(phi2)*sin(phi);
	rotate(3,1)=sin(phi1)*sin(phi);
	rotate(3,2)=-cos(phi1)*sin(phi);
	rotate(3,3)=cos(phi);
    
    rotate=rotate';

