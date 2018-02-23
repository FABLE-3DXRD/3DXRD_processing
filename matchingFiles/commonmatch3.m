
function [idcom,idcomm]=commonmatch3(idmatcha2a1,idmatcha3a1,idmatcha3a2,idrepgm21,idrepgm31,idrepgm32)

idcom=0; %% number of common IDs
idcomm=nan;
for i=1:size(idmatcha3a2,1)
    if (isempty(find(idrepgm32(:) ==i)) ==1 )
        if (isempty(find(idrepgm31(:) ==i)) ==1 )
            id2=idmatcha3a2(i);
            id1=idmatcha3a1(i);       
            if (isempty (id2) ==0 & isempty (id1) ==0 & id2 ~=0 & id1 ~=0) 
                if isempty (find(idrepgm21(:) == id2))
                    if (idmatcha2a1(id2) == id1) %% same grains are assigned!!!
                        idcom=idcom+1;
                        idcomm(idcom,1)=id1;  % ID1
                        idcomm(idcom,2)=id2;  % ID2
                        idcomm(idcom,3)=i;  % ID3
                    end
                end
            end
        end
    end
end