function [idmatchf,idmatch, idnotmttomanym, idnomatchm, dmin, idrepnm, XMEANDIS, COMAssm1, COMAssm2]=mathingtwosteps(assm1,assm2,TOLDistNei, MisOrTol, MeanNeighNom,TOLDistGrainXC )
%% assm1 and assm2 are the outputs from two loading steps
%% mathingtwosteps cross correlates assm2 with assm1
%% idmatchf(i)=j:: ith element of assm2 matches with jth element of assm1 
%% idmatch(i,j):: ith element of assm2 has j=1 number of matched grains, j=2:j1+1 are IDs of matched grains. 
%% idnotmttomanym(i)=j:: jth elemnt of assm2 has many neighbours but non of the neighbours matched! 
%% idnomatchm(i)=j ; jth element of assm2 did not match with any of the grains of assm1 
%% dmin is the minimum distance for ith grain in assm2
%% XMEANDIS is rigid body movement 
%% Paramaters: 
%% TOLDistNei=25; % distance tolerance for finding neighbours! average number of neighbour per grain is 6!
%% MisOrTol=5;
%% MeanNeighNom=6; %mean number of neighbours per grains
%% TOLDistGrainXC=25; % distance tolerance for finding grains accross two layers

%% 1: #  grainno; 2: mean_IA; 3: grainvolume;
%% 4:x;5: y; 6: z; 7:rodx; 8: rody; 9:rodz; 10: U11; 11:U12; 12: U13;
%% 13:U21; 14: U22; 15: U23; 16: U31; 17: U32; 18: U33; 19: eps11;
%% 20: eps22;  21: eps33; 22: ;eps23; 23: eps13; 24: eps12; 25: eps11_s;26: eps22_s
%% 27:eps33_s; 28: eps23_s ; 29:  eps13_s; 30: eps12_s; 
%% 31: sig11; 32: sig22; 33: sig33; 34: sig23; 35: sig13; 36: sig12; 37:sig11_s;
%% 38: sig22_s; 39: sig33_s; 40: sig23_s; 41:sig13_s; 42: sig12_s; 43: eta; 44: eta; 45: previous ID; 46: previous layer. 
%%
%%


%% stratogy: (1) cross corrolate assm1 grains with assem2 using just oerinetations 
%% then look for neighbouring grains. The grains that has more than 10 common neighbouring grains will be recorded. 
%% calculating COM 
COMAssm1=[sum(assm1(:,4).*assm1(:,3))/sum(assm1(:,3)); sum(assm1(:,5).*assm1(:,3))/sum(assm1(:,3));sum(assm1(:,6).*assm1(:,3))/sum(assm1(:,3))]
COMAssm2=[sum(assm2(:,4).*assm2(:,3))/sum(assm2(:,3)); sum(assm2(:,5).*assm2(:,3))/sum(assm2(:,3));sum(assm2(:,6).*assm2(:,3))/sum(assm2(:,3))]
%%
neigha=zeros(1,50); % maximum 50 neighbours assumed!
for i=1:size(assm1,1)
    coma=assm1(i,4:6);
    idn=0;
    for j=1:size(assm1,1)
        if (j ~= i)
            dist=((coma(1)-assm1(j,4))^2+(coma(2)-assm1(j,5))^2+(coma(3)-assm1(j,6))^2)^0.5;
            if (dist<=TOLDistNei)
                idn=idn+1;
                neigha(i,1)=idn;
                neigha(i,1+idn)=j;
            end
        end
    end
end

neighb=zeros(1,50); % maximum 50 neighbours assumed!
for i=1:size(assm2,1)
    comb=assm2(i,4:6);
    idn=0;
    for j=1:size(assm2,1)
        if (j ~= i)
            dist=((comb(1)-assm2(j,4))^2+(comb(2)-assm2(j,5))^2+(comb(3)-assm2(j,6))^2)^0.5;
            if (dist<=TOLDistNei)
                idn=idn+1;
                neighb(i,1)=idn;
                neighb(i,1+idn)=j;
            end
        end
    end
end

disp('neighbours are identified')
disp('fixing rigid body movement')

idm=zeros(size(assm1,1),size(assm2,1));
for i=1:size(assm1,1)
    if (neigha(i,1)>MeanNeighNom)
        if (mod (i,50)==1)
            disp (['working on grain # ', num2str(i)]);
        end
        eula=assm1(i,7:9);
        coma=assm1(i,4:6);
        for j=1:size(assm2,1)
            idtst=samequat(coma,assm2(j,4:6));
            if (idtst==1)
                eulb=assm2(j,7:9);
                misori=totalmisor(eula,eulb);
                if (misori(1)<MisOrTol & misori(2)< MisOrTol)
                    idm(i,j)=0; %number of matched grains.
                    comb=assm2(j,4:6);
                    for k=1:neigha(i,1)
                        idng=neigha(i,1+k);
                        eulna=assm1(idng,7:9);
                        dna=[assm1(idng,4)-coma(1);assm1(idng,5)-coma(2);assm1(idng,6)-coma(3)]';
                        for kk=1:neighb(j,1)
                            idngb=neighb(j,1+kk);
                            eulnb=assm2(idngb,7:9);
                            dnb=[assm2(idngb,4)-comb(1);assm2(idngb,5)-comb(2);assm2(idngb,6)-comb(3)]';
                            idtst1=samequat(dna,dnb);
                            if (idtst1 ==1)
                                misornab=totalmisor(eulna,eulnb);
                                if (misornab(1)<MisOrTol & misornab(2)<MisOrTol)
                                    idm(i,j)=idm(i,j)+1;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

idone=0;
idcen=0;
for i=1:size(assm1,1)
    idf=find(idm(i,:)>0);
    if (size(idf)==1)
        idone=idone+1;
        if (idm(i,idf)>= MeanNeighNom)
            idcen=idcen+1;
            xcdis(idcen,1)=assm1(i,4)-assm2(idf,4);
            xcdis(idcen,2)=assm1(i,5)-assm2(idf,5);
            xcdis(idcen,3)=assm1(i,6)-assm2(idf,6);
            xcdis(idcen,4)=assm1(i,3);
            xcdis(idcen,5)=assm2(idf,3);
        end
    end        
end

XMEANDIS=[sum(xcdis(:,1).*xcdis(:,5))/sum(xcdis(:,5));
          sum(xcdis(:,2).*xcdis(:,5))/sum(xcdis(:,5));
          sum(xcdis(:,3).*xcdis(:,5))/sum(xcdis(:,5))]

assm2rxyz=assm2;  %assm2rxyz is assam2 with COM moved to the COM of assm1 
for i=1:size(assm2,1)
    assm2rxyz(i,4)=assm2rxyz(i,4)+XMEANDIS(1);
    assm2rxyz(i,5)=assm2rxyz(i,5)+XMEANDIS(2);
    assm2rxyz(i,6)=assm2rxyz(i,6)+XMEANDIS(3);
end

%%
%% Rigid body movement is now fixed
%% Cross correlating two layers
%%
%% TOLDistGrainXC
disp('Rigid body movements are fixed')
disp('Cross corrolating two layers')

idmatch=zeros(size(assm2rxyz,1),10); %maxmimum 10 match allowed!
idmatchf=zeros(size(assm2rxyz,1),1); %final matrix that has the matched IDs. 
idnotmttomany=0; %IDs of the grains with too many matches 
idnomatch=0;     %IDs of the grains with no match. 
dmin=nan(size(assm2rxyz,1),1);
idnotmttomanym=NaN;
idnomatchm=NaN;
idrepnm=NaN;

for i=1:size(assm2rxyz,1)
    
    if (mod(i,50) ==1)
        disp(['working on grain # ', num2str(i)])
    end
    eula=assm2rxyz(i,7:9);
    coma=assm2rxyz(i,4:6);
    idxcm=0;
    distab=nan(size(assm1,1),1);
    for j=1:size(assm1,1)
        eulb=assm1(j,7:9);
        misorab=totalmisor(eula,eulb);
        if(misorab(1)<MisOrTol & misorab(2)< MisOrTol)
            distab(j)=((coma(1)-assm1(j,4))^2+(coma(2)-assm1(j,5))^2+(coma(3)-assm1(j,6))^2)^0.5;
            if (distab(j)<= TOLDistGrainXC)
               idxcm=idxcm+1;
               idmatch(i,1)=idxcm; %number of matches
               idmatch(i,1+idxcm)=j;
            end
        end
    end
    dmin(i)=min(distab);
    if (idmatch(i,1) == 1 )
       idmatchf(i)= idmatch(i,2);
    elseif (idmatch(i,1) > 1)  %% if distance and orientation matcheds, then neighbours will be examined. 
        clear nmgn misorab1
        nmgn(1)=0;
        for k=1:idmatch(i,1)
            
            ngam1=idmatch(i,1+k); %mastched grain 
            nmgn(k+1)=0;
            for j=1:neighb(i,1)
                idng=neighb(i,1+j);
                eulnb=assm2rxyz(idng,7:9);
                comnb=assm2rxyz(idng,4:6);
                
                for jj=1:neigha(ngam1,1)
                    idng1=neigha(ngam1,1+jj);
                    eulna=assm1(idng1,7:9);
                    comna=assm1(idng1,4:6);
                    misornab=totalmisor(eulna,eulnb);
                    if (misornab(1)<MisOrTol & misornab(2)<MisOrTol)
                        disnab=((comna(1)-comnb(1))^2+(comna(2)-comnb(2))^2+(comna(3)-comnb(3))^2)^0.5;
                        if (disnab < TOLDistGrainXC)
                            nmgn(k+1)=nmgn(k+1)+1;
                        end                        
                    end
                end
            end
        end
        maxmatch=max(nmgn);
        if (maxmatch > 0)
            idfm=find(nmgn==maxmatch);
            if (size(idfm)==1)
                idmatchf(i)=idmatch(i,idfm);
            else  %%%if there are more than one match with common neighbours
                for jj1=1:max(size(idfm))
                    misorab1(jj1,:)=totalmisor(eula,assm1(idmatch(i,idfm(jj1)),7:9));
                end
                jj1=find(misorab1(:,1)==min(misorab1(:,1)));
                idmatchf(i)=idmatch(i,idfm(jj1));
            end
        else
            idnotmttomany=idnotmttomany+1;
            idnotmttomanym(idnotmttomany)=i;
        end
    elseif (idmatch(i,1)==0)
        idnomatch=idnomatch+1;
        idnomatchm(idnomatch)=i;
    end
end




%%
%%
