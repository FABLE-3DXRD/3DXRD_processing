%% This is to import gff files and convert rod>euler
%% 1: #  grainno; 2: mean_IA; 3: grainvolume;
%% 4:x;5: y; 6: z; 7:rodx; 8: rody; 9:rodz; 10: U11; 11:U12; 12: U13;
%% 13:U21; 14: U22; 15: U23; 16: U31; 17: U32; 18: U33; 19: eps11;
%% 20: eps22;  21: eps33; 22: ;eps23; 23: eps13; 24: eps12; 25: eps11_s;26: eps22_s
%% 27:eps33_s; 28: eps23_s ; 29:  eps13_s; 30: eps12_s; 
%% 31: sig11; 32: sig22; 33: sig33; 34: sig23; 35: sig13; 36: sig12; 37:sig11_s;
%% 38: sig22_s; 39: sig33_s; 40: sig23_s; 41:sig13_s; 42: sig12_s; 43: eta; 44: eta
%%
%%

clear;
close all
fdir='Y:\HA\ESRFDATA\ESRF\id11\Processing_XRD_ZrReal\MeasuredPatterns\XRD_Zr_PlasticityOnset_edna\merge_XRD_Zr_PlasticityOnset_3489_4360__post\merge_XRD_Zr_PlasticityOnset_3489_4360_FitAllB\';
filename = strcat(fdir,'merge_XRD_Zr_PlasticityOnset_3489_4360_FitAllB_final.gff');
filenamerr = strcat(fdir,'merge_XRD_Zr_PlasticityOnset_3489_4360_FitAllB_final_errors.gff');
delimiter = ' ';
startRow = 2;
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);

fileIDerr = fopen(filenamerr,'r');
dataArrayerr = textscan(fileIDerr, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);
fclose(fileIDerr)

%% Create output variable
datafinal = [dataArray{1:end-1}];
dataerrors = [dataArrayerr{1:end-1}];

datafinal(:,45)='';
dataerrors(:,45)='';

datanormal=datafinal;
datanormalerr=dataerrors;
%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans filenamerr fileIDerr dataArrayerr;


%% postprocessing on the data

%%
nog=size(datanormal,1);
id1=0;
for i=1:nog
    test=0;
    %%%% Find the related error
    [idfsbt,tstfsbt]=find(datanormalerr(:,1)==datanormal(i,1));
    
    %%%% stress check
    for j=37:42
        if (abs(datanormal(i,j))>3000000 | abs(datanormalerr(idfsbt,j)>100000))
            test=1;
        end 
    end
    
    
    if (abs(1000*datanormal(i,6)) > 1800)  %%% position check
        test=1;
    end 
    
    if (test == 0)
        [idf,tst]=find(datanormalerr(:,1)==datanormal(i,1));
        if (isnan(idf)==0)
            %if (sum(isnan(datanormalerr(idf,:)))<1)
                id1=id1+1;
                datanormal2(id1,:)=datanormal(i,:);
                [idf,tst]=find(datanormalerr(:,1)==datanormal2(id1,1));
                datanormal2err(id1,:)=datanormalerr(idf,:);
            %end
        end
        
    end
end


nogf=id1;
xdir=datanormal2(:,4)*1000;
ydir=datanormal2(:,5)*1000;
zdir=datanormal2(:,6)*1000;
strsval=datanormal2(:,37:42);

%%
%%

% mmin=min(strsval(:,3));mmax=max(strsval(:,3));
% figure(1)
% scatter3(xdir,ydir,zdir,200,strsval(:,3),'filled')    % draw the scatter plot
% ax = gca;
% view(-31,14)
% xlabel('X- beam direction')
% ylabel('Y- Transverse')
% zlabel('Z- Loading Direction')
% title ('S_{33} Pre-Load Step')
% colormap('jet'); % get current colormap
% 
% gvolume=datanormal2(:,3);
% figure(3)
% scatter3(xdir,ydir,zdir,100*gvolume.^(1/3),strsval(:,3),'filled')    % draw the scatter plot
% ax = gca;
% view(-31,14)
% xlabel('X- beam direction')
% ylabel('Y- Transverse')
% zlabel('Z- Loading Direction')
% title ('S_{33} Pre-Load Step, Size proportional to Volume')
% colormap('jet');   
% caxis([mmin mmax])
% colorbar
% 
% %%
% mmins=-350; mmaxs=350;
% idgrains=1:nogf;
% fig100 = figure(100);
% % fname=strcat(text,'datanormal2.dat');
% % save(fname,'datanormal2','-ascii')
% axes100 = axes('Parent',fig100,'FontSize',12);
% errorbar(idgrains,datanormal2(:,39),datanormal2err(:,39),'bo' )
% % plot(idgrains,datanormal(:,39),'bo' )
% xlabel('grain ID','fontsize',14)
% ylabel('S_{33}: Loaidng direction','fontsize',14)
% grid on
% title ('CPZr Preload','fontsize',14 )  
% axis([-5 nogf+10 mmins mmaxs])
% 
% 
% fig101 = figure(101);
% axes101 = axes('Parent',fig101,'FontSize',12);
% errorbar(idgrains,datanormal2(:,37),datanormal2err(:,37),'bo' )
% xlabel('grain ID','fontsize',14)
% ylabel('S_{11}: Beam Direction','fontsize',14)
% grid on
% title ('CPZr Preload','fontsize',14 )  
% axis([-5 nogf+10 mmins mmaxs])
%%
%%
%%
fname=strcat('datanormal2.dat');
save(fname,'datanormal2','-ascii')
system(['Rod2Euler.exe < ' fname]); %run executable with c
%% Euler angles: 1:grain num; 2: grain num; 3:5 Eulers; 6: tetavsload; 7: volume
graineuler=importdata ('EulersGL.dat');
datanormal3=zeros(size(datanormal2,1),45);
for i=1:size(datanormal2,1)
    [idf,tst]=find(graineuler(:,1)==datanormal2(i,1));
    
    if (isnan(idf)==0)    
        for j=7:9
            datanormal2(i,j)=graineuler(idf,j-4);
            datanormal3(i,1:44)=datanormal2(i,1:44);
            datanormal3(i,45)=graineuler(idf,6);
        end
    else
        disp(['error',i])
    end
end
%%
%%
graineulerLG=importdata ('grain_euler1Vol300PosLoc2Glob.dat');

fname=strcat(fdir,'datanormal2.dat');
save(fname,'datanormal2','-ascii')
fnamerr=strcat(fdir,'datanormal2err.dat');
save(fnamerr,'datanormal2err','-ascii')
fnameul=strcat(fdir,'EulersG2L.dat');
save(fnameul,'graineuler','-ascii')
fnameulg=strcat(fdir,'EulersL2G.dat');
save(fnameulg,'graineulerLG','-ascii')

delete('grain_euler1Vol300PosLoc2Glob.dat', 'grain_euler1Vol300Pos.dat', 'EulersGL.dat','datanormal2.dat')
    
lz=[0 0 1]';
ly=[0 1 0]';
lx=[1 0 0]';
vbasal=[0 0 1]';

xdirn=datanormal3(:,4)*1000;
ydirn=datanormal3(:,5)*1000;
zdirn=datanormal3(:,6)*1000;
figure(4)
scatter3(xdirn,ydirn,zdirn,150*datanormal3(:,3).^(1/3),datanormal3(:,45),'filled')    % draw the scatter plot
ax = gca;
view(-31,14)
xlabel('X- beam direction')
ylabel('Y- Transverse')
zlabel('Z- Loading Direction')
title ('c-axis-load misorientation')
colormap('jet');    
caxis([0 90])
colorbar

totalv=sum(datanormal2(:,3));
for i=1:6
    strs(i)=sum(datanormal2(:,i+36).*datanormal2(:,3))/totalv;
    strserr(i)=sum(datanormal2err(:,i+36).*datanormal2(:,3))/totalv; 
    
    strn(i)=sum(datanormal2(:,i+24).*datanormal2(:,3))/totalv;
    strng(i)=sum(datanormal2(:,i+18).*datanormal2(:,3))/totalv;
    
end
strs
strserr

strn
