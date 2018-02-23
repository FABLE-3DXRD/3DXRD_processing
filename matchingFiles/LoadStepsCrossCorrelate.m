%% This program is the main program for cross corrolating grains from two loading steps 
%%
%% 1: #  grainno; 2: mean_IA; 3: grainvolume;
%% 4:x;5: y; 6: z; 7:rodx; 8: rody; 9:rodz; 10: U11; 11:U12; 12: U13;
%% 13:U21; 14: U22; 15: U23; 16: U31; 17: U32; 18: U33; 19: eps11;
%% 20: eps22;  21: eps33; 22: ;eps23; 23: eps13; 24: eps12; 25: eps11_s;26: eps22_s
%% 27:eps33_s; 28: eps23_s ; 29:  eps13_s; 30: eps12_s; 
%% 31: sig11; 32: sig22; 33: sig33; 34: sig23; 35: sig13; 36: sig12; 37:sig11_s;
%% 38: sig22_s; 39: sig33_s; 40: sig23_s; 41:sig13_s; 42: sig12_s; 43: eta; 44: eta; 45: previous ID; 46: previous layer. 
%%
%%

clear
%% parameters 
TOLDistNei=25; % distance tolerance for finding neighbours! average number of neighbour per grain is 6!
MisOrTol=5;
MeanNeighNom=6; %mean number of neighbours per grains
TOLDistGrainXC=25; % distance tolerance for finding grains accross two layers


assm1=importdata('Assembly_PreLoad.dat');
assm1err=importdata('AssemblyErr_PreLoad.dat');

assm2=importdata('Assembly_Onset.dat');
assm2err=importdata('AssemblyErr_Onset.dat');

assm3=importdata('Assembly_OnePercent.dat');
assm3err=importdata('AssemblyErr_OnePercent.dat');

assm4=importdata('Assembly_Unload.dat');
assm4err=importdata('AssemblyErr_Unload.dat');

[idmatchf_a2a1,idmatch_a2a1, idnotmttomanym_a2a1, idnomatchm_a2a1, dmin_a2a1, idrepnm_a2a1, XMEANDIS_a2a1, COMAssmA1, COMAssmA2]= ...
    mathingtwosteps(assm1,assm2,TOLDistNei, MisOrTol, MeanNeighNom,TOLDistGrainXC );
disp ('assm2 & assm1 are now matched')
save('Matchouput')

[idmatchf_a3a1,idmatch_a3a1, idnotmttomanym_a3a1, idnomatchm_a3a1, dmin_a3a1, idrepnm_a3a1, XMEANDIS_a3a1, COMAssmA1, COMAssmA3]= ...
    mathingtwosteps(assm1,assm3,TOLDistNei, MisOrTol, MeanNeighNom,TOLDistGrainXC );
disp ('assm3 & assm1 are now matched')
save('Matchouput')
 
[idmatchf_a4a1,idmatch_a4a1, idnotmttomanym_a4a1, idnomatchm_a4a1, dmin_a4a1, idrepnm_a4a1, XMEANDIS_a4a1, COMAssmA1, COMAssmA4]= ...
    mathingtwosteps(assm1,assm4,TOLDistNei, MisOrTol, MeanNeighNom,TOLDistGrainXC );
disp ('assm4 & assm1 are now matched')
save('Matchouput')
%% match with 2
[idmatchf_a3a2,idmatch_a3a2, idnotmttomanym_a3a2, idnomatchm_a3a2, dmin_a3a2, idrepnm_a3a2, XMEANDIS_a3a2, COMAssmA2, COMAssmA3]= ... 
    mathingtwosteps(assm2,assm3,TOLDistNei, MisOrTol, MeanNeighNom,TOLDistGrainXC );
disp ('assm3 & assm2 are now matched')
save('Matchouput')
 
[idmatchf_a4a2,idmatch_a4a2, idnotmttomanym_a4a2, idnomatchm_a4a2, dmin_a4a2, idrepnm_a4a2, XMEANDIS_a4a2, COMAssmA2, COMAssmA4]= ...
    mathingtwosteps(assm2,assm4,TOLDistNei, MisOrTol, MeanNeighNom,TOLDistGrainXC );
disp ('assm4 & assm2 are now matched')
save('Matchouput')
%% match with 3
[idmatchf_a4a3,idmatch_a4a3, idnotmttomanym_a4a3, idnomatchm_a4a3, dmin_a4a3, idrepnm_a4a3, XMEANDIS_a4a3, COMAssmA3, COMAssmA4]= ... 
    mathingtwosteps(assm3,assm4,TOLDistNei, MisOrTol, MeanNeighNom,TOLDistGrainXC );
disp ('assm4 & assm3 are now matched')
save('Matchouput')