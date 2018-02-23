%% Import Script for EBSD Data
%
% This script was automatically created by the import wizard. You should
% run the whoole script or parts of it in order to import your data. There
% is no problem in making any changes to this script.

%% Specify Crystal and Specimen Symmetries

% crystal symmetry
CS = symmetry('6/m', [1 1 1.63], 'X||a*', 'Y||b', 'Z||c', 'mineral', 'unkown', 'color', 'light blue');

% specimen symmetry
SS = symmetry('triclinic');

% plotting convention
setMTEXpref('xAxisDirection','north');
setMTEXpref('zAxisDirection','outOfPlane');

%% Specify File Names

% path to files
pname = 'C:\Users\Hamid\Documents\OxfordDocs\CPFE\3dxrd\ForJon\OldVersion\RevisedVoronoi\RevisedVoronoi';

% which files to be imported
fname = [pname '\GRAI001.dat'];

%% Import the Data

% create an EBSD variable containing the data
ebsd = loadEBSD(fname,CS,SS,'interface','generic',...
  'ColumnNames', { 'x' 'y' 'z' 'phi1' 'Phi' 'phi2'}, 'Columns', [2 3 4 5 6 7], 'Bunge');

plot(ebsd)

