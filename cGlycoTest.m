%% This is the test file for cGlyco with examples for Bcells and Monocytes.
%% To run this program, simply change the rootdir depending on the location of the installation directory of cGlyco
%% Then step throught the program to execute individual program modules
%% Define the root path of the cGlyco
rootdir = 'C:\Work\cGlyco';

%% Other dir
Bcelldir      = fullfile(rootdir,'Testcase\Bcell');
Monocytesdir  = fullfile(rootdir,'Testcase\Monocytes');
CompDir       = fullfile(rootdir,'Testcase\BcellMonocytes');
colormapdir   = fullfile(rootdir,'func');
%% Step 1 test case - Peak list generation
% B cell
[MSdata,fitPara]=peakList('Bcell_Human','dir',Bcelldir);
%% Step 2 test case - database building 
% N-linked glycan
antennalib = {'';...                      % define all possible form of antennas here 1
    'GlcNAc(b1-?)';...                             % 2
    'Gal(b1-4)GlcNAc(b1-?)';...                       % 3
    'Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-?)';...                 % 4
    'NeuAc(a2-3)Gal(b1-4)GlcNAc(b1-?)';...                 % 5
    'NeuAc(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-?)';...           % 6
    'Gal(b1-4)GlcNAc(b1-3)Gal(b1-4)GlcNAc(b1-?)';...           % 7
    'Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-4)GlcNAc(b1-?)';...     % 8
    'NeuAc(a2-3)Gal(b1-4)GlcNAc(b1-3)Gal(b1-4)GlcNAc(b1-?)';...     % 9
    'NeuAc(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-4)GlcNAc(b1-?)'}; % 10

coreStruct  ='core4';
BranchType  ='tetra';
withLinkage = 1;
outputdir   = Bcelldir;
highmannose = 1;
GlyType     = 'N-';

[glycanDB,outputname,LibSize] = genGlyLib(antennalib,coreStruct,withLinkage,outputdir,GlyType,BranchType,highmannose);

% O-linked glycan
antennalib = {'';...                        % define all possible form of antennas here 1
    'NeuAc(a2-?)';...                               % 2
    'Gal(b1-4)GlcNAc(b1-?)';...                         % 3
    'Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-?)';...                   % 4
    'NeuAc(a2-3)Gal(b1-4)GlcNAc(b1-?';...                   % 5
    'NeuAc(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-?)';...             % 6
    'Gal(b1-4)GlcNAc(b1-3)Gal(b1-4)GlcNAc(b1-?)';...             % 7
    'Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-4)GlcNAc(b1-?)';...       % 8
    'NeuAc(a2-3)Gal(b1-4)GlcNAc(b1-3)Gal(b1-4)GlcNAc(b1-?)';...       % 9
    'NeuAc(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-4)GlcNAc(b1-?)'};   % 10

coreStruct  = 'core1';
withLinkage = 1;
outputdir   = Bcelldir;
GlyType     = 'O-';

[glycanDB,outputname,LibSize] = genGlyLib(antennalib,coreStruct,withLinkage,outputdir,GlyType);

%% Step 3 test case - MS abundance
% B cell
MSDataFileName  = 'Bcell_Human';
MSdir           = Bcelldir;
glycanDBName    = 'BcellglycanDB';
glycnDBdir      = Bcelldir;
OverSegFilter   = 1.5;
outputExcel     = 'B_cell';
outputdir       = Bcelldir;
[newglycanDB,outputfilename] = glycanAbundance(MSDataFileName,glycanDBName,MSdir,glycnDBdir,OverSegFilter,outputExcel,outputdir);

%% Step 4 test case - MS analysis
% listOfStruct
listOfStruct = {'LeX','sialylLeX','LacNAc','Sia_LacNAc','Bisecting','coreFuc'};
% B cell
MSfilename  = 'Bcell';
glycanDB    = 'Bcell_HumanglycanDB';
glycanDBdir = Bcelldir;
outputdir   = Bcelldir;
newglycanDB = classification(MSfilename,glycanDB,glycanDBdir,outputdir,listOfStruct);

%% Comparison analysis
% listOfStruct
listOfStruct = {'LeX','sialylLeX','LacNAc','Sia_LacNAc','Bisecting','coreFuc'};
% Bcell
glycanDB = fullfile(Bcelldir,'Bcell_HumanglycanDB.mat');
load(glycanDB,'newglycanDB');
ResidueResult = MSResidueAnalysisGNAT('Bcell',newglycanDB);
% Save results to BcellMonocytes folder
storedir = fullfile(CompDir,'BcellResidue.mat');
save(storedir,'ResidueResult');
StructResult  = MSStructAnalysisGNAT('Bcell',newglycanDB,listOfStruct);
% Save results to BcellMonocytes folder
storedir = fullfile(CompDir,'BcellStruct.mat');
save(storedir,'StructResult');
% Monocytes
glycanDB = fullfile(Monocytesdir,'Monocytes_HumanglycanDB.mat');
load(glycanDB,'newglycanDB');
ResidueResult = MSResidueAnalysisGNAT('Monocytes',newglycanDB);
% Save results to BcellMonocytes folder
storedir = fullfile(CompDir,'MonocytesResidue.mat');
save(storedir,'ResidueResult');
StructResult  = MSStructAnalysisGNAT('Monocytes',newglycanDB,listOfStruct);
% Save results to BcellMonocytes folder
storedir = fullfile(CompDir,'MonocytesStruct.mat');
save(storedir,'StructResult');
% Comparison
ResidueMSFileList = {'BcellResidue','MonocytesResidue'};
StructMSFileList  = {'BcellStruct','MonocytesStruct'};
outputfilename    = 'Bcell_Monocytes';
[ResidueTypeExist,StructTypeExist] = Plotheatmap(ResidueMSFileList,StructMSFileList,outputfilename,CompDir,colormapdir,CompDir);
% Structure heatmap
glycanDB1  = 'Bcell_HumanglycanDB';
glycanDB2  = 'Monocytes_HumanglycanDB';
loadpath   = CompDir;
storedir  = CompDir;
PlotheatmapStruct(glycanDB1,glycanDB2,loadpath,colormapdir,storedir)