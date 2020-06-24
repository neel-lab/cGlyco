%% Step 1 test case - Peak list generation
% B cell
loadpath = 'C:\Users\GoEason\Desktop\Paper&dissertation\P0_GDAT\cGlyco\Testcase\Bcell';
[MSdata,fitPara]=peakList('Bcell_Human','msrawdatapathstr',loadpath);
%% Step 2 test case - database building 
% N-linked glycan
antennalib = {'';...                      % define all possible form of antennas here 1
    '{n}';...                             % 2
    '{n{h_b4}}';...                       % 3
    '{n{f_a3}{h_b4}}';...                 % 4
    '{n{h_b4{s_a3}}}';...                 % 5
    '{n{f_a3}{h_b4{s_a3}}}';...           % 6
    '{n{h_b4{n_b3{h_b4}}}}';...           % 7
    '{n{h_b4{n_b3{f_a3}{h_b4}}}}';...     % 8
    '{n{h_b4{n_b3{h_b4{s_a3}}}}}';...     % 9
    '{n{h_b4{n_b3{f_a3}{h_b4{s_a3}}}}}'}; % 10

coreStruct  ='core1';
BranchType  ='Bi';
withLinkage = 1;
storepath   = 'C:\Users\GoEason\Desktop\Paper&dissertation\P0_GDAT\cGlyco\Testcase\';
highmannose = 1;
GlyType     = 'N-';

[glycanDB,outputname,LibSize] = genGlyLib(antennalib,coreStruct,withLinkage,storepath,GlyType,BranchType,highmannose);

% O-linked glycan
antennalib = {'';...                        % define all possible form of antennas here 1
    '{s}';...                               % 2
    '{n{h_b4}}';...                         % 3
    '{n{f_a3}{h_b4}}';...                   % 4
    '{n{h_b4{s_a3}}}';...                   % 5
    '{n{f_a3}{h_b4{s_a3}}}';...             % 6
    '{n{h_b4{n_b3{h_b4}}}}';...             % 7
    '{n{h_b4{n_b3{f_a3}{h_b4}}}}';...       % 8
    '{n{h_b4{n_b3{h_b4{s_a3}}}}}';...       % 9
    '{n{h_b4{n_b3{f_a3}{h_b4{s_a3}}}}}'};   % 10

coreStruct  = 'core1';
withLinkage = 1;
storepath   = 'C:\Users\GoEason\Desktop\Paper&dissertation\P0_GDAT\cGlyco\Testcase\';
GlyType     = 'O-';

[glycanDB,outputname,LibSize] = genGlyLib(antennalib,coreStruct,withLinkage,storepath,GlyType);

%% Step 3 test case - MS abundance
% B cell
MSDataFileName         = 'BcellMSdata';
MSDataLoadpath         = 'C:\Users\GoEason\Desktop\Paper&dissertation\P0_GDAT\cGlyco\Testcase\Bcell';
glycanDBName           = 'BcellglycanDB';
glycanDBLoadpath       = 'C:\Users\GoEason\Desktop\Paper&dissertation\P0_GDAT\cGlyco\Testcase\Bcell';
OverSegmentationFilter = 1.5;
OutputExcelFileName    = 'B_cell';
storepath              = 'C:\Users\GoEason\Desktop\Paper&dissertation\P0_GDAT\cGlyco\Testcase\';
[newglycanDB,outputfilename] = glycanAbundance(MSDataFileName,glycanDBName,MSDataLoadpath,glycanDBLoadpath,OverSegmentationFilter,OutputExcelFileName,storepath);

%% Step 4 test case - MS analysis
% B cell
MSfilename = 'Bcell';
GlycanDB   = 'Bcell_HumanglycanDB';
loadpath   = 'C:\Users\GoEason\Desktop\Paper&dissertation\P0_GDAT\cGlyco\Testcase\Bcell';
storepath  = 'C:\Users\GoEason\Desktop\Paper&dissertation\P0_GDAT\cGlyco\Testcase';
newglycanDB = classification(MSfilename,GlycanDB,loadpath,storepath);
