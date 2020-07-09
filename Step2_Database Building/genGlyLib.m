function [glycanDB,outputname,LibSize] = genGlyLib(antennalib,coreStruct,withLinkage,storedir,GlyType,varargin)
% genGlyLib: generates glycan library by combining user defined glycan
% antenna with likage infomation
%
% Syntax:
% NglyLib = genGlyLib(antennalib)
%
% Input:
%   antennalib: n x 1 cell array of strings containing antenna information in
%   SGP format
%
%   CoreStruct: core1-normal, core2-including core fucosylation, core3-including bisecting,
%   core4-including core fucosylation and bisecting
%
%   withLinkage: Define if the user need linkage information
%
%   For N-glycan only
%   BranchType: Define the branch type the user want to include in the
%   library.'Bi', 'Tri', 'Tetra'
%
%   highmannose: Define if the user need to include high mannose structure
%
% Output:
%   glycanDB: Glycan data structure including 4 fields: i) glycan structure
%   list in SGP format, ii) glycan composition, iii)monoisotopic mass and iv)
%   isotopic distribution
%
% Example1:
% antennalib = {'';...                      % define all possible form of antennas here 1
%     'GlcNAc(b1-?)';...                             % 2
%     'Gal(b1-4)GlcNAc(b1-?)';...                       % 3
%     'Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-?)';...                 % 4
%     'NeuAc(a2-3)Gal(b1-4)GlcNAc(b1-?)';...                 % 5
%     'NeuAc(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-?)';...           % 6
%     'Gal(b1-4)GlcNAc(b1-3)Gal(b1-4)GlcNAc(b1-?)';...           % 7
%     'Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-4)GlcNAc(b1-?)';...     % 8
%     'NeuAc(a2-3)Gal(b1-4)GlcNAc(b1-3)Gal(b1-4)GlcNAc(b1-?)';...     % 9
%     'NeuAc(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-4)GlcNAc(b1-?)'}; % 10
%
%   coreStruct  ='core1';
%   BranchType  ='Bi';
%   withLinkage = 1;
%   storepath   = 'C:\Users\GoEason\Desktop\Paper&dissertation\P0_GDAT\GDAT\Testcase\';
%   highmannose = 1;
%   GlyType     = 'N-';
%
%   [glycanDB,outputname,LibSize] = genGlyLib(antennalib,coreStruct,withLinkage,storepath,GlyType,BranchType,highmannose);
%
% Example2:
% antennalib = {'';...                        % define all possible form of antennas here 1
%     'NeuAc(a2-?)';...                               % 2
%     'Gal(b1-4)GlcNAc(b1-?)';...                         % 3
%     'Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-?)';...                   % 4
%     'NeuAc(a2-3)Gal(b1-4)GlcNAc(b1-?';...                   % 5
%     'NeuAc(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-?)';...             % 6
%     'Gal(b1-4)GlcNAc(b1-3)Gal(b1-4)GlcNAc(b1-?)';...             % 7
%     'Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-4)GlcNAc(b1-?)';...       % 8
%     'NeuAc(a2-3)Gal(b1-4)GlcNAc(b1-3)Gal(b1-4)GlcNAc(b1-?)';...       % 9
%     'NeuAc(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-4)GlcNAc(b1-?)'};   % 10
% 
%   coreStruct  = 'core1';
%   withLinkage = 1;
%   storepath   = 'C:\Users\GoEason\Desktop\Paper&dissertation\P0_GDAT\GDAT\Testcase\';
%   GlyType     = 'O-';
%
%   [glycanDB,outputname,LibSize] = genGlyLib(antennalib,coreStruct,withLinkage,storepath,GlyType);
%
% Author: Yusen Zhou
% Date Lastly Updated: 06/25/2020

% Convert IUPAC to SGP format
numOfantenna  = length(antennalib);
newantennalib = cell(numOfantenna,1);
for i = 1:numOfantenna
    ithantenna = antennalib{i};
    newantennalib{i} = IUPAC2Sgp(ithantenna);
end

if(strcmp(GlyType,'N-'))
    BranchType = varargin{1};
    highmannose = varargin{2};
    [glycanDB,outputname,LibSize] = genNglylib(newantennalib,coreStruct,BranchType,withLinkage,highmannose,storedir);
else
    [glycanDB,outputname,LibSize] = genOglylib(newantennalib,coreStruct,withLinkage,storedir);
end
end