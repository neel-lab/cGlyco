function [glycanDB,outputname,LibSize] = genGlyLib(antennalib,coreStruct,withLinkage,storepath,GlyType,varargin)
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
%   antennalib = {'';...                    % define all possible forms of antennas here 1
%     '{n}';...                             % 2
%     '{n{h_b4}}';...                       % 3
%     '{n{f_a3}{h_b4}}';...                 % 4
%     '{n{h_b4{s_a3}}}';...                 % 5
%     '{n{f_a3}{h_b4{s_a3}}}';...           % 6
%     '{n{h_b4{n_b3{h_b4}}}}';...           % 7
%     '{n{h_b4{n_b3{f_a3}{h_b4}}}}';...     % 8
%     '{n{h_b4{n_b3{h_b4{s_a3}}}}}';...     % 9
%     '{n{h_b4{n_b3{f_a3}{h_b4{s_a3}}}}}'}; % 10
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
%   antennalib = {'';...                        % define all possible form of antennas here 1
%     '{s}';...                               % 2
%     '{n{h_b4}}';...                         % 3
%     '{n{f_a3}{h_b4}}';...                   % 4
%     '{n{h_b4{s_a3}}}';...                   % 5
%     '{n{f_a3}{h_b4{s_a3}}}';...             % 6
%     '{n{h_b4{n_b3{h_b4}}}}';...             % 7
%     '{n{h_b4{n_b3{f_a3}{h_b4}}}}';...       % 8
%     '{n{h_b4{n_b3{h_b4{s_a3}}}}}';...       % 9
%     '{n{h_b4{n_b3{f_a3}{h_b4{s_a3}}}}}'};   % 10
%
%   coreStruct  = 'core1';
%   withLinkage = 1;
%   storepath   = 'C:\Users\GoEason\Desktop\Paper&dissertation\P0_GDAT\GDAT\Testcase\';
%   GlyType     = 'O-';
%
%   [glycanDB,outputname,LibSize] = genGlyLib(antennalib,coreStruct,withLinkage,storepath,GlyType);
%
% Author: Yusen Zhou
% Date Lastly Updated: 05/18/2020

if(strcmp(GlyType,'N-'))
    BranchType = varargin{1};
    highmannose = varargin{2};
    [glycanDB,outputname,LibSize] = genNglylib(antennalib,coreStruct,BranchType,withLinkage,highmannose,storepath);
else
    [glycanDB,outputname,LibSize] = genOglylib(antennalib,coreStruct,withLinkage,storepath);
end
end