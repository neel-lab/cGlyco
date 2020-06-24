function [glycanDB,outputname,LibSize] = genNglylib(antennalib,coreStruct,BranchType,withLinkage,highmannose,storepath )
% genlib: generate N-glycan library by combining user defined glycan
% antenna with likage infomation
%
% Syntax:
% NglyLib = genlib_human(antennalib)
%
% Input:
% antennalib: n x 1 cell array of strings containing antenna information in
% SGP format
%
% CoreStruct: core1-normal, core2-including core fucosylation, core3-including bisecting,
% core4-including core fucosylation and bisecting
%
% BranchType: Define the branch type the user want to include in the
% library.'Bi', 'Tri', 'Tetra'
%
% withLinkage: Define if the user need linkage information
%
% Output:
% NglyLib: m x 1 cell array of strings containing glycan combined
%
% Note:
% If certain restrictions on glycan antennas are required, modifications
% need to be made on variable "tetstrmapnumber"
%
% Children function: UNREPNUMGEN.m redsizebyrule_linkage.m
% unrepnumgen_linkage.m
% Author: Kai Cheng, Yusen Zhou
% Date Lastly Updated: 05/18/2020

%% Generate N_Glycan structure list
fullpath = [storepath 'Nglycanlib'];
if(withLinkage)
    antenna_linkage = addlinkage(antennalib);
    antennanum      = numel(antenna_linkage)/4;
    tetstrmapnumber = unrepnumgen_linkage(antennanum);
    badlist = redsizebyrule_linkage(tetstrmapnumber,antennalib,BranchType);
else
    antennanum = numel(antennalib);
    tetstrmapnumber = unrepnumgen(antennanum);
    badlist = redsizebyrule(tetstrmapnumber,BranchType);
end
goodlist = setdiff(1:size(tetstrmapnumber,1),unique(badlist));
tetstrmapnumber = tetstrmapnumber(goodlist,:);
if(withLinkage)
    NglyLib = buidlibrary_linkage(tetstrmapnumber,coreStruct,antenna_linkage,highmannose);
else
    NglyLib = buidlibrary(tetstrmapnumber,coreStruct,antennalib,highmannose,coreChange);
end

%% Generate glycan structure object
massoption.acetylation = false;
massoption.methylation = true;
massoption.ion = 'Na';

[glycanDB,outputname,LibSize] = reformatStrlib(NglyLib,1,massoption,storepath);
end

function antenna_linkage = addlinkage(antennalib)
% addlinkage: Add linkage information of first residue to the SGP format
%
% Syntax:
% antenna_linkage = addlinkage(antennalib,linkageinfo)
%
% Input:
% antennalib: n x 1 cell array of strings containing antenna information in
%             SGP format
%
% Example:
% antennalib = {'';...                      % define all possible form of antennas here 1
%     '{n{h_b4}}';...                       % 2
%     '{n{f_a3}{h_b4}}';...                 % 3
%     '{n{h_b4{s_a3}}}';...                 % 4
%     '{n{f_a3}{h_b4{s_a3}}}';...           % 5
%     '{n{h_b4{n_b3{h_b4}}}}';...           % 6
%     '{n{h_b4{n_b3{f_a3}{h_b4}}}}';...     % 7
%     '{n{h_b4{n_b3{h_b4{s_a3}}}}}';...     % 8
%     '{n{h_b4{n_b3{f_a3}{h_b4{s_a3}}}}}'}; % 9
%
% antenna_linkage = addlinkage(antennalib);
%
%
% Author: Yusen Zhou 6/9/2015

antenna_linkage = cell(length(antennalib),4);
for i = 1 : 4
    cellnum = 0;
    for j = 1 : length(antennalib)
        jthantenna = antennalib{j};
        if(i==1||i==2)
            jthantenna = regexprep(jthantenna,'{n','{n_b2',1);
        elseif(i==3)
            jthantenna = regexprep(jthantenna,'{n','{n_b4',1);
        elseif(i==4)
            jthantenna = regexprep(jthantenna,'{n','{n_b6',1);
        end
        antenna_linkage{cellnum+1,i} = jthantenna;
        cellnum = cellnum + 1;
    end
end

end

function combi = unrepnumgen_linkage(maxnum)
% UNREPNUMGEN: generate 4 digit combinations for describing glycan
% structure which will not contain repeated motif
%
% Syntax:
% combi = UNREPNUMGEN(maxnum)
%
% Input:
% maxnum: maximum number that will appear in the combination, all numbers
% varies from 1 to maxnum
%
% Output:
% combi: n x 4 matrix containing the 4 digit combination
%
% Note:
% none
% Example:
% combi = unrepnumgen(2)
%
% combi =
%     1     1     1     2     1     1
%     1     1     2     1     1     1
%     1     1     2     2     2     2
%     1     2     1     1    10     1
%     1     2     1     2    11     2
%     1     2     2     1    11     2
%     1     2     2     2    12     3
%     2     1     1     1   100     1
%     2     1     1     2   101     2
%     2     1     2     1   101     2
%     2     1     2     2   102     3
%     2     2     1     1   110     2
%     2     2     1     2   111     3
%     2     2     2     1   111     3
%     2     2     2     2   112     4
%
% Children function: none

% Author: Kai Cheng, Yusen Zhou 6/9/2015

twingrp = zeros(maxnum*maxnum,2);
ind = 1;
for i = 1:maxnum
    for j = 1:maxnum
        twingrp(ind,1) = i;
        twingrp(ind,2) = j;
        ind = ind+1;
    end
end
a = length(twingrp);
quadgrp = zeros(a*a,6);
ind = 1;
for i = 1:a
    for j = 1:a
        quadgrp(ind,1:2) = twingrp(i,:);
        quadgrp(ind,3:4) = twingrp(j,:);
        quadgrp(ind,5)   = 0;
        quadgrp(ind,6)   = 0;
        for k = 1:4
            if(quadgrp(ind,k)==1)
                kthantenna = 0;
            else
                kthantenna = 1;
                quadgrp(ind,6) = quadgrp(ind,6)+1;
            end
            if(k==1)
                quadgrp(ind,5) = quadgrp(ind,5)+kthantenna*100;
            elseif(k==2)
                quadgrp(ind,5) = quadgrp(ind,5)+kthantenna*10;
            else
                quadgrp(ind,5) = quadgrp(ind,5)+kthantenna;
            end
        end
        ind = ind+1;
    end
end
combi = quadgrp;
end

function badlist = redsizebyrule_linkage(tetstrmapnumber,antennalib,BranchType)
badlist = [];
for i = 1:size(tetstrmapnumber,1)
    if(tetstrmapnumber(i,5)~=0)
        % rule1: GlcNAcT-I is the first step of glycosylation.
        if(tetstrmapnumber(i,5)<100)
            badlist = [badlist i];
        else
            % rule2: MGAT5 can only act after MGAT2
            if(tetstrmapnumber(i,4)~=1)
                if(tetstrmapnumber(i,2)==1)
                    badlist = [badlist i];
                end
            end
            % rule3: no more than 1 LacNAc motif is allowed in structure
            % with branch less than 3.
            isexist = 0;
            pattern = 'n';
            pattern2 = 'h';
            for j = 1 : length(antennalib)
                if(length(strfind(antennalib{j},pattern))==2)&&...
                        (length(strfind(antennalib{j},pattern2))==2)
                    border = j;
                    isexist = 1;
                    break
                end
            end
            if(tetstrmapnumber(i,6)<4&&isexist)
                if(sum(tetstrmapnumber(i,1:4)>=border)>0)
                    badlist = [badlist i];
                end
            end
        end
    end
    % rule4: Define the branch type the user want to include in the library
    if(strcmp(BranchType,'Bi'))
        if(tetstrmapnumber(i,6)>2)
            badlist = [badlist i];
        end
    elseif(strcmp(BranchType,'Tri'))
        if(tetstrmapnumber(i,6)>3)
            badlist = [badlist i];
        end
    end
end
end

function combi = unrepnumgen(maxnum)
% UNREPNUMGEN: generate 4 digit combinations for describing glycan
% structure which will not contain repeated motif
%
% Syntax:
% combi = UNREPNUMGEN(maxnum)
%
% Input:
% maxnum: maximum number that will appear in the combination, all numbers
% varies from 1 to maxnum
%
% Output:
% combi: n x 4 matrix containing the 4 digit combination
%
% Note:
% none
% Example:
% combi = unrepnumgen(2)
%
% combi =
%
%      1     1     1     1    0
%      1     1     1     2    1
%      1     1     2     2    2
%      1     2     1     2    2
%      1     2     2     2    3
%      2     2     2     2    4
%
% Children function: none

twingrp = zeros((1+maxnum)*maxnum/2,2);
ind = 1;
for i = 1:maxnum
    for j = i:maxnum
        twingrp(ind,1) = i;
        twingrp(ind,2) = j;
        ind = ind+1;
    end
end
a = length(twingrp);
quadgrp = zeros((a+1)*a/2,5);
ind = 1;
for i = 1:a
    for j = i:a
        quadgrp(ind,1:2) = twingrp(i,:);
        quadgrp(ind,3:4) = twingrp(j,:);
        for k = 1:4
            if(quadgrp(ind,k)~=1)
                quadgrp(ind,5) = quadgrp(ind,5)+1;
            end
        end
        ind = ind+1;
    end
end
combi = quadgrp;
end

function badlist = redsizebyrule(tetstrmapnumber,BranchType)
badlist = [];
for i = 1:size(tetstrmapnumber,1)
    if(strcmp(BranchType,'Bi'))
        if(tetstrmapnumber(i,5)>2)
            badlist = [badlist i];
        end
    elseif(strcmp(BranchType,'Tri'))
        if(tetstrmapnumber(i,5)>3)
            badlist = [badlist i];
        end
    end
end
end

function NglyLib = buidlibrary_linkage(tetstrmapnumber,CoreStruct,antenna_linkage,highmannose)
ncore = {'{n_b{n_b4{h_b4{h_a3','}{h_a6','}}}}'};
tetanta = cell(size(tetstrmapnumber,1),1);
tetantb = tetanta;
tetantc = tetanta;
tetantd = tetanta;
tetante = cell(size(antenna_linkage,1)-1,1);
tetantf = cell(size(antenna_linkage,1)-1,1);
for i = 1:size(tetstrmapnumber,1)
    tetanta{i} = strcat(ncore{1},antenna_linkage{tetstrmapnumber(i,1),1},antenna_linkage{tetstrmapnumber(i,3),3},ncore{2},antenna_linkage{tetstrmapnumber(i,2),2}...
        ,antenna_linkage{tetstrmapnumber(i,4),4},ncore{3});
    if(strcmp(CoreStruct,'core2')||strcmp(CoreStruct,'core4'))
        ncoremodify = '{n_b{f_a6}{n_b4{h_b4{h_a3';
        tetantb{i} = strcat(ncoremodify,antenna_linkage{tetstrmapnumber(i,1),1},antenna_linkage{tetstrmapnumber(i,3),3},ncore{2},antenna_linkage{tetstrmapnumber(i,2),2}...
            ,antenna_linkage{tetstrmapnumber(i,4),4},ncore{3});
    end
    if(strcmp(CoreStruct,'core3')||strcmp(CoreStruct,'core4'))
        ncoremodify = '{n_b{n_b4{h_b4{n_b4}{h_a3';
        tetantc{i} = strcat(ncoremodify,antenna_linkage{tetstrmapnumber(i,1),1},antenna_linkage{tetstrmapnumber(i,3),3},ncore{2},antenna_linkage{tetstrmapnumber(i,2),2}...
            ,antenna_linkage{tetstrmapnumber(i,4),4},ncore{3});
    end
    if(strcmp(CoreStruct,'core4'))
        ncoremodify = '{n_b{f_a6}{n_b4{h_b4{n_b4}{h_a3';
        tetantd{i} = strcat(ncoremodify,antenna_linkage{tetstrmapnumber(i,1),1},antenna_linkage{tetstrmapnumber(i,3),3},ncore{2},antenna_linkage{tetstrmapnumber(i,2),2}...
            ,antenna_linkage{tetstrmapnumber(i,4),4},ncore{3});
    end
end
for i = 2 : size(antenna_linkage)
    tetante{i-1} = strcat(ncore{1},antenna_linkage{i,1},'}{h_a6{h_a3}{h_a6}',ncore{3});
end
if(strcmp(CoreStruct,'core3')||strcmp(CoreStruct,'core4'))
    for i = 2 : size(antenna_linkage)
        ncoremodify = '{n_b{n_b4{h_b4{n_b4}{h_a3';
        tetantf{i-1} = strcat(ncoremodify,antenna_linkage{i,1},'}{h_a6{h_a3}{h_a6}',ncore{3});
    end
end
NglyLib= [tetanta;tetantb;tetantc;tetantd;tetante;tetantf];
NglyLib = NglyLib(~cellfun('isempty',NglyLib));
if(highmannose)
    NglyLib = [NglyLib;{'{n_b{n_b4{h_b4{h_a3{n_b2}}{h_a6{h_a3}}}}}';...
        '{n_b{n_b4{h_b4{h_a3}{h_a6{h_a3}{h_a6}}}}}';...
        '{n_b{n_b4{h_b4{h_a3{h_a2}}{h_a6{h_a3}{h_a6}}}}}';...
        '{n_b{n_b4{h_b4{h_a3}{h_a6{h_a3{h_a2}}{h_a6}}}}}';...
        '{n_b{n_b4{h_b4{h_a3}{h_a6{h_a3}{h_a6{h_a2}}}}}}';...
        '{n_b{n_b4{h_b4{h_a3}{h_a6{h_a3{h_a2}}{h_a6{h_a2}}}}}}';...
        '{n_b{n_b4{h_b4{h_a3{h_a2}{h_a2}}{h_a6{h_a3}{h_a6}}}}}';...
        '{n_b{n_b4{h_b4{h_a3{h_a2}}{h_a6{h_a3{h_a2}}{h_a6}}}}}';...
        '{n_b{n_b4{h_b4{h_a3{h_a2}}{h_a6{h_a3}{h_a6{h_a2}}}}}}';...
        '{n_b{n_b4{h_b4{h_a3{h_a2}}{h_a6{h_a3{h_a2}}{h_a6{h_a2}}}}}}';...
        '{n_b{n_b4{h_b4{h_a3{h_a2{h_a2}}}{h_a6{h_a3{h_a2}}{h_a6}}}}}';...
        '{n_b{n_b4{h_b4{h_a3{h_a2{h_a2}}}{h_a6{h_a3}{h_a6{h_a2}}}}}}';...
        '{n_b{n_b4{h_b4{h_a3{h_a2{h_a2}}}{h_a6{h_a3{h_a2}}{h_a6{h_a2}}}}}}'}];
end
end

function NglyLib = buidlibrary(tetstrmapnumber,CoreStruct,antennalib,highmannose)
ncore = {'{n{n{h{h','}{h','}}}}'};
tetanta = cell(size(tetstrmapnumber,1),1);
tetantb = tetanta;
tetantc = tetanta;
tetantd = tetanta;
tetante = cell(size(antennalib,1)-1,1);
tetantf = cell(size(antennalib,1)-1,1);
for i = 1:size(tetstrmapnumber,1)
    tetanta{i} = strcat(ncore{1},antennalib{tetstrmapnumber(i,1)},antennalib{tetstrmapnumber(i,2)},ncore{2},antennalib{tetstrmapnumber(i,3)}...
        ,antennalib{tetstrmapnumber(i,4)},ncore{3});
    if(strcmp(CoreStruct,'core2')||strcmp(CoreStruct,'core4'))
        ncoremodify = '{n{f}{n{h{h';
        tetantb{i} = strcat(ncoremodify,antennalib{tetstrmapnumber(i,1)},antennalib{tetstrmapnumber(i,2)},ncore{2},antennalib{tetstrmapnumber(i,3)}...
            ,antennalib{tetstrmapnumber(i,4)},ncore{3});
    end
    if(strcmp(CoreStruct,'core3')||strcmp(CoreStruct,'core4'))
        ncoremodify = '{n{n{h{n}{h';
        tetantc{i} = strcat(ncoremodify,antennalib{tetstrmapnumber(i,1)},antennalib{tetstrmapnumber(i,2)},ncore{2},antennalib{tetstrmapnumber(i,3)}...
            ,antennalib{tetstrmapnumber(i,4)},ncore{3});
    end
    if(strcmp(CoreStruct,'core4'))
        ncoremodify = '{n{f}{n{h{n}{h';
        tetantd{i} = strcat(ncoremodify,antennalib{tetstrmapnumber(i,1)},antennalib{tetstrmapnumber(i,2)},ncore{2},antennalib{tetstrmapnumber(i,3)}...
            ,antennalib{tetstrmapnumber(i,4)},ncore{3});
    end
end
for i = 2 : size(antennalib)
    tetante{i-1} = strcat(ncore{1},antennalib{i},'}{h{h}{h}',ncore{3});
end
if(strcmp(CoreStruct,'core3')||strcmp(CoreStruct,'core4'))
    for i = 2 : size(antennalib)
        ncoremodify = '{n{n{h{n}{h';
        tetantf{i-1} = strcat('{n{n{h{n}{h',antennalib{i},'}{h{h}{h}',ncore{3});
    end
end
NglyLib= [tetanta;tetantb;tetantc;tetantd;tetante;tetantf];
NglyLib = NglyLib(~cellfun('isempty',NglyLib));
if(highmannose)
    NglyLib = [NglyLib;{'{n{n{h{h{n}}{h{h}}}}}';...
        '{n{n{h{h}{h{h}{h}}}}}';...
        '{n{n{h{h{h}}{h{h}{h}}}}}';...
        '{n{n{h{h{h}}{h{h{h}}{h}}}}}';...
        '{n{n{h{h{h{h}}}{h{h{h}}{h}}}}}';...
        '{n{n{h{h{h{h}}}{h{h{h}}{h{h}}}}}}'}];
end
end



% function antenna_linkage = addlinkage(antennalib)
% antenna_linkage = cell(length(antennalib),4);
% for i = 1 : 4
%     cellnum = 0;
%     for j = 1 : length(antennalib)
%         jthantenna = antennalib{j};
%         resiudes   = fieldnames(linkageinfo);
%         for k = 1 : length(resiudes)
%             kthresidue = resiudes{k};
%             if(strcmp(kthresidue,'n'))
%                 if(length(strfind(jthantenna,kthresidue))==1)
%                     expression = ['{' kthresidue];
%                     if(i==1||i==2)
%                         jthantenna = regexprep(jthantenna,expression,'{n_b2');
%                     elseif(i==3)
%                         jthantenna = regexprep(jthantenna,expression,'{n_b4');
%                     elseif(i==4)
%                         jthantenna = regexprep(jthantenna,expression,'{n_b6');
%                     end
%                 elseif(length(strfind(jthantenna,kthresidue))==2)
%                     expression = ['{' kthresidue];
%                     replace    = linkageinfo.(kthresidue);
%                     if(i==1||i==2)
%                         jthantenna = regexprep(jthantenna,expression,'{n_b2',1);
%                     elseif(i==3)
%                         jthantenna = regexprep(jthantenna,expression,'{n_b4',1);
%                     elseif(i==4)
%                         jthantenna = regexprep(jthantenna,expression,'{n_b6',1);
%                     end
%                     jthantenna = regexprep(jthantenna,expression,replace,2);
%                 end
%             else
%                 if(strfind(jthantenna,kthresidue))
%                     expression = ['{' kthresidue ];
%                     replace    = linkageinfo.(kthresidue);
%                     jthantenna = regexprep(jthantenna,expression,replace);
%                 end
%             end
%         end
%         antenna_linkage{cellnum+1,i} = jthantenna;
%         cellnum = cellnum + 1;
%     end
% end
% end