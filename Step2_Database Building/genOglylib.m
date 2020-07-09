function [glycanDB,outputname,LibSize] = genOglylib(antennalib,coreStruct,withLinkage,storedir)
% genOglylib: generate O-glycan library by combining user defined glycan
% antenna with likage infomation
%
% Syntax:
% OglyLib = genOglylib(antennalib,linkageinfo,coreStruct)
%
% Input:
% antennalib: n x 1 cell array of strings containing antenna information in
% SGP format
%
% coreStruct: choose the coreStruct
%
% Output:
% OglyLib: m x 1 cell array of strings containing glycan combined
%
% Note:
% If certain restrictions on glycan antennas are required, modifications
% need to be made on variable "tetstrmapnumber"
%
% Children function: UNREPNUMGEN.m
% Author: Kai Cheng, Yusen Zhou 5/18/2020
if(withLinkage)
    if(strcmp(coreStruct,'core1')||strcmp(coreStruct,'core2'))
        ocore = {'{n_a{h_b3','}','}'};
    elseif(strcmp(coreStruct,'core3')||strcmp(coreStruct,'core4')||strcmp(coreStruct,'core5')...
            ||strcmp(coreStruct,'core6')||strcmp(coreStruct,'core7'))
        ocore = {'{n_a','','}'};
    elseif(strcmp(coreStruct,'core8'))
        ocore = {'{n_a{h_a3','}','}'};
    end
    antenna_linkage = crelinkage(antennalib,coreStruct);
    antennanum      = numel(antenna_linkage)/2;
    bistrmapnumber  = unrepnumgen_linkage(antennanum);
    badlist         = redsizebyrule_linkage(bistrmapnumber,antennalib,coreStruct);
    goodlist        = setdiff(1:size(bistrmapnumber,1),unique(badlist));
    bistrmapnumber  = bistrmapnumber(goodlist,:);
    OglyLib         = buildlib_linkage(bistrmapnumber,antenna_linkage,ocore);
else
    if(strcmp(coreStruct,'core1')||strcmp(coreStruct,'core2')||strcmp(coreStruct,'core8'))
        ocore = {'{n{h','}','}'};
    else
        ocore = {'{n','','}'};
    end
    antennanum      = numel(antennalib);
    if(strcmp(coreStruct,'core1'))
        bistrmapnumber  = unrepnumgen_linkage(antennanum);
    else
        bistrmapnumber  = unrepnumgen(antennanum);
    end
    badlist         = redsizebyrule(bistrmapnumber,antennalib,coreStruct);
    goodlist        = setdiff(1:size(bistrmapnumber,1),unique(badlist));
    bistrmapnumber  = bistrmapnumber(goodlist,:);
    OglyLib         = buildlib(bistrmapnumber,antennalib,ocore);
end

%% Generate glycan structure object
massoption.acetylation = false;
massoption.methylation = true;
massoption.ion = 'Na';

[glycanDB,outputname,LibSize] = reformatStrlib(OglyLib,0,massoption,storedir);
end

function antenna_linkage = crelinkage(antennalib,coreStruct)
% crelinkage: Add linkage information of first residue to the SGP format
%
% Syntax:
% antenna_linkage = crelinkage(antennalib,coreStruct)
%
% Input:
% antennalib: n x 1 cell array of strings containing antenna information in
%             SGP format
% coreStruct: O-glycan core Structure
%
%
% Example:
% antennalib = {'';...                        % define all possible form of antennas here 1
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
% coreStruct = 'core1';
%
% antenna_linkage = crelinkage(antennalib,coreStruct);
%
%
% Author: Yusen Zhou 6/10/2015


antenna_linkage = cell(length(antennalib),2);
for i = 1:2
    cellnum = 0;
    if(strcmp(coreStruct,'core1')||strcmp(coreStruct,'core2')||...
            strcmp(coreStruct,'core3')||strcmp(coreStruct,'core4')...
            ||strcmp(coreStruct,'core5')||strcmp(coreStruct,'core8'))
        for j = 1 : length(antennalib)
            jthantenna = antennalib{j};
            if(i==1)
                if(~isempty(regexp(jthantenna,'{s[{}]', 'once')))
                    jthantenna = regexprep(jthantenna,'{s','{s_a3',1);
                elseif(~isempty(regexp(jthantenna,'{n[{}]', 'once'))||...
                        ~isempty(regexp(jthantenna,'{m[{}]', 'once')))
                    if(strcmp(coreStruct,'core5'))
                        jthantenna = regexprep(jthantenna,'{n','{n_a3',1);
                    else
                        jthantenna = regexprep(jthantenna,'{n','{n_b3',1);
                    end
                end
            elseif(i==2)
                if(~isempty(regexp(jthantenna,'{s[{}]', 'once')))
                    jthantenna = regexprep(jthantenna,'{s','{s_a6',1);
                elseif(~isempty(regexp(jthantenna,'{n[{}]', 'once'))||...
                        ~isempty(regexp(jthantenna,'{m[{}]', 'once')))
                    jthantenna = regexprep(jthantenna,'{n','{n_b6',1);
                end
            end
            antenna_linkage{cellnum+1,i} = jthantenna;
            cellnum = cellnum + 1;
        end
    elseif(strcmp(coreStruct,'core6')||strcmp(coreStruct,'core7'))
        for j = 1 : length(antennalib)
            jthantenna = antennalib{j};
            if(i==1)
                if(~isempty(regexp(jthantenna,'{s[{}]', 'once')))
                    jthantenna = regexprep(jthantenna,'{s','{s_a?',1);
                elseif(~isempty(regexp(jthantenna,'{n[{}]', 'once'))||...
                        ~isempty(regexp(jthantenna,'{m[{}]', 'once')))
                    if(strcmp(coreStruct,'core6'))
                        jthantenna = regexprep(jthantenna,'{n','{n_b6',1);
                    else
                        jthantenna = regexprep(jthantenna,'{n','{n_a6',1);
                    end
                end
            elseif(i==2)
                if(~isempty(regexp(jthantenna,'{s[{}]', 'once')))
                    jthantenna = regexprep(jthantenna,'{s','{s_a3',1);
                elseif(~isempty(regexp(jthantenna,'{n[{}]', 'once'))||...
                        ~isempty(regexp(jthantenna,'{m[{}]', 'once')))
                    jthantenna = regexprep(jthantenna,'{n','{n_b3',1);
                end
            end
            antenna_linkage{cellnum+1,i} = jthantenna;
            cellnum = cellnum + 1;
        end
    end
end
end

function combi = unrepnumgen_linkage(maxnum)
% UNREPNUMGEN: generate 2 digit combinations for describing glycan
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
% combi: n x 2 matrix containing the 2 digit combination
%
% Note:
% none
% Example:
% combi = unrepnumgen(2)
%
% combi =
%    1     1
%    1     2
%    2     1
%    2     2
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
combi = twingrp;
end

function combi = unrepnumgen(maxnum)
% UNREPNUMGEN: generate 2 digit combinations for describing glycan
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
% combi: n x 2 matrix containing the 2 digit combination
%
% Note:
% none
% Example:
% combi = unrepnumgen(2)
%
% combi =
%    1     1
%    1     2
%    2     2
%
% Children function: none

% Author: Kai Cheng, Yusen Zhou 6/9/2015

twingrp = zeros((1+maxnum)*maxnum/2,2);
ind = 1;
for i = 1:maxnum
    for j = i:maxnum
        twingrp(ind,1) = i;
        twingrp(ind,2) = j;
        ind = ind+1;
    end
end
combi = twingrp;
end

function OglyLib = buildlib_linkage(bistrmapnumber,antenna_linkage,ocore)
OglyLib = cell(size(bistrmapnumber,1),1);
x = size(bistrmapnumber,1);
wtbar = waitbar(0);
for i = 1:size(bistrmapnumber,1)
    OglyLib{i} = strcat(ocore{1},antenna_linkage{bistrmapnumber(i,1),1},ocore{2},antenna_linkage{bistrmapnumber(i,2),2},ocore{3});
    waitbar(i/x);
end
OglyLib = OglyLib(~cellfun('isempty',OglyLib));
close(wtbar)
end

function OglyLib = buildlib(bistrmapnumber,antennalib,ocore)
OglyLib = cell(size(bistrmapnumber,1),1);
x = size(bistrmapnumber,1);
wtbar = waitbar(0);
for i = 1:size(bistrmapnumber,1)
    OglyLib{i} = strcat(ocore{1},antennalib{bistrmapnumber(i,1)},ocore{2},antennalib{bistrmapnumber(i,2)},ocore{3});
    waitbar(i/x);
end
OglyLib = OglyLib(~cellfun('isempty',OglyLib));
close(wtbar)
end

function badlist = redsizebyrule_linkage(bistrmapnumber,antennalib,coreStruct)
badlist=[];
for i = 1 : length(antennalib)
    pattern = '{n{';
    if(strfind(antennalib{i},pattern))
        nonHexNAcP = i-1;
        break
    end
end
for i = 1:size(bistrmapnumber,1)
    if(strcmp(coreStruct,'core1')||strcmp(coreStruct,'core3')||strcmp(coreStruct,'core5'))
        if(bistrmapnumber(i,2)>nonHexNAcP)
            badlist = [badlist i];
        end
    elseif(strcmp(coreStruct,'core3')||strcmp(coreStruct,'core4')...
            ||strcmp(coreStruct,'core5')||strcmp(coreStruct,'core6')...
            ||strcmp(coreStruct,'core7'))
        if(bistrmapnumber(i,1)<=nonHexNAcP)
            badlist = [badlist i];
        end
        if(strcmp(coreStruct,'core3'))
            if(bistrmapnumber(i,2)>nonHexNAcP)
                badlist = [badlist i];
            end
        end
    elseif(strcmp(coreStruct,'core6')||strcmp(coreStruct,'core7'))
        if(bistrmapnumber(i,1)>nonHexNAcP)
            badlist = [badlist i];
        end
    end
end
end

function badlist = redsizebyrule(bistrmapnumber,antennalib,coreStruct)
badlist=[];
for i = 1 : length(antennalib)
    pattern = '{n{';
    if(strfind(antennalib{i},pattern))
        nonHexNAcP = i-1;
        break
    end
end
for i = 1:size(bistrmapnumber,1)
    if(strcmp(coreStruct,'core1')||strcmp(coreStruct,'core3')||strcmp(coreStruct,'core5')...
            ||strcmp(coreStruct,'core6')||strcmp(coreStruct,'core7'))
        if(bistrmapnumber(i,1)>nonHexNAcP&&bistrmapnumber(i,2)>nonHexNAcP)
            badlist = [badlist i];
        end
    elseif(strcmp(coreStruct,'core3')||strcmp(coreStruct,'core4')...
            ||strcmp(coreStruct,'core5')||strcmp(coreStruct,'core6')...
            ||strcmp(coreStruct,'core7'))
        if(bistrmapnumber(i,1)<=nonHexNAcP&&bistrmapnumber(i,2)<=nonHexNAcP)
            badlist = [badlist i];
        end
        if(strcmp(coreStruct,'core3'))
            if(bistrmapnumber(i,1)>nonHexNAcP&&bistrmapnumber(i,2)>nonHexNAcP)
                badlist = [badlist i];
            end
        end
    end
end
end