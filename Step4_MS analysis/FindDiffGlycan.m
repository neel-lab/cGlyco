function [glycanListdiff, glycanListcommon]= FindDiffGlycan(glycanDB1,glycanDB2,loadpath,varargin)
%FindDiffGlycan is to compare two glycanlist and find the glycans that are
% unique in each glycanlist.
% 
% glycanList = FindDiffGlycan(glycanlist1,glycanlist2) compare the two
%  glycanlist2 to glycanlist1, and find the differences between these two
%  glycanlist.
%
% Input  : Two glycanlist.
% Output : glycanlist with unique glycans.
%
%
%Author: Yusen Zhou
%Date Lastly Updated: 05/18/2020

if(isempty(glycanDB1))||...
        (isempty(glycanDB2))||...
        (isempty(loadpath))
    error(message('IncorrectInputs'));
end
glycanDB1 = [glycanDB1 '.mat'];
fullpath  = fullfile(loadpath, glycanDB1);
load(fullpath, 'newglycanDB');
comparegroup1 = getcompos(newglycanDB);
glycanDB2 = [glycanDB2 '.mat'];
fullpath  = fullfile(loadpath, glycanDB2);
load(fullpath, 'newglycanDB');
comparegroup2 = getcompos(newglycanDB);

glycansincommon1 = [];
glycansindiff1   = [];
count = 0;
num   = 0;
for i = 1 : length(comparegroup1)
    ithcomp    = comparegroup1{i};
    isincommon = sum(strcmpi(comparegroup2,ithcomp));
    if(isincommon)
        glycansincommon1(count+1,1) = i;
        count = count+1;
    else
        glycansindiff1(num+1,1) = i;
        num = num+1;
    end  
end
glycansincommon2 = [];
count = 0;
for i = 1 : length(comparegroup2)
    ithcomp    = comparegroup2{i};
    isincommon  = sum(strcmpi(comparegroup1,ithcomp));
    if(isincommon)
        glycansincommon2(count+1,1) = i;
        count = count+1;
    end
end
glycanListdiff = struct('ListA',[],'ListB',[]);
comparegroupdiff1 = comparegroup1;
comparegroupdiff1(glycansincommon1)= '';
comparegroupdiff2 = comparegroup2;
comparegroupdiff2(glycansincommon2)= '';
glycanListdiff.ListA = comparegroupdiff1;
glycanListdiff.ListB = comparegroupdiff2;
comparegroup1(glycansindiff1)= '';
glycanListcommon = comparegroup1;
savepath1 = fullfile(loadpath, 'glycandiff.mat');
savepath2 = fullfile(loadpath, 'glycancommon.mat');
save(savepath1,'glycanListdiff');
save(savepath2,'glycanListcommon');

if(~isempty(varargin))
    outputfilename = varargin{1};
    outputfilename = [outputfilename '.xlsx'];
    excelfullpath = fullfile(loadpath, outputfilename);
    A1=cellstr('Made by Glycomics Analysis tool');
    A2=cellstr('MS1_unique glycan');
    B2=cellstr('MS2_unique glycan');
    C2=cellstr('Common glycan');
    xlswrite(excelfullpath,A1,1,'A1');
    xlswrite(excelfullpath,A2,1,'A2');
    xlswrite(excelfullpath,B2,1,'B2');
    xlswrite(excelfullpath,C2,1,'C2');
    xlswrite(excelfullpath,comparegroupdiff1,1,'A3');
    xlswrite(excelfullpath,comparegroupdiff2,1,'B3');
    xlswrite(excelfullpath, glycanListcommon,1,'C3');
end
end

function comparegroup = getcompos(newglycanDB)
expecGlycan = newglycanDB.expecGlycan;
comparegroup = {};
count = 0;
for i = 1 : length(expecGlycan)
    comparegroup{count+1,1}   = expecGlycan{i};
    count = count+1;
end
end