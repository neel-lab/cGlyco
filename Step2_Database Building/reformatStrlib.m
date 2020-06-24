function [glycanDB,outputname,LibSize] = reformatStrlib(Strlib,N_Oglycan,massoption,storepath,varargin)
% reformatStrlib: Generate glycan structure data with four fields, i)SGP2.0 glycan structure, ii)glycan composition,
% iii)Monoisotopic mass and iv)isotopic distribution, based on the SGP2.0 list.
%
% Syntax:
% [glycanDB,outputname,LibSize] = reformatStrlib(Strlib,N_Oglycan,massoption,storepath,varargin)
%
% Input:
%   Strlib: glycan structure list in SGP2.0 format
%
%   N_Oglycan: 1 for N-linked glycan library and 0 for O-linked glycan
%   library.
%
%   massoption: Modification for monoshacchrides. For example:
%             massoption.acetylation = true;
%             massoption.methylation = false;
%             massoption.ion = 'Na';
%
%
%   storepath: Directory to store the results.
%
% Output:
%   glycanDB: Glycan data structure including 4 fields: i) glycan structure
%   list in SGP format, ii) glycan composition, iii)monoisotopic mass and iv)
%   isotopic distribution
%
%
%
% Author: Yusen Zhou
% Date Lastly Updated: 05/18/2020

outputfilename = [];
if(~isempty(varargin))
    outputfilename = varargin{1};
end
glycanDB = struct('glycanexpec',[],'expecGlycan',[],'monoisomw',[],'isomwarray',[],'abundance',[]);
if(N_Oglycan)
    outputname = 'NglycanDB';
else
    outputname = 'OglycanDB';
end
Strlibrary       = Strlib;
monoisomw        = '';
glycancomp       = '';
count            = 0;
GlycanMap        = containers.Map;
GlycanstrMap     = containers.Map;
for i = 1 : length(Strlibrary)
    ithStr = Strlibrary{i};
    ithcomposition = Str2Comp(ithStr);
    ithglycanMWt   = glycanMolWt(ithStr,'optionsstr',massoption);
    [isnewcomp,index]    = checkcomp(glycancomp,ithcomposition);
    if(isnewcomp)
        GlycanMap(num2str(ithglycanMWt))   = ithcomposition;
        GlycanstrMap(ithcomposition)       = ithStr;
        glycancomp{count+1,1}              = ithcomposition;
        monoisomw{count+1,1}               = ithglycanMWt;
        count = count+1;
    else
        glystrwithsamecomp = GlycanstrMap(glycancomp{index});
        if(ischar(glystrwithsamecomp))
            Strlist = cell(2,1);
            samecompstr = glystrwithsamecomp;
            Strlist{1,1} = samecompstr;
            Strlist{2,1} = ithStr;
            GlycanstrMap(glycancomp{index}) = Strlist;
        elseif(isa(glystrwithsamecomp,'cell'))
            number = length(glystrwithsamecomp);
            glystrwithsamecomp{number+1,1}  = ithStr;
            GlycanstrMap(glycancomp{index}) = glystrwithsamecomp;
        end
    end
end
[expecGlycan,monoisomw] = sortfromlarge2small(monoisomw,GlycanMap);
[expecGlycan,glycanexpec,monoisomw] = resortbymass(monoisomw,expecGlycan,GlycanstrMap);
glycanDB.glycanexpec    = glycanexpec;
glycanDB.expecGlycan    = expecGlycan;
glycanDB.monoisomw      = monoisomw;
glycanmwarray = Calisomwarray(expecGlycan,massoption);
glycanDB.isomwarray = glycanmwarray;
LibSize   = length(expecGlycan);
savefile1 = [storepath outputname];
save(savefile1,'glycanDB');
% export to excel file
if(~isempty(outputfilename))
    excelfullpath = [storepath outputfilename '.xlsx'];
    A1=cellstr('Composition');
    B1=cellstr('Monoisotopic mass');
    xlswrite(excelfullpath,expecGlycan,1,'A2');
    xlswrite(excelfullpath, monoisomw,1,'B2');
    xlswrite(excelfullpath,A1,1,'A1');
    xlswrite(excelfullpath,B1,1,'B1');
end
end

function Composition = Str2Comp(String)
Composition    = '';

NeuAcposition  = strfind(String,'s');
numofNeuAc     = length(NeuAcposition);
if(numofNeuAc~=0)
    Composition = [Composition 'NeuAc' num2str(numofNeuAc)];
end

NeuGcposition  = strfind(String,'g');
numofNeuGc     = length(NeuGcposition);
if(numofNeuGc~=0)
    Composition = [Composition 'NeuGc' num2str(numofNeuGc)];
end

Fucposition = strfind(String,'f');
numofFuc    = length(Fucposition);
if(numofFuc~=0)
    Composition = [Composition 'Fuc' num2str(numofFuc)];
end

Hexposition    = strfind(String,'h');
numofNeHex     = length(Hexposition);
if(numofNeHex~=0)
    Composition = [Composition 'Hex' num2str(numofNeHex)];
end

HexNAcposition    = strfind(String,'n');
numofNeHexNAc     = length(HexNAcposition);
if(numofNeHexNAc~=0)
    Composition = [Composition 'HexNAc' num2str(numofNeHexNAc)];
end
end

function [isnewcomp,index] = checkcomp(glycancomp,ithcomposition)
isnewcomp = 1;
index     = [];
for i = 1 : length(glycancomp)
    if(strcmp(glycancomp{i},ithcomposition))
        isnewcomp = 0;
        index     = i;
        break;
    end
end
end

function [expecGlycan,monoisomw] = sortfromlarge2small(monoisomw,GlycanMap)
expecGlycan = cell(length(monoisomw),1);
monodouble  = cell2mat(monoisomw);
monodouble  = sort(monodouble);
for i = 1 : length(monodouble)
   expecGlycan{i} = GlycanMap(num2str(monodouble(i)));
   monoisomw{i}   = monodouble(i);
end
end

function [expecGlycan,glycanexpec,monoisomw] = resortbymass(monoisomw,expecGlycan,GlycanstrMap)
for i = 1 : length(monoisomw)
    try
        ithmonoisomw     = monoisomw{i,1};
    catch 
        break
    end
    glycanexpec{i,1} = GlycanstrMap(expecGlycan{i,1});
    num = i; 
    glycanlist = cell(1,1);
    for j = num+1 : length(monoisomw)
        num = num+1;
        if(abs(monoisomw{num,1}-ithmonoisomw)<1)
            if(ischar(glycanexpec{i,1}))
                glycanlist{1,1} = glycanexpec{i,1};
            elseif(isa(glycanexpec{i,1},'cell'))
                glycanlist = glycanexpec{i,1};
            end
            if(ischar(GlycanstrMap(expecGlycan{num,1})))
                glycanlist{2,1} = GlycanstrMap(expecGlycan{num,1});
            elseif(isa(GlycanstrMap(expecGlycan{num,1}),'cell'))
                Length = length(glycanlist);
                for k = 1 : length(GlycanstrMap(expecGlycan{num,1}))
                    Length = Length+1;
                    kthglycanStr = GlycanstrMap(expecGlycan{num,1});
                    glycanlist{Length,1} = kthglycanStr{k};
                end
            end
            glycanexpec{i,1} = glycanlist;
            monoisomw(num)   = '';
            expecGlycan(num) = '';
            num = num-1;
        end
    end
end
end

function glycanmwarray = Calisomwarray(glycancomp,varargin)
if(~isempty(varargin))
    massoption = varargin{1};
end
if(~iscellstr(glycancomp))
    error('WRONG INPUT TYPE')
end
glycanstringarray     = cellfun(@gly1charformat,glycancomp,...
    'UniformOutput', false);
glycanformulaarray    = cellfun(@(x)glycanFormula(x,'optionsstr',massoption),glycanstringarray);
glycanmwarray         = arrayfun(@(x)isotopicdist(x,'SHOWPLOT',false),...
    glycanformulaarray,'UniformOutput', false);
end


