function glycan1letstring=gly1charformat(glycanresiduestring,varargin)
% gly1charformat: Displaying glycan composition using 1 letter 
%
% Syntax:
% glycan1letstring=gly1charformat(glycanresiduestring)
%
% Input:
% glycanresiduestring: glycan composition string
%
% Output:
% glycan1letstring: glycan composition in 1 character format
%
%Example 1:
%      glycan3letterstring='Hex5HexNAc2'; 
%      glycan1letterstring=gly1charformat(glycan3letterstring);
%
% Author: Yusen Zhou
% Date Lastly Updated: 05/18/2020

glycan1letstring =[];

% HexNAC -h
nHexNAC  =  regexp(glycanresiduestring,'(HexNAc\d+)','match');
if(~isempty(nHexNAC))
    numHexNAc     =  str2num(nHexNAC{1,1}(7:end));
    for i = 1 : numHexNAc
        glycan1letstring = strcat(glycan1letstring,'n');
    end
end

% Hex
nHex       =  regexp(glycanresiduestring,'(Hex\d+)','match');
if(~isempty(nHex))
    numHex     =  str2num(nHex{1,1}(4:end));
    for i = 1 : numHex
        glycan1letstring = strcat(glycan1letstring,'h');
    end
end

% Sailic Acid NeuAc
% nNeuAc=findstr('s',gly);
nNeuAc             =  regexp(glycanresiduestring,'(NeuAc\d+)','match');
if(~isempty(nNeuAc))
    numNeuAc       =  str2num(nNeuAc{1,1}(6:end));
    for i = 1 : numNeuAc
        glycan1letstring = strcat(glycan1letstring,'s');
    end
    
end

% Sailic Acid NeuGc
% nNeuAc=findstr('g',gly);
nNeuGc              =  regexp(glycanresiduestring,'(NeuGc\d+)','match') ;
if(~isempty(nNeuGc))
    numNeuGc       =  str2num(nNeuGc{1,1}(6:end));
    for i = 1 : numNeuGc
        glycan1letstring = strcat(glycan1letstring,'g');
    end
end

% Fucose f
nFuc               =  regexp(glycanresiduestring,'(Fuc\d+)','match') ;
if(~isempty(nFuc))
    numFuc         =  str2num(nFuc{1,1}(4:end));
    for i = 1 : numFuc
        glycan1letstring = strcat(glycan1letstring,'f');
    end
end

% pentose
nPen           = regexp(glycanresiduestring,'(Pen\d+)','match') ;
if(~isempty(nPen))
    numPen         =  str2num(nPen{1,1}(4:end));
    for i = 1 : numPen
        glycan1letstring = strcat(glycan1letstring,'x');
    end
end

if(~isempty(varargin))
    name3letter = varargin{1};
    name1letter = varargin{2};
    pattern = ['(' name3letter '\d+)'];
    nMol               =  regexp(glycanresiduestring,pattern,'match') ;
    if(~isempty(nMol))
        numMol         =  str2num(nMol{1,1}(4:end));
        for i = 1 : numMol
            glycan1letstring = strcat(glycan1letstring,name1letter);
        end
    end
end

end

