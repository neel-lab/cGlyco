function Sgpout = IUPAC2Sgp(IUPAC)
%% IUPAC2Sgp: Convert IUPAC to SGP format.
%
%  Syntax:
%     Sgpout = IUPAC2Sgp(IUPAC)
%     Sgpout = IUPAC2Sgp(IUPAC,version)
%
%  Input:
%       IUPAC: A type of glycan nomenclature in linear format.
%  
%  Output:
%     Sgpout: glycan in SGP format
%
%  Example: 
%     M5 = 'Man(a1-3)[Man(a1-6)]Man(a1-6)[Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-r)';
%     Sgpout = IUPAC2Sgp(M5);
%
%     Bibrabch = 'NeuAc(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-2)Man(a1-6)...
%               [NeuAc(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(b1-r)';
%     Sgpout = IUPAC2Sgp(Bibrabch);
%
%  Children function:
%     N/A
%
%  Author: Yusen Zhou
%  Date Lastly Updated: 03/30/16

%% Main function
if(isempty(IUPAC))
    Sgpout = '';
    return
end

% Split the whole sequence into several cells, each cell contains one
%   monosaccharide or bracket
[MonoStr,MonoIndex,BracStr,BracIndex] = Split(IUPAC);

% Convert to SGP format
Sgpout = Convert2Sgp(MonoStr,MonoIndex,BracStr,BracIndex);
end

function [MonoStr,MonoIndex,BracStr,BracIndex] = Split(IUPAC)
% 'Mono(b?-?)' or 'Monob-?' or 'Monob?'
MonoPattern = '[A-Z_a-z_0-9]*\((.*?)\)'; % consider using mass to present the monosaccharide
% '[' or ']'
BracketPattern = '[\[,\]]';
[MonoStr,MonoIndex] = regexp(IUPAC,MonoPattern,'Match','start');
[BracStr,BracIndex] = regexp(IUPAC,BracketPattern,'Match','start');
if(isempty(MonoIndex))
    MonoPattern = '[A-Z_a-z_0-9]*';
    [MonoStr,MonoIndex] = regexp(IUPAC,MonoPattern,'Match','start');
end
end

function Sgpout = Convert2Sgp(MonoStr,MonoIndex,BracStr,BracIndex)
Hex = {'Glc[^AN]','Man[^AN]','Gal[^AN]','Gul[^AN]','Alt[^AN]','All[^AN]',...
            'Tal[^AN]','Ido[^AN]'};
HexNAc = {'GlcNAc','ManNAc','GalNAc','GulNAc','AltNAc','AllNAc',...
            'TalNAc','IdoNAc'};
HexN = {'GlcN[^A]','ManN[^A]','GalN[^A]','GulN[^A]','AltN[^A]',...
            'AllN[^A]','TalN[^A]','IdoN[^A]'};
HexA = {'GlcA','ManA','GalA','GulA','AltA','AllA','TalA','IdoA'};
DHexNAc = {'QuiNAc','RhaNAc','FucNAc'};
DHex = {'Qui[^N]','Rha[^N]','6dAlt','6dTal','Fuc[^N]'};
DiDHex = {'Oli','Tyv','Abe','Par','Dig','Col'};
Pen = {'Ara','Lyx','Xyl','Rib'};
Othersin3 = {'Kdn','NeuGc','Neu5Gc','NeuAc','Neu5Ac','Neu','Bac','LDMan',...
    'Kdo','Dha','DDMan','MurNAc','MurNGc','Mur','Api','Fru',...
    'Tag','Sor','Psi'};
Othersin1 = {'k','g','g','s','s','u','b','L',...
    'o','y','D','A','M','m','i','c','t','r','i'};
Sgpout     = [];
WholeSeq   = sort(cat(2,MonoIndex,BracIndex));
Startpoint = length(WholeSeq);
BracLvl    = 0;
Brackpos   = zeros(length(BracIndex)/2,1);
for i = 1:length(WholeSeq)
    ithIndex = WholeSeq(Startpoint);
    if(~isempty(find(ithIndex-MonoIndex==0,1)))
        ithCellStr = MonoStr{find(ithIndex-MonoIndex==0,1)};
        % Linkage info.
        linkageM = regexp(ithCellStr,'\(.*\)', 'match');
        linkage  = linkageM{1};
        if(strcmp(linkage(5),'r'))
            SgpLink = ['_' linkage(2)]; 
        elseif(strcmp(linkage(5),'?'))
            SgpLink = '';
        else
            SgpLink = ['_' linkage(2) linkage(5)]; 
        end
        if(sum(cellfun(@(x)matchMono(x,ithCellStr),Hex)==1))
            Sgpunit = ['{h' SgpLink];
            Sgpout  = [Sgpout Sgpunit];
        elseif(sum(cellfun(@(x)matchMono(x,ithCellStr),HexNAc)==1))
            Sgpunit = ['{n' SgpLink];
            Sgpout  = [Sgpout Sgpunit];
        elseif(sum(cellfun(@(x)matchMono(x,ithCellStr),HexN)==1))
            Sgpunit = ['{e' SgpLink];
            Sgpout  = [Sgpout Sgpunit];
        elseif(sum(cellfun(@(x)matchMono(x,ithCellStr),HexA)==1))
            Sgpunit = ['{a' SgpLink];
            Sgpout  = [Sgpout Sgpunit];
        elseif(sum(cellfun(@(x)matchMono(x,ithCellStr),DHexNAc)==1))
            Sgpunit = ['{N' SgpLink];
            Sgpout  = [Sgpout Sgpunit];
        elseif(sum(cellfun(@(x)matchMono(x,ithCellStr),DHex)==1))
            Sgpunit = ['{f' SgpLink];
            Sgpout  = [Sgpout Sgpunit];
        elseif(sum(cellfun(@(x)matchMono(x,ithCellStr),DiDHex)==1))
            Sgpunit = ['{d' SgpLink];
            Sgpout  = [Sgpout Sgpunit];
        elseif(sum(cellfun(@(x)matchMono(x,ithCellStr),Pen)==1))
            Sgpunit = ['{x' SgpLink];
            Sgpout  = [Sgpout Sgpunit];
        elseif(sum(cellfun(@(x)matchMono(x,ithCellStr),Othersin3)==1))
            MonoId = find(cellfun(@(x)matchMono(x,ithCellStr),Othersin3));
            SgpSymbol = Othersin1{MonoId};
            Sgpunit   = ['{' SgpSymbol SgpLink];
            Sgpout    = [Sgpout Sgpunit];
        end
    else
        BracPos   = find(ithIndex-BracIndex==0,1);
        BracSym   = BracStr(BracPos);
        if(strcmp(BracSym,']'))
            startpos = find(ithIndex-WholeSeq==0,1);
            BracLvl  = BracLvl+1;
            Brackpos(BracLvl) = length(Sgpout)+1;
        else
            StrLvl       = Sgpout(Brackpos(BracLvl):end);
            MonoStart    = length(strfind(StrLvl,'{'));
            Monoend      = length(strfind(StrLvl,'}'));
            Mononum      = MonoStart - Monoend;
            for j = 1 : Mononum
                Sgpout       = [Sgpout '}'];
            end
            BracLvl      = BracLvl-1;
            if(BracLvl==0)
                Brackpos   = zeros(length(BracIndex)/2,1);
            end
        end
    end
    Startpoint = Startpoint-1;
end
MonoStart    = length(strfind(Sgpout,'{'));
Monoend      = length(strfind(Sgpout,'}'));
for i = 1 : MonoStart-Monoend
    Sgpout  = [Sgpout '}'];
end
end


function isMatch = matchMono(Mono,Monolist)
isMatch = 0;
if(~isempty(regexp(Monolist,Mono,'match')))
    isMatch = 1;
end
end
