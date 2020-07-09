function StructResult = MSStructAnalysis(MSfilename,newglycanDB,listOfStruct,storepath)
% MSResidueAnalysis: Calculate relative abundance of different
% substructures for multi glycan profiles
%
% List of structure option: LeX(LewisX), sialylLeX(sialylLewisx),
% Ley(Lewisy), LacNAc, Sia_LacNAc(SialylLacNAc), Hantigen(H_Antigen),
% Aantigen(A_Antigen), Bantigen(B_Antigen), Bisecting, coreFuc.
%
% Author: Yusen Zhou
% Date Lastly Updated: 05/18/2020
StructResult = struct('MSID',[]);
% calculate 1st MS subStructure relativeabundance(lewis x_a, sialyl-Lewis x, 
%              Lewis y_b, H_Antigen, A_Antigen, B_Antigen, Bisecting, coreFuc)
StructResult.MSID = MSfilename;
for i = 1:length(listOfStruct)
    ithStruct = listOfStruct{i};
    StructResult.(ithStruct) = calculatebyStr(newglycanDB,ithStruct);
end
Storepath = [storepath MSfilename 'Struct.mat'];
save(Storepath,'StructResult');
end

function StructAbundance  = calculatebyStr(newglycanDB,StructStr)
StructAbundance = struct('glycanSpecies',[],'glycanStruct',[]);
glycanStructStr = newglycanDB.glycanexpec;
abundance       = newglycanDB.abundance;
SpAbundance = 0;
StAbundance = 0;
for i = 1 : length(glycanStructStr)
    ithglycan = glycanStructStr{i};
    ismatch = 0;
    if(ischar(ithglycan))
        [ismatch,NumberOfStruct] = calbyDiffStruct(ithglycan,StructStr,ismatch);
        if(ismatch~=0)
            SpAbundance = SpAbundance+abundance(i);
            StAbundance = StAbundance+abundance(i)*NumberOfStruct;
        end
    elseif(isa(ithglycan,'cell'))
        for j = 1 : length(ithglycan)
            str = ithglycan{j};
            [ismatch,NumberOfStruct] = calbyDiffStruct(str,StructStr,ismatch);
            if(ismatch~=0)
                SpAbundance = SpAbundance + ((abundance(i))/length(ithglycan));
                StAbundance = StAbundance + ((abundance(i)*NumberOfStruct)/length(ithglycan));
            end
        end
    end
end
SpAbundance = roundn(SpAbundance*100,-1);
StAbundance = roundn(StAbundance*100,-1);
StructAbundance.glycanSpecies = SpAbundance;
StructAbundance.glycanStruct  = StAbundance;
end

function [ismatch,NumberOfStruct] = calbyDiffStruct(ithglycan,StructStr,ismatch)
NumberOfStruct = 0;
if(strcmp(StructStr,'LeX'))
    if(contains(ithglycan,'{n{f}{h}}')||~isempty(regexp(ithglycan,'{n_b\d{f_a\d}{h_b\d}}','once')))
        ismatch = ismatch+1;
        NumberOfStruct = max(length(strfind(ithglycan,'{n{f}{h}}')),length(regexp(ithglycan,'{n_b\d{f_a\d}{h_b\d}}')));
    end
elseif(strcmp(StructStr,'sialylLeX'))
    if(contains(ithglycan,'{n{f}{h{s}}}')||~isempty(regexp(ithglycan,'{n_b\d{f_a\d}{h_b\d{s_a3}}}','once')))
        ismatch = ismatch+1;
        NumberOfStruct = max(length(strfind(ithglycan,'{n{f}{h{s}}}')),length(regexp(ithglycan,'{n_b\d{f_a\d}{h_b\d{s_a3}}}')));
    end
elseif(strcmp(StructStr,'Ley'))
    if(contains(ithglycan,'{n{f}{h{f}}}')||~isempty(regexp(ithglycan,'{n_b\d{f_a\d}{h_b\d{f_a2}}}','once')))
        ismatch = ismatch+1;
        NumberOfStruct = max(length(strfind(ithglycan,'{n{f}{h{f}}}')),length(regexp(ithglycan,'{n_b\d{f_a\d}{h_b\d{f_a2}}}')));
    end
elseif(strcmp(StructStr,'LacNAc'))
    % avoid bisecting structure.
    if(contains(ithglycan,'{n{h{n{')||contains(ithglycan,'{n{h}')...
            ||~isempty(regexp(ithglycan,'{n_b\d{h_b\d}','once'))||~isempty(regexp(ithglycan,'{n_b\d{h_b\d{n_b\d{','once')))
        ismatch = ismatch+1;
        NumberOfStruct = max(length(strfind(ithglycan,'{n{h{n{'))+length(strfind(ithglycan,'{n{h}'))...
            ,length(regexp(ithglycan,'{n_b\d{h_b\d}'))+length(regexp(ithglycan,'{n_b\d{h_b\d{n_b\d{')));
    end
elseif(strcmp(StructStr,'Sia_LacNAc'))
    if(contains(ithglycan,'{n{h{s}}}')||~isempty(regexp(ithglycan,'{n_b\d{h_b\d{s_a3}}}','once')))
        ismatch = ismatch+1;
        NumberOfStruct = max(length(strfind(ithglycan,'{n{h{s}}}')),length(regexp(ithglycan,'{n_b\d{h_b\d{s_a3}}}')));
    end
elseif(strcmp(StructStr,'Hantigen'))
    if(contains(ithglycan,'{n{h{f}}}')||~isempty(regexp(ithglycan,'{n_b\d{h_b\d{f_a2}}}','once')))
        ismatch = ismatch+1;
        NumberOfStruct = max(length(strfind(ithglycan,'{n{h{f}}}')),length(regexp(ithglycan,'{n_b\d{h_b\d{f_a2}}}')));
    end
elseif(strcmp(StructStr,'Aantigen'))
    if(contains(ithglycan,'{n{h{f}{n}}}')||~isempty(regexp(ithglycan,'{n_b\d{h_b\d{f_a2}{n_a3}}}','once')))
        ismatch = ismatch+1;
        NumberOfStruct = max(length(strfind(ithglycan,'{n{h{f}{n}}}')),length(regexp(ithglycan,'{n_b\d{h_b\d{f_a2}{n_a3}}}')));
    end
elseif(strcmp(StructStr,'Bantigen'))
    if(contains(ithglycan,'{n{h{f}{h}}}')||~isempty(regexp(ithglycan,'{n_b\d{h_b\d{f_a2}{h_a3}}}','once')))
        ismatch = ismatch+1;
        NumberOfStruct = max(length(strfind(ithglycan,'{n{h{f}{h}}}')),length(regexp(ithglycan,'{n_b\d{h_b\d{f_a2}{h_a3}}}')));
    end
elseif(strcmp(StructStr,'Bisecting'))
    if(contains(ithglycan,'{n{n{h{n}{h')||contains(ithglycan,'{n_b{n_b4{h_b4{n_b4}{h_a3')...
            ||contains(ithglycan,'{n{f}{n{h{n}{h')||contains(ithglycan,'{n_b{f_a6}{n_b4{h_b4{n_b4}{h_a3'))
        ismatch = ismatch+1;
        NumberOfStruct = max(length(strfind(ithglycan,'{n{n{h{n}{h'))+length(strfind(ithglycan,'{n{f}{n{h{n}{h'))...
            ,length(regexp(ithglycan,'{n_b{n_b4{h_b4{n_b4}{h_a3'))+length(regexp(ithglycan,'{n_b{f_a6}{n_b4{h_b4{n_b4}{h_a3')));
    end
elseif(strcmp(StructStr,'coreFuc'))
    if(contains(ithglycan,'{n{f}{n{h{h')||contains(ithglycan,'{n_b{f_a6}{n_b4{h_b4{h_a3')...
            ||contains(ithglycan,'{n{f}{n{h{n}{h')||contains(ithglycan,'{n_b{f_a6}{n_b4{h_b4{n_b4}{h_a3'))
        ismatch = ismatch+1;
        NumberOfStruct = max(length(strfind(ithglycan,'{n{f}{n{h{h'))+length(strfind(ithglycan,'{n{f}{n{h{n}{h'))...
            ,length(regexp(ithglycan,'{n_b{f_a6}{n_b4{h_b4{h_a3'))+length(regexp(ithglycan,'{n_b{f_a6}{n_b4{h_b4{n_b4}{h_a3')));
    end
end
end