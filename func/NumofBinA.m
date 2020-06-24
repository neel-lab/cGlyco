function numOfsubst = NumofBinA(structA,structB,MonoIdx,ver,varargin)
%% NumofBinA: To calculate how many structB are formed in structA.
%
%  Syntax:
%     numOfsubst = NumofBinA(structA,structB,MonoIdx,ver)
%
%  Input:
%     structA:  glycan structure string A
%     structB:  Sub-structure string B
%     MonoIdx:  Unique ID number for different monosaccharides.
%     ver    :  '0', normal comparison
%               '1', only used for calculating number of LacNAc extesion in one branch.
%
%  Output:
%     numOfsubst:  Number of structB in structA
%
%  Example:
%     MonoIdx = MonoIDgeneration;
%     structA = '{(b1-r)GlcNAc{(b1-4)GlcNAc{(b1-4)Man[{(a1-3)Man{(b1-2)GlcNAc{(b1-4)Gal{(a2-3)NeuAc}}}}]{(a1-6)Man{(b1-2)GlcNAc[{(a1-3)Fuc}]{(b1-4)Gal{(a2-3)NeuAc}}}}}}}';
%     structB = '{(b1-?)GlcNAc{(b1-4)Gal}}';
%     numOfsubst = NumofBinA(structA,structB,MonoIdx,0)
%
%  Children functions:
%     N/A.
%
%  Author: Yusen Zhou
%  Date Lastly Updated: 11/2/16

%% Main function
SulfoGp = 0;
if(~isempty(varargin))
    SulfoGp = varargin{1};
end

numOfsubst = 0;
RepeatUnits = [];
% Create structure Map
% Branch number
BranchA = numel(regexp(structA,']'))+1;
BranchB = numel(regexp(structB,']'))+1;
numOfMonoA = numel(regexp(structA,'\('));
numOfMonoB = numel(regexp(structB,'\('));
[ALinkMap,AMonoMap,structALvl] = createStructMap(structA,numOfMonoA,MonoIdx);
[BLinkMap,BMonoMap,structBLvl] = createStructMap(structB,numOfMonoB,MonoIdx);
[BranchAMap,LinkAMap,IDAMap] = createBranchMap(BranchA,AMonoMap,ALinkMap,structALvl,SulfoGp);
[BranchBMap,LinkBMap,~] = createBranchMap(BranchB,BMonoMap,BLinkMap,structBLvl,SulfoGp);

% Compare two structure map
StartNum = cell(BranchB,1);
isMatch  = 0;
for i = 1 : BranchB
    ithBBranch = BranchBMap(i,:);
    ithBBranch = ithBBranch(ithBBranch~=0);
    ithBLink   = LinkBMap(i,:);
    ithBLink   = ithBLink(ithBBranch~=0);
    NumofSub   = length(ithBBranch(1,:));
    for j = 1 : BranchA
        jthABranch = BranchAMap(j,:);
        StartLoc   = find(jthABranch==ithBBranch(1));
        MatchedNum = 0;
        for k = 1 : length(StartLoc)
            kthStart = StartLoc(k);
            if(kthStart+NumofSub-1>length(jthABranch(1,:)))
                continue;
            end
            Target   = jthABranch(1,kthStart:kthStart+NumofSub-1);
            isEqual  = isequal(Target,ithBBranch);           
            if(isEqual)
                % Check if linkage matches
                ithALink     = LinkAMap(j,kthStart:kthStart+NumofSub-1);
                unmatched    = ithALink==ithBLink;
                unmatchedID  = ~unmatched;
                unmatchedPos = ithBLink(unmatchedID);
                if(~any(unmatchedPos~=9))
                    StartNum{i,1} = [StartNum{i,1} IDAMap(j,kthStart)];
                    MatchedNum = MatchedNum+1;
                end
            end
        end
        if(ver)
            RepeatUnits = [RepeatUnits MatchedNum];
        end
    end
    if(~isempty(StartNum{i}))
        isMatch = isMatch+1;
    end
end
if(~ver)
    if(isMatch==BranchB)
        if(BranchB==1)
            Samepart = unique(StartNum{1,1});
        else
            Samepart = unique(StartNum{1,1});
            [StartSize,~] = size(StartNum);
            for i = 1 : StartSize
                Samepart = intersect(Samepart,StartNum{i,1});
            end
        end
        numOfsubst = numOfsubst + length(Samepart);
    end
else
    numOfsubst = max(RepeatUnits);
end
end

function [LinkMap,MonoMap,structLvl] = createStructMap(structure,numOfMono,MonoIdx)
[UnitBlock,index]  = regexp(structure,'[{}\[\]]','Match');
LinkMap   = cell(numOfMono,1);
MonoMap   = zeros(numOfMono,1);
structLvl  = zeros(numOfMono,1);
nodelvl    = 0;
MonoID     = 0;
for i = 1 : numel(UnitBlock)
    ithBlock = UnitBlock{i};
    if(ithBlock=='{')
        nodelvl     = nodelvl+1;
        MonoID      = MonoID+1;
        ithIndex    = index(i);
        nextIndex   = index(i+1);
        MononameStr = structure(ithIndex+1:nextIndex-1);
        structLvl(MonoID) = nodelvl;
        LinkEndPos        = strfind(MononameStr,')');
        LinkMap{MonoID,1} = MononameStr(1:LinkEndPos);
        MonoMap(MonoID,1) = MonoIdx.(MononameStr(LinkEndPos+1:end));
    elseif(ithBlock=='}')
        nodelvl = nodelvl-1;
    end
end
end

function [BranchIdx,LinkIdx,IDIdx] = createBranchMap(Branch,MonoMap,LinkMap,structLvl,SulfoGp)
BranchIdx = zeros(Branch,max(structLvl));
LinkIdx = zeros(Branch,max(structLvl));
IDIdx   = zeros(Branch,max(structLvl));
BranchNum = 1;
for i = 1 : numel(structLvl)
    ithMonoName = MonoMap(i);
    ithLink     = LinkMap{i};
    ParentPos   = ithLink(5);
    Config      = ithLink(2);
    if(strcmp(Config,'a'))
        ConfigV = 1;
    else
        ConfigV = 2;
    end
    if(strcmp(ParentPos,'r'))
        Pos = 7;
    elseif(strcmp(ParentPos,'1'))
        Pos = 1;
    elseif(strcmp(ParentPos,'2'))
        Pos = 2;
    elseif(strcmp(ParentPos,'3'))
        Pos = 3;
    elseif(strcmp(ParentPos,'4'))
        Pos = 4;
    elseif(strcmp(ParentPos,'5'))
        Pos = 5;
    elseif(strcmp(ParentPos,'6'))
        Pos = 6;
    elseif(strcmp(ParentPos,'8'))
        Pos = 8;
    elseif(strcmp(ParentPos,'?'))
        Pos = 9;
    end
    % Check sulfation
    SulPos = [];
    if(SulfoGp)
        SStr = regexp(ithLink,'S-\d','match');
        if(~isempty(SStr))
            numofSul = numel(SStr);
            for j = 1:numofSul
                ithSStr   = SStr{j};
                SulLink   = ithSStr(3);
                SulPos(j) = str2double(SulLink);
                SulPos = sort(SulPos);
            end
            for j = 1:numofSul
                Pos = Pos*10+SulPos(j);
            end
        end
    end
    FinalNum  = ithMonoName*10+ConfigV;
    ithLvl = structLvl(i);
    if(i==1)
        BranchIdx(BranchNum,ithLvl) = FinalNum;
        LinkIdx(BranchNum,ithLvl)   = Pos;
        IDIdx(BranchNum,ithLvl)     = i;                 
        continue
    end
    preithLvl = structLvl(i-1);
    if(ithLvl<=preithLvl)
        BranchNum = BranchNum+1;
    end
    BranchIdx(BranchNum,ithLvl) = FinalNum;
    LinkIdx(BranchNum,ithLvl)   = Pos;
    IDIdx(BranchNum,ithLvl)     = i;        
end

for i = 2:Branch
    ithBranch = IDIdx(i,:);
    index     = find(ithBranch, 1 );
    preBranch = IDIdx(i-1,:);
    preBranch(index:end) = 0;
    IDIdx(i,:) = ithBranch+preBranch;
    ithBranch = LinkIdx(i,:);
    preBranch = LinkIdx(i-1,:);
    preBranch(index:end) = 0;
    LinkIdx(i,:) = ithBranch+preBranch;
    ithBranch = BranchIdx(i,:);
    preBranch = BranchIdx(i-1,:);
    preBranch(index:end) = 0;
    BranchIdx(i,:) = ithBranch+preBranch;
end
end
