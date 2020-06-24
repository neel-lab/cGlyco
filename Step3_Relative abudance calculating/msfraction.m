function [newglycanDB,matchedpeakindex] = msfraction(peaklist,PWFH,glycanDB,OverSegmentationFilter,varargin)
% msfraction: peak assignment and relative abundance calculation
% 
% Syntax:
%  [newglycanDB,matchedpeakindex]= msfraction(peaklist,pfwhh,glycanDB,OverSegmentationFilter,varargin)
%
% Input:
%   peaklist: m/z and signal intensity data
%   PWFH: Full width at half maximum data
%   glycanDB: Candidate glycan structure data 
%   OverSegmentationFilter: Maximum distance between two adjacent peaks. 
%
% Output:
%   newglycanDB: glycan structure list including 9 fields: i) glycan structure
%   list in SGP format, ii) glycan composition, iii)monoisotopic mass, iv)
%   isotopic distribution, v) relative abundance, vi) assignment window, 
%   vii)isotopic deviation, viii)centroid mass and ix)residue value.
%
%   matchedpeakindex: Matched peak list in spectra
%
% Author: Yusen Zhou
% Date Lastly Updated: 05/18/2020

%% Input check
if(~isempty(varargin))
    fitOption = 1;
    fitPara   = varargin{1};
    standardV = 0;
else
    fitOption = 0;
    fitPara   = '';
    standardV = '';
end

%% Load data
expecGlycan     = glycanDB.expecGlycan;
glycanmwarray   = glycanDB.isomwarray;
newglycanStruct = glycanDB.glycanexpec;
newglycanDB     = glycanDB;

%% create variables
glycanindex         = 0;
windowindex         = 1;
allpeaktotalarea    = 0;
poitglycanIndexber  = 0;
monoisomw           = [];
deleteindex         = [];
matchedpeakindex    = [];
relativeabundance   = [];
finalAssignmentList = [];
WindowMass          = [];
WindowIDvalue       = [];
Residuemz           = [];

%% Peak assignment and calculation
for i=1:length(expecGlycan)
    ithmddouble = '';
    if(poitglycanIndexber>=i)
        continue
    end
    ithmddouble{1} = glycanmwarray{i,1};
    % Step 1: Define Assignment Window
    [AssignWindow,overlapindexArray,ithmddouble] =...
        extractAssignWindow(i,ithmddouble,expecGlycan,OverSegmentationFilter,glycanmwarray);
    % Step 2: Modify AssignWindow based on the expt peak intensity. Check the glycan existence and match monoisotopic peak
    [isglycaninMS,monopeak,finalAssignWindow,finalithmddouble,finalIndexList,ID] = ...
        matchglycaninMS(ithmddouble,AssignWindow,overlapindexArray,peaklist,OverSegmentationFilter,Residuemz);
    if(~isglycaninMS)
        deleteindex(end+1) = i;
        continue
    end
    % Step 3: Relative abundance calculation
    % Calculate the peak area
    [peakarea,matchedpeakindex,isexist,IDvalue,centroidMass,targetWindow,targetMono,Residuemz] =...
        findoverlaparea(finalithmddouble,monopeak,ID,peaklist,PWFH,matchedpeakindex,...
        finalAssignWindow,finalIndexList,Residuemz);
    if(~isexist)
        deleteindex(end+1) = i;
        continue
    end
    finalAssignmentList{windowindex,1} = targetWindow;
    WindowMass{windowindex,1}          = peaklist(targetMono,1);
    WindowMass{windowindex,2}          = centroidMass;
    WindowIDvalue{windowindex,1}       = IDvalue;
    windowindex = windowindex+1;
    
    % Calcualte total peak area and modify the results
    [allpeaktotalarea,monoisomw,glycanindex,poitglycanIndexber,relativeabundance,...
        deleteindex,standardV] = ResultModifying(i,peakarea,fitOption,relativeabundance,...
        allpeaktotalarea,monoisomw,glycanindex,fitPara,standardV,deleteindex,glycanmwarray);
end

% Delete unmatched glcyans
expecGlycan(deleteindex)     = '';
newglycanStruct(deleteindex) = '';
glycanmwarray(deleteindex)   = '';  

% Calculate relative abundance
for i = 1 : length(expecGlycan)
    relativeabundance(i,1)   = relativeabundance(i,1)/allpeaktotalarea;
end

%% Store Results
newglycanDB.abundance           = relativeabundance;
newglycanDB.monoisomw           = monoisomw;
newglycanDB.glycanexpec         = newglycanStruct;
newglycanDB.isomwarray          = glycanmwarray;
newglycanDB.expecGlycan         = expecGlycan;
newglycanDB.finalAssignmentList = finalAssignmentList;
newglycanDB.WindowMass          = WindowMass;
newglycanDB.WindowIDvalue       = WindowIDvalue;
newglycanDB.Residue             = Residuemz;
end

function [AssignWindow,overlapindexArray,ithmddouble] =...
    extractAssignWindow(startindex,ithmddouble,expecGlycan,OverSegmentationFilter,glycanmwarray)
targetmddouble       = ithmddouble{1};
counter              = 1;
overlapindexArray{1} = 1;
% 1st criteria: the Theoretical relative abundance of last peak in
% assignment window should be > 5%.
lastpeakIndex = find(targetmddouble(:,2)-0.1>0);
lastpeakIndex = 1 : 1 : lastpeakIndex(length(lastpeakIndex));
AssignWindow  = zeros(length(lastpeakIndex),1);
for j = 1 : length(lastpeakIndex)
    AssignWindow(j)  = targetmddouble(lastpeakIndex(j),1);
end
% Check if exists glycan overlapping.
for j = startindex+1 : length(expecGlycan)
    jthmddouble       = glycanmwarray{j,1};
    jthmonoisomass    = jthmddouble(1,1);
    if(jthmonoisomass<=targetmddouble(length(AssignWindow),1)+OverSegmentationFilter)
        ithmddouble{counter+1} = jthmddouble;
        jthlastpeakIndex       = find(jthmddouble(:,2)-0.1>0);
        jthlastpeakIndex       = 1 : 1 : lastpeakIndex(length(jthlastpeakIndex));
        jthAssignWindow        = zeros(length(jthlastpeakIndex),1);
        for k = 1 : length(jthlastpeakIndex)
            jthAssignWindow(k)  = jthmddouble(jthlastpeakIndex(k),1);
        end
        [AssignWindow,overlappeakindex]  = combineAssigWindow(AssignWindow,jthAssignWindow,OverSegmentationFilter);
        overlapindexArray{counter+1} = overlappeakindex;
        counter = counter+1;
    else
        break
    end
end
end

function [AssignWindow,overlappeakindex] = combineAssigWindow(AssignWindow,SecondAssignWindow,OverSegmentationFilter)
startindex                = 1;
[overlappeakindex,~]      = findclosetpeak(SecondAssignWindow(1),AssignWindow,OverSegmentationFilter,startindex);
glycanIndexofAssignWindow = length(SecondAssignWindow)+overlappeakindex-1;
startnum = length(AssignWindow)+1;
for i = startnum:glycanIndexofAssignWindow
    AssignWindow = [AssignWindow; SecondAssignWindow(i-overlappeakindex+1,1)];
end
end

function [isotopicindex,allpeak] = findclosetpeak(isotopicmz,mzlist,OverSegmentationFilter,startIndex)
checkisotopicmz        = abs(mzlist(startIndex:length(mzlist),1)-isotopicmz)<OverSegmentationFilter*1.5;
allpeak                = find(checkisotopicmz)+startIndex-1;
% if multiple peaks are found close to isotopic mass, select the closet peak as the match
if(length(allpeak)>1)
    isotopicmzsdif         = abs(mzlist(allpeak)-isotopicmz);
    isotopicmultipleindex  = allpeak;
    isotopicnewindex       = ~abs(isotopicmzsdif-min(isotopicmzsdif));
    isotopicindex          = isotopicmultipleindex(isotopicnewindex);
else
    isotopicindex = allpeak;
end
end

function [isglycaninMS,monopeak,finalAssignWindow,finalmddouble,finalIndexList,ID]...
    = matchglycaninMS(mddouble,AssignWindow,overlapindexArray,peaklist,OverSegmentationFilter,Residuemz)
isglycaninMS      = 0;
monopeak          = [];
finalAssignWindow = '';
finalIndexList    = '';
finalmddouble     = '';
% Two situations: i) only one glycan; ii) multiple glycan
mz            = mddouble{1}(:,1);
TheoRA        = mddouble{1}(:,2);
firstmz       = mz(1,1);
if(length(overlapindexArray)==1)
    ismultiglycan = 0;
    [monopeakindex,ID] = findpeak(firstmz,peaklist,TheoRA,AssignWindow,ismultiglycan,OverSegmentationFilter,Residuemz);
else
    ismultiglycan = 1;
    [monopeakindex,ID] = findpeak(firstmz,peaklist,TheoRA,AssignWindow,ismultiglycan,OverSegmentationFilter,Residuemz);
end
if(~isempty(monopeakindex))
    [AssignWindowList,mddoubleList,starIndexList] =...
        modifyAssignWindow(AssignWindow,monopeakindex,overlapindexArray,peaklist,mddouble);
    
    % Isotopic mass distribution can not distance too far away.
    monopeak = [];
    for i = 1 : length(monopeakindex)
        ithAssignWindow = AssignWindowList{i};
        peakindex = monopeakindex(i);
        count     = 1;
        for j = 2 : length(ithAssignWindow) 
            peakindex = peakindex+1;
            if((roundn(peaklist(peakindex,1),-1)-roundn(peaklist(peakindex-1,1),-1))<=OverSegmentationFilter*1.2)
                count = count+1;
            else
                break
            end
        end
        if(ismultiglycan)
            if(firstmz<3000)
                Minpeaknumber = 3;
            else
                Minpeaknumber = 4;% MinIsopeak(1);
            end
            if(count>=Minpeaknumber)
                monopeak   = [monopeak monopeakindex(i)];
                finalithmddouble = mddoubleList{i};
                finalmddouble{end+1}     = finalithmddouble;
                finalAssignWindow{end+1} = AssignWindowList{i}(1:count);
                finalIndexList{end+1}    = starIndexList{i};
            end
        else
            if(firstmz<3000)
                Minpeaknumber = 3;
            else
                Minpeaknumber = 4;%length(ithAssignWindow); 
            end
            if(count>=Minpeaknumber)
                monopeak = [monopeak monopeakindex(i)];
                finalmddouble{end+1}     = mddoubleList{i};
                finalAssignWindow{end+1} = AssignWindowList{i};
                finalIndexList{end+1}    = starIndexList{i};
            end
        end
    end
    if(~isempty(monopeak))
        isglycaninMS = 1;
    end
end
end

function [peakindex,ID] = findpeak(peakmz,peaklist,TheoRA,AssignWindow,ismultiglycan,OverSegmentationFilter,Residuemz)
peakindex = [];
ID        = 0;
if(~isempty(Residuemz))
    [col,~]  = size(Residuemz);
    StartIndex = Residuemz(col,1)+1;
else
    StartIndex = 1;
end
[isotopicindex,allpeak] = findclosetpeak(peakmz,peaklist,OverSegmentationFilter,StartIndex);
% Only one candidate monoisotopic peak matched.
if(length(allpeak)==1)
    glycanIndexofPeak    = length(AssignWindow);
    CompareRatio = 0;
    if(isotopicindex+glycanIndexofPeak-1<=length(peaklist))
        startpeak = isotopicindex;
        for i = 1 : glycanIndexofPeak-1
            ithPeakratio = peaklist(startpeak,2)/peaklist(startpeak+1,2);
            ithTheoratio = TheoRA(i,1)/TheoRA(i+1,1);
            Calratio = (log10(ithPeakratio/ithTheoratio)*(TheoRA(i,1)+TheoRA(i+1,1)))^2;
            CompareRatio = CompareRatio+Calratio;
            startpeak = startpeak+1;
        end
        CompareRatio = CompareRatio^0.5;
    else
        CompareRatio = inf;
    end
    if(ismultiglycan)
        peakindex = isotopicindex;
    else
        Tratio = 0;
        for i = 1 : glycanIndexofPeak-1
            Theoratio = (log10(1.5)*(TheoRA(i,1)+TheoRA(i+1,1)))^2;
            Tratio = Tratio+Theoratio;
        end
        Tratio = Tratio^0.5;
        if(CompareRatio<Tratio)
            peakindex = isotopicindex;
            ID        = CompareRatio;
        end
    end
% >1 candidate monoisotopic peaks matched.
elseif(length(allpeak)>1)
    if(~ismultiglycan)
        [firstpeakindex,IDvalue] = choosemono(allpeak,peaklist,TheoRA,AssignWindow);
        peakindex = firstpeakindex;
        ID        = IDvalue;
    else
        peakindex      = allpeak;
    end
end
end

function [firstpeakindex,IDvalue] = choosemono(allpeak,peaklist,TheoRA,AssignWindow)
firstpeakindex    = [];
IDvalue           = 0;
Tratio            = 0; 
CompareRatio      = [];
glycanIndexofPeak = length(AssignWindow);
MaxTheoRA = [];
for i = 1 : glycanIndexofPeak-1
    MaxTheoRA = [MaxTheoRA TheoRA(i,1)+TheoRA(i+1,1)];
end
StdRA  = max(MaxTheoRA);
Tratio = (glycanIndexofPeak-1)*((log10(2))^2);
Tratio = Tratio^0.5;
for i = 1 : length(allpeak)
    ithoeak   = allpeak(i);
    ratio     = 0;
    if(ithoeak+glycanIndexofPeak-1<=length(peaklist))
        startpeak = ithoeak;
        for j = 1 : glycanIndexofPeak-1
            ithPeakratio = peaklist(startpeak,2)/peaklist(startpeak+1,2);
            ithTheoratio = TheoRA(j,1)/TheoRA(j+1,1);
            Calratio  = (log10(ithPeakratio/ithTheoratio)*((TheoRA(j,1)+TheoRA(j+1,1))/StdRA))^2;
            ratio = ratio+Calratio;
            startpeak = startpeak+1;
        end
        CompareRatio(i) = ratio^0.5;
    else
        CompareRatio(i) = inf;
    end
end
if(min(CompareRatio)<Tratio)
    [IDvalue,index]       = min(CompareRatio);
    firstpeakindex = allpeak(index);
end
end

function [AssignWindowList,mddoubleList,starIndexList] =...
    modifyAssignWindow(AssignWindow,monopeak,overlapindexArray,peaklist,mddouble)
AssignWindowList = cell(length(monopeak),1);
mddoubleList     = cell(length(monopeak),1);
starIndexList    = cell(length(monopeak),1);
for i = 1 :length(monopeak)
    peakindex         = monopeak(i);
    maxIntensity      = peaklist(peakindex,2);
    ithAssignWindow   = [];
    for j = 1 : length(AssignWindow)
        if(peakindex<=length(peaklist))
            if((peaklist(peakindex,2)/maxIntensity)>0.15)
                ithAssignWindow(end+1,1) = peaklist(peakindex,1);
                maxIntensity = max(maxIntensity,peaklist(peakindex,2));
                peakindex    = peakindex+1;
            else
                break
            end
        end
    end
    % Recheck if multiple glycan exist in the new Assignment Window
    ithmddouble = mddouble;
    glycanIndexofGlycan = length(ithmddouble);
    startIndex  = [];
    deleteindex = [];
    for j = 1 : glycanIndexofGlycan
        intialstart = overlapindexArray{j};
        if(intialstart>length(ithAssignWindow))
            deleteindex(end+1) = j;
        else
            startIndex = [startIndex intialstart];
        end
    end
    ithmddouble(deleteindex)= '';
    mddoubleList{i}      = ithmddouble;
    AssignWindowList{i}  = ithAssignWindow;
    starIndexList{i}     = startIndex;
end
end

function [allpeaktotalarea,monoisomw,index,poitglycanIndexber,relativeabundance,deleteindex,standardV] = ResultModifying(i,peakarea,fitOption,...
    relativeabundance,allpeaktotalarea,monoisomw,index,fitPara,standardV,deleteindex,glycanmwarray)
glycanindex = i;
for j = 1 : length(peakarea)
    jthpeakarea = peakarea{j};
    if(jthpeakarea==0)
        deleteindex(end+1) = glycanindex{j};
        continue
    end
    if(fitOption)
        if(iscell(fitPara)) % adjust signal intensity by noise.
            if(index == 0)
                standardV = fitPara{1}*exp(fitPara{2}*glycanmwarray{glycanindex,1}(1,1))+...
                    fitPara{3}*exp(fitPara{4}*glycanmwarray{glycanindex,1}(1,1));
            end
            jthV = fitPara{1}*exp(fitPara{2}*glycanmwarray{glycanindex,1}(1,1))+...
                fitPara{3}*exp(fitPara{4}*glycanmwarray{glycanindex,1}(1,1));
        elseif(isstruct(fitPara)) % adjust signal intensity by user-defined function.
            ZoneBoundary = fitPara.Zone;
            Zone1Func = fitPara.Zone1Func;
            Zone2Func = fitPara.Zone2Func;
            Z1Func = str2func(Zone1Func);
            Z2Func = str2func(Zone2Func);
            m2zratio = glycanmwarray{glycanindex,1}(1,1);
            if(index == 0)
                standardV = feval(Z1Func,m2zratio);
            end
            if(m2zratio<=ZoneBoundary)
                jthV = feval(Z1Func,m2zratio);
            elseif(m2zratio>ZoneBoundary)
                jthV = feval(Z2Func,m2zratio);
            end
        end
        jthpeakarea                  = jthpeakarea*(standardV/jthV);
        relativeabundance(index+1,1) = jthpeakarea;
    else
% Linear signal reduction
        jthpeakarea = jthpeakarea*(1+(glycanmwarray{glycanindex,1}(1,1)-1500)/1000);
        relativeabundance(index+1,1) = jthpeakarea;
    end
    allpeaktotalarea      = allpeaktotalarea+jthpeakarea;
    monoisomw{index+1,1}  = glycanmwarray{glycanindex,1}(1,1);
    index                 = index+1;
    glycanindex           = glycanindex+1;
end
poitglycanIndexber = i+length(peakarea)-1;
end
