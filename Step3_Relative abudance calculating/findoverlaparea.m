function  [glycanarea,matchedpeakindex,isexist,IDvalue,centroidMass,targetWindow,targetMono,Residuemz]= findoverlaparea(mdList,...
    monopeak,ID,peaklist,FWHM,matchedpeakindex,finalAssignWindow,starIndexList,Residuemz)
% findoverlaparea: find peak area when there are multi-glycans
%
% See also msfraction.
%
% Author: Yusen Zhou
% Data Lastly Updated: 05/20/2020
%% Initiate the parameters
isexist         = 0;
IDvalue         = [];
centroidMass    = [];
glycanarea      = [];
Fvalue          = [];
TheoAreaList    = cell(length(mdList),1);
exptAreaList    = cell(length(mdList),1);
TotalAreaList   = cell(length(mdList),1);

%% Calculate part
for i = 1 : length(mdList)
    ithmd           = mdList{i};
    ithmonopeak     = monopeak(i);
    ithAssignWindow = finalAssignWindow{i};
    startIndex      = starIndexList{i};
    % Step1 Calcualte experimental area
    [exptPeakarea,Totalarea] = CalExptArea(ithmonopeak,ithAssignWindow,FWHM,peaklist);
    
    % Step2 Calcualte theoretical area
    TheoArea = CalInitialGuess(ithmd,ithmonopeak,FWHM,peaklist,startIndex);
    
    % minimize function
    [TheoArea,FVAL] = fmin(ithmd,TheoArea,exptPeakarea,Totalarea,startIndex);
    TheoAreaList{i}  = TheoArea;
    exptAreaList{i}  = exptPeakarea;
    TotalAreaList{i} = Totalarea;
    Fvalue           = [Fvalue FVAL];
end
% Find minmum FVAL if >1 candidate monoisotopic peak
[~,pos]      = min(Fvalue);
TheoArea     = TheoAreaList{pos};
md           = mdList{pos};
targetMono   = monopeak(pos);
startIndex   = starIndexList{pos};
targetExpt   = exptAreaList{pos};
targetTotal  = TotalAreaList{pos};
targetWindow = finalAssignWindow{pos};
% centroid mass calculation and ID value verification.
[ismatch,glycantotalarea,Centroid,CompareRatio,Residuemz] = ResultsFit(md,TheoArea,targetExpt,...
    startIndex,targetMono,peaklist,ID,targetTotal,Residuemz);
if(ismatch)
    isexist          = 1;
    IDvalue          = CompareRatio;
    centroidMass     = Centroid;
    glycanarea       = glycantotalarea;
    matchedpeakindex = getMatchedpeaks(matchedpeakindex,targetWindow,targetMono);
    
end
end

function [exptPeakarea,Totalarea] = CalExptArea(monopeak,AssignWindow,FWHM,peaklist)
Totalarea    = 0;
startpeak    = monopeak;
exptPeakarea = [];
for i = 1 : length(AssignWindow)
    ithPeakarea     = abs(FWHM(startpeak,1)-FWHM(startpeak,2))*peaklist(startpeak,2);
    Totalarea       = Totalarea+ithPeakarea;
    exptPeakarea(i) = ithPeakarea;
    startpeak       = startpeak+1;
end
end

function [TheoArea,FVAL] = fmin(md,TheoArea,exptPeakarea,Totalarea,startIndex)
AssignWindow = length(exptPeakarea);
glycannum   = length(md);
glycancoeff = zeros(AssignWindow,glycannum);
glycanDis   = zeros(1,glycannum);
for i = 1 : glycannum
    dist = md{1,i};
    for j = startIndex(i):AssignWindow %min(AssignWindow,startIndex(i)+SingleAssignWindow(i)-1)
        glycancoeff(j,i) = dist(j-startIndex(i)+1,2);
    end
    glycanDis(1,i) = sum(glycancoeff(:,i));
end
funString = '@(x)';
for i = 1 : AssignWindow
    funString = [funString 'abs('];
    firstItem = num2str(exptPeakarea(i));
    funString = [funString firstItem];
    for j = 1 : glycannum
        jthcoeff = num2str(glycancoeff(i,j));
        variableIndex = num2str(j);
        jthvariable   = ['-' jthcoeff '*x(' variableIndex ')']; 
        funString     = [funString jthvariable];
    end
    funString = [funString ')'];
    if(i~=AssignWindow)
        funString = [funString '+'];
    end
end
fun = str2func(funString);
for i = 1 : glycannum
    x0(i) = TheoArea{i};
    A(i)  = glycanDis(1,i);%glycanDis(1,i)
    lb(i) = 0;
    ub(i) = Inf;
end
options = optimoptions('fmincon','Display','off');
[x,FVAL] = fmincon(fun,x0,[],[],A,Totalarea,lb,ub,[],options);%fmincon(fun,x0,[],[],A,Totalarea,lb,ub,[],options)
for i = 1 : glycannum
    TheoArea{i} = x(i);
end
end

function TheoArea = CalInitialGuess(md,monopeak,FWHM,peaklist,startIndex)
numofglycans    = length(md);
TheoArea        = cell(numofglycans,1);
% Calculate initial guess by theoretical relative abundance
for i = 1 : numofglycans
    ithglycandist    = md{i}(:,2);
    if(i==1)
        ithpeakarea = abs(FWHM(monopeak,1)-FWHM(monopeak,2))...
            *peaklist(monopeak,2);
        TheoArea{i} = ithpeakarea/ithglycandist(1,1);
    else
        isotopicindex = monopeak+startIndex(i)-1;
        ithpeakarea   = abs(FWHM(isotopicindex,1)-FWHM(isotopicindex,2))...
                      *peaklist(isotopicindex,2);
        for j = 1 :  i-1
            jthglycandist     = md{j}(:,2);
            areafromjthglycan = TheoArea{j}*jthglycandist(startIndex(j+1),1);
            ithpeakarea       = ithpeakarea-areafromjthglycan;
        end
        if(ithpeakarea<=0)
            TheoArea{i} = 0;
        else
            TheoArea{i} = ithpeakarea/ithglycandist(1,1);
        end
    end
end
end

function [ismatch,glycantotalarea,Centroid,CompareRatio,Residuemz] = ResultsFit(md,TheoArea,exptPeakarea,...
    startIndex,monopeak,peaklist,ID,Totalarea,Residuemz)
ismatch = 1;
glycantotalarea = '';
for i = 1 : length(TheoArea)
    glycantotalarea{i} = TheoArea{i};
end
glycanarea = zeros(length(exptPeakarea),1);
theoDist   = zeros(length(exptPeakarea),1);
for i = 1 : length(glycantotalarea)
    ithratio  = glycantotalarea{i}/Totalarea;
    ithglycan = glycantotalarea{i};
    dist = md{1,i};
    ithglycanarea = zeros(length(exptPeakarea),1);
    ithTheoDist   = zeros(length(exptPeakarea),1);
    for j = startIndex(i):length(exptPeakarea)%min(SingleAssignWindow(i)+startIndex(i)-1,length(exptPeakarea))
        ithglycanarea(j) = ithglycan*dist(j-startIndex(i)+1,2);
        ithTheoDist(j)   = dist(j-startIndex(i)+1,2)*ithratio;
    end
    glycanarea = glycanarea+ithglycanarea;
    theoDist   = theoDist+ithTheoDist;
end

% Calculate Centroid value
Centroid = getCentroid(glycanarea,monopeak,peaklist);

% Calculate ID value
if(length(md)>1)
    CompareRatio = getIDvalue(glycanarea,theoDist);
    Theoratio    = ((length(glycanarea)-1)*((log10(1.5))^2))^0.5;
    if(CompareRatio>Theoratio)
        ismatch  = 0;
    end
else
    CompareRatio = ID;
end

% calculate Residue
if(ismatch)
    startpeak = monopeak;
    for i = 1 : length(glycanarea)
        ithCalarea  = glycanarea(i);
        ithExptarea = exptPeakarea(i);
        if(~isempty(Residuemz))
            listnumber  = length(Residuemz(:,1));
        else
            listnumber  = length(Residuemz);
        end
        Residuemz(listnumber+1,1) = startpeak;
        Residuemz(listnumber+1,2) = peaklist(startpeak,1);
        Residuemz(listnumber+1,3) = ithCalarea-ithExptarea; %((ithCalarea-ithExptarea)/ithExptarea)*100;
        startpeak = startpeak+1;
    end
end
end

function CompareRatio = getIDvalue(glycanarea,theoDist)
CompareRatio = 0;
for i = 1 : length(glycanarea)-1
    ithTheoratio = theoDist(i)/theoDist(i+1);
    ithCalratio  = glycanarea(i)/glycanarea(i+1);
    Calratio     = (log10(ithCalratio/ithTheoratio))^2;
    CompareRatio = CompareRatio+Calratio;
end
CompareRatio = CompareRatio^0.5;
end

function Centroid = getCentroid(glycanarea,monopeak,peaklist)
CalX         = 0;
CalY         = 0;
starpeak     = monopeak; 
for i = 1: length(glycanarea)
    CalX     = CalX+peaklist(starpeak,1)*(glycanarea(i));
    CalY     = CalY+(glycanarea(i));
    starpeak = starpeak+1;
end
Centroid     = CalX/CalY;
end

function matchedpeakindex = getMatchedpeaks(matchedpeakindex,AssignWindow,targetMono)
listpeak = targetMono;
for i  = 1 : length(AssignWindow)
    listnumber = length(matchedpeakindex);
    if(isempty(find((matchedpeakindex-listpeak)==0, 1)))
        matchedpeakindex(listnumber+1) = listpeak;
    end
    listpeak = listpeak+1;
end
end
