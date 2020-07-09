function newglycanDB = classification(MSfilename,glycanDBFile,glycanDBdir,outputdir,listOfStruct)
% classification: Monosaccharide and substructure comparison analysis
%
% Input:
%   MSfilename: MS data for glycan profiles.
%     glycanDB: glycan structure data with relative abundance information
%     loadpath: Directory to load glycanDB.
%    storepath: Directory to output results.
%
% Author: Yusen Zhou
% Date Lastly Updated: 05/18/2020


glycanDBname = [glycanDBFile '.mat'];
glycanDBfullpath   = fullfile(glycanDBdir,glycanDBname);
load(glycanDBfullpath,'newglycanDB');
ResidueResult = MSResidueAnalysisGNAT(MSfilename,newglycanDB);
StructResult = MSStructAnalysisGNAT(MSfilename,newglycanDB,listOfStruct);
newglycanDB.MonoAnalysis   = ResidueResult;
newglycanDB.StructAnalysis = StructResult;
MSfilename = [MSfilename 'FinalResult' '.mat'];
savepath = fullfile(outputdir, MSfilename);
save(savepath,'newglycanDB');
end