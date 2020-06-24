function SignalReduceFigure(newglycanDB_Nonfit,newglycanDB)
FitglycanMass     = newglycanDB.monoisomw;
FitglycanMass     = cell2mat(FitglycanMass);
FitglycanAbund    = newglycanDB.abundance.*100;
NonfitglycanAbund = newglycanDB_Nonfit.abundance.*100;

bar(FitglycanMass-0.5,FitglycanAbund,'blue','EdgeColor','blue');
hold on
bar(FitglycanMass+0.5,NonfitglycanAbund,'r','EdgeColor','r');
end