function MonoID = MonoIDgeneration
%% MonoIDgeneration: Assign unique ID to each monosaccharide to reduce the computation loads.

%  Author: Yusen Zhou
%  Date Lastly Updated: 01/08/2020

%% Create empty containers.Map
MonoID = struct();
% Hexose
MonoID.Glc = 1;
MonoID.Man = 2;
MonoID.Gal = 3;
MonoID.Gul = 4;
MonoID.Alt = 5;
MonoID.All = 6;
MonoID.Tal = 7;
MonoID.Ido = 8;

% HexNAc
MonoID.GlcNAc = 9;
MonoID.ManNAc = 10;
MonoID.GalNAc = 11;
MonoID.GulNAc = 12;
MonoID.AltNAc = 13;
MonoID.AllNAc = 14;
MonoID.TalNAc = 15;
MonoID.IdoNAc = 16;

% HexN
MonoID.GlcN = 17;
MonoID.ManN = 18;
MonoID.GalN = 19;
MonoID.GulN = 20;
MonoID.AltN = 21;
MonoID.AllN = 22;
MonoID.TalN = 23;
MonoID.IdoN = 24;

% HexA
MonoID.GlcA = 25;
MonoID.ManA = 26;
MonoID.GalA = 27;
MonoID.GulA = 28;
MonoID.AltA = 29;
MonoID.AllA = 30;
MonoID.TalA = 31;
MonoID.IdoA = 32;

% Deoxyhexose
MonoID.Qui   = 33;
MonoID.Rha   = 34;
MonoID.dGul  = 35;
MonoID.dAlt  = 36;
MonoID.dTal  = 37;
MonoID.Fuc   = 38;

% DeoxyhexNAc
MonoID.QuiNAc   = 39;
MonoID.RhaNAc   = 40;
MonoID.dAltNAc  = 41;
MonoID.dTalNAc  = 42;
MonoID.FucNAc   = 43;

% Di-Deoxyhexose
MonoID.Oli   = 44;
MonoID.Tyv   = 45;
MonoID.Abe   = 46;
MonoID.Par   = 47;
MonoID.Dig   = 48;
MonoID.Col   = 49;

% Pentose
MonoID.Ara   = 50;
MonoID.Lyx   = 51;
MonoID.Xyl   = 52;
MonoID.Rib   = 53;

% Nonulosonate
MonoID.Kdn    = 54;
MonoID.NeuAc  = 55;
MonoID.Neu5Ac = 55;
MonoID.NeuGc  = 56;
MonoID.Neu5Gc = 56;
MonoID.Neu    = 57;
MonoID.Sia    = 58;
MonoID.Pse    = 59;
MonoID.Leg    = 60;
MonoID.Aci    = 61;
MonoID.eLeg   = 62;

% Unknown
MonoID.Bac      = 63;
MonoID.LDManHep = 64;
MonoID.Kdo      = 65;
MonoID.Dha      = 66;
MonoID.DDManHep = 67;
MonoID.MurNAc   = 68;
MonoID.MurNGc   = 69;
MonoID.Mur      = 70;

% Assigned
MonoID.Api = 71;
MonoID.Fru = 72;
MonoID.Tag = 73;
MonoID.Sor = 74;
MonoID.Psi = 75;
end