function [ WUi, WTi, PjumpFi ] = intEnerg( We_U_all_av, We_T_all_av, Fe_U_all_av, Fe_T_all_av, Pe_U_all_av, Pe_T_all_av, intParam )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fracFaces = intParam.fracFaces;
intNodes = intParam.intNodes;

NintP = size(intNodes,1);

WUi = zeros(NintP,1);
WTi = zeros(NintP,1);
PjumpFi = zeros(NintP,1);

for ii=1:NintP

    n1 = intNodes(ii,1);
    n2 = intNodes(ii,2);
    w1 = fracFaces(ii);
    w2 = 1 - fracFaces(ii);

    WU1 = We_U_all_av(n1,:);
    WT1 = We_T_all_av(n1,:);
    WU2 = We_U_all_av(n2,:);
    WT2 = We_T_all_av(n2,:);
    FU1 = Fe_U_all_av(n1,:);
    FT1 = Fe_T_all_av(n1,:);
    FU2 = Fe_U_all_av(n2,:);
    FT2 = Fe_T_all_av(n2,:);
    PU1 = Pe_U_all_av(n1,:);
    PT1 = Pe_T_all_av(n1,:);
    PU2 = Pe_U_all_av(n2,:);
    PT2 = Pe_T_all_av(n2,:);
    
    WUa = WU1*w1 + WU2*w2;
    WTa = WT1*w1 + WT2*w2;
    FUa = FU1*w1 + FU2*w2;
    FTa = FT1*w1 + FT2*w2;
    PUa = PU1*w1 + PU2*w2;
    PTa = PT1*w1 + PT2*w2;
    
    FUr = vec2tens(FUa);
    FTr = vec2tens(FTa);
    PUr = vec2tens(PUa);
    PTr = vec2tens(PTa);

    PjumpF = tensTrace3( ( PTr.' + PUr.' ) * (1/2) * ( FTr - FUr ) );
    
    WUi(ii) = WUa;
    WTi(ii) = WTa;
    PjumpFi(ii) = PjumpF;

end

end

