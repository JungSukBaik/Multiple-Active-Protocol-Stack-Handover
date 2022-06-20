close all;
clear all; clc;

mode = 0;
Ptx_UE   =   23;             % UE transmission power
Ptx_BS   =   33;             % Cell Transmission power
V        =   30;           % Velocity [km/h]                                               

Off1     =   5;              % A3 Offset
Off2     =   -5;
dev_LOS  =   5.8;            % ShadowFading deviation LoS
dev_NLOS =   8.7;            % ShadowFading deviation NLoS 
Rout     =   -100;           
TH_A2    =   -69;            % mMAPS Threshold % 68이었음
TH_A41   =   -80;            % Execution Event
blockage =   -100;            % blockage threshold
TTP      =   0.080;          % TTT [ms]
TTE      =   0.080;          % TTE [ms] % TTT = ms0, ms40, ms64, ms80, ms100, ms128, ms160, ms256, ms320, ms480, ms512, ms640, ms1024, ms1280, ms2560, ms5120
T_MR     =   0.040;          % Measurement Gap Repetition Period 40 or 80
Qin      =   -12;            % [dB]
Qout     =   -12;            % [dB]
Noise_figure =   -1000;
Noise_figure = 10^(Noise_figure/10);
T_recovery = 0.4;            % [s]
T_HIT = 0.05;                % [s]

Num_gNB    = 7;                   % The Number of gNBs
Num_Path   = 12;                   % The Number of Paths
Num_Data   = 3;
num_area = 16;
ccc = zeros(8,Num_Data,Num_Path);    

%% gNB & UE & Path
ISD      = 200;
gNB(1,:) = [0,0];
gNB(2,:) = [ISD,0];  gNB(3,:) = [ISD/2,sqrt(3)*ISD/2];    gNB(4,:) = [-ISD/2,sqrt(3)*ISD/2];
gNB(5,:) = [-ISD,0]; gNB(6,:) = [-ISD/2,-sqrt(3)*ISD/2];  gNB(7,:) = [ISD/2,-sqrt(3)*ISD/2];
for j = 1 : Num_Path
    UE_wayp1(j,:) = [(randi(5)-3)*100,(randi(5)-3)*100];
    UE_end(j,:) = [(randi(5)-3)*100,(randi(5)-3)*100];
end

interval = [2; 2];
gNB_Building = 1+interval(1,1);
gNB = gNB + interval.';

UE_ini(1,:) = [0,0];     UE_ini(2,:) = [0,0];
UE_ini(3,:) = [0,0];     UE_ini(4,:) = [0,0];
UE_ini(5,:) = [0,0];     UE_ini(6,:) = [0,0];
UE_ini(7,:) = [0,0];     UE_ini(8,:) = [0,0];
UE_ini(9,:) = [0,0];     UE_ini(10,:) = [0,0];
UE_ini(11,:) = [0,0];    UE_ini(12,:) = [0,0];
%% Path Information Setting
UE_ini2(1,:) = [-100,0];   UE_ini2(2,:) = [100,0];
UE_ini2(3,:) = [0,-100];   UE_ini2(4,:) = [0,100];
UE_ini2(5,:) = [-100,0];   UE_ini2(6,:) = [100,0];
UE_ini2(7,:) = [0,-100];   UE_ini2(8,:) = [0,100];
UE_ini2(9,:) = [-100,0];   UE_ini2(10,:) = [100,0];
UE_ini2(11,:) = [0,-100];  UE_ini2(12,:) = [0,100];

%% Loop (lambda)
for l = 1 : 8
    lambda       = l / 2 ;
    D_HOsucc_A3 = 0; D_HOsucc_A4 = 0; D_HOsucc_A3_A4 = 0; D_HOsucc_PCHO = 0; D_HOsucc_PCHO_A3_A4 = 0; D_HOsucc_mMAPS = 0;
    D_RLF_A3 = 0; D_RLF_A4 = 0; D_RLF_A3_A4 = 0; D_RLF_PCHO = 0; D_RLF_PCHO_A3_A4 = 0; D_RLF_PCHO_DAPS = 0; D_RLF_mMAPS0 = 0; D_RLF_mMAPS = 0; D_RLF_mMAPS2 = 0; D_RLF_mMAPS3 = 0;
    zone_size    = 100 - gNB_Building - max(interval);  
    %% Monte Carlo Loop
    for n = 1 : Num_Data
        %Parameters for choose target cell to avoid blockage
        HO_ex = ones(Num_Path, Num_gNB);
        HO_ex(:,1) = 0;
        HO_ini = zeros(1,Num_Path);
        RSRP_avg = zeros(Num_Path, Num_gNB);
        HO_min = zeros(1,Num_Path);
        %% Blockage
        Num_Buildings = 0;
        while (Num_Buildings < 1)
            for i = 1 : 16
                Num_Build_area(i) = poissrnd(lambda); % Number of Buildings (/unit_region)
            end
            Num_Buildings = sum(Num_Build_area);
        end        
        temp = rand(2,Num_Buildings);
        for i = 1 : Num_Buildings
            max_length(i) = 2*zone_size*min(temp(1,i), 1-temp(1,i)) * 0.95;
            max_width(i)  = 2*zone_size*min(temp(2,i), 1-temp(2,i)) * 0.95;
        end
stack = zeros(num_area,1);
for i = 1 : Num_Buildings
    if stack(1) < Num_Build_area(1)
        XYCenter(:,i) = zone_size*temp(:,i) + [-200; -200] + gNB_Building;
        stack(1) = stack(1) + 1;
    elseif stack(2) < Num_Build_area(2)
        XYCenter(:,i) = zone_size*temp(:,i) + [-200; -100] + gNB_Building;
        stack(2) = stack(2) + 1;
    elseif stack(3) < Num_Build_area(3)
        XYCenter(:,i) = zone_size*temp(:,i) + [-200; 0] + gNB_Building;
        stack(3) = stack(3) + 1;
    elseif stack(4) < Num_Build_area(4)
        XYCenter(:,i) = zone_size*temp(:,i) + [-200; 100] + gNB_Building;
        stack(4) = stack(4) + 1;
    elseif stack(5) < Num_Build_area(5)
        XYCenter(:,i) = zone_size*temp(:,i) + [-100; -200] + gNB_Building;
        stack(5) = stack(5) + 1;
    elseif stack(6) < Num_Build_area(6)
        XYCenter(:,i) = zone_size*temp(:,i) + [-100; -100] + gNB_Building;
        stack(6) = stack(6) + 1;
    elseif stack(7) < Num_Build_area(7)
        XYCenter(:,i) = zone_size*temp(:,i) + [-100; 0] + gNB_Building;
        stack(7) = stack(7) + 1;
    elseif stack(8) < Num_Build_area(8)
        XYCenter(:,i) = zone_size*temp(:,i) + [-100; 100] + gNB_Building;
        stack(8) = stack(8) + 1;
    elseif stack(9) < Num_Build_area(9)
        XYCenter(:,i) = zone_size*temp(:,i) + [0; -200] + gNB_Building;
        stack(9) = stack(9) + 1;
    elseif stack(10) < Num_Build_area(10)
        XYCenter(:,i) = zone_size*temp(:,i) + [0; -100] + gNB_Building;
        stack(10) = stack(10) + 1;
    elseif stack(11) < Num_Build_area(11)
        XYCenter(:,i) = zone_size*temp(:,i) + [0; 0] + gNB_Building;
        stack(11) = stack(11) + 1;
    elseif stack(12) < Num_Build_area(12)
        XYCenter(:,i) = zone_size*temp(:,i) + [0; 100] + gNB_Building;
        stack(12) = stack(12) + 1;
    elseif stack(13) < Num_Build_area(13)
        XYCenter(:,i) = zone_size*temp(:,i) + [100; -200] + gNB_Building;
        stack(13) = stack(13) + 1;
    elseif stack(14) < Num_Build_area(14)
        XYCenter(:,i) = zone_size*temp(:,i) + [100; -100] + gNB_Building;
        stack(14) = stack(14) + 1;
    elseif stack(15) < Num_Build_area(15)
        XYCenter(:,i) = zone_size*temp(:,i) + [100; 0] + gNB_Building;
        stack(15) = stack(15) + 1;
    elseif stack(16) < Num_Build_area(16)
        XYCenter(:,i) = zone_size*temp(:,i) + [100; 100] + gNB_Building;
        stack(16) = stack(16) + 1;
    end
end
clear temp stack;

        B_length = max_length.*rand(1,Num_Buildings);
        B_width  = max_width.*rand(1,Num_Buildings);
        for i = 1 : Num_Buildings
            XYCoord(i,:,:) = [XYCenter(1,i)-B_length(1,i)/2 XYCenter(2,i)-B_width(1,i)/2; ...
                              XYCenter(1,i)+B_length(1,i)/2 XYCenter(2,i)-B_width(1,i)/2; ...
                              XYCenter(1,i)+B_length(1,i)/2 XYCenter(2,i)+B_width(1,i)/2; ...
                              XYCenter(1,i)-B_length(1,i)/2 XYCenter(2,i)+B_width(1,i)/2; ...
                              XYCenter(1,i)-B_length(1,i)/2 XYCenter(2,i)-B_width(1,i)/2];
        end
        klast = 5 * ones(1,size(XYCoord,1));    % (Number of Vertex) + 1
        for i = 1 : Num_Buildings                 % Maximum distance from center
            MD(i) = dist2(XYCenter(:,i),[XYCoord(i,1,1),XYCoord(i,1,2)]); 
        end
        clear B_length B_width Num_Build_area;

        %%
     for p = 1 : Num_Path    
            ini_ini2 = XYdist(UE_ini(p,:),UE_ini2(p,:));
        ini2_wayp1 = XYdist(UE_ini2(p,:),UE_wayp1(p,:));
        wayp1_end = XYdist(UE_wayp1(p,:),UE_end(p,:));
            N_MR = floor(abs(ini_ini2(1))/(V*T_MR))+floor(abs(ini_ini2(2))/(V*T_MR))+floor(abs(ini2_wayp1(1))/(V*T_MR))+floor(abs(ini2_wayp1(2))/(V*T_MR))+floor(abs(wayp1_end(1))/(V*T_MR))+floor(abs(wayp1_end(2))/(V*T_MR));
            N_TTP = ceil(TTP/T_MR);            %%% # of Measure during TTT
            N_TTE = ceil(TTE/T_MR);            %%% # of Measure during TTE
            UE = zeros(N_MR,2);
        tmp = ini_ini2;
        %%% RSRP Generation
        for i = 1 : floor(abs(ini_ini2(1))/(V*T_MR))     
                if tmp(1) > 0
                    tmp(1) = tmp(1) - V*T_MR;
                    if i ~= 1
                        UE(i,:) = UE(i-1,:) + [V*T_MR,0];
                    else
                        UE(i,:) = UE_ini(p,:) + [V*T_MR,0];
                    end
                else
                    tmp(1) = tmp(1) + V*T_MR;
                    if i ~= 1
                        UE(i,:) = UE(i-1,:) - [V*T_MR,0];
                    else
                        UE(i,:) = UE_ini(p,:) - [V*T_MR,0];
                    end
                end
                
            for j = 1 : Num_gNB
                testpoint(:,1) = UE(i,:);
                testpoint(:,2) = gNB(j,:);
                d(i,j)         = dist2(testpoint(:,1),testpoint(:,2));
                flag(p,i,j)      = LOS_NLOS(testpoint,XYCoord,klast,Num_Buildings,XYCenter,max(MD));
                RSRP_DL(i,j)     = RSS(Ptx_BS,d(i,j),flag(p,i,j),mode);
                RSRP_DL1(i,j)     = RSS(Ptx_BS,d(i,j),flag(p,i,j),1);
                RSRP_UL(i,j)     = RSS(Ptx_UE,d(i,j),flag(p,i,j),mode);
                RSRP_linear(i,j) = 10^(RSRP_DL(i,j)/10);
                RSRP_linear1(i,j) = 10^(RSRP_DL1(i,j)/10);
            end
            for j = 1 : Num_gNB
                S            = RSRP_linear(i,j);
                I            = sum(RSRP_linear(i,:)) - S + Noise_figure;   % Noise figure 추가 필요 %
                S1            = RSRP_linear1(i,j);
                I1            = sum(RSRP_linear1(i,:)) - S1 + Noise_figure;   % Noise figure 추가 필요 %
                SINR(i,j)    = S/I;
                SINR1(i,j)   = S1/I1;
                SINR_dB(i,j) = 10*log10(SINR(i,j));     
                SINR_dB1(i,j) = 10*log10(SINR1(i,j));
            end        
        end
        UE(i,1) = UE_ini2(p,1);
        
       for i = floor(abs(ini_ini2(1))/(V*T_MR))+1:floor(abs(ini_ini2(1))/(V*T_MR))+floor(abs(ini_ini2(2))/(V*T_MR))      
                if tmp(2) > 0
                    tmp(2) = tmp(2) - V*T_MR;
                    if i ~= 1
                        UE(i,:) = UE(i-1,:) + [0, V*T_MR];
                    else
                        UE(i,:) = UE_ini(p,:) + [0, V*T_MR];
                    end
                else
                    tmp(2) = tmp(2) + V*T_MR;
                    if i ~= 1
                        UE(i,:) = UE(i-1,:) - [0, V*T_MR];
                    else
                        UE(i,:) = UE_ini(p,:) - [0, V*T_MR];
                    end
                
                end
            
            for j = 1 : Num_gNB
                testpoint(:,1) = UE(i,:);
                testpoint(:,2) = gNB(j,:);
                d(i,j)         = dist2(testpoint(:,1),testpoint(:,2));
                flag(p,i,j)      = LOS_NLOS(testpoint,XYCoord,klast,Num_Buildings,XYCenter,max(MD));
                RSRP_DL(i,j)     = RSS(Ptx_BS,d(i,j),flag(p,i,j),mode);
                RSRP_DL1(i,j)     = RSS(Ptx_BS,d(i,j),flag(p,i,j),1);
                RSRP_UL(i,j)     = RSS(Ptx_UE,d(i,j),flag(p,i,j),mode);
                RSRP_linear(i,j) = 10^(RSRP_DL(i,j)/10);
                RSRP_linear1(i,j) = 10^(RSRP_DL1(i,j)/10);
            end
            for j = 1 : Num_gNB
                S            = RSRP_linear(i,j);
                I            = sum(RSRP_linear(i,:)) - S + Noise_figure;   % Noise figure 추가 필요 %
                S1            = RSRP_linear1(i,j);
                I1            = sum(RSRP_linear1(i,:)) - S1 + Noise_figure;   % Noise figure 추가 필요 %
                SINR(i,j)    = S/I;
                SINR1(i,j)   = S1/I1;
                SINR_dB(i,j) = 10*log10(SINR(i,j));      
                SINR_dB1(i,j) = 10*log10(SINR1(i,j));
            end
        end
        UE(i,2) = UE_ini2(p,2);
        
        tmp = ini2_wayp1;
        if floor(abs(ini_ini2(1))/(V*T_MR))+floor(abs(ini_ini2(2))/(V*T_MR))+1 <= N_MR
            for i = floor(abs(ini_ini2(1))/(V*T_MR))+floor(abs(ini_ini2(2))/(V*T_MR))+1 : floor(abs(ini_ini2(1))/(V*T_MR))+floor(abs(ini_ini2(2))/(V*T_MR))+floor(abs(ini2_wayp1(1))/(V*T_MR))
                    if tmp(1) > 0
                        tmp(1) = tmp(1) - V*T_MR;
                        UE(i,:) = UE(i-1,:) + [V*T_MR,0];              
                    else
                        tmp(1) = tmp(1) + V*T_MR;
                        UE(i,:) = UE(i-1,:) - [V*T_MR,0];
                    end
                    
            for j = 1 : Num_gNB
                testpoint(:,1) = UE(i,:);
                testpoint(:,2) = gNB(j,:);
                d(i,j)         = dist2(testpoint(:,1),testpoint(:,2));
                flag(p,i,j)      = LOS_NLOS(testpoint,XYCoord,klast,Num_Buildings,XYCenter,max(MD));
                RSRP_DL(i,j)     = RSS(Ptx_BS,d(i,j),flag(p,i,j),mode);
                RSRP_DL1(i,j)     = RSS(Ptx_BS,d(i,j),flag(p,i,j),1);
                RSRP_UL(i,j)     = RSS(Ptx_UE,d(i,j),flag(p,i,j),mode);
                RSRP_linear(i,j) = 10^(RSRP_DL(i,j)/10);
                RSRP_linear1(i,j) = 10^(RSRP_DL1(i,j)/10);
            end
            for j = 1 : Num_gNB
                S            = RSRP_linear(i,j);
                I            = sum(RSRP_linear(i,:)) - S + Noise_figure;   % Noise figure 추가 필요 %
                S1            = RSRP_linear1(i,j);
                I1            = sum(RSRP_linear1(i,:)) - S1 + Noise_figure;   % Noise figure 추가 필요 %
                SINR(i,j)    = S/I;
                SINR1(i,j)   = S1/I1;
                SINR_dB(i,j) = 10*log10(SINR(i,j));       
                SINR_dB1(i,j) = 10*log10(SINR1(i,j));
            end
            end
            UE(i,1) = UE_wayp1(p,1);
            for  i = floor(abs(ini_ini2(1))/(V*T_MR))+floor(abs(ini_ini2(2))/(V*T_MR))+floor(abs(ini2_wayp1(1))/(V*T_MR))+1 : floor(abs(ini_ini2(1))/(V*T_MR))+floor(abs(ini_ini2(2))/(V*T_MR))+floor(abs(ini2_wayp1(1))/(V*T_MR))+floor(abs(ini2_wayp1(2))/(V*T_MR))  
               
                    if tmp(2) > 0
                        tmp(2) = tmp(2) - V*T_MR;
                        UE(i,:) = UE(i-1,:) + [0, V*T_MR];
                    else
                        tmp(2) = tmp(2) + V*T_MR;
                        UE(i,:) = UE(i-1,:) - [0, V*T_MR];
                    end

            for j = 1 : Num_gNB
                testpoint(:,1) = UE(i,:);
                testpoint(:,2) = gNB(j,:);
                d(i,j)         = dist2(testpoint(:,1),testpoint(:,2));
                flag(p,i,j)      = LOS_NLOS(testpoint,XYCoord,klast,Num_Buildings,XYCenter,max(MD));
                RSRP_DL(i,j)     = RSS(Ptx_BS,d(i,j),flag(p,i,j),mode);
                RSRP_DL1(i,j)     = RSS(Ptx_BS,d(i,j),flag(p,i,j),1);
                RSRP_UL(i,j)     = RSS(Ptx_UE,d(i,j),flag(p,i,j),mode);
                RSRP_linear(i,j) = 10^(RSRP_DL(i,j)/10);
                RSRP_linear1(i,j) = 10^(RSRP_DL1(i,j)/10);
            end
            for j = 1 : Num_gNB
                S            = RSRP_linear(i,j);
                I            = sum(RSRP_linear(i,:)) - S + Noise_figure;   % Noise figure 추가 필요 %
                S1            = RSRP_linear1(i,j);
                I1            = sum(RSRP_linear1(i,:)) - S1 + Noise_figure;   % Noise figure 추가 필요 %
                SINR(i,j)    = S/I;
                SINR1(i,j)   = S1/I1;
                SINR_dB(i,j) = 10*log10(SINR(i,j));    
                SINR_dB1(i,j) = 10*log10(SINR1(i,j));
            end
            end
        end
        UE(i,2) = UE_wayp1(p,2);
        
        tmp = wayp1_end;
        if floor(abs(ini2_wayp1(1))/(V*T_MR))+floor(abs(ini2_wayp1(2))/(V*T_MR))+1 <= N_MR
            for i = floor(abs(ini_ini2(1))/(V*T_MR))+floor(abs(ini_ini2(2))/(V*T_MR))+floor(abs(ini2_wayp1(1))/(V*T_MR))+floor(abs(ini2_wayp1(2))/(V*T_MR))+1 : floor(abs(ini_ini2(1))/(V*T_MR))+floor(abs(ini_ini2(2))/(V*T_MR))+floor(abs(ini2_wayp1(1))/(V*T_MR))+floor(abs(ini2_wayp1(2))/(V*T_MR))+floor(abs(wayp1_end(1))/(V*T_MR))
            
                    if tmp(1) > 0
                        tmp(1) = tmp(1) - V*T_MR;
                        UE(i,:) = UE(i-1,:) + [V*T_MR,0];              
                    else
                        tmp(1) = tmp(1) + V*T_MR;
                        UE(i,:) = UE(i-1,:) - [V*T_MR,0];
                    end
                    for j = 1 : Num_gNB
                testpoint(:,1) = UE(i,:);
                testpoint(:,2) = gNB(j,:);
                d(i,j)         = dist2(testpoint(:,1),testpoint(:,2));
                flag(p,i,j)      = LOS_NLOS(testpoint,XYCoord,klast,Num_Buildings,XYCenter,max(MD));
                RSRP_DL(i,j)     = RSS(Ptx_BS,d(i,j),flag(p,i,j),mode);
                RSRP_DL1(i,j)     = RSS(Ptx_BS,d(i,j),flag(p,i,j),1);
                RSRP_UL(i,j)     = RSS(Ptx_UE,d(i,j),flag(p,i,j),mode);
                RSRP_linear(i,j) = 10^(RSRP_DL(i,j)/10);
                RSRP_linear1(i,j) = 10^(RSRP_DL1(i,j)/10);
                    end
            for j = 1 : Num_gNB
                S            = RSRP_linear(i,j);
                I            = sum(RSRP_linear(i,:)) - S + Noise_figure;   % Noise figure 추가 필요 %
                S1            = RSRP_linear1(i,j);
                I1            = sum(RSRP_linear1(i,:)) - S1 + Noise_figure;   % Noise figure 추가 필요 %
                SINR(i,j)    = S/I;
                SINR1(i,j)   = S1/I1;
                SINR_dB(i,j) = 10*log10(SINR(i,j));     
                SINR_dB1(i,j) = 10*log10(SINR1(i,j));
            end
            end
            UE(i,1) = UE_end(p,1);
            
             for i = floor(abs(ini_ini2(1))/(V*T_MR))+floor(abs(ini_ini2(2))/(V*T_MR))+floor(abs(ini2_wayp1(1))/(V*T_MR))+floor(abs(ini2_wayp1(2))/(V*T_MR))+floor(abs(wayp1_end(1))/(V*T_MR))+1 : floor(abs(ini_ini2(1))/(V*T_MR))+floor(abs(ini_ini2(2))/(V*T_MR))+floor(abs(ini2_wayp1(1))/(V*T_MR))+floor(abs(ini2_wayp1(2))/(V*T_MR))+floor(abs(wayp1_end(1))/(V*T_MR))+floor(abs(wayp1_end(2))/(V*T_MR))
                    if tmp(2) > 0
                        tmp(2) = tmp(2) - V*T_MR;
                        UE(i,:) = UE(i-1,:) + [0, V*T_MR];
                    else
                        tmp(2) = tmp(2) + V*T_MR;
                        UE(i,:) = UE(i-1,:) - [0, V*T_MR];
                    end
                 
            for j = 1 : Num_gNB
                testpoint(:,1) = UE(i,:);
                testpoint(:,2) = gNB(j,:);
                d(i,j)         = dist2(testpoint(:,1),testpoint(:,2));
                flag(p,i,j)      = LOS_NLOS(testpoint,XYCoord,klast,Num_Buildings,XYCenter,max(MD));
                RSRP_DL(i,j)     = RSS(Ptx_BS,d(i,j),flag(p,i,j),mode);
                RSRP_DL1(i,j)     = RSS(Ptx_BS,d(i,j),flag(p,i,j),1);
                RSRP_UL(i,j)     = RSS(Ptx_UE,d(i,j),flag(p,i,j),mode);
                RSRP_linear(i,j) = 10^(RSRP_DL(i,j)/10);
                RSRP_linear1(i,j) = 10^(RSRP_DL1(i,j)/10);
            end
            for j = 1 : Num_gNB
                S            = RSRP_linear(i,j);
                I            = sum(RSRP_linear(i,:)) - S + Noise_figure;   % Noise figure 추가 필요 %
                S1            = RSRP_linear1(i,j);
                I1            = sum(RSRP_linear1(i,:)) - S1 + Noise_figure;   % Noise figure 추가 필요 %
                SINR(i,j)    = S/I;
                SINR1(i,j)   = S1/I1;
                SINR_dB(i,j) = 10*log10(SINR(i,j));    
                SINR_dB1(i,j) = 10*log10(SINR1(i,j)); 
            end
            end
        end
            UE(i,2) = UE_end(p,2);
            clear S I;

            %%% LOS/NLOS를 구분하여 각 위치별 Event 및 HO_Msg의 확률을 계산 %%%
            for i = 1 : N_MR      %%% Path 길이가 동일하지 않으면 N_MR이 달라짐
                for j = 1 : Num_gNB
                    if j == 1
                        P_A3Off1(i,j)  = 0;
                        P_A3Off2(i,j)  = 0;
                        P_A41(i,j)     = 0;
                        P_CFMsucc(i,j) = 0;
                        if flag(p,i,1) == 0
                            P_A2(i)      = qfunc((RSRP_DL(i,1)-TH_A2)/(dev_LOS));   %%% A2가 일어날 확률
                            P_MRsucc(i)  = qfunc((Rout-RSRP_UL(i,1))/(dev_LOS));    %%% MR 전송이 성공할 확률
                            P_CMDsucc(i) = qfunc((Qout-SINR_dB(i,1))/(dev_LOS));    %%% HO_CMD 전송이 성공할 확률
                        else
                            P_A2(i)      = qfunc((RSRP_DL(i,1)-TH_A2)/(dev_NLOS));
                            P_MRsucc(i)  = qfunc((Rout-RSRP_UL(i,1))/(dev_NLOS));
                            P_CMDsucc(i) = qfunc((Qout-SINR_dB(i,1))/(dev_NLOS));
                        end
                    else
                        if flag(p,i,1) == 0 && flag(p,i,j) == 0
                            P_A3Off1(i,j) = qfunc((RSRP_DL(i,1)+Off1-RSRP_DL(i,j))/(dev_LOS));  %%% A3(Off1)가 일어날 확률
                            P_A3Off2(i,j) = qfunc((RSRP_DL(i,1)+Off2-RSRP_DL(i,j))/(dev_LOS));  %%% A3(Off2)가 일어날 확률
                        else
                            P_A3Off1(i,j) = qfunc((RSRP_DL(i,1)+Off1-RSRP_DL(i,j))/(dev_NLOS));
                            P_A3Off2(i,j) = qfunc((RSRP_DL(i,1)+Off2-RSRP_DL(i,j))/(dev_NLOS));
                        end
                        if flag(p,i,j) == 0
                            P_A41(i,j)     = qfunc((TH_A41-RSRP_DL(i,j))/(dev_LOS));     %%% A4가 일어날 확률(EEVENT)
                            P_CFMsucc(i,j) = qfunc((Rout-RSRP_UL(i,j))/(dev_LOS));
                            P_CMPsucc(i,j) = qfunc((Rout-RSRP_UL(i,j))/(dev_LOS));
                        else
                            P_A41(i,j)     = qfunc((TH_A41-RSRP_DL(i,j))/(dev_NLOS));
                            P_CFMsucc(i,j) = qfunc((Rout-RSRP_UL(i,j))/(dev_NLOS));
                            P_CMPsucc(i,j) = qfunc((Rout-RSRP_UL(i,j))/(dev_NLOS));
                        end            
                    end        
                end 
            end


            %% NR CHO
             for k = 2 : 7
                 for i = 1 : N_MR-(N_TTP+N_TTE+1)    %%% HO Preparation Triggering %%%
                     if i == 1
                          Pst_mMAPS(k,1,1)  = 1;
                     end  
                     for j = 1 : N_TTP
                          Pst_mMAPS(k,i,j+1)  = Pst_mMAPS(k,i,j) * P_A2(i+j);
                     end    
                     P_HST1_mMAPS(k,i)  = 0;
                     for j = 1 : min(i,N_TTP)
                         P_HST1_mMAPS(k,i)  = P_HST1_mMAPS(k,i) + Pst_mMAPS(k,i+1-j, j);
                     end
                     Pst_mMAPS(k,i+1,1)  = P_HST1_mMAPS(k,i) * (1-P_A2(i+1));
                     P_PTrig_mMAPS(i+N_TTP,k)  = Pst_mMAPS(k,i,N_TTP+1);
                 end
                 P_PTrig_mMAPS(N_MR,k)  = 0;
                 for i = (N_TTP+2) : (N_MR-N_TTE)    %%% HO Execution Triggering %%%
                     if i == (N_TTP+2)
                         Pst_mMAPS(k,1,N_TTP+2)  = Pst_mMAPS(k,1,N_TTP+1);
                     end    
                     for j = 1 : N_TTE
                         Pst_mMAPS(k,i-N_TTP-1,j+N_TTP+2)  = Pst_mMAPS(k,i-N_TTP-1,j+N_TTP+1) * P_A3Off1(i+j,k)*P_A41(i+j,k);
                     end
                     P_HST2_mMAPS(k,i)  = 0;
                     for j = 1 : min(i-N_TTP-1,N_TTE)
                         P_HST2_mMAPS(k,i)  = P_HST2_mMAPS(k,i) + Pst_mMAPS(k, i-N_TTP-j, j+N_TTP+1);
                     end
                     Pst_mMAPS(k,i-N_TTP,N_TTP+2)  = Pst_mMAPS(k,i-N_TTP,N_TTP+1) + P_HST2_mMAPS(k,i) * (1-P_A3Off1(i+1,k)*P_A41(i+1,k));                
                     P_ETrig_mMAPS(i+N_TTE,k)  = Pst_mMAPS(k,i-N_TTP-1,N_TTP+N_TTE+2);    
                 end
             end

            %% Event가 Trigger 되었을 때, Average Triggering Point
            dist_path = (1:N_MR) * T_MR * V;
            for j = 2 : Num_gNB
                sw=0;
                cnt=0;
                end_p=1;
                for i = 1:N_MR
                    if sw==0 && P_ETrig_mMAPS(i,j)>=10^-3
                        sw=1;
                    elseif sw==1 && P_ETrig_mMAPS(i,j)<=10^-5
                        cnt=cnt+1;
                    end
                    if cnt>=4
                        end_p=i;
                        break;
                    end
                    if i == N_MR
                        end_p=i;
                    end
                end
                
                 Avg_ETP_mMAPS(j)  = dist_path(1:end_p) * P_ETrig_mMAPS(1:end_p,j) / sum(P_ETrig_mMAPS(1:end_p,j));
                 [temp,Avg_ETP_mMAPS_idx(j)]  = min(abs(Avg_ETP_mMAPS(j)-dist_path(:)));
                 
                  Avg_PTP_mMAPS(j)  = dist_path(1:end_p) * P_PTrig_mMAPS(1:end_p,j) / sum(P_PTrig_mMAPS(1:end_p,j));
                 [temp,Avg_PTP_mMAPS_idx(j)]  = min(abs(Avg_PTP_mMAPS(j)-dist_path(:)));
                if j==2
                    PTP_start(p)=Avg_PTP_mMAPS(j);
                    ETP_start(p)=Avg_ETP_mMAPS(j);
                else
                    if PTP_start(p)>Avg_PTP_mMAPS(j)
                        PTP_start(p)=Avg_PTP_mMAPS(j);
                    end
                    if ETP_start(p)>Avg_ETP_mMAPS(j)
                        ETP_start(p)=Avg_ETP_mMAPS(j);
                    end
                end
            end
            
            
            min_path(p) = min(horzcat(ETP_start(p),PTP_start(p)));
            max_path = ETP_start;
             [temp,temp_idx] = min(Avg_ETP_mMAPS([2:Num_gNB]));    TgNB_mMAPS  = temp_idx + 1;
%% Throughput 측정
            for o = 2 : Num_gNB
                if o==2
                    max_gNB = 2;
                else
                    if mean(SINR_dB1(floor(PTP_start(p)/(T_MR*V)):floor(ETP_start(p)/(T_MR*V)),o))>mean(SINR_dB1(floor(PTP_start(p)/(T_MR*V)):floor(ETP_start(p)/(T_MR*V)),max_gNB))
                        max_gNB=o;
                    end
                end
            end
            
            for o = 2 : Num_gNB
                if o==2 && max_gNB ~= o
                    max_gNB2 = 2;
                elseif max_gNB ~= o
                    if mean(SINR_dB1(floor(PTP_start(p)/(T_MR*V)):floor(ETP_start(p)/(T_MR*V)),o))>mean(SINR_dB1(floor(PTP_start(p)/(T_MR*V)):floor(ETP_start(p)/(T_MR*V)),max_gNB2))
                        max_gNB2=o;
                    end
                end
            end
            
            for o = 2 : Num_gNB
                if o==2 && max_gNB ~= o && max_gNB2 ~= o
                    max_gNB3 = 2;
                elseif max_gNB ~= o && max_gNB2 ~= o
                    if mean(SINR_dB1(floor(PTP_start(p)/(T_MR*V)):floor(ETP_start(p)/(T_MR*V)),o))>mean(SINR_dB1(floor(PTP_start(p)/(T_MR*V)):floor(ETP_start(p)/(T_MR*V)),max_gNB2))
                        max_gNB3=o;
                    end
                end
            end
            
            for o = 2 : Num_gNB
                if o==2
                    max_gNB1 = 2;
                else
                    if mean(SINR_dB1(floor(min_path(p)/(T_MR*V)):floor(ETP_start(p)/(T_MR*V)),o))>mean(SINR_dB1(floor(min_path(p)/(T_MR*V)):floor(ETP_start(p)/(T_MR*V)),max_gNB1))
                        max_gNB1=o;
                    end
                end
            end
   
%% Outage probability
   if min(SINR_dB1(floor(min_path(p)/(T_MR*V)):floor(max_path(p)/(T_MR*V)),1)) < Qout
            
                if ETP_start(p)<=PTP_start(p)
                SINaR_mMAPS0(p) = sum(SINR_dB1(floor(ETP_start(p)/(T_MR*V)),1)<Qout);
                else
                SINaR_mMAPS0(p) = sum(SINR_dB1(floor(PTP_start(p)/(T_MR*V)):floor(ETP_start(p)/(T_MR*V)),1)<Qout)/(floor(ETP_start(p)/(T_MR*V))-floor(PTP_start(p)/(T_MR*V))+1);
                end
                
                if ETP_start(p)<=PTP_start(p)
                SINaR_mMAPS(p) = sum(SINR_dB1(floor(ETP_start(p)/(T_MR*V)),1)<Qout).*sum(SINR_dB1(floor(ETP_start(p)/(T_MR*V)),max_gNB)<Qout);
                else
                SINaR_mMAPS(p) = sum(SINR_dB1(floor(PTP_start(p)/(T_MR*V)):floor(ETP_start(p)/(T_MR*V)),1)<Qout)/(floor(ETP_start(p)/(T_MR*V))-floor(PTP_start(p)/(T_MR*V))+1).*sum(SINR_dB1(floor(PTP_start(p)/(T_MR*V)):floor(ETP_start(p)/(T_MR*V)),max_gNB)<Qout)/(floor(ETP_start(p)/(T_MR*V))-floor(PTP_start(p)/(T_MR*V))+1);
                end
                
                if ETP_start(p)<=PTP_start(p)
                SINaR_mMAPS2(p) = sum(SINR_dB1(floor(ETP_start(p)/(T_MR*V)),1)<Qout).*sum(SINR_dB1(floor(ETP_start(p)/(T_MR*V)),max_gNB)<Qout).*sum(SINR_dB1(floor(ETP_start(p)/(T_MR*V)),max_gNB2)<Qout);
                else
                SINaR_mMAPS2(p) = sum(SINR_dB1(floor(PTP_start(p)/(T_MR*V)):floor(ETP_start(p)/(T_MR*V)),1)<Qout)/(floor(ETP_start(p)/(T_MR*V))-floor(PTP_start(p)/(T_MR*V))+1).*sum(SINR_dB1(floor(PTP_start(p)/(T_MR*V)):floor(ETP_start(p)/(T_MR*V)),max_gNB)<Qout)/(floor(ETP_start(p)/(T_MR*V))-floor(PTP_start(p)/(T_MR*V))+1).*sum(SINR_dB1(floor(PTP_start(p)/(T_MR*V)):floor(ETP_start(p)/(T_MR*V)),max_gNB2)<Qout)/(floor(ETP_start(p)/(T_MR*V))-floor(PTP_start(p)/(T_MR*V))+1);
                end
                
                if ETP_start(p)<=PTP_start(p)
                SINaR_mMAPS3(p) = sum(SINR_dB1(floor(ETP_start(p)/(T_MR*V)),1)<Qout).*sum(SINR_dB1(floor(ETP_start(p)/(T_MR*V)),max_gNB)<Qout).*sum(SINR_dB1(floor(ETP_start(p)/(T_MR*V)),max_gNB2)<Qout).*sum(SINR_dB1(floor(ETP_start(p)/(T_MR*V)),max_gNB3)<Qout);
                else
                SINaR_mMAPS3(p) = sum(SINR_dB1(floor(PTP_start(p)/(T_MR*V)):floor(ETP_start(p)/(T_MR*V)),1)<Qout)/(floor(ETP_start(p)/(T_MR*V))-floor(PTP_start(p)/(T_MR*V))+1).*sum(SINR_dB1(floor(PTP_start(p)/(T_MR*V)):floor(ETP_start(p)/(T_MR*V)),max_gNB)<Qout)/(floor(ETP_start(p)/(T_MR*V))-floor(PTP_start(p)/(T_MR*V))+1).*sum(SINR_dB1(floor(PTP_start(p)/(T_MR*V)):floor(ETP_start(p)/(T_MR*V)),max_gNB2)<Qout)/(floor(ETP_start(p)/(T_MR*V))-floor(PTP_start(p)/(T_MR*V))+1).*sum(SINR_dB1(floor(PTP_start(p)/(T_MR*V)):floor(ETP_start(p)/(T_MR*V)),max_gNB3)<Qout)/(floor(ETP_start(p)/(T_MR*V))-floor(PTP_start(p)/(T_MR*V))+1);
                end
                
                clear temp temp_idx;
        D_RLF_mMAPS = horzcat(D_RLF_mMAPS,  SINaR_mMAPS); 
        D_RLF_mMAPS0 = horzcat(D_RLF_mMAPS0,  SINaR_mMAPS0); 
        D_RLF_mMAPS2 = horzcat(D_RLF_mMAPS2,  SINaR_mMAPS2); 
        D_RLF_mMAPS3 = horzcat(D_RLF_mMAPS3,  SINaR_mMAPS3); 
   end
            %%  [Performance 2] Success probability of Early Preparation
            for j = 2 : Num_gNB
                P_HOsucc_mMAPS(j)  = P_MRsucc(1,Avg_PTP_mMAPS_idx(j))*P_CMDsucc(1,Avg_PTP_mMAPS_idx(j))*P_CMPsucc(Avg_ETP_mMAPS_idx(j),j)+(1-P_MRsucc(1,Avg_PTP_mMAPS_idx(j))*P_CMDsucc(1,Avg_PTP_mMAPS_idx(j)))*P_MRsucc(1,Avg_ETP_mMAPS_idx(j))*P_CMDsucc(1,Avg_ETP_mMAPS_idx(j))*P_CMPsucc(Avg_ETP_mMAPS_idx(j),j);
            end
            Result_HOsucc_mMAPS(p)  = P_HOsucc_mMAPS(TgNB_mMAPS);
            clear -regexp ^N_ ^P_ ^Pst_ ^T_Hold ^TgNB_ ^RSRP_ ^SINR;
            clear dist_path i j k;
     end
        for p = 1 : Num_Path
            D_HOsucc_mMAPS = vertcat(D_HOsucc_mMAPS,Result_HOsucc_mMAPS(p));
        end
        clear d flag klast max_length max_width MD Num_Build_area Num_Buildings UE;
        clear -regexp ^XY ^TP_ ^Result_;
    end
    Average_HOsucc_mMAPS(l) = sum(D_HOsucc_mMAPS(2:end))/size(D_HOsucc_mMAPS(2:end),1);
    Average_RLF_mMAPS(l) = mean(D_RLF_mMAPS);
    Average_RLF_mMAPS0(l) = mean(D_RLF_mMAPS0);
    Average_RLF_mMAPS2(l) = mean(D_RLF_mMAPS2);
    Average_RLF_mMAPS3(l) = mean(D_RLF_mMAPS3);
    clear -regexp ^D_
    disp(l);
end
Average_MIT_mMAPS = T_recovery*(1-Average_HOsucc_mMAPS);

x = (1:8)/2;
close all;
figure('position', [100, 100, 400, 400]);
plot(x, 1-Average_HOsucc_mMAPS,'-d','markersize',9,'linewidth',1);hold on; grid on;
legend('MAPS','location','southwest','fontsize',9);
xlabel('Blockage density \lambda','fontsize',12);
ylabel('Average Handover Failure rate','fontsize',12);

figure('position', [500, 100, 400, 400]);
plot(x, Average_MIT_mMAPS,'-d','markersize',9,'linewidth',1);hold on; grid on;
legend('MAPS','location','northwest','fontsize',9);
xlabel('Blockage density \lambda','fontsize',12);
ylabel('Average Mobility Interuption Time (MIT) [s]','fontsize',12);

figure('position', [900, 100, 400, 400]); 
plot(x, Average_RLF_mMAPS,'-d','markersize',9,'linewidth',1); hold on; grid on;
legend('MAPS','location','northwest','fontsize',9);
xlabel('Blockage density \lambda','fontsize',12);
ylabel('Outage probability','fontsize',12);