clear all; clc;
close all;

%%% 0: RSRP Pattern 생성 모드       %%%
%%% 1: Training/Test Data 생성 모드 %%%

mode  = 1; 

P_tx   = 33;         % gNB Transmission Power [dBm]
V      = 4;          % User Velocity [m/s]
T_MR   = 0.040;      % Measurement Interval [ms]
A2_Th  = -69;        % Event A2 Threshold [dBm] (-66부터 가능할 듯)
A3_Off = 5;          % Event A3 Offset [dBm]
blockage_th = -100;   % serving의 blockage threshold [dBm]
A4_Th = -80;         % Event A4 Threshold [dBm]

V = 4;

Filter_k = 1;                   % L3 Filter Coefficient
Filter_a = (1/2)^(Filter_k/4);  % The Forgetting Factor
N_MR     = 20;                  % The Length of Time-series MR
N_gNB    = 7;                   % The Number of gNBs
N_Path   = 8;                   % The Number of Paths

if mode == 0
    N_Data = 1;
elseif mode == 1
    N_Data = 1;              % The Number of training/test data
end


%% gNBs Position Setting
ISD      = 200;
gNB(1,:) = [0,0];
gNB(2,:) = [ISD,0];  gNB(3,:) = [ISD/2,sqrt(3)*ISD/2];    gNB(4,:) = [-ISD/2,sqrt(3)*ISD/2];
gNB(5,:) = [-ISD,0]; gNB(6,:) = [-ISD/2,-sqrt(3)*ISD/2];  gNB(7,:) = [ISD/2,-sqrt(3)*ISD/2];


interval = [2; 2];
gNB_Building = 1+interval(1,1);
gNB = gNB + interval.';

%% Building Position Setting
zone_size    = 100 - gNB_Building - max(interval);                     % Building의 중심이 생기는 범위
lambda       = 1;                       % PPP Density
num_area = 16;                          % 건물 구역 수 
for i = 1 : num_area
    N_Build_area(i) = poissrnd(lambda); % Number of Buildings (/unit_region)
end
N_Buildings = sum(N_Build_area);
temp        = rand(2,N_Buildings);
for i = 1 : N_Buildings
    max_length(i) = 2*zone_size*min(temp(1,i), 1-temp(1,i)) ;
    max_width(i)  = 2*zone_size*min(temp(2,i), 1-temp(2,i)) ;
end
stack = zeros(num_area,1);
for i = 1 : N_Buildings
    if stack(1) < N_Build_area(1)
        XYCenter(:,i) = zone_size*temp(:,i) + [-200; -200] + gNB_Building;
        stack(1) = stack(1) + 1;
    elseif stack(2) < N_Build_area(2)
        XYCenter(:,i) = zone_size*temp(:,i) + [-200; -100] + gNB_Building;
        stack(2) = stack(2) + 1;
    elseif stack(3) < N_Build_area(3)
        XYCenter(:,i) = zone_size*temp(:,i) + [-200; 0] + gNB_Building;
        stack(3) = stack(3) + 1;
    elseif stack(4) < N_Build_area(4)
        XYCenter(:,i) = zone_size*temp(:,i) + [-200; 100] + gNB_Building;
        stack(4) = stack(4) + 1;
    elseif stack(5) < N_Build_area(5)
        XYCenter(:,i) = zone_size*temp(:,i) + [-100; -200] + gNB_Building;
        stack(5) = stack(5) + 1;
    elseif stack(6) < N_Build_area(6)
        XYCenter(:,i) = zone_size*temp(:,i) + [-100; -100] + gNB_Building;
        stack(6) = stack(6) + 1;
    elseif stack(7) < N_Build_area(7)
        XYCenter(:,i) = zone_size*temp(:,i) + [-100; 0] + gNB_Building;
        stack(7) = stack(7) + 1;
    elseif stack(8) < N_Build_area(8)
        XYCenter(:,i) = zone_size*temp(:,i) + [-100; 100] + gNB_Building;
        stack(8) = stack(8) + 1;
    elseif stack(9) < N_Build_area(9)
        XYCenter(:,i) = zone_size*temp(:,i) + [0; -200] + gNB_Building;
        stack(9) = stack(9) + 1;
    elseif stack(10) < N_Build_area(10)
        XYCenter(:,i) = zone_size*temp(:,i) + [0; -100] + gNB_Building;
        stack(10) = stack(10) + 1;
    elseif stack(11) < N_Build_area(11)
        XYCenter(:,i) = zone_size*temp(:,i) + [0; 0] + gNB_Building;
        stack(11) = stack(11) + 1;
    elseif stack(12) < N_Build_area(12)
        XYCenter(:,i) = zone_size*temp(:,i) + [0; 100] + gNB_Building;
        stack(12) = stack(12) + 1;
    elseif stack(13) < N_Build_area(13)
        XYCenter(:,i) = zone_size*temp(:,i) + [100; -200] + gNB_Building;
        stack(13) = stack(13) + 1;
    elseif stack(14) < N_Build_area(14)
        XYCenter(:,i) = zone_size*temp(:,i) + [100; -100] + gNB_Building;
        stack(14) = stack(14) + 1;
    elseif stack(15) < N_Build_area(15)
        XYCenter(:,i) = zone_size*temp(:,i) + [100; 0] + gNB_Building;
        stack(15) = stack(15) + 1;
    elseif stack(16) < N_Build_area(16)
        XYCenter(:,i) = zone_size*temp(:,i) + [100; 100] + gNB_Building;
        stack(16) = stack(16) + 1;
    end
end
clear temp stack;
B_length = max_length.*rand(1,N_Buildings);
B_width  = max_width.*rand(1,N_Buildings);
for i = 1 : N_Buildings
    XYCoord(i,:,:) = [XYCenter(1,i)-B_length(1,i)/2 XYCenter(2,i)-B_width(1,i)/2; ...
                      XYCenter(1,i)+B_length(1,i)/2 XYCenter(2,i)-B_width(1,i)/2; ...
                      XYCenter(1,i)+B_length(1,i)/2 XYCenter(2,i)+B_width(1,i)/2; ...
                      XYCenter(1,i)-B_length(1,i)/2 XYCenter(2,i)+B_width(1,i)/2; ...
                      XYCenter(1,i)-B_length(1,i)/2 XYCenter(2,i)-B_width(1,i)/2];
end
% Building 회전 줄 때 필요한 코드
% max_theta = 2*pi;                                  % Maximum Building Orientation
% B_theta   = max_theta * rand(N_buildings,1);       % Building's Orientation
% for i = 1 : Number_build
%     XYCoord(i,:,1) = XYCoord(i,:,1) - XYCenter(1,i);
%     XYCoord(i,:,2) = XYCoord(i,:,2) - XYCenter(2,i);    
%     temp(i,:,1)    = cos(B_theta(i)) .* XYCoord(i,:,1) - sin(B_theta(i)) .* XYCoord(i,:,2);
%     temp(i,:,2)    = sin(B_theta(i)) .* XYCoord(i,:,1) + cos(B_theta(i)) .* XYCoord(i,:,2);
%     XYCoord(i,:,1) = temp(i,:,1) + XYCenter(1,i);
%     XYCoord(i,:,2) = temp(i,:,2) + XYCenter(2,i);
% end
% clear temp;
klast = 5 * ones(1,size(XYCoord,1));    % (Number of Vertex) + 1
for i = 1 : N_Buildings                 % Maximum distance from center
    MD(i) = dist2(XYCenter(:,i),[XYCoord(i,1,1),XYCoord(i,1,2)]); 
end
clear B_length B_width N_build_area zone_size;
%% UE Initial Position Setting
for j = 1 : N_Path
    UE_wayp1(j,:) = [(randi(5)-3)*100,(randi(5)-3)*100];
    UE_end(j,:) = [(randi(5)-3)*100,(randi(5)-3)*100];
end

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
for i = 1:N_Path
figure(i); voronoi(gNB(:,1),gNB(:,2),'--'); hold on; set(gcf,'Position',[0 500 500 500]); grid on; 
set(gca,'XTick',[-200:100:200]); set(gca,'YTick',[-200:100:200]); set(gca,'GridAlpha',1); hold on;
end

% for i = 1 : N_Path
%     Path_Nmeas(i)    = floor(dist2(UE_init(i,:),Path_end(i,:))/(V*T_MR));
% end
%% Main Function

for n = 1 : N_Data
    %Parameters for choose target cell to avoid blockage
    HO_ex = ones(N_Path, N_gNB);
    HO_ex(:,1) = 0;
    HO_ini = zeros(1,N_Path);
    RSRP_avg = zeros(N_Path, N_gNB);
    HO_min = zeros(1,N_Path);
   V = gamrnd(8,1);     % User Random Velocity [m/s]
   Vh(n)=V  *3.6;
    for p = 1 : N_Path
        ini_ini2 = XYdist(UE_ini(p,:),UE_ini2(p,:));
        ini2_wayp1 = XYdist(UE_ini2(p,:),UE_wayp1(p,:));
        wayp1_end = XYdist(UE_wayp1(p,:),UE_end(p,:));
        
             %Path_Nmeas = floor(dist2(UE_init(p,:),Path_end(p,:))/(V*T_MR));
             % Path_Nmeas = floor((dist2(UE_init(p,:),UE_wayp(p,:))+dist2(UE_wayp(p,:),UE_init2(p,:))+dist2(UE_init2(p,:),Path_end(p,:)))/(V*T_MR));
            Path_Nmeas = floor(abs(ini_ini2(1))/(V*T_MR))+floor(abs(ini_ini2(2))/(V*T_MR))+floor(abs(ini2_wayp1(1))/(V*T_MR))+floor(abs(ini2_wayp1(2))/(V*T_MR))+floor(abs(wayp1_end(1))/(V*T_MR))+floor(abs(wayp1_end(2))/(V*T_MR));
             %    Path_Nmeas = 1000;
        
        UE = zeros(Path_Nmeas,2);
        %UE(:,1) = UE_end(p,1);
        %UE(:,2) = UE_end(p,2);
        
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
                
            for j = 1 : N_gNB
                testpoint(:,1) = UE(i,:);
                testpoint(:,2) = gNB(j,:);
                d(i,j)         = dist2(testpoint(:,1),testpoint(:,2));
                flag(p,i,j)      = LOS_NLOS(testpoint,XYCoord,klast,N_Buildings,XYCenter,max(MD));
                RSRP(i,j)      = RSS(P_tx,d(i,j),flag(p,i,j),mode);        
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
            
            for j = 1 : N_gNB
                testpoint(:,1) = UE(i,:);
                testpoint(:,2) = gNB(j,:);
                d(i,j)         = dist2(testpoint(:,1),testpoint(:,2));
                flag(p,i,j)      = LOS_NLOS(testpoint,XYCoord,klast,N_Buildings,XYCenter,max(MD));
                RSRP(i,j)      = RSS(P_tx,d(i,j),flag(p,i,j),mode);        
            end
        end
        UE(i,2) = UE_ini2(p,2);
        
        tmp = ini2_wayp1;
        if floor(abs(ini_ini2(1))/(V*T_MR))+floor(abs(ini_ini2(2))/(V*T_MR))+1 <= Path_Nmeas
            for i = floor(abs(ini_ini2(1))/(V*T_MR))+floor(abs(ini_ini2(2))/(V*T_MR))+1 : floor(abs(ini_ini2(1))/(V*T_MR))+floor(abs(ini_ini2(2))/(V*T_MR))+floor(abs(ini2_wayp1(1))/(V*T_MR))
                    if tmp(1) > 0
                        tmp(1) = tmp(1) - V*T_MR;
                        UE(i,:) = UE(i-1,:) + [V*T_MR,0];              
                    else
                        tmp(1) = tmp(1) + V*T_MR;
                        UE(i,:) = UE(i-1,:) - [V*T_MR,0];
                    end
                    
            for j = 1 : N_gNB
                testpoint(:,1) = UE(i,:);
                testpoint(:,2) = gNB(j,:);
                d(i,j)         = dist2(testpoint(:,1),testpoint(:,2));
                flag(p,i,j)      = LOS_NLOS(testpoint,XYCoord,klast,N_Buildings,XYCenter,max(MD));
                RSRP(i,j)      = RSS(P_tx,d(i,j),flag(p,i,j),mode);        
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

            for j = 1 : N_gNB
                testpoint(:,1) = UE(i,:);
                testpoint(:,2) = gNB(j,:);
                d(i,j)         = dist2(testpoint(:,1),testpoint(:,2));
                flag(p,i,j)      = LOS_NLOS(testpoint,XYCoord,klast,N_Buildings,XYCenter,max(MD));
                RSRP(i,j)      = RSS(P_tx,d(i,j),flag(p,i,j),mode);        
            end
            end
        end
        UE(i,2) = UE_wayp1(p,2);
        
        tmp = wayp1_end;
        if floor(abs(ini2_wayp1(1))/(V*T_MR))+floor(abs(ini2_wayp1(2))/(V*T_MR))+1 <= Path_Nmeas
            for i = floor(abs(ini_ini2(1))/(V*T_MR))+floor(abs(ini_ini2(2))/(V*T_MR))+floor(abs(ini2_wayp1(1))/(V*T_MR))+floor(abs(ini2_wayp1(2))/(V*T_MR))+1 : floor(abs(ini_ini2(1))/(V*T_MR))+floor(abs(ini_ini2(2))/(V*T_MR))+floor(abs(ini2_wayp1(1))/(V*T_MR))+floor(abs(ini2_wayp1(2))/(V*T_MR))+floor(abs(wayp1_end(1))/(V*T_MR))
            
                    if tmp(1) > 0
                        tmp(1) = tmp(1) - V*T_MR;
                        UE(i,:) = UE(i-1,:) + [V*T_MR,0];              
                    else
                        tmp(1) = tmp(1) + V*T_MR;
                        UE(i,:) = UE(i-1,:) - [V*T_MR,0];
                    end
                    for j = 1 : N_gNB
                        testpoint(:,1) = UE(i,:);
                        testpoint(:,2) = gNB(j,:);
                        d(i,j)         = dist2(testpoint(:,1),testpoint(:,2));
                        flag(p,i,j)      = LOS_NLOS(testpoint,XYCoord,klast,N_Buildings,XYCenter,max(MD));
                        RSRP(i,j)      = RSS(P_tx,d(i,j),flag(p,i,j),mode);        
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
                 
            for j = 1 : N_gNB
                testpoint(:,1) = UE(i,:);
                testpoint(:,2) = gNB(j,:);
                d(i,j)         = dist2(testpoint(:,1),testpoint(:,2));
                flag(p,i,j)      = LOS_NLOS(testpoint,XYCoord,klast,N_Buildings,XYCenter,max(MD));
                RSRP(i,j)      = RSS(P_tx,d(i,j),flag(p,i,j),mode);        
            end
            end
            UE(i,2) = UE_end(p,2);
        end
        %%% RSRP Filtering
        [row, col] = size(RSRP);
        for i = 1 : row - 4
            for j = 1 : N_gNB
                RSRP_L1(i,j) = (RSRP(i,j)+RSRP(i+1,j)+RSRP(i+2,j)+RSRP(i+3,j)+RSRP(i+4,j))/5;
                if i == 1
                    RSRP_L3(i,j) = RSRP_L1(i,j);
                else
                    RSRP_L3(i,j) = (1-Filter_a)*RSRP_L3(i-1,j) + Filter_a*RSRP_L1(i,j);
                end
            end
        end
        
        %%% RSRP Plot
        if n == 1
         figure; hold on; grid on;
         for i = 1 : N_gNB
             k = 1 : length(RSRP_L3(:,:));
             plot(k.*T_MR*V, RSRP_L3(:,i),'linewidth',1);
         end
         xlabel('Distance [m]','fontsize',10,'fontweight','b');
         ylabel('RSRP [dBm]','fontsize',10,'fontweight','b'); 
         legend('gNB #1','gNB #2','gNB #3','gNB #4','gNB #5','gNB #6','gNB #7');
        end
        %%% (N X K) Time-series MR 
        Prep_idx = 1;
        while RSRP_L3(Prep_idx) > A2_Th
            Prep_idx = Prep_idx + 1;
        end
        if Prep_idx >= N_MR
            RSRP_MR = round(RSRP_L3((Prep_idx-N_MR+1):(Prep_idx),:));
        else
            disp('Error: A2 Event가 N이전에 발생!');
        end
        if HO_ini(p) ==0
            HO_ini(p) = Prep_idx;
        elseif HO_ini(p) > Prep_idx
            HO_ini(p) = Prep_idx;
        end
        %%% Results
        if mode == 0
            dlmwrite('pattern.csv',RSRP_MR,'precision','%.4f');
            dlmwrite('dist.csv',d,'precision','%.4f');
            dlmwrite('flag.csv',flag,'precision','%.4f');
        elseif mode == 1
            for i = 1 : N_gNB
                X(p,(i-1)*N_MR+(1:N_MR)) = RSRP_MR(:,i);    % X: (RSRP1_1,RSRP1_2,...,RSRPK_N)
            end        
            Exec_idx = ones(1,N_gNB)*size(RSRP_L3,1);
            for i = 2 : N_gNB
                for j = (Prep_idx+1) : size(RSRP_L3,1)
                    if RSRP_L3(j,i)-RSRP_L3(j,1)-A3_Off >= 0 && RSRP_L3(j,i) > A4_Th
                        Exec_idx(i) = j;
                        if HO_min(p) == 0 
                            HO_min(p) = j;
                        elseif HO_min(p) > j
                            HO_min(p) = j;
                        end
                        break
                    end
                end
            end
            RSRP_L3_srv = ones(1,N_Path);
            for k = HO_ini(p) : HO_min(p)
                if RSRP_L3(k,1) <= blockage_th
                    RSRP_L3_srv(p) = 0;
                end
            end
                    
            [temp,target] = min(Exec_idx);        
            X(p,N_MR*N_gNB+1) = target;                     % Y: Event A3가 가장 먼저 일어나는 gNB
            clear temp target Exec_idx;            % (동시에 일어나는 경우는 가장 idx가 빠른 gNB, RSRP 고려 X)
        end
%         loopchecker = [p,n]; disp(loopchecker);
    loopchecker = [p,n]; disp(loopchecker);
    
    if n==1 
        figure(p); plot(UE(:,1),UE(:,2),'linewidth',2); hold on;
    end
   for k = HO_ini(p) : HO_min(p)
    
        for h = 2 : N_gNB
            if flag(p,k,h) == 1
                HO_ex(p,h) = 0;
            end
        end
    end
    for k = 1 : N_gNB
        RSRP_L3(:,k) = RSRP_L3(:,k) .* HO_ex(p,k);    % LOS=1, NLOS=0 값을 곱해 RSRP 가공
    end
    for k = 1 : N_gNB
        RSRP_avg(p,k) = mean(RSRP_L3(HO_ini(p):HO_min(p), k));
    end
    end
    
       
       Y2 = zeros(1,N_Path);
       for q = 1 : N_Path
           c = 0;
           a = -1000;
           if RSRP_L3_srv(q) == 0
               Y2(q)=-1;
           else
           for j = 1 : N_gNB
               if RSRP_avg(q,j) ~= 0 && a < RSRP_avg(q,j)
                   a = RSRP_avg(q,j);
                   c = j;
               end 
           end
           Y2(q) = c; %save
           end
       end
       
       X = horzcat(X, Y2.');
       X_label = X(:,N_gNB*N_MR+1:N_gNB*N_MR+2);
       
       if mode == 1
         for p = 1 : N_Path
            if min(X_label(p,:)) ==1
                label = prod(X_label(p,:));
            elseif min(X_label(p,:)) == 2
                label = max(X_label(p,:))+6;
            elseif min(X_label(p,:)) == 3
                label = max(X_label(p,:))+11;
            elseif min(X_label(p,:)) == 4
                label = max(X_label(p,:))+15;
            elseif min(X_label(p,:)) == 5
                label = max(X_label(p,:))+18;
            elseif min(X_label(p,:)) == 6
                label = max(X_label(p,:))+20;
            elseif min(X_label(p,:)) == 7
                label = 28;
            elseif min(X_label(p,:)) == 0
                label = max(X_label(p,:))+28;
            else    
                label = max(X_label(p,:))+35;
            end
            Dataframe(n+N_Data*(p-1),:) = horzcat(X(p,:),label);   % Save
         end
       end
      % if mode == 1
       %    for p = 1 : N_Path
        %    [t,r]=size(Dataframe);
         %   Dataframe(n+N_Data*(p-1),r+1) = Y2(p); %save
          % end
       %end


    clear d testpoint RSRP RSRP_L1 RSRP_MR X X_label label;
end
%Save Dataframe
if mode == 1
    dlmwrite('Dataframe_l1_special.csv',Dataframe,'precision','%.4f');
end
% Simulation Settings Plot
% figure; voronoi(gNB(:,1),gNB(:,2)); hold on; set(gcf,'Position',[0 500 500 500]); grid on;
% for i = 1 : N_Path
%     plot([UE_init(i,1), Path_end(i,1)],[UE_init(i,2), Path_end(i,2)],'linewidth',2);
% end

for p = 1:N_Path
for i = 1 : size(gNB,1)
    figure(p);text(gNB(i,1),gNB(i,2),['gNB #',num2str(i)]);
end
end
for p = 1:N_Path
for i = 1 : N_Buildings
    figure(p);plot([XYCoord(i,1,1),XYCoord(i,2,1)],[XYCoord(i,1,2),XYCoord(i,2,2)],'k','linewidth',1); hold on;
    plot([XYCoord(i,2,1),XYCoord(i,3,1)],[XYCoord(i,2,2),XYCoord(i,3,2)],'k','linewidth',1);
    plot([XYCoord(i,3,1),XYCoord(i,4,1)],[XYCoord(i,3,2),XYCoord(i,4,2)],'k','linewidth',1);
    plot([XYCoord(i,4,1),XYCoord(i,1,1)],[XYCoord(i,4,2),XYCoord(i,1,2)],'k','linewidth',1);
end
end
set(gcf,'Position',[100 100 500 500]); xlim([-200 200]); ylim([-200 200]);
% plot(XYCenter(1,:),XYCenter(2,:),'r.')
clear i j k n p;
