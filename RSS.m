function [ Result ] = RSS(P_tx,d,flag,type)
    % [Reference]
    % M. R. Akdeniz, et. al, "Millimeter wave channel modeling and
    %                         cellular capacity evalution"
    % [Parameters]
    % P_tx : Transmission power of gNB [dBm]
    % d    : Distance between UE and gNB [m]
    % flag : 0 (LOS), 1(NLOS)
    % type : 0 (No fading), 1(Small-scale fading)

    if flag == 0
        if type == 0
            PL = 61.4 + 20*log10(d);
        elseif type == 1
            PL = 61.4 + 20*log10(d) + sqrt(5.8)*randn; % dev: 5.8 [dB]
        end
        
    elseif flag == 1
        if type == 0
            PL = 72.0 + 29.2*log10(d);
        elseif type == 1
            PL = 72.0 + 29.2*log10(d) + sqrt(8.7)*randn; % dev: 8.7 [dB]
        end        
    end

    Result = P_tx - PL;
end