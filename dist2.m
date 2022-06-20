function [ Result ] = dist2(a,b)
    % a = [x_a, y_a]
    % b = [x_b, y_b]
    Result = sqrt((a(1)-b(1))^2+(a(2)-b(2))^2);
end
