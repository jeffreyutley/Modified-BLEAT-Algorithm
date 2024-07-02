%
% This file is the Modified BLEAT Algorithm, including error checking at
% the tips:

% Receiving input:
%
% Driving function:
fprintf('Enter the driving function in t, beginning with "@(t)". \n');
function_input = input('Driving function: ', 's');
f = str2func(function_input);

% Upper bound of interval:
fprintf('Enter a value T > 0 such that the input function is continuous on the closed interval [0, T]. \n');
T = input('Your value: ');

% Number of subdivisions:
fprintf('Enter an initial number of subdivisions. \n');
N = input('Your number of subdivisions: ');

% Choice to refine:
fprintf('Do you want to refine the program at the tips? (N will increase)\n');
refine_choice = input('Y/N: ', 's');

% Choice to error check:
fprintf('Do you want to check for error points? (Some points will be removed from the hull)\n');
error_check = input('Y/N: ', 's');

% Evenly partitions [0,T] in an array t:
t = linspace(0, T, N+1);

% If the user selects to refine at the tip, adds more points in this area:
if refine_choice == 'Y'
    % Declares a temporary variable to hold the initial number of subdivisions:
    J = N;

    % Skips to the end of the time interval:
    i = floor(0.99*J);

    % Iterates through the driving function values at the end of the
    % interval:
    while i <= N+1
        % Adds new subdividisions until the difference of consecutive
        % driving function values is lower than a given tolerance: 1/3J
        while abs(f(t(i)) - f(t(i-1))) >= 1/(3*J)
            N = N + 1;
            for x = N+1:-1:i+1
                t(x) = t(x-1);
            end
            t(i) = (t(i) + t(i-1))/2;
        end
        i = i + 1;
    end
    
    % Prints the new number of sudvisions:
    fprintf('Final number of subdivisions: %d\n', N);
end

% Gets the list of values of the driving function:
d = arrayfun(f, t);

% Creates lists for tips in upper and lower half planes
upper_list = zeros(1, N+1);
lower_list = zeros(1, N+1);

% Main loop: finds inverse comformal map for each subinterval and applies 
% them to all necessary points
for j = 1:N
    % First calculates R and s values
    R = d(N + 2 - j) - d(N + 1 - j);
    s = t(N + 2 - j) - t(N + 1 - j);

    % Calculates a parameter depending on these values:
    Y = sqrt(R.^2 + 16*s);
            
    % Finds formula for comformal map h
    h = @(z) (0.5)*(2*z + R + Y)*((2*z + R -Y)./(2*z + R + Y)).^(0.5-(R./(2*Y)));
    
    %Finds c value
    c = (-1)*(R./sqrt(s));

    % Calculates a parameter based on c
    sc = sqrt(c.^2 + 16);
     
    %Finds a, b, alpha, and beta values
    a = c./(2*sc) - 0.5;
    b = (-1)*c./(2*sc) - 0.5;
    alpha = (c + sc)./(2);
    beta = (c - sc)./(2);
        
    % Applies map for this iteration to all points previously added to the
    % hull (which are none if j = 1):
    if j > 1
        for k = 1:j-1
            upper_list(k) = h(upper_list(k));
            lower_list(k) = h(lower_list(k));
        end
    end
                
    % Calculates new tips:
    [upper_list(j), lower_list(j)]  = tip_value(a, b, alpha, beta, s, R);
end

% Shifts the entire hull over by this value (translation property)
% (this is needed since the algorithm only approximates the change in value
% that the driving function takes on)
upper_list = upper_list + f(0);
lower_list = lower_list + f(0);

% Removes points that have been calculated incorrectly due to approximation
% error if prompted by the user:
if (error_check == 'Y')
    % Skips to the 20th index:
    x = 20;

    % Decreases index until the difference between consecutive points of
    % the upper and lower hulls is too large:
    while ((x > 0) && (abs(upper_list(x) - upper_list(x+1)) < 400/N) && (abs(lower_list(x) - lower_list(x+1)) < 400/N))
        x = x - 1;
    end
    
    % If there is a jump in values, removes this portion of the hull:
    if (x ~= 0)
        for a = 1:x
            upper_list(1) = [];
        end
    end
end

% Colors for plots:
str = '#216AB3';
color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

str0 = '#B32020';
color0 = sscanf(str0(2:end),'%2x%2x%2x',[1 3])/255;

%Creates plots, setting the window to have our desired dimensions:
plot(real(upper_list), imag(upper_list), '.', 'Color', color);
hold on;
plot(real(lower_list), imag(lower_list), '.', 'Color', color0);

% Hard-coded formula to find the upper and lower tips of hulls, modified
% for choices of branch cut (0,2*pi) and (-2*pi,0), respectively:
function [u_tip, l_tip] = tip_value(a, b, alpha, beta, s, R)
    % Hard-codes the modulus and argument of the upper tip:
    tip_modulus = -real(b)*log(abs(alpha)) + imag(b)*angle(alpha) - real(a)*log(abs(-beta)) + imag(a)*angle(-beta);
    tip_argument = -imag(b)*log(abs(alpha)) - real(b)*angle(alpha) - imag(a)*log(abs(-beta)) - real(a)*angle(-beta);
    
    % Computes the upper tip:
    u_tip = sqrt(s)*exp(tip_modulus + 1j*(tip_argument - (a*pi))) + R;
    l_tip = sqrt(s)*exp(tip_modulus + 1j*(tip_argument + (a*pi))) + R;
end