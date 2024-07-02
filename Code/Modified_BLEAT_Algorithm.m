%
% Modified_BLEAT_Algorithm.m
%
% Author: Jeffrey Utley
% Date: May 10, 2022
% Title: Modified BLEAT Algorithm for Approximating Complex-Driven Loewner
% Hulls
%
% 
% This algorithm, titled the Modified BLEAT Algorithm, extends a slightly
% modified version the algorithm found in the paper BLEAT: The
% Curve-Creating Black box by Sarah Claiborne and Catherine P. Simpson 
% (Section 4.2). The key difference is that the BLEAT algorithm 
% approximates Loewner hulls in the upper-half complex plane, with 
% real-valued driving functions. The Modified BLEAT Algorithm approximates
% Loewner hulls in the entire complex plane and with complex-valued driving
% functions.
%
% Note: When a real-valued driving function is used in this program, the
% result is the Loewner hull and its reflection about the real line (as in
% the Schwarz Refleciton principle used to analytically extend the solution
% at time t).
%

% Using this program, you may either generate an image or a video of a
% specific Loewner hull. The program calls for the user to input:
% 1. A driving function in t
% 2. A final time T at which the hull is generated
% 3. A number of subdivisions N
% 4. A choice to refine the approximation at the end of the interval [0,T]
% 5. The choice to generate an image or video of the hull forming
% If a video is chosen, the user is then prompted to specify a total number 
% of frames for the video as well as possible dimensions of the plot to 
% generate for each frame.

% Receiving Input:
%
% Driving function:
fprintf('Enter the driving function in t, beginning with "@(t)".\n');
function_input = input('Driving function: ', 's');
f = str2func(function_input);

% Upper bound of interval:
fprintf(['Enter the final time T > 0 for the hull to be approximated to. ' ...
    '\n']);
T = input('Your value: ');

% Number of subdivisions:
fprintf(['Enter an initial number of subdivisions N for the interval [0,T]. ' ...
    '\n']);
N = input('Your number of subdivisions: ');

% Choice to refine:
fprintf(['Do you want to refine the program at the end of the interval? ' ...
    '(N will increase)\n']);
refine_choice = input('Y/N: ', 's');

% Choice to generate an image or static video:
fprintf(['Would you like to generate an image of the time-T hull (1) or a ' ...
    'static video of the hull growth (2)?\n']);
image_video_choice = input(['Enter "1" for an image; enter "2" for a video: ' ...
    ''], 's');

% If a video is selected, the program gathers necessary information:
if image_video_choice == '2'
    % Total number of frames for the video, which turns out to be the 
    % number of plots to create:
    fprintf(['Enter a number of frames for the movie (Recommended: N/10). ' ...
        '\n']);
    frames = input('Your number of frames: ');

    % Choice to specify the dimensions for each plot in the video:
    fprintf('Would you like to set x and y dimensions for each plot? \n');
    dim_choice = input('Y/N: ', 's');

    % If the user chooses to specify dimensions, the program takes in the 
    % lower and upper limits for the x and y axes as input:
    if dim_choice == 'Y'
        x_min = input('Enter the lower limit on the x-axis: ');
        x_max = input('Enter the upper limit on the x-axis: ');
        y_min = input('Enter the lower limit on the y-axis: ');
        y_max = input('Enter the upper limit on the y-axis: ');
    end

    % Creates a VideoWriter object to make and export a video:
    v = VideoWriter('MBA_static_video.avi');
    open(v);
end

% We next want to subdivide the interval [0,T] N times, giving us a
% partition consisting of evenly spaced intervals. The linspace function
% gives an array t of N+1 evenly spaced points from 0 to T, thus giving us
% our desired subdivision:
t = linspace(0, T, N+1);

% If the user selects to refine at the end of the interval, the program 
% adds more subdivisions near the final time T to give us a better
% approximation:
if refine_choice == 'Y'
    % The program first declares a temporary variable to hold the initial 
    % number of subdivisions:
    J = N;

    % The program then skips to the end of the time interval [0,T]:
    i = floor(0.99*J);

    % The program finally iterates through the driving function values at 
    % the end of the interval:
    while i <= N+1
        % Adds new subdividisions until the difference of consecutive
        % driving function values is lower than a given tolerance: 1/3J
        while abs(f(t(i)) - f(t(i-1))) >= 1/(3*J)
            % Increases the number of subdivisions:
            N = N + 1;

            % To add a new subdivision, the program needs to shift the
            % values at the i+1,...,N+1 indices to the i+2,...,N+2 indices.
            % It does this by shifting each value t(j) from the j'th index 
            % to the (j+1)'th index, for j = i+1, ..., N+1. This leaves the 
            % i'th index free, so that a new subdivision may be added. 
            for x = N+1:-1:i+1
                t(x) = t(x-1);
            end
            
            % Lets the new time subdivision correspond to the average its
            % two adjacent time-values:
            t(i) = (t(i) + t(i-1))/2;
        end

        % Increases i for the next interval:
        i = i + 1;
    end
    
    % Prints the new number of sudvisions for the user:
    fprintf('Final number of subdivisions: %d\n', N);
end

% Uses the arrayfun MATLAB function to generate an array d of values
% d(i) = f(t(i)). These values correspond to the driving function at each 
% subdivision (or equivalently, each endpoint of the partition):
d = arrayfun(f, t);

% Creates array lists for values of the upper and lower hulls:
upper_list = zeros(1, N+1);
lower_list = zeros(1, N+1);

% Main loop:
% 
% The program finds an approximation for the user-input driving function on
% each interval in the partition (defined by the subdivision), and the 
% inverse comformal map associated to each partition. This is done using
% functions defined at the end of this m-file.
%
% Note that this loop iterates through each interval in the partition 
% backwards, thus allowing the inverse conformal maps for the j-th interval 
% to be applied to the 1, ... j-1 values in the upper_list and lower_list
% arrays ("shifting back an interval"). Lastly, the new values
% upper_list(j), lower_list(j) are calculated using another function at the
% end of this m-file.
%
for j = 1:N
    % First, calculates paramters used in BLEAT (Sec. 4.2) to approximate
    % the driving function on each interval:
    R = d(N+2-j) - d(N+1-j);
    s = t(N+2-j) - t(N+1-j);
    Y = sqrt(R.^2 + 16*s);
    
    % Next, calculates values used to generate the new "tips" of the hull 
    % (to be stored in upper_list(j) and lower_list(j)):
    c = (-1)*(R./sqrt(s));
    sqrt_c = sqrt(c.^2 + 16);

    % Calculates the parameters (depending on c) which will be used in the
    % computation of the hull tips for the current interval:
    delta = 0.5 - c./(2*sqrt_c);
    eps = 0.5 + c./(2*sqrt_c);
    D = (c + sqrt_c)./(2);
    E = (c - sqrt_c)./(2);
        
    % Uses the inverse_conformal_map function to apply the inverse 
    % conformal map associated to the current interval to all points added
    % to the upper and lower lists on previous intervals:
    if j > 1
        for k = 1:j-1
            upper_list(k) = h(upper_list(k), R, Y);
            lower_list(k) = h(lower_list(k), R, Y);
        end
    end
                
    % Calculates new "tips" for the hull and stores them in upper_list(j)
    % and lower_list(j):
    [upper_list(j), lower_list(j)]  = tip_value(delta, eps, D, E, s, R);
end

% Shifts the entire hull over by the value of the driving function at time
% t = 0:
upper_list = upper_list + f(0);
lower_list = lower_list + f(0);

% Creates code for plot colors:
str = '#216AB3';
color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

str0 = '#B32020';
color0 = sscanf(str0(2:end),'%2x%2x%2x',[1 3])/255;

% If the user specified to create a static video, creates a video:
if image_video_choice == '2'
    % Iterates through each frame in the video (the number of which was 
    % specified by the user) and creates the associate plot:
    for i = 1:frames
        % Upper bound for indices to use in the approximation of the hull
        % for this frame:
        k = N+1 - floor(N*(i/frames));

        % Makes a plot of the intermediate hull using the MATLAB plot 
        % function, with the previously-defined colors:
        plot(real(upper_list(N+1:-1:k)), imag(upper_list(N+1:-1:k)), '.', ...
            'Color', color);
        hold on;
        plot(real(lower_list(N+1:-1:k)), imag(lower_list(N+1:-1:k)), '.', ...
            'Color', color0);
        
        % Uses the dimensions specified by the user:
        if dim_choice == 'Y'
            xlim([x_min x_max]);
            ylim([y_min y_max]);
        end

        % Adds the plot to the video using the writeVideo function:
        writeVideo(v, getframe);
    end
    
    % Closes the file written to:
    close(v);
    
else
    % If the user_specifies an image of the hull at time T, uses the plot
    % function to create an image of the hull, with previously-defined 
    % colors:
    plot(real(upper_list), imag(upper_list), '.', 'Color', color);
    hold on;
    plot(real(lower_list), imag(lower_list), '.', 'Color', color0);
    
end

% Functions defined for this program:
% 
% Hard-coded formula to find the upper and lower tips of hulls on each
% interval (that is, a computation of the upper and lower tips of the
% approximating driving function):
function [u_tip, l_tip] = tip_value(delta, eps, D, E, s, R)
    % Hard-codes the modulus and argument of the upper tip with MATLAB
    % functions log and angle:
    tip_modulus = real(eps)*log(abs(D)) - imag(eps)*angle(D) + real(delta)*log(abs(E)) - imag(delta)*angle(-E);
    tip_argument = imag(eps)*log(abs(D)) + real(eps)*angle(D) + imag(delta)*log(abs(E)) + real(delta)*angle(-E);

    % Computes the upper and lower tips with the exp MATLAB function:
    u_tip = sqrt(s)*exp(tip_modulus + 1j*(tip_argument + (delta*pi))) + R;
    l_tip = sqrt(s)*exp(tip_modulus + 1j*(tip_argument - (delta*pi))) + R;
    
    %
    % In addition to the above code, there is a second computation which 
    % will provide the same result:
    %
    % tip_modulus = real(eps)*log(abs(D)) - imag(eps)*angle(-D) + real(delta)*log(abs(E)) - imag(delta)*angle(E);
    % tip_argument = imag(eps)*log(abs(D)) + real(eps)*angle(-D) + imag(delta)*log(abs(E)) + real(delta)*angle(E);
    % u_tip = sqrt(s)*exp(tip_modulus + 1j*(tip_argument + (eps*pi))) + R;
    % l_tip = sqrt(s)*exp(tip_modulus + 1j*(tip_argument - (eps*pi))) + R;
end

% Inverse conformal map formula for approximating driving functions used on
% each interval:
function inverse_conformal_map = h(z, R, Y)
    factor = 2*z + R + Y;
    inverse_conformal_map = (0.5)*factor*((2*z + R - Y)./factor).^(0.5-(R./(2*Y)));
end
