% 
% Vertical_Slit_Algorithm.m
%
% Author: Jeffrey Utley
% Date: May 10, 2022
% Title: Approximation of Complex-Driven Loewner Hulls by Vertical Slits
% 
%
% As the Modified BLEAT Algorithm builds off of the algorithm described in
% the Section 4.2 of the paper BLEAT: The Curve-Creating Black Box to 
% approximate Loewner hulls will complex-valued driving functions, the 
% Vertical Slit Algorithm builds off of the algorithm described in Section 
% 4 of the paper "Numerical computations for the Schramm-Loewner Evolution"
% by Tom Kennedy to approximate the same object.
%
% The algorithm described in Tom Kennedy's paper approximates the Loewner 
% hull (restricted to the upper halfplane) with real-valued driving
% functions. This is done by approximating the driving function on small 
% intervals by constant driving functions, for which the solution to the 
% Loewner Equation can be explicitly computed. This algorithm analytically 
% extends the solution of the LE at each time to the entire complex plane 
% by Schwarz reflection; to do this, the computation is hard-coded in a 
% function at the end of this file. Thus, users can enter both real-valued 
% and complex-valued driving functions and get a hull after some specified
% time.
%

% Using this program, you may either generate an image or a video of a
% specific Loewner hull. The program calls for the user to input:
% 1. A driving function in t
% 2. A final time T at which the hull is generated
% 3. A number of subdivisions N
% 4. The choice to generate an image or video of the hull forming
% If a video is chosen, the user then specifies a total number of frames 
% for the video as well as possible dimensions of the plot to generate for 
% each frame.

% User input:
%
% Driving function:
fprintf('Enter the driving function in t, beginning with "@(t)".\n');
function_input = input('Driving function: ', 's');
f = str2func(function_input);

% Upper bound of interval:
fprintf(['Enter the final time T > 0 for the hull to be approximated to.' ...
    '\n']);
T = input('Your value: ');

% Number of subdivisions:
fprintf('Enter a number of subdivisions N for the interval [0,T].\n');
N = input('Your number of subdivisions: ');

% Choice to generate an image or static video:
fprintf(['Would you like to generate an image of the time-T hull (1) or a ' ...
    'static video of the hull growth (2)?\n']);
image_video_choice = input(['Enter "1" for an image; enter "2" for a video: ' ...
    ''], 's');

% If a video is selected, the program gathers necessary information:
if image_video_choice == '2'
    % Total number of frames for the video which turns out to be the 
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
    v = VideoWriter('VSA_static_video.avi');
    open(v);
end

% We next want to subdivide the interval [0,T] N times, giving us a
% partition consisting of evenly spaced intervals. The linspace function
% gives an array t of N+1 evenly spaced points from 0 to T, thus giving us
% our desired subdivision:
t = linspace(0, T, N+1);

% Computes the length s of each interval in the partition defined by t:
s = T/N;

% Computes a paramater based on s:
sq = 2i*sqrt(s);

% Uses the arrayfun MATLAB function to generate an array d of values
% d(i) = f(t(i)). These values correspond to the driving function at each 
% subdivision (or equivalently, each endpoint of the partition):
d = arrayfun(f, t);

% Creates array lists for values of the upper and lower hulls:
upper_list = zeros(1, N+1);
lower_list = zeros(1, N+1);

% Main loop:
% 
% The program approximates the user-input driving function on each interval 
% in the partition (defined by the subdivision) by a constant driving 
% function at the lower endpoint. The program can then explicitly compute
% the inverse conformal map associated to each interval. This is done using
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
    % Sets the constant value for the approximating driving function on the 
    % current interval [t(N+1-j),t(N+2-j)] to the value of the user-input
    % driving function f(t) at the lower endpoint t(N+1-j):
    c = d(N+1-j);
    
    % Uses the inverse conformal map associated to this interval on all
    % values previously added to upper_list and lower_list:
    if j > 1
        % Since the solution at time t to the Loewner Equation with 
        % constant driving function at c is different for the halfplanes 
        % above and below the line {z|Im(z)=Im(c)}, the program checks on 
        % the position of each value currently in the upper and lower hulls
        % to decide which inverse conformal map to use. The functions g and
        % g_l, located at the end of this m-file, hard-code the maps
        % corresponding to the halfplanes above and below c, respectively:
        for i = 1:j-1
            % If the upper hull value is in the halfplane above or below c:
            if imag(upper_list(i)) > imag(c)
                upper_list(i) = g(upper_list(i), c, s);
            else
                upper_list(i) = g_l(upper_list(i), c, s);
            end

            % If the lower hull value is in the halfplane above or below c:
            if imag(lower_list(i)) < imag(c)
                lower_list(i) = g_l(lower_list(i), c, s);
            else
                lower_list(i) = g(lower_list(i), c, s);
            end
        end
    end
        
    % Computes the new "tips" of the hull on this interval:
    upper_list(j) = c + sq;
    lower_list(j) = c - sq;
end

% Sets the final value in the upper and lower hulls to be the point f(0),
% the user-input driving function at time 0:
upper_list(N+1) = f(0);
lower_list(N+1) = f(0);

% Creates colors for the plot:
str = '#216AB3';
color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

str0 = '#B32020';
color0 = sscanf(str0(2:end),'%2x%2x%2x',[1 3])/255;

% If the user specified to create a video make plots, creates a video:
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
% This function hard-codes the conformal map associated to the halfplane
% above c, which corresponds to the branch cut (0, 2pi) being chosen in the
% logarithm used in the formula for the conformal map:
function upper_value = g(z, c, s)
    % Uses the MATLAB functions log and abs:
    modulus = 0.5*log(abs((z-c).^2 - 4*s));

    % Uses the MATLAB functions mod and angle to set the branch cut to 
    % (0,2pi):
    argument = mod(angle((z-c).^2 - 4*s),2*pi);

    % Finds the upper tip of the hull for the current interval:
    upper_value = c + exp(modulus + 0.5j*argument);
end

% This function is the formula for the conformal map associated to the 
% halfplane below c, which corresponds to the branch cut [-2pi,0) being 
% taken in the logarithm:
function lower_value = g_l(z, c, s)
    % Uses the MATLAB functions log and abs:
    modulus = 0.5*log(abs((z-c).^2 - 4*s));

    % Uses the MATLAB functions mod and angle to (temporarily) set the 
    % branch cut to (0,2pi):
    argument = mod(angle((z-c).^2 - 4*s),2*pi);

    % Finds the lower tip of the hull for the current interval by 
    % subtracting the argument by 2pi to set the branch cut to (-2pi,0):
    lower_value = c + exp(modulus + 0.5j*(argument - 2*pi));
end