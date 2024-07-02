% 
% VSA_Dynamic_Videos.m
%
% Author: Jeffrey Utley
% Date: May 10, 2022
% Title: Dynamic Videos Using the Vertical Slit Algorithm
%
%
% This program uses the Vertical Slit Algorithm to create dynamic videos
% of the Loewner hull forming. The "dynamic video" is a sequence of evenly 
% distributed plots which show the right hull "morphing" into the left 
% hull. The program makes use the Duality Property which allows us to 
% calculate the right hull for a driving function as the rotation of a left
% hull for a slightly altered driving function.
%
% This program takes the following as input:
% 1. A driving function in t
% 2. A final time T at which the hull is generated
% 3. A number of subdivisions N
% 4. A choice to view a frame for each interval in the partition
% 5. A choice of lower and upper bounds for the dimensions of the plot
%
% This algorithm uses the same main loop of the Vertical Slit Algorithm
% (VSA), but stores the result after each iteration to be used in the movie
% generated. This is because the dynamic formation of the hull is reflected
% by the way that VSA "builds up" the hull. After generating both the right
% and left hulls while storing each "intermediate" right and left hull
% along the way, this program creates plots of these intermediate hulls.

% User input:
%
% Driving function:
fprintf('Enter the driving function in t, beginning with "@(t)". \n');
function_input = input('Driving function: ', 's');
driving_func = str2func(function_input);

% Upper bound of interval:
fprintf(['Enter the final time T > 0 for the hull to be approximated to.' ...
    '\n']);
T = input('Your value: ');

% Number of subdivisions:
fprintf(['Enter an initial number of subdivisions N for the interval [0,T]' ...
    '.\n']);
N = input('Your number of subdivisions: ');

% Choice to add a frame for every intermediate hull (that is, after every
% interval of the main loop) or not:
fprintf(['Do you want to view every frame? (if not, the movie consits of ' ...
    'every 50th frame)\n']);
frame_choice = input('Y/N: ', 's');

% Window size:
fprintf(['Enter a positive number for the dimensions of the square plot:' ...
    '\n']);
x1 = input('Enter the lower limit on the x-axis: ');
x2 = input('Enter the upper limit on the x-axis: ');
y1 = input('Enter the lower limit on the y-axis: ');
y2 = input('Enter the upper limit on the y-axis: ');

% The value in the variable frames is the number of indices to skip in
% between plots used in the movie:
if frame_choice == 'Y'
    % If the user wanted to view every frame, then every index contributes 
    % a plot to the movie:
    frames = 1;
else
    % If the user does not choose to view every frame, the program adds the
    % plot from every 50th index to the video:
    frames = floor(N/50);
end

% Creates the "modified driving function" which generates (a rotation of)
% the right hull as its left hull:
mod_driving_func = @(t) -1j*driving_func(T-t);

% We next want to subdivide the interval [0,T] N times, giving us a
% partition consisting of evenly spaced intervals. The linspace function
% gives an array t of N+1 evenly spaced points from 0 to T, thus giving us
% our desired subdivision:
t = linspace(0, T, N+1);

% Finds the length of each interval in the partition:
s = T/N;

% Computes a paramater based on s:
sq = 2i*sqrt(s);

% Uses the arrayfun MATLAB function to generate an array d of values
% d(i) = f(t(i)). These values correspond to the driving function at each 
% subdivision (or equivalently, each endpoint of the partition):
d = arrayfun(driving_func, t);

% Similarly, fills in values of the modified driving function at the
% endpoint of each partition:
d_mod = arrayfun(mod_driving_func, t);

% Creates an array for the list of tips of the hull (for the modified 
% driving function) at every time interval in the partition:
%
% Note: This three-dimensional array stores the hull generated after the 
% j-th interval in the following loop to the set mod_data(:,:,j). The upper
% hull is stored in mod_data(:,1,j) and the lower hull in mod_data(:,2,j).
% The k-th point in the upper hull added on the j-th interval will then be
% given by mod_data(k,1,j).
mod_data = zeros(N+1, 2, N);

% Main loop from VSA used by approximate the left hull driven by the
% modified driving function:
for j = 1:N
    % Sets the constant value for the approximating driving function on the 
    % current interval [t(N+1-j),t(N+2-j)] to the value of the modified
    % driving function mod_driving_func(t) at the lower endpoint t(N+1-j):
    c_mod = d_mod(N+1-j);
    
    % Uses the inverse conformal map associated to this interval on all
    % values previously added to mod_data:
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
            if imag(mod_data(i, 1, j)) > imag(c_mod)
                mod_data(i, 1, (j:N)) = g(mod_data(i, 1, j), c_mod, s);
            else
                mod_data(i, 1, (j:N)) = g_l(mod_data(i, 1, j), c_mod, s);
            end

            % If the lower hull value is in the halfplane above or below c:
            if imag(mod_data(i, 2, j)) < imag(c_mod)
                mod_data(i, 2, (j:N)) = g_l(mod_data(i, 2, j), c_mod, s);
            else
                mod_data(i, 2, (j:N)) = g(mod_data(i, 2, j), c_mod, s);
            end
        end
    end

    % Computes the new "tips" of the hull on this interval:
    mod_data(j, 1, (j:N)) = c_mod + sq;
    mod_data(j, 2, (j:N)) = c_mod - sq;
end

% Sets the value of the final point of both the upper and lower hulls on
% the last iteration to the value of the modified driving function at time
% 0:
mod_data(N+1, (1:2), N) = mod_driving_func(0);

% Multiplies the "modified" left hull by i to store the right hull (from 
% each iteration) into the array mod_data:
mod_data = 1j*mod_data;

% Creates an array for the input driving function, storing the hull after
% each interval the same as in mod_data:
driv_data = zeros(N+1, 2, N);

% Main loop from VSA used to approximate the Loewner hull driven by the
% original driving function:
for j = 1:N
    % Sets the constant value for the approximating driving function on the 
    % current interval [t(N+1-j),t(N+2-j)] to the value of the user-input
    % driving function driving_func(t) at the lower endpoint t(N+1-j):
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
            if imag(driv_data(i, 1, j)) > imag(c)
                driv_data(i, 1, (j:N)) = g(driv_data(i, 1, j), c, s);
            else
                driv_data(i, 1, (j:N)) = g_l(driv_data(i, 1, j), c, s);
            end

            % If the lower hull value is in the halfplane above or below c:
            if imag(driv_data(i, 2, j)) < imag(c)
                driv_data(i, 2, (j:N)) = g_l(driv_data(i, 2, j), c, s);
            else
                driv_data(i, 2, (j:N)) = g(driv_data(i, 2, j), c, s);
            end
        end
    end

    % Computes the new "tips" of the left hull on this interval:
    driv_data(j, 1, (j:N)) = c + sq;
    driv_data(j, 2, (j:N)) = c - sq;
end

% Sets the value of the final point of both the upper and lower hulls on
% the final iteration to the value of the (original) driving function at 
% time 0:
driv_data(N+1, (1:2), N) = driving_func(0);

% Creates a VideoWriter object to make and export a video:
v = VideoWriter('VSA_dynamic_video.avi');
open(v);

% Colors for plots:
str = '#216AB3';
color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

str0 = '#B32020';
color0 = sscanf(str0(2:end),'%2x%2x%2x',[1 3])/255;

str1 = '#115F0C';
color1 = sscanf(str1(2:end),'%2x%2x%2x',[1 3])/255;

str2 = '#524301';
color2 = sscanf(str2(2:end),'%2x%2x%2x',[1 3])/255;

% Final loop:
% 
% This loop goes through and adds plots to the video for each index to be
% included. That is, on the i-th interval of this loop, this program plots
% the i-th intermediate left hull (the approximation of the left hull after 
% the i-th loop) and the N+1-i -th intermediate right hull. By going
% backwards with the right hull and forward with the left hull, we see the
% right hull "deform" as the left grows.
for i = 1:N
    % On the first interval, just plots the full right hull:
    if i == 1
        % Plots the right hull using the MATLAB plot function:
        plot(real(mod_data(:,1,N)), imag(mod_data(:,1,N)), '.', 'Color', ...
            color0);
        hold on;
        plot(real(mod_data(:,2,N)), imag(mod_data(:,2,N)), '.', 'Color', ...
            color2);
        xlim([x1 x2]);
        ylim([y1 y2]);
        
        % Uses the WriteVideo function to add this plot as a frame to the 
        % VideoWriter object v:
        writeVideo(v, getframe);
            
        % Clears the plot:
        hold off; 
        
    % Only creates a plot if i is a multiple of "frames" (the variable
    % storing the number of indices to skip between creating each plot):
    elseif mod(i, frames) == 0
        % Plots the left hull:
        plot(real(driv_data((1:i-1),1,i-1)), ...
            imag(driv_data((1:i-1),1,i-1)),'.', 'Color', color);
        hold on;
        plot(real(driv_data((1:i-1),2,i-1)), ...
            imag(driv_data((1:i-1),2,i-1)), '.', 'Color', color1);
        xlim([x1 x2]);
        ylim([y1 y2]);

        % Keeps the left hull on the plot:
        hold on;
        
        % Plots the right hull:
        plot(real(mod_data((1:N+1-i),1,N+1-i)), ...
            imag(mod_data((1:N+1-i),1,N+1-i)), '.', 'Color', color0);
        hold on;
        plot(real(mod_data((1:N+1-i),2,N+1-i)), ...
            imag(mod_data((1:N+1-i),2,N+1-i)), '.', 'Color', color2);
        xlim([x1 x2]);
        ylim([y1 y2]);
        
        % Uses the WriteVideo function to add this plot as a frame to the 
        % VideoWriter object v:
        writeVideo(v, getframe);
        
        % Clears the frame:
        hold off;   
    end
end

% Plots the final-time left hull:
plot(real(driv_data(:,1,N)), imag(driv_data(:,1,N)), '.', 'Color', color);
hold on;
plot(real(driv_data(:,2,N)), imag(driv_data(:,2,N)), '.', 'Color', color1);
xlim([x1 x2]);
ylim([y1 y2]);

% Uses the WriteVideo function to add this plot as a frame to the 
% VideoWriter object v:
writeVideo(v, getframe);

% Closes the VideoWriter object in order to export the video:
close(v);

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

% This function is the formula for the conformal map for the halfplane
% below c, which corresponds to the branch cut (-2pi,0) being taken in the
% logarithm:
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