%
% MBA_Dynamic_Videos.m
%
% Author: Jeffrey Utley
% Date: May 10, 2022
% Title: Dynamic Videos Using the Modified BLEAT Algorithm
%
%
% This program uses the Modified BLEAT Algorithm to create dynamic videos
% of the Loewner hull forming. The "dynamic video" is a sequence of evenly 
% distributed plots which show the right hull "morphing" into the left 
% hull. The program makes use the Duality Property, which allows us to 
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
% This algorithm uses the same main loop of the Modified BLEAT Algorithm 
% (MBA), but stores the result after each iteration to be used in the movie
% generated. This is because the dynamic formation of the hull is reflected
% by the way that MBA "builds up" the hull. After generating both the right
% and left hulls while storing each "intermediate" right and left hull
% along the way, this program creates plots of these intermediate hulls.

% Receiving input:
%
% Driving function:
fprintf('Enter the driving function in t, beginning with "@(t)".\n');
driving_func_input = input('Driving function: ', 's');
driving_func = str2func(driving_func_input);

% Upper bound of interval:
fprintf('Enter the final time T > 0 for the hull to be approximated to.\n');
T = input('Your value: ');

% Number of subdivisions:
fprintf(['Enter an initial number of subdivisions N for the interval [0,T]' ...
    '.\n']);
N = input('Your value: ');

% Choice to add a frame for every intermediate hull (that is, after every
% interval of the main loop) or not:
fprintf(['Do you want to view every frame? (if not, the movie consits of ' ...
    'every 50th frame)\n']);
frame_choice = input('Y/N: ', 's');

% Window size:
fprintf('Enter values for the dimensions of the square plot:\n');
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

% Main loop from MBA used by approximate the left hull driven by the
% modified driving function:
for j = 1:N
    % First, calculates paramters used in BLEAT (Sec. 4.2) to approximate
    % the driving function on each interval:
    R_mod = d_mod(N+2-j) - d_mod(N+1-j);
    Y_mod = sqrt(R_mod.^2 + 16*s);

    % Next, calculates values to generate the new "tips" of the hull (to be
    % stored in upper_list(j) and lower_list(j)):
    c_mod = (-1)*(R_mod./sqrt(s));
    sqrt_cmod = sqrt(c_mod.^2 + 16);
    
    % Calculates the parameters (depending on c) which will be used in the
    % computation of the hull tips for the current interval:
    delta_mod = 0.5 - c_mod./(2*sqrt_cmod);
    eps_mod = 0.5 + c_mod./(2*sqrt_cmod);
    D_mod = (c_mod + sqrt_cmod)./(2);
    E_mod = (c_mod - sqrt_cmod)./(2);

    % Uses the inverse_conformal_map function to apply the inverse 
    % conformal map associated to the current interval to all points added
    % to the upper and lower lists on previous intervals:
    if j > 1
        for k = 1:j-1
            mod_data(k, 1, (j:N)) = h(mod_data(k, 1, j), R_mod, Y_mod);
            mod_data(k, 2, (j:N)) = h(mod_data(k, 2, j), R_mod, Y_mod);
        end
    end
        
    % Calculates new "tips" for the hull and stores them in 
    % mod_data(j, 1, (j:N)) and mod_data(j, 2, (j:N)):
    [mod_data(j, 1, (j:N)), mod_data(j, 2, (j:N))] = tip_value(delta_mod, ...
        eps_mod, D_mod, E_mod, s, R_mod);
end

% Multiplies the "modified" left hull by i to store the right hull (from 
% each iteration) into the array mod_data:
mod_data = 1j*mod_data;

% Creates an array for the input driving function, storing the hull after
% each interval the same as in mod_data:
driv_data = zeros(N+1, 2, N);

% Main loop from MBA used to approximate the Loewner hull driven by the
% original driving function:
for j = 1:N
    % First, calculates paramters used in BLEAT (Sec. 4.2) to approximate
    % the driving function on each interval:
    R = d(N+2-j) - d(N+1-j);
    Y = sqrt(R.^2 + 16*s);

    % Next, calculates values to generate the new "tips" of the hull (to be
    % stored in upper_list(j) and lower_list(j)):
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
            driv_data(k, 1, (j:N)) = h(driv_data(k, 1, j), R, Y);
            driv_data(k, 2, (j:N)) = h(driv_data(k, 2, j), R, Y);
        end
    end
        
    % Calculates new "tips" for the hull and stores them in 
    % driv_data(j, 1, (j:N)) and driv_data(j, 2, (j:N)):
    [driv_data(j, 1, (j:N)), driv_data(j, 2, (j:N))] = tip_value(delta, ...
        eps, D, E, s, R);
end

% Creates a VideoWriter object to make and export a video:
v = VideoWriter('MBA_dynamic_video.avi');
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
% included. Since the Modified BLEAT Algorithm approximates each hull up to
% translated, the left and right hulls must first be shifted by the value of
% the driving function at the left endpoint of each interval in the
% partition. Then, on the i-th interval of this loop, this program plots the
% the i-th intermediate left hull (the approximation of the left hull after 
% the i-th loop) and the N+1-i -th intermediate right hull. By going
% backwards with the right hull and forward with the left hull, we see the 
% right hull "deform" as the left grows.
for i = 1:N
    % On the first interval, just plots the full right hull:
    if i == 1
        % Shifts the right hull:
        mod_data(:,:,N) = mod_data(:,:,N) + driving_func(t(N));
        
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
        % Shifts the left and right hulls by the value of the driving 
        % function correpsonding to the current interval:
        driv_data(:,:,i-1) = driv_data(:,:,i-1) + driving_func(t(N+1-i));
        mod_data(:,:,N+1-i) = mod_data(:,:,N+1-i) + driving_func(t(N+1-i));
            
        % Plots the left hull:
        plot(real(driv_data(:,1,i-1)), imag(driv_data(:,1,i-1)), '.', ...
            'Color', color);
        hold on;
        plot(real(driv_data(:,2,i-1)), imag(driv_data(:,2,i-1)), '.', ...
            'Color', color1);
        xlim([x1 x2]);
        ylim([y1 y2]);

        % Keeps the left hull on the plot:
        hold on;
        
        % Plots the right hull:
        plot(real(mod_data(:,1,N+1-i)), imag(mod_data(:,1,N+1-i)), '.', ...
            'Color', color0);
        hold on;
        plot(real(mod_data(:,2,N+1-i)), imag(mod_data(:,2,N+1-i)), '.', ...
            'Color', color2);
        xlim([x1 x2]);
        ylim([y1 y2]);
        
        % Uses the WriteVideo function to add this plot as a frame to the 
        % VideoWriter object v:
        writeVideo(v, getframe);
        
        % Clears the frame:
        hold off;   
    end
end

% Shifts the final-time left hull to its correct place:
driv_data(:,:,N) = driv_data(:,:,N) + driving_func(t(1));

% Plots the final-time left hull:
plot(real(driv_data(:,1,N)), imag(driv_data(:,1,N)), '.', 'Color', color);
hold on;
plot(real(driv_data(:,2,N)), imag(driv_data(:,2,N)), '.', 'Color', color1);
xlim([x1 x2]);
ylim([y1 y2]);

% Adds the plot to the video:
writeVideo(v, getframe);

% Closes the VideoWriter object in order to export the video:
close(v);

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

% Inverse conformal map formula from BLEAT for approximating driving 
% functions used on each interval:
function inverse_conformal_map = h(z, R, Y)
    factor = 2.*z + R + Y;
    inverse_conformal_map = (0.5).*factor.*((2.*z + R - Y)./factor).^(0.5-(R./(2.*Y)));
end
