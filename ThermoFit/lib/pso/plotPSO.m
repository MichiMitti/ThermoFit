clc;
clear all;
close all;

%% Plot file for PSO

%% Creates animated Plots and Videos for Benchmark Functions
% list of all csv files
FileList = struct2cell(dir('*.csv'));

%loops over all csv files
% for i =  1:length(FileList(1,:))
%     filename = string(FileList{1,i});
%     createPlots(filename);
% end
createPlots("log_rastrigin_without_rl.csv");

%% Function to create plot from filename
function [returnf] = createPlots(pathToCSV)

[filepath,name,ext] = fileparts(pathToCSV);
videoName = name;

pause_video = 0.1;

%% Initialize video
delete(videoName + '.mp4')
logVideo = VideoWriter(videoName,'MPEG-4'); %open video file
logVideo.FrameRate = 5;  %can be adjusted is speed of video
open(logVideo)

%% Read in "log.csv"
results_table = readtable(pathToCSV);
results_array = table2array(results_table);

dimensions = table2array(readtable(pathToCSV,'Range','B9:B9'));
boundaries = table2array(readtable(pathToCSV,'Range','B10:C10'));
swarm_size = table2array(readtable(pathToCSV,'Range','B11:B11'));

header = readcell(pathToCSV,"Range",[1,1,13,2*dimensions+2]);
problem_description = header(10,1);


%% Selection benchmark problems
resolution = 1000;

x = linspace(boundaries(1),boundaries(2),resolution);
y = linspace(boundaries(1),boundaries(2),resolution);
[X,Y] = meshgrid(x,y);

if problem_description == "=== Benchmark Sphere ==="
    Z = reshape(sphere([X(:)'; Y(:)']), resolution, resolution);
    root = [0,0];
elseif  problem_description == "=== Benchmark Rosenbrock ==="
    null = 1;
    Z = reshape(rosenbrock([X(:)'; Y(:)'],null), resolution, resolution);
    root = [null,null^2];
elseif problem_description == "=== Benchmark Rastrigin ==="
    Z = reshape(rastrigin([X(:)'; Y(:)']), resolution, resolution);
    root = [0,0];
end


%% Plotting

% creates figure
fig = figure('Name', "Plot of Data");
set(gcf,'color','w');
hold on

% sets position and size of plot
x0=10;
y0=10;
width=1000;
height=1000;
fig.Position = [x0,y0,width,height];

% calculates scaling for colorbar
max_fval = max(max(Z));;
power10 = ceil(log10(max_fval));
scale10 = logspace(0, power10, 10);

% creates contour f polt for benchmark problem with root
benchmark_contour = contourf(X,Y,Z, scale10);
colormap(flipud(hot));
colorbar;
clim([1 10^power10]);
set(gca,'ColorScale','log');
bench_root = scatter(root(1),root(2), 200,'+','MarkerEdgeColor','b','LineWidth',2);

% sets title for benchmark problem
tit = title(problem_description(1,1));
tit.FontSize = 20;

% shows swarm size
text_size = annotation('textbox', [0, 0.75, 0, 0], 'String', ...
     "Swarm size: " + swarm_size);
text_size.FontSize = 14;

% setps to next iteration
l = 1;
lines_to_next_iter = swarm_size + 1;

iterations = (length(results_array(:,1))-4)/(swarm_size+1);
if iterations < 200
    loopsize = iterations;
else
    loopsize = 200;
end

scatter_particle_size = 70;

% loop over all iterations
for i = 1:loopsize -1
    
    % current position
    px = results_array(l:l+swarm_size,2);
    py = results_array(l:l+swarm_size,3);
    
    % best found swarm position
    bestx = results_array(l,4);
    besty = results_array(l,5);
    
    % shows iteration and function value (of best paritcle)
    text_iter = annotation('textbox', [0, 0.6, 0, 0], 'String', ...
        "Iteration: " + results_array(l,2));
    text_iter.FontSize = 14;
    text_fval = annotation('textbox', [0, 0.5, 0, 0], 'String', ...
            " Function value: " + results_array(l,4+dimensions));
    text_fval.FontSize = 14;

    % goes to next iteration
    l = l + lines_to_next_iter;
    
    % current velocity
    vx = results_array(l:l+swarm_size,2+dimensions);
    vy = results_array(l:l+swarm_size,3+dimensions);

    % scatter plot for all particles with quiver (velocity)
    position = scatter(px,py,scatter_particle_size,'MarkerEdgeColor','k','MarkerFaceColor','k');
    velocity = quiver(px, py, vx, vy, 0,'Color','k');

    % scatter plot for best position
    swarm_best = scatter(bestx,besty, 200,'x','MarkerEdgeColor','k','LineWidth',2);

    % write frame to video
    pause(pause_video)              % wait for next frame
    frame = getframe(gcf);          % get frame
    writeVideo(logVideo, frame);     % write to video

    % deletes old plots
    delete(position);
    delete(velocity);
    delete(swarm_best);
    delete(text_iter);
    delete(text_fval);
end


% shows iteration and function value (of best paritcle)
text_iter = annotation('textbox', [0, 0.6, 0, 0], 'String', ...
    "Iteration: " + results_array(l,2));
text_iter.FontSize = 14;
text_fval = annotation('textbox', [0, 0.5, 0, 0], 'String', ...
     " Function value: " + results_array(l,4+dimensions));
text_fval.FontSize = 14;
text_feval = annotation('textbox', [0, 0.35, 0, 0], 'String', ...
     " Function calls: " + results_array(end - 1, 2));
text_feval.FontSize = 14;


scatter(bestx,besty, 200,'x','MarkerEdgeColor','k','LineWidth',3)

pause(pause_video)                      % wait for next frame

frame = getframe(gcf); %get frame
writeVideo(logVideo, frame);

close(logVideo)

close all;

returnf = 0;
end


%% Test Functions
function [y] = sphere(x)
    y = sum(x.^2);
end

function [y] = rastrigin(x)
    l = length(x);
    A = 10; 
    y = zeros(1, l);

    for i = 1:l
        y(i) = 2*A + ( x(1,i).^2 - A * cos(2*pi*x(1,i)) + (x(2,i).^2 - A * cos(2*pi*x(2,i))) );
    end

end

function [y] = rosenbrock(x, null)
% zero at b^2, boundaries ]-inf, inf[
  
    a = 100;
    b = null;
    l = length(x);
    y = zeros(1, l);

    for i = 1:l
        y(i) = a * sum((x(1,i)-x(2,i).^2).^2 + (b-x(1,i)).^2);
    end
end
