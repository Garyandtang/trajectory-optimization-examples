function [] = draw_result(soln, trajhandle)
xSoln = [soln.qSoln; soln.dqSoln];


params = sys_params;



%% **************************** FIGURES *****************************

h_fig = figure;
sz = [790 607]; % figure size
screensize = get(0,'ScreenSize');
xpos = ceil((screensize(3)-sz(1))/2); % center the figure on the screen horizontally
ypos = ceil((screensize(4)-sz(2))/2); % center the figure on the screen vertically
set(h_fig, 'Position', [xpos ypos sz])

h_3d = subplot(3,3,[1,2,4,5,7,8]);
axis equal
grid on
view(90,0);
ylabel('y [m]'); zlabel('z [m]');

% h_2d = subplot(1,2,2);
% plot_2d = plot(h_2d, 0, 0);
% grid on;
% xlabel('t [s]'); ylabel('z [m]');

quadcolors = lines(1);

set(gcf,'Renderer','OpenGL')

%% *********************** INITIAL CONDITIONS ***********************
t_total  = 5;             % Total simulated time
tstep    = 0.01;          % this determines the time step at which the solution is given
cstep    = 0.05;          % image capture time interval
max_iter = t_total/cstep; % max iteration
nstep    = cstep/tstep;
time     = 0; % current time
err = []; % runtime errors
% Get start and stop position
des_start = trajhandle(0,[]);
des_stop  = trajhandle(inf,[]);

% Get boundary
d_state = nan(max_iter,2);
for iter = 1:max_iter
    dd = trajhandle(cstep*iter,[]);
    d_state(iter,:) = dd.pos(1:2)';
end
y_lim = [min(d_state(:,1)) - 0.1, max(d_state(:,1)) + 0.1];
z_lim = [min(d_state(:,2)) - 0.1, max(d_state(:,2)) + 0.1];
if(4*(z_lim(2) - z_lim(1)) < y_lim(2) - y_lim(1))
    z_lim(1) = z_lim(1) - (y_lim(2) - y_lim(1))/8;
    z_lim(2) = z_lim(2) + (y_lim(2) - y_lim(1))/8;
end
stop_pos = des_stop.pos;
x0        = [des_start.pos; 0; des_start.vel; 0];
xtraj     = nan(max_iter*nstep, length(x0));
ttraj     = nan(max_iter*nstep, 1);

x         = x0;        % state

pos_tol = 0.01;
vel_tol = 0.03;
ang_tol = 0.05;

%% ************************* RUN SIMULATION *************************
disp('Simulation Running....')
% Main loop
for iter = 1:size(soln.qSoln,2)

  timeint = time:tstep:time+cstep;

  tic;
  % Initialize quad plot
  if iter == 1
    subplot(3,3,[1,2,4,5,7,8]);
    quad_state = simStateToQuadState(xSoln(:, 1));
    QP = QuadPlot(1, quad_state, params.arm_length, 0.05, quadcolors(1,:), max_iter, h_3d);
    ylim(y_lim); zlim(z_lim);
    quad_state = simStateToQuadState(x);
    QP.UpdateQuadPlot(quad_state, time);
    h_title = title(h_3d, sprintf('iteration: %d, time: %4.2f', iter, time));
  end
  
  quad_state = simStateToQuadState(xSoln(:, iter));
  time = soln.tSoln(iter);
  QP.UpdateQuadPlot(quad_state, time);
    
  pause(0.1)
end

end
