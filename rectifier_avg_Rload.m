%% run_rectifier_averaged_17state.m
% Averaged nonlinear model of the single-phase rectifier stage.
%
% State order:
% x(1)  = v_alpha
% x(2)  = v_beta
% x(3)  = e2
% x(4)  = theta
% x(5)  = i_alpha
% x(6)  = i_beta
% x(7)  = e3
% x(8)  = e4
% x(9)  = e5
% x(10) = iL_d
% x(11) = iL_q
% x(12) = vC_d
% x(13) = vC_q
% x(14) = i_dc
% x(15) = vC_dc
% x(16) = ig_d
% x(17) = ig_q
%
% Sign conventions:
%   iL   > 0 : converter -> AC terminal/grid
%   ig   > 0 : AC capacitor node -> grid
%   i_dc > 0 : DC load/output node -> DC-link node/converter
%   i_b  > 0 : DC-link node -> H-bridge
%
% IMPORTANT:
%   The grid current is modeled dynamically through Lg + Rg.
%   vg_d and vg_q are stiff-grid dq voltage inputs specified by the user.
%   A shunt load resistor Rload closes the DC output after Ldc + RLdc.

clear;
clc;
% close all;

%% 1) USER PARAMETERS
p.tStop = 1.5;
p.dtOut = 1e-5;
p.maxStep = 1e-4;

p.Kad = 1;  % virtual resistance for the LCL capacitor current feedback. Set to 0 to disable.

p.f0 = 60;
p.w0 = 2*pi*p.f0;
p.kSOGI = 1;            % SOGI damping factor

p.k1 = 100;               % PLL proportional
p.k2 = 1;               % PLL integral
p.k3 = 0.1;               % DC-voltage proportional
p.k4 = 5;               % DC-voltage integral
p.k5 = 0.5;               % d-current integral
p.k6 = 20;               % d-current proportional
p.k7 = 0.5;               % q-current integral
p.k8 = 20;               % q-current proportional

p.L = 800e-6;                % AC converter-side inductance [H]
p.RL = 0;                 % AC inductor resistance [ohm]
p.C = 100e-6;                % AC shunt capacitance [F]

p.Lg = 800e-6;               % grid-side inductance [H]
p.Rg = 0.2;               % grid-side inductor ESR [ohm]
p.vg_d = 170;             % stiff-grid d-axis voltage [V]
p.vg_q = 0;             % stiff-grid q-axis voltage [V]

p.Ldc  = 1e-6;             % DC series inductance [H]
p.RLdc = 0.1;             % DC inductor ESR [ohm]
p.Cdc  = 400e-6;             % DC-link capacitance [F]
p.RCdc = 0.01;             % DC capacitor ESR [ohm]
p.Rsh  = 1e6;             % DC shunt resistance [ohm]
p.Rload = 88;            % DC output/load resistance [ohm]

p.vdc_ref = 300;
p.iq_ref  = 0;
p.igq_ref = 0;
p.vd_ff = 170;            % constant feedforward d voltage [V]
p.vq_ff = 0;            % constant feedforward q voltage [V]

% the code iteratively solves for a consistent vdc at every ODE evaluation
p.vdcFloor = 1;       % prevent devision by 0 when calculating mod-index
p.algTol = 1e-10;       % convergence tol for algebraic-loop iteration
p.algMaxIter = 50;      % maximum number of iterations for solving algebraic loop at one ODE eval
p.algRelax = 0.7;       % relaxation factor to prevent oscillatory convergence

%% 2) NETWORK CLOSURE
% Grid-side current is generated dynamically by Lg + Rg.
% The stiff-grid dq voltages are specified above as p.vg_d and p.vg_q.
% The DC output is closed by p.Rload, connected from the node after
% Ldc + RLdc to the DC return. No external DC-side voltage input is used.

%% 3) INITIAL CONDITIONS
x0 = zeros(17,1);
x0(4)  = 0;               % theta [rad]
x0(15) = p.vdc_ref;       % ideal DC capacitor voltage [V]
x0(16) = 0;               % grid-side d-axis current [A]
x0(17) = 0;               % grid-side q-axis current [A]

%% 4) VERIFY REQUIRED PARAMETERS
required = {'kSOGI','k1','k2','k3','k4','k5','k6','k7','k8', ...
            'L','C','Lg','Rg','vg_d','vg_q','Ldc','RLdc','Cdc','RCdc','Rsh','Rload','vd_ff','vq_ff'};
missing = {};
for k = 1:numel(required)
    val = p.(required{k});
    if ~(isscalar(val) && isfinite(val))
        missing{end+1} = required{k}; %#ok<SAGROW>
    end
end
if ~isempty(missing)
    error(['Enter finite numeric values for these parameters before running:' newline ...
           strjoin(missing, ', ')]);
end

%% 5) RUN AVERAGED MODEL
tEval = (0:p.dtOut:p.tStop).';
opts = odeset('RelTol',1e-3,'AbsTol',1e-8,'MaxStep',p.maxStep);

fprintf('Running averaged rectifier model from 0 to %.3f s...\n',p.tStop);
[t,x] = ode15s(@(t,x) rectifierODE(t,x,p),tEval,x0,opts);
fprintf('Averaged-model integration complete.\n');

%% 6) RECONSTRUCT ALL SIGNALS
avg = reconstructSignals(t,x,p);
avg.x = x;
save('rectifier_averaged_results.mat','avg','p');

%% 7) BASIC PLOTS
fprintf('Plotting the terminal signals...');
figure('Name','Averaged Rectifier - Terminal Signals');
tiledlayout(3,2);

nexttile;
plot(avg.t,avg.v_grid_alpha);
grid on;
xlabel('Time [s]'); ylabel('v_{AC} [V]');
title('AC terminal voltage');

nexttile;
plot(avg.t,avg.i_ac_terminal);
grid on;
xlabel('Time [s]'); ylabel('i_{AC} [A]');
title('AC terminal current');

nexttile;
plot(avg.t,avg.v_ac_terminal);
grid on;
xlabel('Time [s]'); ylabel('Voltage [V]');
title('Capacitor voltage');

nexttile;
plot(avg.t,avg.iL_alpha);
grid on;
xlabel('Time [s]'); ylabel('Current [A]');
title('Converter-side inductor current');

nexttile;
plot(avg.t,avg.v_dc);
grid on;
xlabel('Time [s]'); ylabel('v_{dc} [V]');
title('Measured DC-link voltage');

nexttile;
plot(avg.t,avg.i_dc);
grid on;
xlabel('Time [s]'); ylabel('i_{dc} [A]');
title('DC-link inductor current');


fprintf('Plotting the controller signals...');
figure('Name','Averaged Rectifier - Controller Signals');
tiledlayout(3,2);

nexttile;
plot(avg.t,avg.theta);
grid on;
xlabel('Time [s]'); ylabel('\theta [rad]');
title('PLL angle');

nexttile;
plot(avg.t,avg.vhat_q);
grid on;
xlabel('Time [s]'); ylabel('\hat{v}_q [V]');
title('PLL q-axis error');

nexttile;
plot(avg.t,avg.id_ref,avg.t,avg.ihat_d);
grid on;
xlabel('Time [s]'); ylabel('Current [A]');
legend('i_d^*','\hat{i}_d');
title('d-axis current loop');

nexttile;
plot(avg.t,avg.iq_ref,avg.t,avg.ihat_q);
grid on;
xlabel('Time [s]'); ylabel('Current [A]');
legend('i_q^*','\hat{i}_q');
title('q-axis current loop');

nexttile;
plot(avg.t,avg.m);
grid on;
xlabel('Time [s]'); ylabel('m');
title('Modulation index');

nexttile;
plot(avg.t,avg.vconv_alpha_ref,avg.t,avg.vconv_alpha);
grid on;
xlabel('Time [s]'); ylabel('Voltage [V]');
legend('v_{conv,\alpha}^*','v_{conv,\alpha}');
title('Converter voltage command vs average');

figure('Name','LCL dq Current Diagnostics');
tiledlayout(2,1);

nexttile;
plot(avg.t,avg.iL_d,'LineWidth',1.2); hold on;
plot(avg.t,avg.ig_d,'LineWidth',1.2);
grid on;
xlabel('Time [s]');
ylabel('Current [A]');
legend('i_{L,d}','i_{g,d}');
title('d-axis currents');
xlim([0 p.tStop]);

nexttile;
plot(avg.t,avg.iL_q,'LineWidth',1.2); hold on;
plot(avg.t,avg.ig_q,'LineWidth',1.2);
grid on;
xlabel('Time [s]');
ylabel('Current [A]');
legend('i_{L,q}','i_{g,q}');
title('q-axis currents');
xlim([0 p.tStop]);

% figure;
% plot(avg.t,avg.mstar);
% hold on;
% plot(avg.t,avg.m);
% yline(1,'--');
% yline(-1,'--');
% grid on;
% legend('m^*','m');

% PLL estimated angular frequency
avg.omega_hat = p.w0 ...
              + p.k1*avg.vhat_q ...
              + p.k2*avg.e2;

figure('Name','PLL Frequency');
plot(avg.t,avg.omega_hat,'LineWidth',1.2);
hold on;
yline(p.w0,'--','\omega_0');
grid on;

xlabel('Time [s]');
ylabel('Angular frequency [rad/s]');
legend('\hat{\omega}','\omega_0');
title('PLL Estimated Frequency vs Nominal Frequency');
xlim([0 p.tStop]);


figure('Name','Grid Voltage Reference Check');
plot(avg.t,avg.v_grid_alpha,'LineWidth',1.2);
hold on;
plot(avg.t,avg.v_grid_alpha_theta,'--','LineWidth',1.2);
grid on;

xlabel('Time [s]');
ylabel('Voltage [V]');
legend('\phi = \omega_0t','\theta');
title('Grid Voltage Reconstruction');



%% 8) OPTIONAL SIMULINK COMPARISON
% Set cmp.enabled = true and replace the model/signal names.
% The Simulink signals must be logged to logsout.
cmp.enabled = false;
cmp.modelName = 'your_switched_model';
cmp.vacLogName = 'v_ac_terminal';
cmp.iacLogName = 'i_ac_terminal';
cmp.vdcLogName = 'v_dc';
cmp.idcLogName = 'i_dc';

if cmp.enabled
    fprintf('Running explicitly-switched Simulink model...\n');
    simIn = Simulink.SimulationInput(cmp.modelName);
    simIn = simIn.setModelParameter('StopTime',num2str(p.tStop));
    swOut = sim(simIn);
    logs = swOut.logsout;

    sw.vac = getLoggedTimeseries(logs,cmp.vacLogName);
    sw.iac = getLoggedTimeseries(logs,cmp.iacLogName);
    sw.vdc = getLoggedTimeseries(logs,cmp.vdcLogName);
    sw.idc = getLoggedTimeseries(logs,cmp.idcLogName);

    figure('Name','Averaged vs Explicitly Switched');
    tiledlayout(2,2);

    nexttile;
    plot(sw.vac.Time,sw.vac.Data); hold on;
    plot(avg.t,avg.v_ac_terminal,'LineWidth',1.2);
    grid on; xlabel('Time [s]'); ylabel('Voltage [V]');
    legend('Switched','Averaged'); title('AC terminal voltage');

    nexttile;
    plot(sw.iac.Time,sw.iac.Data); hold on;
    plot(avg.t,avg.i_ac_terminal,'LineWidth',1.2);
    grid on; xlabel('Time [s]'); ylabel('Current [A]');
    legend('Switched','Averaged'); title('AC terminal current');

    nexttile;
    plot(sw.vdc.Time,sw.vdc.Data); hold on;
    plot(avg.t,avg.v_dc,'LineWidth',1.2);
    grid on; xlabel('Time [s]'); ylabel('Voltage [V]');
    legend('Switched','Averaged'); title('DC-link voltage');

    nexttile;
    plot(sw.idc.Time,sw.idc.Data); hold on;
    plot(avg.t,avg.i_dc,'LineWidth',1.2);
    grid on; xlabel('Time [s]'); ylabel('Current [A]');
    legend('Switched','Averaged'); title('DC terminal current');
end

%% LOCAL FUNCTIONS
function dx = rectifierODE(t,x,p)
    a = rectifierAlgebraic(t,x,p);

    v_alpha = x(1);
    v_beta  = x(2);
    i_alpha = x(5);
    i_beta  = x(6);
    iL_d    = x(10);
    iL_q    = x(11);
    vC_d    = x(12);
    vC_q    = x(13);
    i_dc    = x(14);
    vC_dc   = x(15);
    ig_d     = x(16);
    ig_q     = x(17);

    dx = zeros(17,1);

    dx(1) = p.kSOGI*p.w0*(a.v_ac_terminal-v_alpha)-p.w0*v_beta;
    dx(2) = p.w0*v_alpha;
    dx(3) = a.vhat_q;
    dx(4) = p.w0+p.k1*a.vhat_q+p.k2*x(3);

    dx(5) = p.kSOGI*p.w0*(a.iL_alpha-i_alpha)-p.w0*i_beta;
    dx(6) = p.w0*i_alpha;

    dx(7) = p.vdc_ref-a.v_dc;
    dx(8) = a.id_ref-a.ihat_d;
    dx(9) = p.iq_ref-a.ihat_q;

    dx(10) = (a.vconv_d_phi-vC_d-p.RL*iL_d)/p.L+p.w0*iL_q;
    dx(11) = (a.vconv_q_phi-vC_q-p.RL*iL_q)/p.L-p.w0*iL_d;

    dx(12) = (iL_d-ig_d)/p.C+p.w0*vC_q;
    dx(13) = (iL_q-ig_q)/p.C-p.w0*vC_d;

    dx(14) = (-a.v_dc-(p.RLdc+p.Rload)*i_dc)/p.Ldc;
    dx(15) = (p.Rsh*(i_dc-a.i_b)-vC_dc)/(p.Cdc*(p.Rsh+p.RCdc));

    % Grid-side inductor dynamics: capacitor node -> stiff grid
    dx(16) = (vC_d-p.vg_d-p.Rg*ig_d)/p.Lg+p.w0*ig_q;
    dx(17) = (vC_q-p.vg_q-p.Rg*ig_q)/p.Lg-p.w0*ig_d;
end

function a = rectifierAlgebraic(t,x,p)
    v_alpha = x(1);
    v_beta  = x(2);
    theta   = x(4);
    i_alpha = x(5);
    i_beta  = x(6);
    e3      = x(7);
    e4      = x(8);
    e5      = x(9);
    iL_d    = x(10);
    iL_q    = x(11);
    vC_d    = x(12);
    vC_q    = x(13);
    i_dc    = x(14);
    vC_dc   = x(15);
    ig_d     = x(16);
    ig_q     = x(17);

    s = sin(theta);
    c = cos(theta);

    % Stiff-grid dq quantities
    a.vg_d = p.vg_d;
    a.vg_q = p.vg_q;

    phi = p.w0*t;
    sp = sin(phi);
    cp = cos(phi);  

    % Physical alpha-axis terminal/grid quantities
    a.v_ac_terminal = vC_d*sp+vC_q*cp;
    a.iL_alpha      = iL_d*sp+iL_q*cp;
    a.i_ac_terminal = ig_d*sp+ig_q*cp;
    a.v_grid_alpha  = p.vg_d*sp+p.vg_q*cp;
    a.v_grid_alpha_theta  = p.vg_d*s+p.vg_q*c;

    % DC output/load quantities. Because i_dc is defined positive from
    % the load node toward the DC link, normal rectification has i_dc < 0.
    a.i_load = -i_dc;
    a.vo = p.Rload*a.i_load;       % equivalently, -p.Rload*i_dc

    a.vhat_d = v_alpha*s-v_beta*c;
    a.vhat_q = v_alpha*c+v_beta*s;
    a.ihat_d = i_alpha*s-i_beta*c;
    a.ihat_q = i_alpha*c+i_beta*s;

    % Scalar fixed-point solve for the ESR/modulator/controller algebraic loop.
    vdc = max(vC_dc,p.vdcFloor);
    A = p.Rsh/(p.Rsh+p.RCdc);

    % % version 3
        iC_d = iL_d - ig_d;
        iC_q = iL_q - ig_q;
        % vconv_d_ref = vPI_d - p.w0*p.L*a.ihat_q + p.vd_ff- p.Kad*iC_d;
        % vconv_q_ref = vPI_q + p.w0*p.L*a.ihat_d + p.vq_ff - p.Kad*iC_q;

        iC_alpha = iC_d*sp + iC_q*cp;
        iC_beta  = -iC_d*cp + iC_q*sp;

        % rotate from plant frame (phi) to PLL frame (theta) before entering the controller equations
        iC_d = ...
            iC_alpha*s ...
            - iC_beta*c;

        iC_q = ...
            iC_alpha*c ...
            + iC_beta*s;

    for k = 1:p.algMaxIter
        ev = p.vdc_ref-vdc;
        id_ref = -p.k3*ev - p.k4*e3;

        vPI_d = p.k5*e4+p.k6*(id_ref-a.ihat_d);
        vPI_q = p.k7*e5+p.k8*(p.iq_ref-a.ihat_q);

        vconv_d_ref = vPI_d - p.w0*p.L*a.ihat_q + p.vd_ff- p.Kad*iC_d;
        vconv_q_ref = vPI_q + p.w0*p.L*a.ihat_d + p.vq_ff - p.Kad*iC_q;

        vconv_alpha_ref = vconv_d_ref*s+vconv_q_ref*c;

        vdc_for_mod = max(vdc,p.vdcFloor);
        mstar = vconv_alpha_ref/vdc_for_mod;
        m = min(max(mstar,-1),1);

        i_b = m*a.iL_alpha;
        vdc_new = A*(vC_dc+p.RCdc*(i_dc-i_b));

        if abs(vdc_new-vdc) <= p.algTol*max(1,abs(vdc))
            vdc = vdc_new;
            break;
        end
        vdc = (1-p.algRelax)*vdc+p.algRelax*vdc_new;
    end

    a.v_dc = vdc;
    a.ev_dc = p.vdc_ref-a.v_dc;
    a.id_ref = -p.k3*a.ev_dc-p.k4*e3;
    a.iq_ref = p.iq_ref;

    a.vPI_d = p.k5*e4+p.k6*(a.id_ref-a.ihat_d);
    a.vPI_q = p.k7*e5+p.k8*(p.iq_ref-a.ihat_q);

    a.vconv_d_ref = a.vPI_d - p.w0*p.L*a.ihat_q + p.vd_ff- p.Kad*iC_d;
    a.vconv_q_ref = a.vPI_q + p.w0*p.L*a.ihat_d + p.vq_ff- p.Kad*iC_q;

    a.vconv_alpha_ref = a.vconv_d_ref*s+a.vconv_q_ref*c;
    a.vconv_beta_ref  = -a.vconv_d_ref*c+a.vconv_q_ref*s;

    % rotate from controller frame (theta) to phi (plant frame) before entering the plant equations
    a.vconv_d_phi = ...
        a.vconv_alpha_ref*sp ...
        - a.vconv_beta_ref*cp;
    a.vconv_q_phi = ...
        a.vconv_alpha_ref*cp ...
        + a.vconv_beta_ref*sp;


    vdc_for_mod = max(a.v_dc,p.vdcFloor);
    a.mstar = a.vconv_alpha_ref/vdc_for_mod;
    a.m = min(max(a.mstar,-1),1);

    if abs(a.mstar) > 1e-12
        a.gamma = a.m/a.mstar;
    else
        a.gamma = 1;
    end

    a.vconv_d = a.gamma*a.vconv_d_ref;
    a.vconv_q = a.gamma*a.vconv_q_ref;
    a.vconv_alpha = a.m*a.v_dc;

    a.i_b = a.m*a.iL_alpha;
    a.iC_dc = i_dc-a.i_b-a.v_dc/p.Rsh;
    a.i_sh = a.v_dc/p.Rsh;

    vdc_check = A*(vC_dc+p.RCdc*(i_dc-a.i_b));
    a.algResidual = a.v_dc-vdc_check;
end

function avg = reconstructSignals(t,x,p)
    N = numel(t);
    avg.t = t;

    avg.v_alpha = x(:,1);
    avg.v_beta  = x(:,2);
    avg.e2      = x(:,3);
    avg.theta   = mod(x(:,4),2*pi);
    avg.i_alpha = x(:,5);
    avg.i_beta  = x(:,6);
    avg.e3      = x(:,7);
    avg.e4      = x(:,8);
    avg.e5      = x(:,9);
    avg.iL_d    = x(:,10);
    avg.iL_q    = x(:,11);
    avg.vC_d    = x(:,12);
    avg.vC_q    = x(:,13);
    avg.i_dc    = x(:,14);
    avg.i_load  = -x(:,14);
    avg.vC_dc   = x(:,15);
    avg.ig_d     = x(:,16);
    avg.ig_q     = x(:,17);

    names = {'vg_d','vg_q','vo','i_load','v_ac_terminal','i_ac_terminal','v_grid_alpha','v_grid_alpha_theta','iL_alpha', ...
             'vhat_d','vhat_q','ihat_d','ihat_q','ev_dc','id_ref','iq_ref', ...
             'vPI_d','vPI_q','vconv_d_ref','vconv_q_ref','vconv_alpha_ref', ...
             'vconv_beta_ref','mstar','m','gamma','vconv_d','vconv_q', ...
             'vconv_alpha','i_b','v_dc','iC_dc','i_sh','algResidual'};

    for k = 1:numel(names)
        avg.(names{k}) = zeros(N,1);
    end

    for n = 1:N
        a = rectifierAlgebraic(t(n),x(n,:).',p);
        for k = 1:numel(names)
            avg.(names{k})(n) = a.(names{k});
        end
    end
end

function ts = getLoggedTimeseries(logs,signalName)
    el = logs.getElement(signalName);
    if isempty(el)
        error('Signal "%s" was not found in logsout.',signalName);
    end
    ts = el.Values;
    if ~isa(ts,'timeseries')
        error('Logged signal "%s" is not a timeseries.',signalName);
    end
end
