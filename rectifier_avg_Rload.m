%% rectifier_avg_Rload.m
% Averaged nonlinear model of the single-phase rectifier stage.
%
% State order:
% x(1)  = v_alpha       % voltage SOGI state
% x(2)  = v_beta        % voltage SOGI state
% x(3)  = e2            
% x(4)  = theta
% x(5)  = i_alpha       % grid-side current SOGI state
% x(6)  = i_beta        % grid-side current SOGI state
% x(7)  = e3            % DC voltage controller PI integrator state
% x(8)  = e4            % current controller d-axis PI integrator state
% x(9)  = e5            % current controller q-axis PI integrator state
% x(10) = iL_d          % converter-side inductor d-axis current in plant frame (phi)
% x(11) = iL_q          % converter-side inductor q-axis current in plant frame (phi)
% x(12) = vC_d          % shunt capacitor d-axis voltage in plant frame (phi)
% x(13) = vC_q          % shunt capacitor q-axis voltage in plant frame (phi)
% x(14) = i_dc          % DC-link inductor current
% x(15) = vC_dc         % DC-link capacitor voltage
% x(16) = ig_d          % grid-side inductor d-axis current in plant frame (phi)
% x(17) = ig_q          % grid-side inductor q-axis current in plant frame (phi)

% Sign conventions:
%   iL   > 0 : converter -> AC terminal/grid
%   ig   > 0 : AC capacitor node -> grid
%   i_dc > 0 : DC load/output node -> DC-link node/converter
%   i_b  > 0 : DC-link node -> H-bridge
%
% IMPORTANT:
%   The grid current is modeled dynamically through Lg + Rg.
%   A shunt load resistor Rload closes the DC output after Ldc + RLdc.

clear;
clc;
close all;

%% 1) USER PARAMETERS
p.tStop = 1;
p.dtOut = 1e-5;
p.maxStep = 1e-4;

p.Kad = 0;  % virtual resistance for the LCL capacitor current feedback. Set to 0 to disable.

p.f0 = 60;
p.w0 = 2*pi*p.f0;
p.kSOGI = 1;            % SOGI damping factor

p.k1 = 100;             % PLL proportional
p.k2 = 1;               % PLL integral
p.k3 = 0.1;             % DC-voltage proportional
p.k4 = 5;               % DC-voltage integral
p.k5 = 0.5;             % d-current integral
p.k6 = 20;              % d-current proportional
p.k7 = 0.5;             % q-current integral
p.k8 = 20;              % q-current proportional

p.L = 800e-6;           % AC converter-side inductance [H]
p.RL = 0;               % AC inductor resistance [ohm]
p.C = 100e-6;           % AC shunt capacitance [F]

p.Lg = 800e-6;          % grid-side inductance [H]
p.Rg = 0.2;             % grid-side inductor ESR [ohm]
p.vg_d = 170;           % stiff-grid d-axis voltage [V]
p.vg_q = 0;             % stiff-grid q-axis voltage [V]

p.Ldc  = 1e-6;          % DC series inductance [H]
p.RLdc = 0.1;           % DC inductor ESR [ohm]
p.Cdc  = 400e-6;        % DC-link capacitance [F]
p.RCdc = 0.01;          % DC capacitor ESR [ohm]
p.Rsh  = 1e6;           % DC shunt resistance [ohm]
p.Rload = 88;           % DC output/load resistance [ohm]

p.vdc_ref = 300;        % DC voltage reference [V]
p.iq_ref  = 0;          % q-axis current reference [A]
% p.igq_ref = 0;          % grid q-axis current reference [A]
p.vd_ff = 170;          % constant feedforward d voltage [V]
p.vq_ff = 0;            % constant feedforward q voltage [V]

% the code iteratively solves for a consistent vdc at every ODE evaluation
p.vdcFloor = 1;        % prevent devision by 0 when calculating mod-index
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
% save('rectifier_averaged_results.mat','avg','p');

f1 = plotGridStates(t,p,avg);
idx = avg.t > 0.8;
checkGridStateAvg(idx, avg);
checkControllerSetpoints(idx, avg);

%% 7) BASIC PLOTS
fprintf('Plotting the terminal signals...');
figure('Name','Averaged Rectifier - Terminal Signals');
tiledlayout(3,2);

nexttile;
plot(avg.t,avg.vg);
grid on;
xlabel('Time [s]'); ylabel('v_{AC} [V]');
title('AC terminal voltage');

nexttile;
plot(avg.t,avg.ig);
grid on;
xlabel('Time [s]'); ylabel('i_{AC} [A]');
title('AC terminal current');

nexttile;
plot(avg.t,avg.vC);
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
plot(avg.t,avg.id_ref,avg.t,avg.ighat_d);
grid on;
xlabel('Time [s]'); ylabel('Current [A]');
legend('i_d^*','\hat{i}_d');
title('d-axis current loop');

nexttile;
plot(avg.t,avg.iq_ref,avg.t,avg.ighat_q);
grid on;
xlabel('Time [s]'); ylabel('Current [A]');
legend('i_q^*','\hat{i}_q');
title('q-axis current loop');

% nexttile;
% plot(avg.t,avg.m);
% grid on;
% xlabel('Time [s]'); ylabel('m');
% title('Modulation index');

% nexttile;
% plot(avg.t,avg.vconv_alpha_ref,avg.t,avg.vconv_alpha);
% grid on;
% xlabel('Time [s]'); ylabel('Voltage [V]');
% legend('v_{conv,\alpha}^*','v_{conv,\alpha}');
% title('Converter voltage command vs average');

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

% % figure;
% % plot(avg.t,avg.mstar);
% % hold on;
% % plot(avg.t,avg.m);
% % yline(1,'--');
% % yline(-1,'--');
% % grid on;
% % legend('m^*','m');

% % PLL estimated angular frequency
% avg.omega_hat = p.w0 ...
%               + p.k1*avg.vhat_q ...
%               + p.k2*avg.e2;

% figure('Name','PLL Frequency');
% plot(avg.t,avg.omega_hat,'LineWidth',1.2);
% hold on;
% yline(p.w0,'--','\omega_0');
% grid on;

% xlabel('Time [s]');
% ylabel('Angular frequency [rad/s]');
% legend('\hat{\omega}','\omega_0');
% title('PLL Estimated Frequency vs Nominal Frequency');
% xlim([0 p.tStop]);


figure('Name','Grid Voltage Reference Check');
plot(avg.t,avg.vg,'LineWidth',1.2);
hold on;
plot(avg.t,avg.vg_theta,'--','LineWidth',1.2);
grid on;

xlabel('Time [s]');
ylabel('Voltage [V]');
legend('\phi = \omega_0t','\theta');
title('Grid Voltage Reconstruction');

%% DIAGNOSTIC PLOTS
function fig = plotGridStates(t, p, avg)
    % ============================================================
    % Current dq quantities in plant frame (phi) vs controller
    % frame (theta)
    % ============================================================
    phi = p.w0*t;
    sp = sin(phi);
    cp = cos(phi);
    st = sin(avg.theta);
    ct = cos(avg.theta);
    % ------------------------------------------------------------
    % Converter-side inductor current
    %
    % Existing states iL_d and iL_q are dq quantities in the
    % plant frame rotating at phi.
    %
    % First reconstruct alpha-beta quantities from phi-frame dq.
    % ------------------------------------------------------------
    avg.iL_alpha = avg.iL_d.*sp + avg.iL_q.*cp;
    avg.iL_beta  = -avg.iL_d.*cp + avg.iL_q.*sp;

    % Rotate the same alpha-beta quantities using theta
    avg.iL_d_theta = avg.iL_alpha.*st - avg.iL_beta.*ct;
    avg.iL_q_theta = avg.iL_alpha.*ct + avg.iL_beta.*st;

    % ------------------------------------------------------------
    % Grid-side inductor current
    % ------------------------------------------------------------
    avg.ig_alpha = avg.ig_d.*sp + avg.ig_q.*cp;
    avg.ig_beta  = -avg.ig_d.*cp + avg.ig_q.*sp;

    % Rotate using theta
    avg.ig_d_theta = avg.ig_alpha.*st - avg.ig_beta.*ct;
    avg.ig_q_theta = avg.ig_alpha.*ct + avg.ig_beta.*st;

    fprintf('Plotting AC-side states.\n');
    fig = figure('Name','AC-side states');
    tiledlayout(2,1);
    % ================================================================
    % 1. Grid voltage and capacitor voltage
    % ================================================================
    ax1 = nexttile;
    plot(avg.t,avg.vg,'LineWidth',1.2);
    hold on;
    plot(avg.t,avg.vC,'LineWidth',1.2);
    hold on;
    % plot(avg.t,avg.v_alpha,'LineWidth',1.2);
    % hold on;
    % plot(avg.t,avg.v_beta,'LineWidth',1.2);

    grid on;
    xlabel('Time [s]');
    ylabel('Voltage [V]');
    legend('v_g','v_C','v_\alpha','v_\beta','Location','best');
    title('Grid and Capacitor Voltages');
    xlim([0 p.tStop]);

    % ================================================================
    % 2. Grid-side and converter-side currents
    % ================================================================
    ax2 = nexttile;
    plot(avg.t,avg.ig,'LineWidth',1.2);
    hold on;
    plot(avg.t,avg.iL_alpha,'LineWidth',1.2);
    % hold on;
    % plot(avg.t,avg.i_alpha,'LineWidth',1.2);
    % hold on;
    % plot(avg.t,avg.i_beta,'LineWidth',1.2);

    grid on;
    xlabel('Time [s]');
    ylabel('Current [A]');
    legend('i_g','i_L','i_\alpha','i_\beta','Location','best');
    title('Grid-side and Converter-side Currents');
    xlim([0 p.tStop]);

    linkaxes([ax1 ax2],"x")


end

function checkGridStateAvg(idx, avg)
    fprintf('\nChecking average grid-side and converter-side currents...\n');
    
    igd = mean(avg.ig_d(idx));
    igq = mean(avg.ig_q(idx));
    delta_ig = rad2deg(atan2(igq,igd));

    iLd = mean(avg.iL_d(idx));
    iLq = mean(avg.iL_q(idx));
    delta_iL = rad2deg(atan2(iLq,iLd));

    fprintf('Average ig: %.3f + j%.3f A (angle = %.2f deg)\n',igd,igq,delta_ig);
    fprintf('Average iL: %.3f + j%.3f A (angle = %.2f deg)\n',iLd,iLq,delta_iL);
end

function checkControllerSetpoints(idx, avg)
    fprintf('\nChecking controller setpoints...\n');
    id_ref = mean(avg.id_ref(idx));
    id_hat = mean(avg.ighat_d(idx));

    iq_ref = mean(avg.iq_ref(idx));
    iq_hat = mean(avg.ighat_q(idx));

    vconv_d_ref = mean(avg.vconvhat_d_ref(idx));
    vconv_d = mean(avg.vconv_d(idx));   

    vconv_q_ref = mean(avg.vconvhat_q_ref(idx));
    vconv_q = mean(avg.vconv_q(idx));

    fprintf('Average id_ref: %.3f A, Average i_hat_d: %.3f A\n',id_ref,id_hat);
    fprintf('Average iq_ref: %.3f A, Average i_hat_q: %.3f A\n',iq_ref,iq_hat);
    fprintf('Average vconv_d_ref: %.3f V, Average vconv_d: %.3f V\n',vconv_d_ref,vconv_d);
    fprintf('Average vconv_q_ref: %.3f V, Average vconv_q: %.3f V\n',vconv_q_ref,vconv_q);
end 


%% 8) OPTIONAL SIMULINK COMPARISON
% Set cmp.enabled = true and replace the model/signal names.
% The Simulink signals must be logged to logsout.
cmp.enabled = false;
cmp.modelName = 'your_switched_model';
cmp.vacLogName = 'vC';
cmp.iacLogName = 'ig';
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
    plot(avg.t,avg.vC,'LineWidth',1.2);
    grid on; xlabel('Time [s]'); ylabel('Voltage [V]');
    legend('Switched','Averaged'); title('AC terminal voltage');

    nexttile;
    plot(sw.iac.Time,sw.iac.Data); hold on;
    plot(avg.t,avg.ig,'LineWidth',1.2);
    grid on; xlabel('Time [s]'); ylabel('Current [A]');
    legend('Switched','Averaged'); title('AC terminal current');

    nexttile;
    plot(sw.vdc.Time , sw.vdc.Data); hold on;
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

    dx(1) = p.kSOGI*p.w0*(a.vC-v_alpha) - p.w0*v_beta;
    dx(2) = p.w0*v_alpha;
    dx(3) = a.vhat_q;
    dx(4) = p.w0+p.k1*a.vhat_q+p.k2*x(3);

    dx(5) = p.kSOGI*p.w0*(a.ig-i_alpha)-p.w0*i_beta;
    dx(6) = p.w0*i_alpha;

    dx(7) = p.vdc_ref-a.v_dc;   % = 0;
    % dx(7) = 0;
    dx(8) = a.id_ref-a.ighat_d;
    dx(9) = p.iq_ref-a.ighat_q;

    dx(10) = (a.vconv_d - vC_d - p.RL*iL_d)/p.L + p.w0*iL_q;
    dx(11) = (a.vconv_q - vC_q - p.RL*iL_q)/p.L - p.w0*iL_d;

    dx(12) = (iL_d-ig_d)/p.C+p.w0*vC_q;
    dx(13) = (iL_q-ig_q)/p.C-p.w0*vC_d;

    dx(14) = (-a.v_dc-(p.RLdc+p.Rload)*i_dc)/p.Ldc;
    dx(15) = (p.Rsh*(i_dc-a.i_b)-vC_dc)/(p.Cdc*(p.Rsh+p.RCdc));

    % Grid-side inductor dynamics: capacitor node -> stiff grid
    dx(16) = (vC_d - p.vg_d - p.Rg*ig_d)/p.Lg + p.w0*ig_q;
    dx(17) = (vC_q - p.vg_q - p.Rg*ig_q)/p.Lg - p.w0*ig_d;
end

function a = rectifierAlgebraic(t,x,p)
%RECTIFIERALGEBRAIC Calculates all quantities that are not independent dynamic states but are needed by the ODEs.
%   x(t) --> frame transformations --> controller eqns --> vconv_ref --> m --> i_b --> vdc
%   for-loop solves an instantaneous algebraic loop for vdc, which is used in the modulator and controller equations.
%   vdc --> ic_ref --> vconv_ref --> m --> i_b --> vdc

    % unpack states to be used in the algebraic equations
    v_alpha = x(1);
    v_beta  = x(2);
            % e2 = x(3) is not used in the algebraic equations
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

    % Stiff-grid dq quantities
    a.vg_d = p.vg_d;
    a.vg_q = p.vg_q;

    % plant frame angle (phi) is independent of PLL angle output (theta)
    phi = p.w0*t;
    sp = sin(phi);
    cp = cos(phi);  
    % Physical grid keeps rotating perfectly at the nominal frequency
    a.ig = ig_d*sp + ig_q*cp;
    a.vg = p.vg_d*sp + p.vg_q*cp;
    a.vC = vC_d*sp + vC_q*cp;
    a.iL_alpha = iL_d*sp + iL_q*cp;
    iL_beta  = -iL_d*cp + iL_q*sp;

    % DC output/load quantities. Because i_dc is defined positive from
    % the load node toward the DC link, normal rectification has i_dc < 0.
    a.i_load = -i_dc;
    a.vo = p.Rload*a.i_load;       % equivalently, -p.Rload*i_dc

    % controller frame angle (theta) = the PLL-estimated angle
    s = sin(theta);
    c = cos(theta);
            % !!!! for DIAGNOSTIC ONLY: does theta lock with phi?
            a.vg_theta  = p.vg_d*s + p.vg_q*c;

    % controller uses hatted-dq quantities = plant-frame alpha-beta quantities rotated by theta
    a.vhat_d = v_alpha*s - v_beta*c;
    a.vhat_q = v_alpha*c + v_beta*s;
    a.ighat_d = i_alpha*s - i_beta*c;     % using grid current
    a.ighat_q = i_alpha*c + i_beta*s;     % using grid current

    % Scalar fixed-point solve for the ESR/modulator/controller algebraic loop.
    vdc = max(vC_dc,p.vdcFloor);
    A = p.Rsh/(p.Rsh+p.RCdc);

    % iL in controller theta frame
    a.iLhat_d = a.iL_alpha*s - iL_beta*c;
    a.iLhat_q = a.iL_alpha*c + iL_beta*s;

    % capacitor current in plant frame (phi)
    iC_d = iL_d - ig_d;
    iC_q = iL_q - ig_q;

    % transform to stationary frame
    iC_alpha = iC_d*sp + iC_q*cp;
    iC_beta  = -iC_d*cp + iC_q*sp;

    % rotate by theta before entering the controller algebraic equations
    iChat_d = iC_alpha*s - iC_beta*c;
    iChat_q = iC_alpha*c + iC_beta*s;

    % Iteratively solve for a consistent vdc at every ODE evaluation
    % the controller uses vdc to determine current controller voltage reference output
    % which is used to determine the modulation
    % which dictates the bridge DC current
    % which affects the DC capacitor ESR voltage drop 
    % which affects the DC capacitor voltage
    for k = 1:p.algMaxIter
        ev = p.vdc_ref-vdc;
        id_ref = -p.k3*ev - p.k4*e3;   % = -12;

        % Current controller outputs dq voltage references
        vPI_d = p.k5*e4 + p.k6*(id_ref - a.ighat_d);
        vPI_q = p.k7*e5 + p.k8*(p.iq_ref - a.ighat_q);

        % add decoupling terms and feedforward terms to get the final dq voltage references
        % last term: virtual resistance for the LCL capacitor current feedback. Set to 0 to disable.
        vconvhat_d_ref = vPI_d - p.w0*p.L*a.ighat_q + p.vd_ff - p.Kad*iChat_d;
        vconvhat_q_ref = vPI_q + p.w0*p.L*a.ighat_d + p.vq_ff - p.Kad*iChat_q;

        % undo rotation by theta to get the stationary-frame voltages
        vconv_alpha_ref = vconvhat_d_ref*s + vconvhat_q_ref*c;
                % vconv_beta_ref  = -vconvhat_d_ref*c + vconvhat_q_ref*s;     % not used 

        % determine modulation index
        vdc_for_mod = max(vdc,p.vdcFloor);
        mstar = vconv_alpha_ref/vdc_for_mod;    % definition of modulation index for averaged bridge model
        m = min(max(mstar,-1),1);   % enforce saturation limits

        % calculate bridge current; this preserves p(t) for ideal averaged bridge model
        i_b = m*a.ig;           
        % calculate new vdc based on the bridge current and the DC capacitor ESR drop
        vdc_new = A*(vC_dc+p.RCdc*(i_dc-i_b));  

        if abs(vdc_new-vdc) <= p.algTol*max(1,abs(vdc)) % converged to a consistent vdc value
            vdc = vdc_new;
            break;
        end
        % relax the iteration by going just 0.7 of the way to the new value to prevent oscillatory convergence
        vdc = (1-p.algRelax)*vdc+p.algRelax*vdc_new;
    end

    % recalculate everthing with the final vdc value to return to the ODE function
    a.v_dc = vdc;
    a.ev_dc = p.vdc_ref-a.v_dc;     % error signal for the DC voltage controller
    a.id_ref = -p.k3*a.ev_dc - p.k4*e3;    
    % a.id_ref = -12; 
    a.iq_ref = p.iq_ref;

    a.vPI_d = p.k5*e4 + p.k6*(a.id_ref - a.ighat_d);
    a.vPI_q = p.k7*e5 + p.k8*(p.iq_ref - a.ighat_q);

    a.vconvhat_d_ref = a.vPI_d - p.w0*p.L*a.ighat_q + p.vd_ff - p.Kad*iChat_d;
    a.vconvhat_q_ref = a.vPI_q + p.w0*p.L*a.ighat_d + p.vq_ff - p.Kad*iChat_q;

    % undo rotation by theta to get the stationary-frame voltages
    a.vconv_alpha_ref = a.vconvhat_d_ref*s+a.vconvhat_q_ref*c;
    a.vconv_beta_ref  = -a.vconvhat_d_ref*c+a.vconvhat_q_ref*s;

    % modulation index
    vdc_for_mod = max(a.v_dc,p.vdcFloor);
    a.mstar = a.vconv_alpha_ref/vdc_for_mod;
    a.m = min(max(a.mstar,-1),1);

    % m=m* == gamma=1 when the modulation is NOT saturated. Otherwise, gamma < 1.
    if abs(a.mstar) > 1e-12
        a.gamma = a.m/a.mstar;
    else
        a.gamma = 1;
    end

    % a.gamma = a.m/a.mstar;  % this is always well-defined because vdc_for_mod >= p.vdcFloor > 0

    % back to physical quantities 
    % actual voltage applied to the AC side of the rectifier, by definition of the averaged bridge model
    a.vconv_alpha = a.m*a.v_dc;
    a.vconv_beta= a.gamma*a.vconv_beta_ref;

    % rotate by phi (plant frame) before entering the plant equations
    a.vconv_d = a.vconv_alpha*sp - a.vconv_beta*cp;
    a.vconv_q = a.vconv_alpha*cp + a.vconv_beta*sp;


    % DC-side KCL
    a.i_b = a.m*a.ig;  
    a.iC_dc = i_dc-a.i_b-a.v_dc/p.Rsh;
    a.i_sh = a.v_dc/p.Rsh;

    % check whether fixed-point solution converged to a consistent vdc value
    vdc_check = A*(vC_dc+p.RCdc*(i_dc-a.i_b));
    a.algResidual = a.v_dc-vdc_check;     % ideally this should be 0; may not be if fixed-point iteration did not converge
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
    avg.ig_d    = x(:,16);
    avg.ig_q    = x(:,17);

    names = {'vg_d','vg_q','vo','i_load','vC','ig','vg','vg_theta','iL_alpha', ...
             'vhat_d','vhat_q','ighat_d','ighat_q','iLhat_d','iLhat_q', ...
             'ev_dc','id_ref','iq_ref', ...
             'vPI_d','vPI_q','vconvhat_d_ref','vconvhat_q_ref','vconv_alpha_ref', ...
             'vconv_beta_ref','mstar','m','gamma','vconv_d','vconv_q', ...
             'vconv_alpha','vconv_beta','i_b','v_dc','iC_dc','i_sh', ...
             'algResidual'};

    for n = 1:N
        a = rectifierAlgebraic(t(n),x(n,:).',p);

        for k = 1:numel(names)
            name = names{k};

            % Skip this name if it does not exist in a
            if ~isfield(a,name)
                continue
            end

            % Initialize the field the first time it is encountered
            if ~isfield(avg,name)
                avg.(name) = zeros(N,1);
            end

            avg.(name)(n) = a.(name);
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