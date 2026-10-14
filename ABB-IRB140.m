% Longitudes (ABB IRB-140)
L0 = 0.070;   
L1 = 0.3521;   
L2 = 0.360;  
L3 = 0.380;    
L4 = 0.065;   
 
% Parametros DH;
L(1) = Link([0,  L1, L0,  pi/2,  0,  0]);
L(2) = Link([0,  0,  L2,  0,     0,  0]);
L(3) = Link([0,  0,  0,   pi/2,  0,  pi/2]);  % offset +pi/2
L(4) = Link([0,  L3, 0,  -pi/2,  0,  0]);
L(5) = Link([0,  0,  0,   pi/2,  0,  0]);
L(6) = Link([0,  L4, 0,   0,     0,  0]);
 
% Limites articulares
L(1).qlim = deg2rad([-180,  180]);
L(2).qlim = deg2rad([-90,   110]);
L(3).qlim = deg2rad([-230,   50]);
L(4).qlim = deg2rad([-200,  200]);
L(5).qlim = deg2rad([-115,  115]);
L(6).qlim = deg2rad([-400,  400]);
 
robot = SerialLink(L, 'name', 'ABB IRB-140');
robot.display();

%---------------------------------------------------------------------------------
% Definir dh_matrix como funcion anonima (equivalente a la Toolbox)
dh_matrix = @(theta, d, a, alpha) ...
    [cos(theta), -sin(theta)*cos(alpha),  sin(theta)*sin(alpha), a*cos(theta);
     sin(theta),  cos(theta)*cos(alpha), -cos(theta)*sin(alpha), a*sin(theta);
     0,           sin(alpha),             cos(alpha),            d;
     0,           0,                      0,                     1];
 
% Variables simbolicas
syms q1 q2 q3 q4 q5 q6
q = [q1 q2 q3 q4 q5 q6];
 
% Calculo de las matrices i-1_T_i mediante bucle
for i = 1:robot.n
    Ti(:,:,i) = dh_matrix(q(i) + robot.offset(i), ...
                           robot.d(i), ...
                           robot.a(i), ...
                           robot.alpha(i));
end
 
% Mostrar cada matriz con 4 cifras significativas
% Los valores del tipo 6.123e-17 son cero numerico 
vpa(Ti, 4)
 
%---------------------------------------------------------------------------------
% Cinematica directa completa
% Forma manual (producto de matrices):
T06 = Ti(:,:,1)*Ti(:,:,2)*Ti(:,:,3)*Ti(:,:,4)*Ti(:,:,5)*Ti(:,:,6);

% Forma con Toolbox:
T06 = vpa(robot.fkine(q), 4);

% 10 puntos articulares [deg] -> radianes
Q_deg = [
    0,    0,    0,    0,    0,    0;
   30,  -20,   30,    0,   45,    0;
   60,  -45,   60,   90,   30,   45;
   90,   10,  -30,  -45,   60,   30;
  120,  -60,   20,   60,  -30,   90;
  -30,   20,  -60,   30,   15,  -45;
  150,  -10,   40,  -90,   45,   60;
  -60,   40,  -20,   45,  -60,  -30;
  -90,  -30,   50,   60,   30,  -90;
  180,   50,  -10,  -30,  -45,   45];

Q = deg2rad(Q_deg);   

%---------------------------------------------------------------------------------
% Cinematica inversa
tam_movimientos = size(Q, 1);

Tfkine = zeros(4, 4, tam_movimientos);
Tikunc = zeros(4, 4, tam_movimientos);
qinv   = zeros(tam_movimientos, robot.n);

for i = 1:tam_movimientos
    %Toma la configuracion articular del punto i y calcula CD: T deseada.

    Tfkine(:,:,i) = robot.fkine(Q(i,:))

    %Con la T deseada, se calcula la cinematica inversa
    qinv(i,:)

    %Con los angulos que devolvio ikunc, se calcula de nuevo la cinematica directa.
    %Si la cinematica inversa es correcta, esta T tiene que ser igual a Tfkine(:,:,i)
    Tikunc(:,:,i) = robot.fkine(qinv(i,:));
end

%---------------------------------------------------------------------------------
% Jacobiano y deteccion de singularidades
J_matrices = zeros(6, robot.n, tam_movimientos);
J_det      = zeros(tam_movimientos, 1);

for i = 1:tam_movimientos
    % Matriz jacobiana en cada configuracion
    J_matrices(:,:,i) = robot.jacob0(Q(i,:));

    % Determinante (indica singularidad si = 0)
    J_det(i) = det(J_matrices(:,:,i));
end

% Mostrar resultados
fprintf('%-6s | %-20s | %-12s\n', 'Punto', 'det(J)', 'Estado');
fprintf('%s\n', repmat('-', 1, 44));

for i = 1:tam_movimientos
    if abs(J_det(i)) < 1e-4
        estado = 'SINGULAR';
    else
        estado = 'No singular';
    end
    fprintf('  %2d   | %20.6f  | %s\n', i, J_det(i), estado);
end

%%
% ---------------------------------------------------------------------------------
% PARÁMETROS DINÁMICOS: 
%Fuente: "Kinematics, Dynamics and Evaluation of Energy Consumption
%          for ABB IRB-140" - valores obtenidos con SOLID EDGE ST4 (Baquero
%          & Ramírez, Universidad Nacional Colombia)

% Masas [kg]
mass = [35, 25, 18, 6.5, 1.5, 0.5]; % la herramienta se ha estimado su peso
% Centros de gravedad [m] respecto al SR DH de cada eslabón
pos_cdg = [ 0.06491  -0.13089   0.02913;   % link 1
         -0.00973  -0.16171   0.09243;   % link 2
         -0.00089  -0.00867  -0.13649;   % link 3
         -0.04077   0         0;         % link 4
         -0.0000602 -0.0004162 -0.04867; % link 5
          0         0         0.032];    % link 6

% Matrices de inercia [kg·m²]
I1 = [ 0.765  -0.3    0.067;
      -0.3     0.31  -0.14;
       0.067  -0.14   0.78 ];

I2 = [ 1.38    0.045  -0.02;
       0.045   0.364  -0.31;
      -0.02   -0.31    1.09 ];

I3 = [ 0.349   0       0.002;
       0       0.338   0.021;
       0.002   0.021   0.014];

I4 = [ 0.005   0       0;
       0       0.011   0;
       0       0       0.016];

I5 = [ 0.005   0       0;
       0       0.005   4e-5;
       0       4e-5    0.0006];

I6 =  [ 1.5e-4   0   0;
        0   1.5e-4   0
        0   0        0.4e-4];

I_Tot = {I1, I2, I3, I4, I5, I6}

% Introducir en el objeto robot
for i = 1:6
    robot.links(i).m = mass(i);
    robot.links(i).r = pos_cdg(i,:);
    robot.links(i).I = I_Tot{i};
    robot.links(i).B = 0.05;
    robot.links(i).Tc = [0.01; -0.01];
    robot.links(i).Jm = 0;
    robot.links(i).G  = 1;
end

%ahora calculo Tau, G, M, C
qd  = zeros(1,6);   % velocidad
qdd = zeros(1,6);   %  aceleración

% Par total para los 10 puntos (qd=0, qdd=0)
tau_total = zeros(tam_movimientos, 6);
G = zeros(tam_movimientos, 6);
M = zeros(tam_movimientos, 6);
C = zeros(tam_movimientos, 6);

for i = 1:tam_movimientos
    % Par total -> τ = M(q)·q̈ + C(q,q̇)·q̇ + G(q)
    tau_total(i,:) = robot.rne(Q(i,:),qd, qdd);
    
    % Gravedad ->  G(q)
    G(i,:) = robot.gravload(Q(i,:));
    
    % Inercia (con aceleración ejemplo = 1 rad/s²) -> M(q)·q̈
    M(i,:) = robot.itorque(Q(i,:), qdd);
    
    % Coriolis (con velocidad ejemplo = 0.5 rad/s) -> τ_coriolis = C(q, qd) · qd
    Cmat = robot.coriolis(Q(i,:), qd);
    C(i,:) = Cmat * qd';  
    
    % Verificación
    tau = G(i,:) + M(i,:) + C(i,:);
end

disp(tau_total)


% -------------------------------------------------------------------------
%  3. DEFINICIÓN DE LOS PUNTOS DE LA TAREA
% -------------------------------------------------------------------------
%  Escenario: el robot parte de home, se aproxima a una pieza y la recoge,
%  la transporta a una posición intermedia, la deposita y vuelve a home.
%
%  Posiciones articulares [deg]:
%   q_home   → posición de reposo
%   q_aprox  → aproximación sobre la pieza
%   q_agarre → agarre de la pieza (robot toca la pieza)
%   q_trans  → posición intermedia de transporte
%   q_deposit→ posición de depósito
%   q_retiro → retiro tras depositar
%
%  Las trayectorias cartesianas conectan poses donde la herramienta debe
%  moverse en línea recta (aproximación, transporte y retiro).

q_home    = deg2rad([  0,    0,    0,   0,   0,   0]);
q_aprox   = deg2rad([ 45,  -30,   60,   0,  30,   0]);
q_agarre  = deg2rad([ 45,  -15,   45,   0,  60,   0]);
q_trans   = deg2rad([ 90,  -45,   45,   0,  45,   0]);
q_deposit = deg2rad([120,  -30,   30,  30,  60,   0]);
q_retiro  = deg2rad([120,  -60,   45,   0,  30,   0]);

% Número de muestras temporales por trayectoria
N = 100;
t = linspace(0, 1, N);

% -------------------------------------------------------------------------
%  4. TRAYECTORIAS ARTICULARES
% -------------------------------------------------------------------------
fprintf('Calculando trayectorias articulares...\n');

% Trayectoria A1: home → aproximación
[qA1, qdA1, qddA1] = jtraj(q_home, q_aprox, N);

% Trayectoria A2: agarre → transporte(sin carga)
[qA2, qdA2, qddA2] = jtraj(q_agarre, q_trans, N);

% Trayectoria A3: depósito → home
[qA3, qdA3, qddA3] = jtraj(q_deposit, q_home, N);

% Trayectoria articular completa para animación
q_art_total  = [qA1; qA2; qA3];
qd_art_total = [qdA1; qdA2; qddA3];

% -------------------------------------------------------------------------
%  5. TRAYECTORIAS CARTESIANAS
% -------------------------------------------------------------------------
fprintf('Calculando trayectorias cartesianas...\n');

% Poses cartesianas obtenidas de la cinemática directa
T_aprox   = robot.fkine(q_aprox);
T_agarre  = robot.fkine(q_agarre);
T_trans   = robot.fkine(q_trans);
T_deposit = robot.fkine(q_deposit);
T_retiro  = robot.fkine(q_retiro);
T_home    = robot.fkine(q_home);

% Trayectoria C1: aproximación → agarre  (descenso hacia la pieza)
Ts_C1 = ctraj(T_aprox,  T_agarre,  N);

% Trayectoria C2: transporte → depósito  (movimiento con carga)
Ts_C2 = ctraj(T_trans,  T_deposit, N);

% Trayectoria C3: depósito → retiro  (ascenso tras depositar)
Ts_C3 = ctraj(T_deposit, T_retiro,  N);

% Opciones para ikine (cinemática inversa numérica sin Optimization Toolbox)
ik_opt = {'ilimit', 1000, 'tol', 1e-6, 'quiet'};

% Cinemática inversa para cada trayectoria cartesiana
fprintf('  Calculando IK trayectoria C1 (aprox -> agarre)...\n');
q_C1 = zeros(N, 6);
q0   = q_aprox;
for k = 1:N
    q_sol = robot.ikine(Ts_C1(k), q0, ik_opt{:});
    if ~isempty(q_sol)
        q_C1(k,:) = q_sol;
        q0 = q_sol;
    else
        q_C1(k,:) = q0;  % mantener último valor si no converge
    end
end

fprintf('  Calculando IK trayectoria C2 (trans -> deposit) CON CARGA...\n');
% --- Añadir carga antes de la trayectoria de transporte ---
robot.payload(2.0, [0, 0, 0.08]);   % pieza de 2 kg, a 8 cm del flange
q_C2 = zeros(N, 6);
q0   = q_trans;
for k = 1:N
    q_sol = robot.ikine(Ts_C2(k), q0, ik_opt{:});
    if ~isempty(q_sol)
        q_C2(k,:) = q_sol;
        q0 = q_sol;
    else
        q_C2(k,:) = q0;
    end
end

fprintf('  Calculando IK trayectoria C3 (deposit -> retiro)...\n');
% --- Quitar carga al depositar ---
robot.payload(0, [0, 0, 0.08]);
q_C3 = zeros(N, 6);
q0   = q_deposit;
for k = 1:N
    q_sol = robot.ikine(Ts_C3(k), q0, ik_opt{:});
    if ~isempty(q_sol)
        q_C3(k,:) = q_sol;
        q0 = q_sol;
    else
        q_C3(k,:) = q0;
    end
end

% -------------------------------------------------------------------------
%  6. CÁLCULO DINÁMICO EN TODAS LAS TRAYECTORIAS
% -------------------------------------------------------------------------
fprintf('Calculando dinámica...\n');

function [tau_tot, tau_G, tau_I, tau_C, manip] = calc_dynamics(rob, Q, Qd, Qdd, has_load, load_mass, load_pos)
    N = size(Q,1);
    tau_tot = zeros(N,6);
    tau_G   = zeros(N,6);
    tau_I   = zeros(N,6);
    tau_C   = zeros(N,6);
    manip   = zeros(N,1);

    if has_load
        rob.payload(load_mass, load_pos);
    else
        rob.payload(0, load_pos);
    end

    for k = 1:N
        q   = Q(k,:);
        qd  = Qd(k,:);
        qdd = Qdd(k,:);

        tau_tot(k,:) = rob.rne(q, qd, qdd);
        tau_G(k,:)   = rob.gravload(q);
        tau_I(k,:)   = rob.itorque(q, qdd);
        C_mat        = rob.coriolis(q, qd);
        tau_C(k,:)   = (C_mat * qd')';

        J = rob.jacob0(q);
        manip(k) = abs(det(J));
    end
end

qd_zero  = zeros(N,6);
qdd_zero = zeros(N,6);

% Velocidades y aceleraciones para trayectorias cartesianas
% (calculamos diferencias finitas)
dt = 1/(N-1);
qdC1  = [diff(q_C1,1,1)/dt;  zeros(1,6)];
qddC1 = [diff(qdC1,1,1)/dt;  zeros(1,6)];
qdC2  = [diff(q_C2,1,1)/dt;  zeros(1,6)];
qddC2 = [diff(qdC2,1,1)/dt;  zeros(1,6)];
qdC3  = [diff(q_C3,1,1)/dt;  zeros(1,6)];
qddC3 = [diff(qdC3,1,1)/dt;  zeros(1,6)];

load_pos = [0, 0, 0.08];

% Segmento A1: home → aprox  (sin carga)
[tA1, gA1, iA1, cA1, mA1] = calc_dynamics(robot, qA1, qdA1, qddA1, false, 0, load_pos);

% Segmento C1: aprox → agarre  (sin carga)
[tC1, gC1, iC1, cC1, mC1] = calc_dynamics(robot, q_C1, qdC1, qddC1, false, 0, load_pos);

% Segmento A2: agarre → trans  (CON carga)
[tA2, gA2, iA2, cA2, mA2] = calc_dynamics(robot, qA2, qdA2, qddA2, true, 2.0, load_pos);

% Segmento C2: trans → deposit  (CON carga)
[tC2, gC2, iC2, cC2, mC2] = calc_dynamics(robot, q_C2, qdC2, qddC2, true, 2.0, load_pos);

% Segmento C3: deposit → retiro  (sin carga)
[tC3, gC3, iC3, cC3, mC3] = calc_dynamics(robot, q_C3, qdC3, qddC3, false, 0, load_pos);

% Segmento A3: deposit → home  (sin carga)
[tA3, gA3, iA3, cA3, mA3] = calc_dynamics(robot, qA3, qdA3, qddA3, false, 0, load_pos);

Q_total   = [qA1; q_C1; qA2; q_C2; q_C3; qA3];
Qd_total  = [qdA1; qdC1; qdA2; qdC2; qdC3; qdA3];
tau_total = [tA1; tC1; tA2; tC2; tC3; tA3];
tau_G_tot = [gA1; gC1; gA2; gC2; gC3; gA3];
tau_I_tot = [iA1; iC1; iA2; iC2; iC3; iA3];
tau_C_tot = [cA1; cC1; cA2; cC2; cC3; cA3];
manip_tot = [mA1; mC1; mA2; mC2; mC3; mA3];

nTot = size(Q_total, 1);
t_total = linspace(0, 6, nTot);

seg_bounds = [1, N+1, 2*N+1, 3*N+1, 4*N+1, 5*N+1, 6*N];
seg_labels = {'A1: home→aprox', 'C1: aprox→agarre', ...
              'A2: agarre→trans (CARGA)', 'C2: trans→deposit (CARGA)', ...
              'C3: deposit→retiro', 'A3: deposit→home'};

% -------------------------------------------------------------------------
%  7. VISUALIZACIÓN DE POSICIÓN ARTICULAR A LO LARGO DE LA TAREA
% -------------------------------------------------------------------------
fig1 = figure('Name','Posicion articular','NumberTitle','off','Position',[50 50 1200 700]);
joint_names = {'q_1','q_2','q_3','q_4','q_5','q_6'};
colors = lines(6);
for j = 1:6
    subplot(2,3,j);
    plot(t_total, rad2deg(Q_total(:,j)), 'Color', colors(j,:), 'LineWidth', 1.5);
    hold on;
    for s = 1:length(seg_bounds)-1
        xline(t_total(seg_bounds(s)), '--k', 'Alpha', 0.4);
    end
    xlabel('Tiempo (s)'); ylabel('Ángulo (°)');
    title(['Posición articular ' joint_names{j}]);
    grid on; box on;
end
sgtitle('Posiciones articulares — Trayectoria completa', 'FontSize', 13, 'FontWeight', 'bold');

% -------------------------------------------------------------------------
%  8. VISUALIZACIÓN DE VELOCIDAD ARTICULAR
% -------------------------------------------------------------------------
fig2 = figure('Name','Velocidad articular','NumberTitle','off','Position',[50 50 1200 700]);
for j = 1:6
    subplot(2,3,j);
    plot(t_total, rad2deg(Qd_total(:,j)), 'Color', colors(j,:), 'LineWidth', 1.5);
    hold on;
    for s = 1:length(seg_bounds)-1
        xline(t_total(seg_bounds(s)), '--k', 'Alpha', 0.4);
    end
    xlabel('Tiempo (s)'); ylabel('Velocidad (°/s)');
    title(['Velocidad articular ' joint_names{j}]);
    grid on; box on;
end
sgtitle('Velocidades articulares — Trayectoria completa', 'FontSize', 13, 'FontWeight', 'bold');

% -------------------------------------------------------------------------
%  9. PARES ARTICULARES TOTALES τ
% -------------------------------------------------------------------------
fig3 = figure('Name','Par articular total','NumberTitle','off','Position',[50 50 1200 700]);
for j = 1:6
    subplot(2,3,j);
    plot(t_total, tau_total(:,j), 'Color', colors(j,:), 'LineWidth', 1.5);
    hold on;
    for s = 1:length(seg_bounds)-1
        xline(t_total(seg_bounds(s)), '--k', 'Alpha', 0.4);
    end
    xlabel('Tiempo (s)'); ylabel('Par (N·m)');
    title(['Par total ' joint_names{j}]);
    grid on; box on;
end
sgtitle('Pares articulares totales τ — Trayectoria completa', 'FontSize', 13, 'FontWeight', 'bold');

% -------------------------------------------------------------------------
%  10. COMPONENTES DINÁMICOS: GRAVEDAD, INERCIA, CORIOLIS
% -------------------------------------------------------------------------
% --- Par gravitacional ---
fig4 = figure('Name','Par gravitacional','NumberTitle','off','Position',[50 50 1200 700]);
for j = 1:6
    subplot(2,3,j);
    plot(t_total, tau_G_tot(:,j), 'Color', colors(j,:), 'LineWidth', 1.5);
    hold on;
    for s = 1:length(seg_bounds)-1
        xline(t_total(seg_bounds(s)), '--k', 'Alpha', 0.4);
    end
    xlabel('Tiempo (s)'); ylabel('Par (N·m)');
    title(['Par gravedad G(q) — ' joint_names{j}]);
    grid on; box on;
end
sgtitle('Componente gravitacional G(q) — Trayectoria completa', 'FontSize', 13, 'FontWeight', 'bold');

% --- Par inercial ---
fig5 = figure('Name','Par inercial','NumberTitle','off','Position',[50 50 1200 700]);
for j = 1:6
    subplot(2,3,j);
    plot(t_total, tau_I_tot(:,j), 'Color', colors(j,:), 'LineWidth', 1.5);
    hold on;
    for s = 1:length(seg_bounds)-1
        xline(t_total(seg_bounds(s)), '--k', 'Alpha', 0.4);
    end
    xlabel('Tiempo (s)'); ylabel('Par (N·m)');
    title(['Par inercia M(q)q̈ — ' joint_names{j}]);
    grid on; box on;
end
sgtitle('Componente inercial M(q)q̈ — Trayectoria completa', 'FontSize', 13, 'FontWeight', 'bold');

% --- Par de Coriolis ---
fig6 = figure('Name','Par Coriolis','NumberTitle','off','Position',[50 50 1200 700]);
for j = 1:6
    subplot(2,3,j);
    plot(t_total, tau_C_tot(:,j), 'Color', colors(j,:), 'LineWidth', 1.5);
    hold on;
    for s = 1:length(seg_bounds)-1
        xline(t_total(seg_bounds(s)), '--k', 'Alpha', 0.4);
    end
    xlabel('Tiempo (s)'); ylabel('Par (N·m)');
    title(['Par Coriolis C(q,q̇)q̇ — ' joint_names{j}]);
    grid on; box on;
end
sgtitle('Componente de Coriolis C(q,q̇)q̇ — Trayectoria completa', 'FontSize', 13, 'FontWeight', 'bold');

% -------------------------------------------------------------------------
%  11. POSICIÓN CARTESIANA DEL EXTREMO A LO LARGO DE LA TAREA
% -------------------------------------------------------------------------
px = zeros(nTot,1); py = zeros(nTot,1); pz = zeros(nTot,1);
for k = 1:nTot
    T = robot.fkine(Q_total(k,:));
    px(k) = T.t(1);
    py(k) = T.t(2);
    pz(k) = T.t(3);
end

fig7 = figure('Name','Posicion cartesiana','NumberTitle','off','Position',[50 50 1000 700]);
subplot(3,1,1);
plot(t_total, px*1000, 'b', 'LineWidth', 1.5); hold on;
for s = 1:length(seg_bounds)-1, xline(t_total(seg_bounds(s)),'--k','Alpha',0.4); end
ylabel('p_x (mm)'); title('Posición cartesiana del extremo'); grid on; box on;

subplot(3,1,2);
plot(t_total, py*1000, 'r', 'LineWidth', 1.5); hold on;
for s = 1:length(seg_bounds)-1, xline(t_total(seg_bounds(s)),'--k','Alpha',0.4); end
ylabel('p_y (mm)'); grid on; box on;

subplot(3,1,3);
plot(t_total, pz*1000, 'g', 'LineWidth', 1.5); hold on;
for s = 1:length(seg_bounds)-1, xline(t_total(seg_bounds(s)),'--k','Alpha',0.4); end
ylabel('p_z (mm)'); xlabel('Tiempo (s)'); grid on; box on;

sgtitle('Posición cartesiana X, Y, Z del efector final', 'FontSize', 13, 'FontWeight', 'bold');

% Trayectoria 3D del extremo
fig8 = figure('Name','Trayectoria 3D','NumberTitle','off','Position',[50 50 800 600]);
plot3(px*1000, py*1000, pz*1000, 'b-', 'LineWidth', 1.5);
hold on;
% Marcar inicio y fin de cada segmento
seg_starts = seg_bounds(1:end-1);
scatter3(px(seg_starts)*1000, py(seg_starts)*1000, pz(seg_starts)*1000, ...
         60, 'r', 'filled');
for s = 1:length(seg_starts)
    text(px(seg_starts(s))*1000, py(seg_starts(s))*1000, pz(seg_starts(s))*1000, ...
         ['  S' num2str(s)], 'FontSize', 8);
end
xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');
title('Trayectoria 3D del efector final', 'FontSize', 13, 'FontWeight', 'bold');
grid on; box on; view(45, 30);
legend('Trayectoria', 'Inicio segmento', 'Location', 'best');

% -------------------------------------------------------------------------
%  12. MANIPULABILIDAD |det(J)| A LO LARGO DE TODA LA TRAYECTORIA
% -------------------------------------------------------------------------
fig9 = figure('Name','Manipulabilidad','NumberTitle','off','Position',[50 50 1000 400]);
semilogy(t_total, manip_tot, 'm-', 'LineWidth', 1.5);
hold on;
for s = 1:length(seg_bounds)-1
    xline(t_total(seg_bounds(s)), '--k', 'Alpha', 0.4, ...
          'Label', ['S' num2str(s)], 'LabelVerticalAlignment', 'top', 'FontSize', 8);
end
yline(1e-4, ':r', 'LineWidth', 1.2, 'Label', 'Umbral singularidad');
xlabel('Tiempo (s)');
ylabel('|det(J)|  (escala log)');
title('Manipulabilidad — |det(J)| a lo largo de la trayectoria completa', ...
      'FontSize', 13, 'FontWeight', 'bold');
grid on; box on;

% -------------------------------------------------------------------------
%  13. COMPARATIVA: PAR TOTAL CON Y SIN CARGA (articulación 2)
% -------------------------------------------------------------------------
% Recalcular A2 sin carga para comparar
robot.payload(0, load_pos);
[tA2_nc, ~, ~, ~, ~] = calc_dynamics(robot, qA2, qdA2, qddA2, false, 0, load_pos);
robot.payload(2.0, load_pos);
[tA2_c,  ~, ~, ~, ~] = calc_dynamics(robot, qA2, qdA2, qddA2, true,  2.0, load_pos);
robot.payload(0, load_pos);

t_A2 = linspace(0, 1, N);
fig10 = figure('Name','Efecto de la carga','NumberTitle','off','Position',[50 50 900 400]);
plot(t_A2, tA2_nc(:,2), 'b--', 'LineWidth', 1.5, 'DisplayName', 'Sin carga');
hold on;
plot(t_A2, tA2_c(:,2), 'r-',  'LineWidth', 1.5, 'DisplayName', 'Con carga (2 kg)');
xlabel('Tiempo normalizado (s)');
ylabel('Par τ_2 (N·m)');
title('Efecto de la carga en el par de la articulación 2 (tramo agarre→trans)', ...
      'FontSize', 12, 'FontWeight', 'bold');
legend('Location','best'); grid on; box on;

fprintf('\n=============================================================\n');
fprintf('   RESUMEN — PLANIFICACIÓN DE TRAYECTORIAS ABB IRB-140\n');
fprintf('=============================================================\n');
fprintf('\nTrayectorias articulares (jtraj, perfil polinomial 5.º orden):\n');
fprintf('  A1: home           → aproximación    (%d muestras)\n', N);
fprintf('  A2: agarre         → transporte       (%d muestras) [CON CARGA]\n', N);
fprintf('  A3: depósito       → home             (%d muestras)\n', N);

fprintf('\nTrayectorias cartesianas (ctraj, perfil trapezoidal):\n');
fprintf('  C1: aproximación   → agarre           (%d muestras)\n', N);
fprintf('  C2: transporte     → depósito         (%d muestras) [CON CARGA]\n', N);
fprintf('  C3: depósito       → retiro           (%d muestras)\n', N);

fprintf('\nCarga externa: 2.0 kg a 8 cm del flange (eje Z)\n');
fprintf('Activa durante segmentos: A2, C2\n');

fprintf('\nManipulabilidad mínima a lo largo de la tarea: %.4e\n', min(manip_tot));
fprintf('Manipulabilidad máxima a lo largo de la tarea: %.4e\n', max(manip_tot));
n_sing = sum(manip_tot < 1e-4);
if n_sing > 0
    fprintf('AVISO: %d muestras próximas a singularidad (|det(J)| < 1e-4)\n', n_sing);
else
    fprintf('La trayectoria se mantiene alejada de singularidades.\n');
end

fprintf('\nPar máximo articulación 2 (sin carga): %.2f N·m\n', max(abs(tA2_nc(:,2))));
fprintf('Par máximo articulación 2 (con carga): %.2f N·m\n', max(abs(tA2_c(:,2))));
fprintf('Incremento de par por la carga: %.1f %%\n', ...
    100*(max(abs(tA2_c(:,2))) - max(abs(tA2_nc(:,2)))) / max(abs(tA2_nc(:,2))));
fprintf('=============================================================\n');

% -------------------------------------------------------------------------
%  14. ANIMACIÓN DEL ROBOT
% -------------------------------------------------------------------------
 figure('Name','Animacion ABB IRB-140','Position',[50 50 900 700]);
 robot.plot(Q_total, 'trail', 'b-', 'fps', 30);