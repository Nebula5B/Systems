% Luke Palmer % Lucas Hall % Rahul Boddupalli % Jason Lange 
%Mech 3140-002
%Dr. Martin
%MWF 9:00-9:50
%Midterm Project 
clc; clf;

%% Section 1
%% Part A:
% Develop the model that includes the motor inductance.
% 
% Referance pgs. 1-2 of handwritten work

%% Part B:
% Develop a model that neglects the motor inductance. What are J eff and b eff?
% 
% Referance pgs. 3-4 of handwritten work
% Jeff = J*R/Kt
% beff = (R*b + Kt*Kb)/Kt

%% Part C:
% Look up the following values for the 24-volt KN 16 M4 LR motor 
% (specification sheet on Canvas). Note any discrepancies in the published values.
% 
% Referance Pg. 5 of handwritten work

% k_b = 5.8434 * 10^-3;                 % Back EMF Constant; V-s/rad
k_b = .06685;            % Back EMF Constant; V-s/rad
R = .75;                 % Armature Resistance; Ohms
k_t = .067;              % Motor Torque Constant; N-m/A
M_0 = 2.36;              % Small Torque; N*m
omega_nls = 628.31853;   % No Load Speed ; rad/s
J = 5.94 * 10^-4;        % Rotor MMOI; kg-m^
% b = 6.493 * 10^-6;                    % Viscous Damping Constant; N-m-s/rad
b = 0.000592;            % Viscous Damping Constant; N-m-s/rad
tau_motor = 3.31 * 10^-3;% Motor Time Constant ; s
L_a = 2.5 * 10^-3;       % Armature Inductance; H
Vi = 12 ;                % Input Voltage; V
r = .015*2*pi;           %m
I_m = .05;               %Nm
m =  .499;               %kgs
k = 1.27;                %N/m
C = 2.5*10^-2;           %N/m/s

%% Part E: 
% What are the eigenvalue(s) for the system (for both models)? What are the corresponding 
% time constant(s)?
%
% Referance pgs. 1-4 of handwritten work
% calculate exact values for each eiganvalue and time constant and write
% in !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% With Inductance
s_1 = (-J*R - L_a*b + sqrt((J*R + L_a*b)^2 - 4 * J*L_a * (R*b + k_b*k_t)))/(2*J*L_a);
s_2 = (-J*R - L_a*b - sqrt((J*R + L_a*b)^2 - 4 * J*L_a * (R*b + k_b*k_t)))/(2*J*L_a);
tf_L = tf([k_t],[J*L_a, J*R + L_a*b, R*b + k_b*k_t]);
% Without Inductance
s_3 = (-R*b - k_t*k_b)/(J*R);
tf_i = tf([k_t],[J*R, R*b + k_t*k_b]);
%% Part F:
% Simulate the step response of the motor to an input voltage of 12 Volts. Compare the
% model that includes motor inductance and the model without the motor inductance (plot
% the speed output on the same graph). Explain the comparison.

figure(1)
hold on
step(tf_L,'b+')
step(tf_i,'ro')
grid on
xlabel('Time')
ylabel('Motor Speed (rpm)')
legend('System With Inductance','System Without Inductance','location','southeast')

%%
% Figure 1 shows that the angular speed output of both the system without
% inductance and with inductance are incrediblly similiar.
%% Part H:
% Determine and plot the frequency response (gain and phase) for the 1st and 2nd order
% models on the same plot for comparison and discuss the results. Pick a frequency to verify
% your analytical frequency response using the simulation provided
% Bode plot for both circuits

% plot both of the bode plots
figure(2)
hold on
bode(tf_L)
bode(tf_i)
grid on
legend('System With Inductance','System Without Inductance','location','southwest')

t = [0:.0001:100];
A = Vi ;                          %Volts 

% pulled frequencies from Bode plot at various locations 
omega = [.1 100 1000];            %rad/s

% Frequency response of circuit w/ and w/o Inductance
z = 3;                            %figure number
n = 1;
while z <= 5
Gain = k_t/(sqrt((R*b + k_b*k_t - J*L_a*omega(n)^2)^2 + (J*R*omega(n) + L_a*b*omega(n)*omega(n))^2));
Phase = -atan((J*R*omega(n) + L_a*b*omega(n))/(R*b + k_b*k_t - J*L_a*omega(n)^2));
Vp = A*(Gain)*sin(omega(n)*t + Phase);
Gain_2 = k_t / (sqrt((J*R*omega(n))^2 + (R*b + k_t*k_b)^2));
Phase_2 = -atan((J*R*omega(n))/(R*b + k_t*k_b));
Vp_2 = A*(Gain_2)*sin(omega(n)*t + Phase_2);
figure(z)
hold on
grid on
yyaxis right
plot(t,Vp)
xlabel('time (s)')
ylabel('Volts (V) W/ Inductance')
yyaxis left
plot(t,Vp_2)
ylabel('Volts (V) W/o inductance')
if z == 3
    xlim([0 100])
    title('Frequency Response (omega = .1 rad/s)')
elseif z == 4
    xlim([0 .1])
    title('Frequency Response (omega = 100 rad/s)')
else
    xlim([0 .01])
    title('Frequency Response (omega = 1000 rad/s)')
end
n = n + 1;
z = z + 1;
end
hold off

%%
% Figure 3 above shows that both outputs, w/ and w/o inductance, where frequency is .1
% rad/s are the same. If you refer to the body plot on figure 2 you can tell
% that the gain and phase at this frequency is equal in both cases as well.
% 
% Figure 4 shows both outputs, w/ and w/o inductance, where frequency is
% 100 rad/s. If you refer to the Bode plot once again you can see that the
% phase for both systems is off quite a bit and the gain is still
% relativaly similiar. this results in the outputs having quite a large
% phase change and almost unnoticable change in amplitude.
% 
% Figure 5 shows both outputs, w/ and w/o inductance, where frequency is
% 1000 rad/s. the Bode plot at this frequency shows large differences
% between the two systems in both gain and phase, this translates to the
% output in terms of amplitude and phase change, respectively.

%% Part I:
% Discuss when is it reasonable (or not reasonable) to neglect the motor inductance?
% 
% It is reasonable to neglect motor inductance when your system is
% operating at a low enough frequency because the out put will be the same
% either way. This was observed in Part H and is due to both systems having
% the same gain and phase at these low frequencies.

%% Section 2:
% Consider the rack and pinion (schematic below) where the tire force 
% during turning is modeled as a spring (this is the force that provides 
% what is known as the aligning moment in your car which causes the 
% steering wheel to “straighten” if you remove your hand from the steering 
% wheel in a turn). Assume the DC motor from Section I of the project is 
% attached to the steering wheel as shown in the picture. The steering 
% wheel (not shown in the schematic) has inertia Is. Values for the 
% inertia, damping and spring coefficients are given the table below 
% (note Im in the schematic represents the motor inertia from Section I)

%% Part A:
% Derive the model for the motor-rack and pinion assembly assuming the 
% inductance of the motor is negligible.
% 
% Referance pgs. 6-8 of the hand written work.

%% Part B:
% Using the provided model parameter, develop a MATLAB simulation of the 
% steering assembly with input voltage and output steer angle. Plot the 
% step response of the steering assembly. Is the steering assembly over or 
% under damped? What is the DC gain?
% 
% Reference pg. 9 of handwritten work.

figure(6)
tf_new = tf(1,[(R/k_t)*(I_m + m*r^2), (R/k_t)*(C*r^2 + b+ k_b*k_t/R), (R/k_t)*k*r^2]);
info_tf = stepinfo(tf_new)
step(tf_new)
grid on
xlabel('Time')
ylabel('Theta (deg)')

%%
% The steering assembly appears to be under damped, and DC gain is !!!!!!!!

%%
% Effective values
Jeff = (I_m*R/k_t) + (m*r^2*R/k_t);
Beff = R/k_t*((C*r^2) + (b)+ (k_b*k_t/R));
Keff = R/k_t*k*r^2;

%% Part C:
% An outside company has designed a proportion-derivative controller for the steering system
% (i.e. controller output is the voltage supplied to motor and reference is desired steer angle).
% MATLAB code that simulates the closed loop steering dynamics is available on the class
% website (run_car.p). Use the MATLAB function to determine the proportional and
% derivative control gain and to reconstruct the complete control law
%% Part C.a:
% Simulate the closed loop system and plot the desired steer angle and actual steer
% angle on one graph
% 
% To use run_car.p clear all in command window and Run the program twice

% Runcar Simulation
steer =10; %angle desired
i=1;
t = 0:.01:90;
tire_angle(i) = 0;
steer_act(i) = 0;
for i= 1:length(t)
[GPS,yaw_gyro,del_meas,WP]= run_car(steer*pi/180, 0.001);
tire_angle(i+1) = del_meas;
steer_act(i+1) = tire_angle(i)*12;
end
fclose('all')
figure(7)
plot(t,tire_angle(1:length(t))*180/pi,'k-')
xlabel('Time')
ylabel('Theta (deg)')
hold on
step(tf_new)
grid on
xlabel('Time')
ylabel('Theta (deg)')
legend('Simulation Runcar','Closed Loop','location','southeast')

%% Part C.b:
% Find the closed loop eigenvalues and plot them in the complex plane (s-plane)

figure(8)
pzmap(tf_new)

%% Part C.c:
% Find the proportional and derivative control gains
% 
% Reference Pgs. 10-11 of handwritten work

% Constants from runcar
wd = 16.535;
zeta = .4994;
wn = 19.0849;
kp = wn^2*Jeff - Keff;
kd = 2*zeta*wn*Jeff - Beff;

%% Part C.d:
% Provide the complete control law
tf_control = tf([kd, (kp+Keff)],[Jeff, (Beff+kd), (Keff +kp)]);

%% Part D:
% Verify the controller parameters using the model that you develop in part b. Simulate the
% closed loop response with your simulation. Plot the output from run_car.p and from
% your simulation code on the same graph.

t_ = 0:.01:2;
U = zeros(length(t_),1)+10;
figure(9)
plot(t,tire_angle(1:length(t))*180/pi,'k-')
hold on
lsim(tf_control,U,t_)
xlabel('Time')
ylabel('Theta (deg)')
legend('Simulation Runcar','Closed Loop  with Controller','location','southeast')
