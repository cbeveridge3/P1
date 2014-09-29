

%first, initialize time step
dt = .01;
t = 0:dt:100;
ts = 5/dt;

%No stimulus
I = 0;

%set up constants--------------------------------------------------
%These are the maximum conductances for each ion in (mS/cm^2)
gbK = 36;
gbNa = 120;
gbL = 0.3;

%Membrane capacitence is standard (microF/cm^2)
Cm = 1;

%these are the Nernst Potentials in mV
EK = -12;
ENa = 115;
EL = 10.6;

%this is the resting potential in mV
Vrest = -70;

%Initialize displacements
Vm = Vrest; % Vm is potential inside minus potential putside
V = abs(Vm)-abs(Vrest); % difference from Vm and resting value
VK = EK-abs(Vrest);
VNa = ENa-abs(Vrest);
VL = EL-abs(Vrest);

%initialize gating variables (rate constants)
alphaM = .1 * ((25 - V) / (exp((25 - V) / 10) - 1));
betaM = 4*exp(-V/18);
alphaN = .01 * ((10 - V) / (exp((10-V) / 10) - 1));
betaN = .125*exp(-V/80);
alphaH = .07*exp(-V/20);
betaH = 1 / (exp((30 - V) / 10) + 1);


%initialize particle probability 
n = alphaN / (alphaN + betaN);
m = alphaM / (alphaM + betaM);
h = alphaH / (alphaH + betaH);

%initialize conductance
gK = [n^4 * gbK];
gNa = [m^3 * gbNa * h];

%Run a for loop to caluclate Vm at each time step
for i = 1:length(t)-1
    
    % calculate gating variables
    alphaM = .1 * ((25 - V(i)) / (exp((25 - V(i)) / 10) - 1));
    betaM = 4*exp(-V(i)/18);
    alphaN = .01 * ((10 - V(i)) / (exp((10-V(i)) / 10) - 1));
    betaN = .125*exp(-V(i)/80);
    alphaH = .07*exp(-V(i)/20);
    betaH = 1 / (exp((30 - V(i)) / 10) + 1);

    
    dn = alphaN * (1-n) - betaN*n;
    dm = alphaM * (1-m) - betaM*m;
    dh = alphaH * (1-h) - betaH*h;

    %Calculate probabilities
    n = n + dn*dt;
    m = m + dm*dt;
    h = h + dh*dt;
    
    %Determine if stimulation is on.
%     if i <= (ts+150) && i >= 150
%         I = 50;
%     else
%         I = 0;
%     end

    
    % calculate conductances
    gNa(i+1) = m^3*gbNa*h;
    gK(i+1) = n^4*gbK;
    
    % calculate currents
    INa = gNa(i+1)*(Vm(i)-VNa);
    IK = gK(i+1)*(Vm(i)-VK);
    IL = gbL*(Vm(i)-VL);
    
    
    %calculation total current
    Iion = (I-IK-INa-IL);
    
    % calculate membrane potential
    dm = Iion/Cm;
    V(i+1) = V(i)+dm*dt;
    Vm(i+1) = Vrest+V(i+1);
end

% plot Vm
figure(1)
plot(t,Vm)
xlabel('Time (ms)')
ylabel('Membrane Potential (mV)')
axis([0 100 -100 50]);

% plot conductances
figure(2)
plot(t,gNa,t,gK)
legend('g_{Na}','g_K')
xlabel('Time (ms)')
ylabel('Conductance (mS/cm^2)')
