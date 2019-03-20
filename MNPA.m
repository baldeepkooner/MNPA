%{
G = [G1 -G1 0 0 0;
     -G1 G1+G2 0 0 0;
     0 0 G3 0 0;
     0 0 0 G4 -G4;
     0 0 0 -G4 G0+G4]

B = [1 0;
     0 0;
     0 0;
     0 1;
     0 0]

C = [0 0 0 0 0 0 0;
     -C C 0 0 0 0 0; 
     0 0 L 0 0 0 0; 
     0 0 0 0 0 0 0; 
     0 0 0 0 0 0 0; 
     0 0 0 0 0 0 0; 
     0 0 0 0 0 0 0] 

G = 1/R

N1:    i1 + i2 = 0
       i1 + [(V1-V2)/R1 + C(d(V1-V2)/dt)] = 0 
       V1 = Vin (1)
N2:   i2 + i3 + i4 = 0
      [(Vin-V2)/R1 + C(d(Vin-V2)/dt)] + V2/R2 + iL (2)
N3:   iL + i3 = 0
      iL + V3/R3 = 0 (3)
N4:   alpha*i3 + i4 = 0
      alpha*i3 + (V4-V5)/R4 = 0 (4)
      V4 = alpha*i3 (5)
N5:   i4 + i0 = 0 
      (V5-V4)/R4 + V0/R0 = 0 (6)

N1:
      i1 + [(V1-V2)G1 + (V1-V2)sC] = 0 (1)
      V1 = Vin (2)
N2: 
      [(V2-V1)G1 + (V2-V1)sC] + V2G2 + (V2-V3)sL (3)
N3: 
      (V3-V2)sL + V3G3 = 0 (4)
N4:   
      alpha*i3 + (V4-V5)G4 = 0 (5)
      V4 = alpha*i3 (6)
N5:   
      (V5-V4)G4 + V0G0 = 0 (7)

Where s = jw

    V1,     V2, IL,   V3,       V4,     Vo

G:
      1,     0,  0,     0,        0,     0 
    -G1, G1+G2,  1,     0,        0,     0 
      0,    -1,  0,     1,        0,     0 
      0,     0, -1,    G3,        0,     0 
      0,     0,  0,  a*G3,        1,     0  
      0,     0,  0,     0,      -G4, G4+GO 
    

C:
      0,       0,  0,  0,     0,   0, 
     -C,       C,  0,  0,     0,   0, 
      0,       0,  L,  0,     0,   0, 
      0,       0,  0,  0,     0,   0, 
      0,       0,  0,  0,     0,   0, 
      0,       0,  0,  0,     0,   0, 
     
%}

clear all
close all

G = zeros(6, 6); 

%Resistances:
R1 = 1;
R2 = 2;
R3 = 10;
R4 = 0.1; 
R0 = 1000; 

%Conductances:
G1 = 1/R1;
G2 = 1/R2;
G3 = 1/R3;
G4 = 1/R4;
G0 = 1/R0;

%Additional Parameters:
alpha = 100;
Cval = 0.25;
L = 0.2;
vin = zeros(1, 20);
vo = zeros(1, 20);
v3 = zeros(1, 20);

G(1, 1) = 1;                                    % 1
G(2, 1) = -G1; G(2, 2) = G1 + G2; G(2, 3) = 1;  % 2
G(3 ,2) = -1; G(3, 4) = 1;                      % iL
G(4, 3) = -1; G(4, 4) = G3;                     % 3
G(5, 5) = 1; G(5, 4) = alpha*G3;                % 4
G(6, 6) = G4 + G0; G(6, 5) = -G4;               % 5

C = zeros(6, 6);

C(2, 1) = -Cval; C(2, 2) = Cval;
C(3, 3) = L;

F = zeros(1, 6);
v = -10;

for i = 1:21
    vin(i) = v;
    F(1) = vin(i);
    
    Vm = G\F';
    
    vo(i) = Vm(6);
    v3(i) = Vm(4);
    v = v + 1;
end


figure(1)
plot(vin, vo);
title('VO for DC Sweep (Vin): -10 V to 10 V');
xlabel('Vin (V)')
ylabel('Vo (V)')

figure(2)
plot(vin, v3)
title('V3 for DC Sweep (Vin): -10 V to 10 V')
xlabel('Vin (V)')
ylabel('V3 (V)')

%{
N1:
      i1 + [(V1-V2)G1 + (V1-V2)/sC] = 0 (1)
      V1 = Vin (2)
N2: 
      [(V2-V1)G1 + (V2-V1)/sC] + V2G2 + (V2-V3)sL (3)
N3: 
      (V3-V2)sL + V3G3 = 0 (4)
N4:   
      alpha*i3 + (V4-V5)G4 = 0 (5)
      V4 = alpha*i3 (6)
N5:   
      (V5-V4)G4 + V0G0 = 0 (7)

Where s = jw
%}

F(1) = 1;
vo2 = zeros(1, 1000); 
freq = linspace(0, 1000, 1000); % note: in radians
Av = zeros(1, 1000);
Avlog = zeros(1, 1000);

for i = 1:1000
    Vm2 = (G+1i*freq(i)*C)\F';
    vo2(i) = Vm2(6);
    Av(i) = vo2(i)/F(1);
    Avlog(i) = log10(Av(i));
end 
    
figure(3)
semilogx(freq, Avlog)
xlim([0 1000])
title('Vo(w) dB (part C)')
xlabel('w (rad)')
ylabel('Av (dB)')
    
w = pi;
Av2 = zeros(15, 1);
Cper = zeros(15, 1);
vo3 = zeros(1, 15);

for i = 1:1000
    C(2, 1) = normrnd(-Cval, 0.05); 
    C(2, 2) = normrnd(Cval, 0.05);
    C(3, 3) = normrnd(L, 0.05);
    Vm3 = (G+1i*w*C)\F';
    vo3(i) = Vm3(6);
    Av2(i) = vo3(i)/F(1);
end

figure(4)
hist(real(Av2), 25)
title('Gain Distribution for Normal Perturbations in C matrix')
xlabel('Gain at w = pi')


    