%Projeto especializado - Sistema com estimaticao offline
clear all; close all; clc;

%% Informacoes do sistema
[nump, denp] = acoustic_path_iir('primary');
[nums, dens] = acoustic_path_iir('secondary');

%% Entrada x(n)
freq_senoide = 400;
amplitude = 1;
ts = 0.0001;
fs = 1/ts;
tmax = 5;
t = 0:ts:tmax;

%% Ruido branco de perturbacao
p = sqrt(10e-6)*randn(1,length(t));

%% Tipos de entradas
senoidal = amplitude*sin(2*pi*freq_senoide*t);
senoidal_h = amplitude*sin(2*pi*freq_senoide*t)+(amplitude/3)*sin(2*pi*3*freq_senoide*t)+(amplitude/5)*sin(2*pi*5*freq_senoide*t);

ruido_branco = sqrt(2)*randn(1,length(t));
%analisa_espectro(ruido_branco,fs);
ruido_colorido = lowpass(ruido_branco,150,fs);
%analisa_espectro(ruido_colorido,fs);

%% Entrada aplicada
x = senoidal;

%% Inicializacao dos parametros

%Espectros para apresentacao da evolucao dos sinais 
b_espectro_e = zeros(1,256);
b_espectro_d = zeros(1,256);
b_espectro_yd = zeros(1,256);

%Caminho primario
d = zeros(1,length(x));
b_entrada_p = zeros(length(nump),1);
b_saida_p = zeros(length(denp),1);

%Caminho secundario estimado
xd = zeros(1,length(x));
b_entrada_se = zeros(length(nums),1);
b_saida_se = zeros(length(dens),1);

%Filtro adaptativo
N = 100;
mi = 0.0001;
y = zeros(length(x),1);
b_entrada_w = zeros(N,1);
w = zeros(N,1);
b_entrada_lms = zeros(N,1);

%Caminho secundario
yd = zeros(1,length(x));
b_entrada_s = zeros(length(nums),1);
b_saida_s = zeros(length(dens),1);

%Equacao a diferencas
% figure;
% pause(0.01);
for n = 1:length(x)

    %Caminho primario - atualiza buffer de entrada
    b_entrada_p(2:end) = b_entrada_p(1:end-1);
    b_entrada_p(1) = x(n);
    
    %Caminho primario - calcula saida
    d(n) = nump*b_entrada_p - denp(2:end)*b_saida_p(1:end-1);
   
    %Caminho primario - atualiza buffer de saida
    b_saida_p(2:end) = b_saida_p(1:end-1);
    b_saida_p(1) = d(n);
    
    %-----------------------------------------
    
    %Caminho secundario estimado - atualiza buffer de entrada
    b_entrada_se(2:end) = b_entrada_se(1:end-1);
    b_entrada_se(1) = x(n);
    
    %Caminho secundario estimado - calcula saida
    xd(n) = nums*b_entrada_p - dens(2:end)*b_saida_s(1:end-1);
    
    %Caminho secundario estimado - atualiza buffer de saida
    b_saida_se(2:end) = b_saida_se(1:end-1);
    b_saida_se(1) = xd(n);
    
    %-----------------------------------------
    
    %Filtro adaptativo - atualiza buffer de entrada
    b_entrada_w(2:end) = b_entrada_w(1:end-1);
    b_entrada_w(1) = x(n);
    
    %Filtro adaptativo - calcula saida
    y(n) = b_entrada_w'*w;
    
    %-----------------------------------------
    
    %Caminho secundario - atualiza buffer de entrada com a saida do filtro adaptativo
    b_entrada_s(2:end) = b_entrada_s(1:end-1);
    b_entrada_s(1) = y(n);
    
    %Simula a variancia do caminho secundario no tempo
    if(n < 20000)
        %Caminho secundario - calcula saida
        yd(n) = nums*b_entrada_s - dens(2:end)*b_saida_s(1:end-1);
    else
        %Caminho secundario - calcula saida
         yd(n) = 0.0000002*nums*b_entrada_s - 0.000008*dens(2:end)*b_saida_s(1:end-1);
    end
    
    %Caminho secundario - atualiza buffer de saida
    b_saida_s(2:end) = b_saida_s(1:end-1);
    b_saida_s(1) = yd(n);
    
    %-----------------------------------------
    
    %Calcula erro com ruido branco
    e(n) = d(n) - yd(n) + p(n);
    
    %-----------------------------------------
    
    %LMS - atualiza buffer de entrada 
    b_entrada_lms(2:end) = b_entrada_lms(1:end-1);
    b_entrada_lms(1) = xd(n);
    
    %LMS - calcula e atualiza saida
    w = w + mi*b_entrada_lms*e(n);
    
    %-----------------------------------------
    
    %Apresentacao dos espectros
    b_espectro_e(2:end) = b_espectro_e(1:end-1);
    b_espectro_e(1) = e(n);
    
    b_espectro_d(2:end) = b_espectro_d(1:end-1);
    b_espectro_d(1) = d(n);
    
    b_espectro_yd(2:end) = b_espectro_yd(1:end-1);
    b_espectro_yd(1) = yd(n);
    
    % Animação das curvas do Espectro
%     if n >= 255;
%          cla;
% %       title('Espectro dos sinais de saída','FontSize',18);
%         subplot(1,2,1);
%         plot_espectro(b_espectro_d(1,(end-255:end)),fs, 'b'); %Saida do caminho primario
%         hold on
%         plot_espectro(b_espectro_e(1,(end-255:end)),fs, 'r'); %Erro
%         ylim([0 1.2])
%         legend('d(n)', 'e(n)','FontSize',14);
%         hold off
% %         
% %       title('Espectro dos sinais de entrada','FontSize',18);
%         subplot(1,2,2);
%         plot_espectro(b_espectro_d(1,(end-255:end)),fs, 'r'); % Espectro saida primario
%         hold on
%         plot_espectro(b_espectro_yd(1,(end-255:end)),fs, 'r'); % Espectro saida secundario
%         ylim([0 1.2])
%         legend('d(n)','yd(n)','FontSize',14);
%         sgtitle('Espectro dos sinais para entrada','FontSize',20)
%         pause(0.01);
%     end
    
end

% Figuras
figure
subplot(3,1,1)
plot(t,d);
grid
title('Saída do caminho primário','FontSize',12)
xlabel('Tempo (s)')
ylabel('Amplitude')
legend('d(n)','FontSize',11)
subplot(3,1,2)
plot(t,-yd)
grid
title('Saída do caminho secundário','FontSize',12)
xlabel('Tempo (s)')
ylabel('Amplitude')
legend('-yd(n)','FontSize',11)
subplot(3,1,3)
plot(t,e)
hold on
plot(t,p)
grid
title('Sinal de erro com perturbação','FontSize',12)
xlabel('Tempo (s)')
ylabel('Amplitude')
legend('e(n)','p(n)','FontSize',11)
sgtitle('Resposta do sistema','FontSize',16)

% % Geracao dos audios
%  entrada = audioplayer(e,fs);
%  playblocking(entrada)
%  audiowrite('e_400Hz.wav',e,fs)