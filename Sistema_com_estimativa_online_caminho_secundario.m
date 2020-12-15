%Projeto especializado - Sistema com estimativa online do caminho secundario
clear all; close all; clc;

%% Informacoes do sistema
[nump, denp] = acoustic_path_iir('primary');
[nums, dens] = acoustic_path_iir('secondary');

%% Parametros de simulacao
freq_senoide = 400;
amplitude = 3;
ts = 0.00001;
fs = 1/ts;
tmax = 5;
t = 0:ts:tmax;

%% Parametros das interacoes
J = zeros(1,length(t));
R = 1;

%% Tipos de entradas
senoidal = amplitude*sin(2*pi*freq_senoide*t);
senoidal_h = amplitude*sin(2*pi*freq_senoide*t)+(amplitude/3)*sin(2*pi*3*freq_senoide*t)+(amplitude/5)*sin(2*pi*5*freq_senoide*t);

ruido_branco = sqrt(10)*randn(1,length(t));
%analisa_espectro(ruido_branco,fs);
ruido_colorido = lowpass(ruido_branco,150,fs);
% analisa_espectro(ruido_colorido,fs);

for k = 1:R
    %% Entrada x(n)
    x = senoidal;
    
    %% Ruido branco p(n)
    p = sqrt(10e-6)*randn(1,length(t));
    
    %% Ruido auxiliar v(n)
    v = sqrt(5e-5)*randn(1,length(t));
    
    %% Inicializacao dos parametros
    
    %Espectros
    b_espectro_d = zeros(1,256);
    b_espectro_e = zeros(1,256);
    b_espectro_ud = zeros(1,256);
    b_espectro_vd = zeros(1,256);
    
    %Caminho primario
    d = zeros(1,length(x));
    b_entrada_p = zeros(length(nump),1);
    b_saida_p = zeros(length(denp),1);
    
    %Caminho secundario
    ud = zeros(1,length(x));
    b_entrada_s = zeros(length(nums),1);
    b_saida_s = zeros(length(dens),1);
    
    %Caminho secundario estimado
    N_se = 150;
    mi_se = 0.003;
    vd = zeros(length(v),1);
    b_entrada_se = zeros(N_se,1);
    Se = zeros(N_se,1);
    b_entrada_lms_se = zeros(N_se,1);
    
    %Copia do caminho secundario estimado
    xd = zeros(1,length(x));
    b_entrada_cse = zeros(N_se,1);
    
    %Filtro adaptativo W
    N_w = 100;
    mi_w = 5e-4;
    y = zeros(length(x),1);
    b_entrada_w = zeros(N_w,1);
    W = zeros(N_w,1);
    b_entrada_lms_w = zeros(N_w,1);
    
    %% Execucao do algoritmo ANC
%     figure;
%     pause(0.01);
    for n = 1:length(x)
        
        %Caminho primario - atualiza buffer de entrada
        b_entrada_p(2:end) = b_entrada_p(1:end-1);
        b_entrada_p(1) = x(n);
        
        %Caminho primario - calcula saida
        d(n) = nump*b_entrada_p - denp(2:end)*b_saida_p(1:end-1);
        
        %Caminho primario - atualiza buffer de saida
        b_saida_p(2:end) = b_saida_p(1:end-1);
        b_saida_p(1) = d(n);
        
        %---------------------------
        
        %Filtro adaptativo W(z) - atualiza buffer de entrada
        b_entrada_w(2:end) = b_entrada_w(1:end-1);
        b_entrada_w(1) = x(n);
        
        %Filtro adaptativo - calcula saida
        y(n) = b_entrada_w'*W;
        
        %---------------------------
        
        %Sinal u(n)
        u(n) = y(n) + v(n);
        
        %---------------------------
        
        %Caminho secundario - atualiza buffer de entrada
        b_entrada_s(2:end) = b_entrada_s(1:end-1);
        b_entrada_s(1) = u(n);
        
        %Simulacao da alteracao no sistema
        if n < 250000
            ud(n) = nums*b_entrada_s - dens(2:end)*b_saida_s(1:end-1);
        else
            ud(n) = 1.25.*nums*b_entrada_s - 1.1.*dens(2:end)*b_saida_s(1:end-1);
        end
        
        %Caminho secundario - atualiza buffer de saida
        b_saida_s(2:end) = b_saida_s(1:end-1);
        b_saida_s(1) = ud(n);
        
        %----------------------------
        
        %Calculo do erro e(n)
        e(n) = d(n) - ud(n);
        
        %----------------------------
        
        %Caminho secundario estimado - atualiza buffer de entrada
        b_entrada_se(2:end) = b_entrada_se(1:end-1);
        b_entrada_se(1) = v(n);
        
        %Caminho secundario estimado - calcula saida
        vd(n) = Se'*b_entrada_se;
        %----------------------------
        
        %Calculo do erro f(n)
        f(n) = -vd(n) - e(n);
        
        %----------------------------
        
        %LMS Se - atualiza buffer de entrada
        b_entrada_lms_se(2:end) = b_entrada_lms_se(1:end-1);
        b_entrada_lms_se(1) = v(n);
        
        %LMS Se - calcula e atualiza saida
        Se = Se + mi_se*b_entrada_lms_se*f(n);
        %-----------------------------
        
        %Copia do caminho secundario
        b_entrada_cse(2:end) = b_entrada_cse(1:end-1);
        b_entrada_cse(1) = x(n);
        
        %Caminho primario - calcula saida
        xd(n) = b_entrada_cse'*Se;
        
        %----------------------------
        
        %LMS W - atualiza buffer de entrada
        b_entrada_lms_w(2:end) = b_entrada_lms_w(1:end-1);
        b_entrada_lms_w(1) = xd(n);
        
        %LMS - calcula e atualiza saida
        W = W + mi_w*b_entrada_lms_w*e(n);
        
%         %Apresentacao dos espectros
%         b_espectro_d(2:end) = b_espectro_d(1:end-1);
%         b_espectro_d(1) = d(n);
% 
%         b_espectro_e(2:end) = b_espectro_e(1:end-1);
%         b_espectro_e(1) = e(n);
%         
%         b_espectro_ud(2:end) = b_espectro_ud(1:end-1);
%         b_espectro_ud(1) = ud(n);
%         
%         b_espectro_vd(2:end) = b_espectro_vd(1:end-1);
%         b_espectro_vd(1) = vd(n);
%         
%         
% %         Animação das curvas do Espectro
%             if n >= 255;
%                  cla;
%         %         title('Espectro dos sinais de saída','FontSize',18);
%                 subplot(1,2,1);
%                 plot_espectro(b_espectro_d(1,(end-255:end)),fs, 'b'); % Espectro da entrada Ruído
%                 hold on
%                 plot_espectro(b_espectro_e(1,(end-255:end)),fs, 'r'); % Espectro do erro
%                 xlim([0 5000]);
%                 ylim([0 1.2]);
%                 legend('d(n)', 'e(n)','FontSize',14);
%                 hold off
%         %
%         %         title('Espectro dos sinais de entrada','FontSize',18);
%                 subplot(1,2,2);
%                 plot_espectro(b_espectro_ud(1,(end-255:end)),fs, 'b'); % Espectro saida primario
%                 hold on
%                 plot_espectro(b_espectro_vd(1,(end-255:end)),fs, 'r'); % Espectro saida secundario
%                 xlim([0 5000]);
%                 ylim([0 0.5]);
%                 legend('ud(n)','vd(n)','FontSize',14);
%                 sgtitle('Espectro dos sinais para entrada senoidal de 400 Hz','FontSize',20)
%                 pause(0.01);
%             end
    end
%     J = J + e.^2;
end

%Calcula e apresenta a media do erro
% J = J/R;
% figure
% plot(t,J);
% grid;
% title('Média do erro das interações');

% Figuras
figure
subplot(3,1,1)
plot(t,d)
grid
title('Saída do caminho primário','FontSize',12)
xlabel('Tempo (s)')
ylabel('Amplitude')
legend('d(n)','FontSize',11)
subplot(3,1,2)
plot(t,-ud)
grid
title('Saída do caminho secundário','FontSize',12)
xlabel('Tempo (s)')
ylabel('Amplitude')
legend('-ud(n)','FontSize',11)
subplot(3,1,3)
plot(t,e)
hold on
plot(t,p)
grid
title('Sinal de erro com perturbação','FontSize',12)
xlabel('Tempo (s)')
ylabel('Amplitude')
legend('e(n)','p(n)','FontSize',11)
sgtitle('Resposta do sistema a uma entrada','FontSize',16)

%Comparacao do modelo estimado com o teorico
figure;
freqz(Se,1,256,fs);
hold on;
freqz(nums,dens,256,fs);
lines = findall(gcf,'type','line');
lines(1).Color = 'red';
lines(2).Color = 'black';
lines(3).Color = 'black';
legend('$$\hat{S}$$(z)','S(z)','FontSize',11,'Interpreter','Latex');
sgtitle('Resposta em frequência do caminho secundário','FontSize',16);

% %Resposta em frequencia (FIR) da planta
% figure;
% freqz(nump,denp,256,fs);
% xlabel('Frequência (Hz)');
% ylabel('Magnitude (dB)');
% sgtitle('Resposta em frequência do caminho primário identificado','FontSize',16);
% 
% figure;
% freqz(nums,dens,256,fs);
% xlabel('Frequência (Hz)');
% ylabel('Magnitude (dB)');
% sgtitle('Resposta em frequência do caminho secundário identificado','FontSize',16);

% % Audio
%  entrada = audioplayer(e,fs);
%  playblocking(entrada)
%  audiowrite('e_colorido.wav',e,100000);