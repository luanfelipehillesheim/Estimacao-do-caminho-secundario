%Funcao para calculo da FFT
function [sinal,freq] = analisa_espectro(entrada,freq_amostragem)

tam = length(entrada);                      %Tamanho do sinal de entrada
k = 0:tam-1;                                %Vetor baseado na entrada
periodo = tam/freq_amostragem;              %Periodo
freq = k/periodo;                           %Eixo das abscissas
sinal = fftn(entrada)/tam;                  %FFT normalizada sobre o sinal de entrada
fc = ceil(tam/2);                           %Frequencia de corte para ajustar o vetor
sinal = sinal(1:fc);                        %Atualiza o sinal de saida

figure();
plot(freq(1:fc),abs(sinal));        
xlabel('Frequência (Hz)','FontSize',16);
ylabel('Amplitude','FontSize',16); 
title('Espectro do sinal','FontSize',18)
end