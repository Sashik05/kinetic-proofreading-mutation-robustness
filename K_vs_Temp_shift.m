%%%% K vs T plot
%%% K_old = 1; K_new = K*K_old;
close all; clear all; clc;
T = [1:50];
K = (1/293).*((T+293).*exp(-8.4225e03.*((1./(T+293))-(1/293))));
figure
plot(T,K, '--')
xlabel('\Delta T = T_{new} - T_{old}')
ylabel('K_{x, new}/K_{x, old}')