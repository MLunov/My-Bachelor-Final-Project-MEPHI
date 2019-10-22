function y = koeff(M2,M1)
if M2/M1 >= 100
    y = 1;
else
    y = 0.01*(M2/M1);
end
end