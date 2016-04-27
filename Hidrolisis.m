%Modelo de la cinética de primer orden 
function B= Hidrolisis(Var, t) %Var representas las variables del modelo
Bo=Var(1);
Kh=Var(2);
B= Bo*(1-exp(-Kh*t));
end