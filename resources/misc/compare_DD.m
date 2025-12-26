

%%%%% DD -Acoustic Wave Graphs %%%%%%%

iter_Lumped = [9,9,9,10,11,29,34];
iter_NN = [6,6,7,7,7,13,16];
iter_Twolevel = [6,6,6,8,8,8,9];

Nsub = [4,8,16,32,64,128,256];

plot(Nsub,iter_Lumped,'-*',Nsub,iter_NN,'-o',Nsub,iter_Twolevel,'-+')