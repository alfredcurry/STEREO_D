#!/bin/bash
FOLDER=And11_Litasov/Booth_new #Inward #And11_ZH94/Booth_Grid #
mkdir -p Results/${FOLDER}/AVG
Minit=0.15
Finit=6e5
#grid size , mass , Tbound type , Pbound type , Guess temperature , Psurf , # , ftime , tmax , Teq OR initial surface flux , folder , nabtype , MLtype , Massloss type , BLtype , core mass frac , Tss

./MixStruct 800 ${Minit} l f 1700 5e8 1 0.001 2e10 ${Finit} ${FOLDER} c n 0 0 0.3 2320
