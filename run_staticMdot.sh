#!/bin/bash
FOLDER=And11_Litasov/Static_Mdot #Inward #And11_ZH94/Booth_Grid #
mkdir -p Results/${FOLDER}/AVG
Minit=1.0
Finit=6e5
#grid size , mass , Tbound type , Pbound type , surface temp OR L, Psurf , # , ftime , tmax , Teq . folder , nabtype , MLtype , Massloss , BLtype , corefrac

./MixStruct 800 ${Minit} f f 1200 1 1 1 2e10 ${Finit} ${FOLDER} a n K 0 0.3 2600
