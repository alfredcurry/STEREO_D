#!/bin/bash
FOLDER=time/IrrEnv
mkdir -p Results/${FOLDER}/AVG

./IrrEnvStruct2 500 1 1e10 0.01 0.66 300 0
