#!/bin/bash

tessgrd -r70/140/0/55 -b141/111 -z260e03 | \
tesspot dem-tess.txt -lZ-GP-V.log -v | \
tessgx dem-tess.txt -lZ-GV-Vx.log -v | \
tessgy dem-tess.txt -lZ-GV-Vy.log -v | \
tessgz dem-tess.txt -lZ-GV-Vz.log -v | \
tessgxx dem-tess.txt -lZ-GGT-Vxx.log -v | \
tessgxy dem-tess.txt -lZ-GGT-Vxy.log -v | \
tessgxz dem-tess.txt -lZ-GGT-Vxz.log -v | \
tessgyy dem-tess.txt -lZ-GGT-Vyy.log -v | \
tessgyz dem-tess.txt -lZ-GGT-Vyz.log -v | \
tessgzz dem-tess.txt -lZ-GGT-Vzz.log -v | \
tessgxxx dem-tess.txt -lZ-GC-Vxxx.log -v | \
tessgxxy dem-tess.txt -lZ-GC-Vxxy.log -v | \
tessgxxz dem-tess.txt -lZ-GC-Vxxz.log -v | \
tessgxyz dem-tess.txt -lZ-GC-Vxyz.log -v | \
tessgyyx dem-tess.txt -lZ-GC-Vyyx.log -v | \
tessgyyy dem-tess.txt -lZ-GC-Vyyy.log -v | \
tessgyyz dem-tess.txt -lZ-GC-Vyyz.log -v | \
tessgzzx dem-tess.txt -lZ-GC-Vzzx.log -v | \
tessgzzy dem-tess.txt -lZ-GC-Vzzy.log -v | \
tessgzzz dem-tess.txt -lZ-GC-Vzzz.log -v > Result_China_30fen30fen.txt