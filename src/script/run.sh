#!/bin/bash

./ARCS.py -s ../../../d12.cfg -K 25 -d ../../result/d12/25mer/ -p 20 |tee ../../result/d12/25mer/log
#./ARCS.py -s ../../../d12.cfg -K 27 -d ../../result/d12/27mer/ -p 16
./ARCS.py -s ../../../d12.cfg -K 45 -d ../../result/d12/45mer/ -p 20 |tee ../../result/d12/45mer/log
./ARCS.py -s ../../../d12.cfg -K 29 -d ../../result/d12/29mer/ -p 20 |tee ../../result/d12/29mer/log
./ARCS.py -s ../../../d12.cfg -K 31 -d ../../result/d12/31mer/ -p 20 |tee ../../result/d12/31mer/log
./ARCS.py -s ../../../d12.cfg -K 33 -d ../../result/d12/33mer/ -p 20 |tee ../../result/d12/33mer/log
./ARCS.py -s ../../../d12.cfg -K 35 -d ../../result/d12/35mer/ -p 20 |tee ../../result/d12/35mer/log
