#!/bin/bash

echo ""
echo "****************************"
echo "***    V = 0, K = 0      ***"
echo "****************************"

for i in `ls ../traces/`; do
    echo ""
    echo "--- $i ---"
	./cachesim -v 0 -k 0 -i ../traces/$i
done;

