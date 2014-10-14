#!/bin/bash

echo ""
echo "****************************"
echo "***         K = 0        ***"
echo "****************************"

for i in `ls ../traces/`; do
    echo ""
    echo "--- $i ---"
	./cachesim -k 0 -i ../traces/$i
done;

