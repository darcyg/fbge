#!/bin/bash
gcc -o run -O2 -I/usr/include/ncursesw -std=gnu99 -ggdb *.c -lm -lncursesw
