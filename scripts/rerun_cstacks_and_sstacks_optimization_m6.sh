#!/bin/bash

echo 'Running cstacks to build catalog'
/programs/stacks-1.48/bin/cstacks -b 3 -p 4 -o data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/ -g -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/CXD-10 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/CXD-3 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/CXD-5 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/CXD-7 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/CXD-9 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/D10 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/D20 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/D22 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/JJ0030 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/JJC44 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/108 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/JJ-109 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/JJ-35 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/JJ-71 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/JJ-92 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/JJ-103 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/JJ-52 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/JJ-97 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/JJC-73 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/LCG09-2 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/LCG09-3 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/LCG09-4 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/LCG09-5 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/117 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/XH-13 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/XHDC-2 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/XHDC-5 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/XHDC-6 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/XHDC-7 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/FQ-6 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/FQ-10 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/FQ-14 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/FQ-26 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/FQ-27 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/CX-29 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/26 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/27 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/28 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/31 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/DF-5 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/DFDC-11 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/DFDC-13 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/DFDC-14 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/DFDC-19 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/Q-6 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/Q-14 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/Q-19 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/Q-20 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/Q-30 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/H-22 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/H-25 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/H-26 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/HDC-3 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/HDC-4 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/HM2 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/HM3 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/HM4 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/HM5 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/HM6 2>&1

echo 'Running sstacks to rematch to catalog'
/programs/stacks-1.48/bin/sstacks -b 3 -p 4 -o data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/ -g -c data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/batch_3 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/CXD-10 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/CXD-3 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/CXD-5 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/CXD-7 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/CXD-9 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/D10 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/D20 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/D22 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/JJ0030 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/JJC44 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/108 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/JJ-109 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/JJ-35 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/JJ-71 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/JJ-92 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/JJ-103 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/JJ-52 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/JJ-97 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/JJC-73 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/LCG09-2 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/LCG09-3 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/LCG09-4 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/LCG09-5 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/117 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/XH-13 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/XHDC-2 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/XHDC-5 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/XHDC-6 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/XHDC-7 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/FQ-6 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/FQ-10 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/FQ-14 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/FQ-26 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/FQ-27 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/CX-29 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/26 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/27 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/28 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/31 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/DF-5 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/DFDC-11 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/DFDC-13 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/DFDC-14 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/DFDC-19 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/Q-6 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/Q-14 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/Q-19 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/Q-20 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/Q-30 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/H-22 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/H-25 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/H-26 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/HDC-3 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/HDC-4 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/HM2 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/HM3 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/HM4 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/HM5 -s data/STACKS_processed/4_depth_optimization/m6/rxstacks_corrected/HM6 2>&1
