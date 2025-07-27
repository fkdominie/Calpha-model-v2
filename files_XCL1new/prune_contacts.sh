#!/bin/bash

### PDB numbering (1-93) indicated with "r" in file names (e.g r1-93). 

### *** 1j8i (monomer state) ***
### Allow contacts within region 9-68. "Mirror" contacts ij in chain 1 to chain 2 (i+93,j+93).
### Exclude disulfide bonds 11,48 and 104,141.

awk '{ if ( ($2<9) || ($2>68) || ($4<9) || ($4>68) || ($2==11 && $4==48) ) {} else print $2-1 " " $4-1; }' \
    1j8i.2359.pdb.contacts > smog_1j8i_r9-68_mirror
awk '{ if ( ($2<9) || ($2>68) || ($4<9) || ($4>68) || ($2==11 && $4==48) ) {} else print $2-1+93 " " $4-1+93; }' \
    1j8i.2359.pdb.contacts >> smog_1j8i_r9-68_mirror

### *** 2jp1 (dimer state) ***
### Allow contacts within and between regions 8-52 (chain 1) and res 101-145 (chain 2).
### Remove disulfide bonds 11,48 and 104,141. 

awk '{ if ( ($2<8) || ($2>52) || ($4<8) || ($4>52) || ($2==11 && $4==48) ) {} else print $2-1 " " $4-1; }' \
    2jp1_dimer.2363.pdb.contacts > smog_2jp1_r8-52
awk '{ if ( ($2<101) || ($2>145) || ($4<101) || ($4>145) || ($2==104 && $4==141) ) {} else print $2-1 " " $4-1; }' \
    2jp1_dimer.2363.pdb.contacts >> smog_2jp1_r8-52
awk '{ if ( ($2>=8 && $2<=52 && $4>=101 && $4<=145) ) print $2-1 " " $4-1; }' \
    2jp1_dimer.2363.pdb.contacts >> smog_2jp1_r8-52

### "Symmetrize" intra-chain contacts of 2jp1: 
### Include contacts ij in chain 1 if the corresponding contact (i+93,j+93) 
### is present in chain 2, and vice versa for chain 2 contacts. The following
### contacts should therefore be removed  (see contacts_to_exclude_2jp1.c):
### (13,17); (13,40); (28,36); (38,46); (100,141); (106,139); (110,137); (126,144); (132,138);

awk '{ if ( ($1==13 && $2==17) || ($1==13 && $2==40) || ($1==28 && $2==36) || ($1==38 && $2==46) \
         || ($1==100 && $2==141) || ($1==106 && $2==139) || ($1==110 && $2==137) \
	 || ($1==126 && $2==144) || ($1==132 && $2==138) ) {} else print $0;}' \
	 smog_2jp1_r8-52 > smog_2jp1_r8-52_sym



### Swap chains in the dimer state

#awk '{if ($1>=93) {print $1-93 " " $2 " " $3 " " $4;} else {print $1+93 " " $2 " " $3 " " $4;} }' \
#    native_2jp1 > native_2jp1_swap

### Indicate disordered chain regions (column 5):
### 1j8i_9-67 
### 2jp1_9-51 

#awk '{if ( ($1>=8 && $1<=66) || ($1>=8+93 && $1<=66+93) ) {print $0 " 0";} else {print $0 " 1";} }' \
#    native_1j8i > native_1j8i_disreg

#awk '{if ( ($1>=8 && $1<=50) || ($1>=8+93 && $1<=50+93) ) {print $0 " 0";} else {print $0 " 1";} }' \
#    native_2jp1 > native_2jp1_disreg

#awk '{if ( ($1>=8 && $1<=50) || ($1>=8+93 && $1<=50+93) ) {print $0 " 0";} else {print $0 " 1";} }' \
#    native_2jp1_swap > native_2jp1_swap_disreg


### Special: Segment 9-51 of 2jp1 (ordered region).

#awk '{if ($1>=8 && $1<=50) {print a++ " " $2 " " $3 " " $4;}}' native_2jp1 > native_2jp1_9-51_ordreg1
#awk '{if ($1>=8 && $1<=50) {print a++ " " $2 " " $3 " " $4;}}' native_2jp1_swap > native_2jp1_9-51_ordreg2

#awk '{print $1-8 " " $2-8;}' smog_2jp1_9-51_symmetrized > smog_2jp1_9-51_ordreg
#awk '{print $1-8 " " $2-8;}' smog_2jp1_9-51_sym > smog_2jp1_9-51_sym_ordreg

