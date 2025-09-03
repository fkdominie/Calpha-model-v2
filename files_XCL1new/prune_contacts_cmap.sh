#!/bin/bash

### *** 1j8i (monomer state) ***
### Allow contacts within region 9-68. "Mirror" contacts ij in chain 1 to chain 2 (i+93,j+93).
### Exclude disulfide bonds 11,48 and 104,141.

awk '{ if ( ($1<8) || ($1>67) || ($2<8) || ($2>67) || ($1==10 && $2==47) ) {} else print $1 " " $2; }' \
    1j8i_pdb2cont.contacts > cmap_1j8i_r9-68_lc3_rc4.5_mirror
awk '{ if ( ($1<8) || ($1>67) || ($2<8) || ($2>67) || ($1==10 && $2==47) ) {} else print $1+93 " " $2+93; }' \
    1j8i_pdb2cont.contacts >> cmap_1j8i_r9-68_lc3_rc4.5_mirror

### *** 2jp1 (dimer state) ***
### Allow contacts within and between regions 8-52 (chain 1) and res 101-145 (chain 2).
### Remove disulfide bonds 11,48 and 104,141. 

awk '{ if ( ($1<7) || ($1>51) || ($2<7) || ($2>51) || ($1==10 && $2==47) ) {} else print $1 " " $2; }' \
    2jp1_pdb2cont.contacts > cmap_2jp1_r8-52_lc3_rc4.5
awk '{ if ( ($1<100) || ($1>144) || ($2<100) || ($2>144) || ($1==103 && $2==140) ) {} else print $1 " " $2; }' \
    2jp1_pdb2cont.contacts >> cmap_2jp1_r8-52_lc3_rc4.5
awk '{ if ( ($1>=7 && $1<=51 && $2>=100 && $2<=144) ) print $1 " " $2; }' \
    2jp1_pdb2cont.contacts >> cmap_2jp1_r8-52_lc3_rc4.5

### "Symmetrize" intra-chain contacts of 2jp1: 
### Include contacts ij in chain 1 if the corresponding contact (i+93,j+93) 
### is present in chain 2, and vice versa for chain 2 contacts. The following
### contacts should therefore be removed  (see contacts_to_exclude_2jp1.c):
### CMAP CONTACTS: (17,42); (23,38); (32,49); (37,46); (40,44); (108,137); (110,132); (120,127); (127,143); 

awk '{ if ( ($1==17 && $2==42) || ($1==23 && $2==38) || ($1==32 && $2==49) || ($1==37 && $2==46) \
         || ($1==40 && $2==44) || ($1==108 && $2==137) || ($1==110 && $2==132) \
	 || ($1==120 && $2==127) || ($1==127 && $2==143) ) {} else print $0;}' \
    cmap_2jp1_r8-52_lc3_rc4.5 > cmap_2jp1_r8-52_lc3_rc4.5_sym


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

