#!/bin/bash

grep ' CA ' 1j8i_model1.pdb | awk '{print a++ " " $7 " " $8 " " $9 }' - > native_1j8i_model1
grep ' CA ' 1j8i_model1.pdb | awk '{print 93+a++ " " $7+40 " " $8+40 " " $9 }' - >> native_1j8i_model1
grep ' CA ' 2jp1_modeller_8-54_fullchain.pdb | awk '{print a++ " " $7 " " $8 " " $9 }' - > native_2jp1_fullchain

### 1j8i r9-68 ordered
awk '{if ( ($1<8 || $1>67) && ($1<=92) ) print $1 " 1";}' native_1j8i_model1 > tmp1
awk '{ print $1+93 " 1";}' tmp1 > tmp2
cat tmp1 tmp2 > dis_regions_1j8i

### 2jp1 r8-52 ordered
awk '{if ( ($1<7 || $1>51) && ($1<=92) ) print $1 " 1";}' native_2jp1_fullchain > tmp3
awk '{ print $1+93 " 1";}' tmp3 > tmp4
cat tmp3 tmp4 > dis_regions_2jp1 
