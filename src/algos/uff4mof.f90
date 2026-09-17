!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2023 Philipp Pracht
!
! crest is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! crest is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with crest.  If not, see <https://www.gnu.org/licenses/>.
!================================================================================!

module uff4mof_data
!*******************************************************************
!* UFF4MOF force field parameters, transcribed from the UFF4MOF_DATA
!* table of lammps_interface (P. Boyd):
!* https://github.com/peteboyd/lammps_interface/blob/master/lammps_interface/uff4mof.py
!*
!* Columns per atom type, as in the original table:
!*   r1, theta0, x1, D1, zeta, Z1, Vi, Uj, Xi, hard, radius
!*
!* This module only tabulates the parameters and provides a lookup
!* of which atom types exist for a given element (e.g. "Fe" ->
!* "Fe3+2","Fe6+2","Fe6+3","Fe4+2"). It does not yet feed into any
!* calculator setup.
!*******************************************************************
  use crest_parameters,only:wp
  implicit none
  private

  public :: uff4mof_entry
  public :: uff4mof_table
  public :: uff4mof_ntypes
  public :: uff4mof_element_symbol
  public :: uff4mof_options
  public :: uff4mof_label_cn
  public :: uff4mof_autodetect
  public :: uff4mof_get_entry
  public :: uff4mof_candidates_at_cn

  type :: uff4mof_entry
    character(len=8) :: label
    real(wp) :: r1,theta0,x1,D1,zeta,Z1,Vi,Uj,Xi,hard,radius
  end type uff4mof_entry

  integer,parameter :: uff4mof_ntypes = 222

!&<
  type(uff4mof_entry),parameter :: uff4mof_table(uff4mof_ntypes) = [ &
    uff4mof_entry('Du      ',0.01_wp,180_wp,0.4_wp,5000_wp,12_wp,10_wp,0_wp,0_wp,9.66_wp,14.92_wp,0.7_wp),&
    uff4mof_entry('H_      ',0.354_wp,180_wp,2.886_wp,0.044_wp,12_wp,0.712_wp,0_wp,0_wp,4.528_wp,6.9452_wp,0.371_wp),&
    uff4mof_entry('H_b     ',0.46_wp,83.5_wp,2.886_wp,0.044_wp,12_wp,0.712_wp,0_wp,0_wp,4.528_wp,6.9452_wp,0.371_wp),&
    uff4mof_entry('He4+4   ',0.849_wp,90_wp,2.362_wp,0.056_wp,15.24_wp,0.098_wp,0_wp,0_wp,9.66_wp,14.92_wp,1.3_wp),&
    uff4mof_entry('Li      ',1.336_wp,180_wp,2.451_wp,0.025_wp,12_wp,1.026_wp,0_wp,2_wp,3.006_wp,2.386_wp,1.557_wp),&
    uff4mof_entry('Be3+2   ',1.074_wp,109.47_wp,2.745_wp,0.085_wp,12_wp,1.565_wp,0_wp,2_wp,4.877_wp,4.443_wp,1.24_wp),&
    uff4mof_entry('B_3     ',0.838_wp,109.47_wp,4.083_wp,0.18_wp,12.052_wp,1.755_wp,0_wp,2_wp,5.11_wp,4.75_wp,0.822_wp),&
    uff4mof_entry('B_2     ',0.828_wp,120_wp,4.083_wp,0.18_wp,12.052_wp,1.755_wp,0_wp,2_wp,5.11_wp,4.75_wp,0.822_wp),&
    uff4mof_entry('C_3     ',0.757_wp,109.47_wp,3.851_wp,0.105_wp,12.73_wp,1.912_wp,2.119_wp,2_wp,5.343_wp,5.063_wp,0.759_wp),&
    uff4mof_entry('C_R     ',0.729_wp,120_wp,3.851_wp,0.105_wp,12.73_wp,1.912_wp,0_wp,2_wp,5.343_wp,5.063_wp,0.759_wp),&
    uff4mof_entry('C_2     ',0.732_wp,120_wp,3.851_wp,0.105_wp,12.73_wp,1.912_wp,0_wp,2_wp,5.343_wp,5.063_wp,0.759_wp),&
    uff4mof_entry('C_1     ',0.706_wp,180_wp,3.851_wp,0.105_wp,12.73_wp,1.912_wp,0_wp,2_wp,5.343_wp,5.063_wp,0.759_wp),&
    uff4mof_entry('N_3     ',0.7_wp,106.7_wp,3.66_wp,0.069_wp,13.407_wp,2.544_wp,0.45_wp,2_wp,6.899_wp,5.88_wp,0.715_wp),&
    uff4mof_entry('N_R     ',0.699_wp,120_wp,3.66_wp,0.069_wp,13.407_wp,2.544_wp,0_wp,2_wp,6.899_wp,5.88_wp,0.715_wp),&
    uff4mof_entry('N_2     ',0.685_wp,111.2_wp,3.66_wp,0.069_wp,13.407_wp,2.544_wp,0_wp,2_wp,6.899_wp,5.88_wp,0.715_wp),&
    uff4mof_entry('N_1     ',0.656_wp,180_wp,3.66_wp,0.069_wp,13.407_wp,2.544_wp,0_wp,2_wp,6.899_wp,5.88_wp,0.715_wp),&
    uff4mof_entry('O_3     ',0.658_wp,104.51_wp,3.5_wp,0.06_wp,14.085_wp,2.3_wp,0.018_wp,2_wp,8.741_wp,6.682_wp,0.669_wp),&
    uff4mof_entry('O_3_z   ',0.528_wp,146_wp,3.5_wp,0.06_wp,14.085_wp,2.3_wp,0.018_wp,2_wp,8.741_wp,6.682_wp,0.669_wp),&
    uff4mof_entry('O_3_M   ',0.658_wp,109.47_wp,3.5_wp,0.06_wp,14.085_wp,2.3_wp,0.018_wp,2_wp,8.741_wp,6.682_wp,0.669_wp),&
    uff4mof_entry('O_R     ',0.68_wp,110_wp,3.5_wp,0.06_wp,14.085_wp,2.3_wp,0_wp,2_wp,8.741_wp,6.682_wp,0.669_wp),&
    uff4mof_entry('O_2     ',0.634_wp,120_wp,3.5_wp,0.06_wp,14.085_wp,2.3_wp,0_wp,2_wp,8.741_wp,6.682_wp,0.669_wp),&
    uff4mof_entry('O_1     ',0.639_wp,180_wp,3.5_wp,0.06_wp,14.085_wp,2.3_wp,0_wp,2_wp,8.741_wp,6.682_wp,0.669_wp),&
    uff4mof_entry('F_      ',0.668_wp,180_wp,3.364_wp,0.05_wp,14.762_wp,1.735_wp,0_wp,2_wp,10.874_wp,7.474_wp,0.706_wp),&
    uff4mof_entry('Ne4+4   ',0.92_wp,90_wp,3.243_wp,0.042_wp,15.44_wp,0.194_wp,0_wp,2_wp,11.04_wp,10.55_wp,1.768_wp),&
    uff4mof_entry('Na      ',1.539_wp,180_wp,2.983_wp,0.03_wp,12_wp,1.081_wp,0_wp,1.25_wp,2.843_wp,2.296_wp,2.085_wp),&
    uff4mof_entry('Mg3+2   ',1.421_wp,109.47_wp,3.021_wp,0.111_wp,12_wp,1.787_wp,0_wp,1.25_wp,3.951_wp,3.693_wp,1.5_wp),&
    uff4mof_entry('Mg6     ',1.421_wp,90_wp,3.021_wp,0.111_wp,12_wp,1.787_wp,0_wp,1.25_wp,3.951_wp,3.693_wp,1.5_wp),&
    uff4mof_entry('Al3     ',1.244_wp,109.47_wp,4.499_wp,0.505_wp,11.278_wp,1.792_wp,0_wp,1.25_wp,4.06_wp,3.59_wp,1.201_wp),&
    uff4mof_entry('Si3     ',1.117_wp,109.47_wp,4.295_wp,0.402_wp,12.175_wp,2.323_wp,1.225_wp,1.25_wp,4.168_wp,3.487_wp,1.176_wp),&
    uff4mof_entry('P_3+3   ',1.101_wp,93.8_wp,4.147_wp,0.305_wp,13.072_wp,2.863_wp,2.4_wp,1.25_wp,5.463_wp,4_wp,1.102_wp),&
    uff4mof_entry('P_3+5   ',1.056_wp,109.47_wp,4.147_wp,0.305_wp,13.072_wp,2.863_wp,2.4_wp,1.25_wp,5.463_wp,4_wp,1.102_wp),&
    uff4mof_entry('P_3+q   ',1.056_wp,109.47_wp,4.147_wp,0.305_wp,13.072_wp,2.863_wp,2.4_wp,1.25_wp,5.463_wp,4_wp,1.102_wp),&
    uff4mof_entry('S_3+2   ',1.064_wp,92.1_wp,4.035_wp,0.274_wp,13.969_wp,2.703_wp,0.484_wp,1.25_wp,6.928_wp,4.486_wp,1.047_wp),&
    uff4mof_entry('S_3+4   ',1.049_wp,103.2_wp,4.035_wp,0.274_wp,13.969_wp,2.703_wp,0.484_wp,1.25_wp,6.928_wp,4.486_wp,1.047_wp),&
    uff4mof_entry('S_3+6   ',1.027_wp,109.47_wp,4.035_wp,0.274_wp,13.969_wp,2.703_wp,0.484_wp,1.25_wp,6.928_wp,4.486_wp,1.047_wp),&
    uff4mof_entry('S_R     ',1.077_wp,92.2_wp,4.035_wp,0.274_wp,13.969_wp,2.703_wp,0_wp,1.25_wp,6.928_wp,4.486_wp,1.047_wp),&
    uff4mof_entry('S_2     ',0.854_wp,120_wp,4.035_wp,0.274_wp,13.969_wp,2.703_wp,0_wp,1.25_wp,6.928_wp,4.486_wp,1.047_wp),&
    uff4mof_entry('Cl      ',1.044_wp,180_wp,3.947_wp,0.227_wp,14.866_wp,2.348_wp,0_wp,1.25_wp,8.564_wp,4.946_wp,0.994_wp),&
    uff4mof_entry('Ar4+4   ',1.032_wp,90_wp,3.868_wp,0.185_wp,15.763_wp,0.3_wp,0_wp,1.25_wp,9.465_wp,6.355_wp,2.108_wp),&
    uff4mof_entry('K_      ',1.953_wp,180_wp,3.812_wp,0.035_wp,12_wp,1.165_wp,0_wp,0.7_wp,2.421_wp,1.92_wp,2.586_wp),&
    uff4mof_entry('Ca6+2   ',1.761_wp,90_wp,3.399_wp,0.238_wp,12_wp,2.141_wp,0_wp,0.7_wp,3.231_wp,2.88_wp,2_wp),&
    uff4mof_entry('Sc3+3   ',1.513_wp,109.47_wp,3.295_wp,0.019_wp,12_wp,2.592_wp,0_wp,0.7_wp,3.395_wp,3.08_wp,1.75_wp),&
    uff4mof_entry('Ti3+4   ',1.412_wp,109.47_wp,3.175_wp,0.017_wp,12_wp,2.659_wp,0_wp,0.7_wp,3.47_wp,3.38_wp,1.607_wp),&
    uff4mof_entry('Ti6+4   ',1.412_wp,90_wp,3.175_wp,0.017_wp,12_wp,2.659_wp,0_wp,0.7_wp,3.47_wp,3.38_wp,1.607_wp),&
    uff4mof_entry('V_3+5   ',1.402_wp,109.47_wp,3.144_wp,0.016_wp,12_wp,2.679_wp,0_wp,0.7_wp,3.65_wp,3.41_wp,1.47_wp),&
    uff4mof_entry('Cr6+3   ',1.345_wp,90_wp,3.023_wp,0.015_wp,12_wp,2.463_wp,0_wp,0.7_wp,3.415_wp,3.865_wp,1.402_wp),&
    uff4mof_entry('Mn6+2   ',1.382_wp,90_wp,2.961_wp,0.013_wp,12_wp,2.43_wp,0_wp,0.7_wp,3.325_wp,4.105_wp,1.533_wp),&
    uff4mof_entry('Fe3+2   ',1.27_wp,109.47_wp,2.912_wp,0.013_wp,12_wp,2.43_wp,0_wp,0.7_wp,3.76_wp,4.14_wp,1.393_wp),&
    uff4mof_entry('Fe6+2   ',1.335_wp,90_wp,2.912_wp,0.013_wp,12_wp,2.43_wp,0_wp,0.7_wp,3.76_wp,4.14_wp,1.393_wp),&
    uff4mof_entry('Co6+3   ',1.241_wp,90_wp,2.872_wp,0.014_wp,12_wp,2.43_wp,0_wp,0.7_wp,4.105_wp,4.175_wp,1.406_wp),&
    uff4mof_entry('Ni4+2   ',1.164_wp,90_wp,2.834_wp,0.015_wp,12_wp,2.43_wp,0_wp,0.7_wp,4.465_wp,4.205_wp,1.398_wp),&
    uff4mof_entry('Cu3+1   ',1.302_wp,109.47_wp,3.495_wp,0.005_wp,12_wp,1.756_wp,0_wp,0.7_wp,4.2_wp,4.22_wp,1.434_wp),&
    uff4mof_entry('Zn3+2   ',1.193_wp,109.47_wp,2.763_wp,0.124_wp,12_wp,1.308_wp,0_wp,0.7_wp,5.106_wp,4.285_wp,1.4_wp),&
    uff4mof_entry('Ga3+3   ',1.26_wp,109.47_wp,4.383_wp,0.415_wp,11_wp,1.821_wp,0_wp,0.7_wp,3.641_wp,3.16_wp,1.211_wp),&
    uff4mof_entry('Ge3     ',1.197_wp,109.47_wp,4.28_wp,0.379_wp,12_wp,2.789_wp,0.701_wp,0.7_wp,4.051_wp,3.438_wp,1.189_wp),&
    uff4mof_entry('As3+3   ',1.211_wp,92.1_wp,4.23_wp,0.309_wp,13_wp,2.864_wp,1.5_wp,0.7_wp,5.188_wp,3.809_wp,1.204_wp),&
    uff4mof_entry('Se3+2   ',1.19_wp,90.6_wp,4.205_wp,0.291_wp,14_wp,2.764_wp,0.335_wp,0.7_wp,6.428_wp,4.131_wp,1.224_wp),&
    uff4mof_entry('Br      ',1.192_wp,180_wp,4.189_wp,0.251_wp,15_wp,2.519_wp,0_wp,0.7_wp,7.79_wp,4.425_wp,1.141_wp),&
    uff4mof_entry('Kr4+4   ',1.147_wp,90_wp,4.141_wp,0.22_wp,16_wp,0.452_wp,0_wp,0.7_wp,8.505_wp,5.715_wp,2.27_wp),&
    uff4mof_entry('Rb      ',2.26_wp,180_wp,4.114_wp,0.04_wp,12_wp,1.592_wp,0_wp,0.2_wp,2.331_wp,1.846_wp,2.77_wp),&
    uff4mof_entry('Sr6+2   ',2.052_wp,90_wp,3.641_wp,0.235_wp,12_wp,2.449_wp,0_wp,0.2_wp,3.024_wp,2.44_wp,2.415_wp),&
    uff4mof_entry('Y_3+3   ',1.698_wp,109.47_wp,3.345_wp,0.072_wp,12_wp,3.257_wp,0_wp,0.2_wp,3.83_wp,2.81_wp,1.998_wp),&
    uff4mof_entry('Zr3+4   ',1.564_wp,109.47_wp,3.124_wp,0.069_wp,12_wp,3.667_wp,0_wp,0.2_wp,3.4_wp,3.55_wp,1.758_wp),&
    uff4mof_entry('Nb3+5   ',1.473_wp,109.47_wp,3.165_wp,0.059_wp,12_wp,3.618_wp,0_wp,0.2_wp,3.55_wp,3.38_wp,1.603_wp),&
    uff4mof_entry('Mo6+6   ',1.467_wp,90_wp,3.052_wp,0.056_wp,12_wp,3.4_wp,0_wp,0.2_wp,3.465_wp,3.755_wp,1.53_wp),&
    uff4mof_entry('Mo3+6   ',1.484_wp,109.47_wp,3.052_wp,0.056_wp,12_wp,3.4_wp,0_wp,0.2_wp,3.465_wp,3.755_wp,1.53_wp),&
    uff4mof_entry('Tc6+5   ',1.322_wp,90_wp,2.998_wp,0.048_wp,12_wp,3.4_wp,0_wp,0.2_wp,3.29_wp,3.99_wp,1.5_wp),&
    uff4mof_entry('Ru6+2   ',1.478_wp,90_wp,2.963_wp,0.056_wp,12_wp,3.4_wp,0_wp,0.2_wp,3.575_wp,4.015_wp,1.5_wp),&
    uff4mof_entry('Rh6+3   ',1.332_wp,90_wp,2.929_wp,0.053_wp,12_wp,3.5_wp,0_wp,0.2_wp,3.975_wp,4.005_wp,1.509_wp),&
    uff4mof_entry('Pd4+2   ',1.338_wp,90_wp,2.899_wp,0.048_wp,12_wp,3.21_wp,0_wp,0.2_wp,4.32_wp,4_wp,1.544_wp),&
    uff4mof_entry('Ag1+1   ',1.386_wp,180_wp,3.148_wp,0.036_wp,12_wp,1.956_wp,0_wp,0.2_wp,4.436_wp,3.134_wp,1.622_wp),&
    uff4mof_entry('Cd3+2   ',1.403_wp,109.47_wp,2.848_wp,0.228_wp,12_wp,1.65_wp,0_wp,0.2_wp,5.034_wp,3.957_wp,1.6_wp),&
    uff4mof_entry('In3+3   ',1.459_wp,109.47_wp,4.463_wp,0.599_wp,11_wp,2.07_wp,0_wp,0.2_wp,3.506_wp,2.896_wp,1.404_wp),&
    uff4mof_entry('Sn3     ',1.398_wp,109.47_wp,4.392_wp,0.567_wp,12_wp,2.961_wp,0.199_wp,0.2_wp,3.987_wp,3.124_wp,1.354_wp),&
    uff4mof_entry('Sb3+3   ',1.407_wp,91.6_wp,4.42_wp,0.449_wp,13_wp,2.704_wp,1.1_wp,0.2_wp,4.899_wp,3.342_wp,1.404_wp),&
    uff4mof_entry('Te3+2   ',1.386_wp,90.25_wp,4.47_wp,0.398_wp,14_wp,2.882_wp,0.3_wp,0.2_wp,5.816_wp,3.526_wp,1.38_wp),&
    uff4mof_entry('I_      ',1.382_wp,180_wp,4.5_wp,0.339_wp,15_wp,2.65_wp,0_wp,0.2_wp,6.822_wp,3.762_wp,1.333_wp),&
    uff4mof_entry('Xe4+4   ',1.267_wp,90_wp,4.404_wp,0.332_wp,12_wp,0.556_wp,0_wp,0.2_wp,7.595_wp,4.975_wp,2.459_wp),&
    uff4mof_entry('Cs      ',2.57_wp,180_wp,4.517_wp,0.045_wp,12_wp,1.573_wp,0_wp,0.1_wp,2.183_wp,1.711_wp,2.984_wp),&
    uff4mof_entry('Ba6+2   ',2.277_wp,90_wp,3.703_wp,0.364_wp,12_wp,2.727_wp,0_wp,0.1_wp,2.814_wp,2.396_wp,2.442_wp),&
    uff4mof_entry('La3+3   ',1.943_wp,109.47_wp,3.522_wp,0.017_wp,12_wp,3.3_wp,0_wp,0.1_wp,2.8355_wp,2.7415_wp,2.071_wp),&
    uff4mof_entry('Ce6+3   ',1.841_wp,90_wp,3.556_wp,0.013_wp,12_wp,3.3_wp,0_wp,0.1_wp,2.774_wp,2.692_wp,1.925_wp),&
    uff4mof_entry('Pr6+3   ',1.823_wp,90_wp,3.606_wp,0.01_wp,12_wp,3.3_wp,0_wp,0.1_wp,2.858_wp,2.564_wp,2.007_wp),&
    uff4mof_entry('Nd6+3   ',1.816_wp,90_wp,3.575_wp,0.01_wp,12_wp,3.3_wp,0_wp,0.1_wp,2.8685_wp,2.6205_wp,2.007_wp),&
    uff4mof_entry('Pm6+3   ',1.801_wp,90_wp,3.547_wp,0.009_wp,12_wp,3.3_wp,0_wp,0.1_wp,2.881_wp,2.673_wp,2_wp),&
    uff4mof_entry('Sm6+3   ',1.78_wp,90_wp,3.52_wp,0.008_wp,12_wp,3.3_wp,0_wp,0.1_wp,2.9115_wp,2.7195_wp,1.978_wp),&
    uff4mof_entry('Eu6+3   ',1.771_wp,90_wp,3.493_wp,0.008_wp,12_wp,3.3_wp,0_wp,0.1_wp,2.8785_wp,2.7875_wp,2.227_wp),&
    uff4mof_entry('Gd6+3   ',1.735_wp,90_wp,3.368_wp,0.009_wp,12_wp,3.3_wp,0_wp,0.1_wp,3.1665_wp,2.9745_wp,1.968_wp),&
    uff4mof_entry('Tb6+3   ',1.732_wp,90_wp,3.451_wp,0.007_wp,12_wp,3.3_wp,0_wp,0.1_wp,3.018_wp,2.834_wp,1.954_wp),&
    uff4mof_entry('Dy6+3   ',1.71_wp,90_wp,3.428_wp,0.007_wp,12_wp,3.3_wp,0_wp,0.1_wp,3.0555_wp,2.8715_wp,1.934_wp),&
    uff4mof_entry('Ho6+3   ',1.696_wp,90_wp,3.409_wp,0.007_wp,12_wp,3.416_wp,0_wp,0.1_wp,3.127_wp,2.891_wp,1.925_wp),&
    uff4mof_entry('Er6+3   ',1.673_wp,90_wp,3.391_wp,0.007_wp,12_wp,3.3_wp,0_wp,0.1_wp,3.1865_wp,2.9145_wp,1.915_wp),&
    uff4mof_entry('Tm6+3   ',1.66_wp,90_wp,3.374_wp,0.006_wp,12_wp,3.3_wp,0_wp,0.1_wp,3.2514_wp,2.9329_wp,2_wp),&
    uff4mof_entry('Yb6+3   ',1.637_wp,90_wp,3.355_wp,0.228_wp,12_wp,2.618_wp,0_wp,0.1_wp,3.2889_wp,2.965_wp,2.158_wp),&
    uff4mof_entry('Lu6+3   ',1.671_wp,90_wp,3.64_wp,0.041_wp,12_wp,3.271_wp,0_wp,0.1_wp,2.9629_wp,2.4629_wp,1.896_wp),&
    uff4mof_entry('Hf3+4   ',1.611_wp,109.47_wp,3.141_wp,0.072_wp,12_wp,3.921_wp,0_wp,0.1_wp,3.7_wp,3.4_wp,1.759_wp),&
    uff4mof_entry('Ta3+5   ',1.511_wp,109.47_wp,3.17_wp,0.081_wp,12_wp,4.075_wp,0_wp,0.1_wp,5.1_wp,2.85_wp,1.605_wp),&
    uff4mof_entry('W_6+6   ',1.392_wp,90_wp,3.069_wp,0.067_wp,12_wp,3.7_wp,0_wp,0.1_wp,4.63_wp,3.31_wp,1.538_wp),&
    uff4mof_entry('W_3+4   ',1.526_wp,109.47_wp,3.069_wp,0.067_wp,12_wp,3.7_wp,0_wp,0.1_wp,4.63_wp,3.31_wp,1.538_wp),&
    uff4mof_entry('W_3+6   ',1.38_wp,109.47_wp,3.069_wp,0.067_wp,12_wp,3.7_wp,0_wp,0.1_wp,4.63_wp,3.31_wp,1.538_wp),&
    uff4mof_entry('Re6+5   ',1.372_wp,90_wp,2.954_wp,0.066_wp,12_wp,3.7_wp,0_wp,0.1_wp,3.96_wp,3.92_wp,1.6_wp),&
    uff4mof_entry('Re3+7   ',1.314_wp,109.47_wp,2.954_wp,0.066_wp,12_wp,3.7_wp,0_wp,0.1_wp,3.96_wp,3.92_wp,1.6_wp),&
    uff4mof_entry('Os6+6   ',1.372_wp,90_wp,3.12_wp,0.037_wp,12_wp,3.7_wp,0_wp,0.1_wp,5.14_wp,3.63_wp,1.7_wp),&
    uff4mof_entry('Ir6+3   ',1.371_wp,90_wp,2.84_wp,0.073_wp,12_wp,3.731_wp,0_wp,0.1_wp,5_wp,4_wp,1.866_wp),&
    uff4mof_entry('Pt4+2   ',1.125_wp,90_wp,2.754_wp,0.08_wp,12_wp,3.382_wp,0_wp,0.1_wp,4.79_wp,4.43_wp,1.557_wp),&
    uff4mof_entry('Au4+3   ',1.262_wp,90_wp,3.293_wp,0.039_wp,12_wp,2.625_wp,0_wp,0.1_wp,4.894_wp,2.586_wp,1.618_wp),&
    uff4mof_entry('Hg1+2   ',1.34_wp,180_wp,2.705_wp,0.385_wp,12_wp,1.75_wp,0_wp,0.1_wp,6.27_wp,4.16_wp,1.6_wp),&
    uff4mof_entry('Tl3+3   ',1.518_wp,120_wp,4.347_wp,0.68_wp,11_wp,2.068_wp,0_wp,0.1_wp,3.2_wp,2.9_wp,1.53_wp),&
    uff4mof_entry('Pb3     ',1.459_wp,109.47_wp,4.297_wp,0.663_wp,12_wp,2.846_wp,0.1_wp,0.1_wp,3.9_wp,3.53_wp,1.444_wp),&
    uff4mof_entry('Bi3+3   ',1.512_wp,90_wp,4.37_wp,0.518_wp,13_wp,2.47_wp,1_wp,0.1_wp,4.69_wp,3.74_wp,1.514_wp),&
    uff4mof_entry('Po3+2   ',1.5_wp,90_wp,4.709_wp,0.325_wp,14_wp,2.33_wp,0.3_wp,0.1_wp,4.21_wp,4.21_wp,1.48_wp),&
    uff4mof_entry('At      ',1.545_wp,180_wp,4.75_wp,0.284_wp,15_wp,2.24_wp,0_wp,0.1_wp,4.75_wp,4.75_wp,1.47_wp),&
    uff4mof_entry('Rn4+4   ',1.42_wp,90_wp,4.765_wp,0.248_wp,16_wp,0.583_wp,0_wp,0.1_wp,5.37_wp,5.37_wp,2.2_wp),&
    uff4mof_entry('Fr      ',2.88_wp,180_wp,4.9_wp,0.05_wp,12_wp,1.847_wp,0_wp,0_wp,2_wp,2_wp,2.3_wp),&
    uff4mof_entry('Ra6+2   ',2.512_wp,90_wp,3.677_wp,0.404_wp,12_wp,2.92_wp,0_wp,0_wp,2.843_wp,2.434_wp,2.2_wp),&
    uff4mof_entry('Ac6+3   ',1.983_wp,90_wp,3.478_wp,0.033_wp,12_wp,3.9_wp,0_wp,0_wp,2.835_wp,2.835_wp,2.108_wp),&
    uff4mof_entry('Th6+4   ',1.721_wp,90_wp,3.396_wp,0.026_wp,12_wp,4.202_wp,0_wp,0_wp,3.175_wp,2.905_wp,2.018_wp),&
    uff4mof_entry('Pa6+4   ',1.711_wp,90_wp,3.424_wp,0.022_wp,12_wp,3.9_wp,0_wp,0_wp,2.985_wp,2.905_wp,1.8_wp),&
    uff4mof_entry('U_6+4   ',1.684_wp,90_wp,3.395_wp,0.022_wp,12_wp,3.9_wp,0_wp,0_wp,3.341_wp,2.853_wp,1.713_wp),&
    uff4mof_entry('Np6+4   ',1.666_wp,90_wp,3.424_wp,0.019_wp,12_wp,3.9_wp,0_wp,0_wp,3.549_wp,2.717_wp,1.8_wp),&
    uff4mof_entry('Pu6+4   ',1.657_wp,90_wp,3.424_wp,0.016_wp,12_wp,3.9_wp,0_wp,0_wp,3.243_wp,2.819_wp,1.84_wp),&
    uff4mof_entry('Am6+4   ',1.66_wp,90_wp,3.381_wp,0.014_wp,12_wp,3.9_wp,0_wp,0_wp,2.9895_wp,3.0035_wp,1.942_wp),&
    uff4mof_entry('Cm6+3   ',1.801_wp,90_wp,3.326_wp,0.013_wp,12_wp,3.9_wp,0_wp,0_wp,2.8315_wp,3.1895_wp,1.9_wp),&
    uff4mof_entry('Bk6+3   ',1.761_wp,90_wp,3.339_wp,0.013_wp,12_wp,3.9_wp,0_wp,0_wp,3.1935_wp,3.0355_wp,1.9_wp),&
    uff4mof_entry('Cf6+3   ',1.75_wp,90_wp,3.313_wp,0.013_wp,12_wp,3.9_wp,0_wp,0_wp,3.197_wp,3.101_wp,1.9_wp),&
    uff4mof_entry('Es6+3   ',1.724_wp,90_wp,3.299_wp,0.012_wp,12_wp,3.9_wp,0_wp,0_wp,3.333_wp,3.089_wp,1.9_wp),&
    uff4mof_entry('Fm6+3   ',1.712_wp,90_wp,3.286_wp,0.012_wp,12_wp,3.9_wp,0_wp,0_wp,3.4_wp,3.1_wp,1.9_wp),&
    uff4mof_entry('Md6+3   ',1.689_wp,90_wp,3.274_wp,0.011_wp,12_wp,3.9_wp,0_wp,0_wp,3.47_wp,3.11_wp,1.9_wp),&
    uff4mof_entry('No6+3   ',1.679_wp,90_wp,3.248_wp,0.011_wp,12_wp,3.9_wp,0_wp,0_wp,3.475_wp,3.175_wp,1.9_wp),&
    uff4mof_entry('Lw6+3   ',1.698_wp,90_wp,3.236_wp,0.011_wp,12_wp,3.9_wp,0_wp,0_wp,3.5_wp,3.2_wp,1.9_wp),&
    uff4mof_entry('O_3_f   ',0.634_wp,109.47_wp,3.5_wp,0.06_wp,14.085_wp,2.3_wp,0.018_wp,2_wp,8.741_wp,6.682_wp,0.669_wp),&
    uff4mof_entry('O_2_z   ',0.528_wp,120.0_wp,3.5_wp,0.06_wp,14.085_wp,2.3_wp,0_wp,2_wp,8.741_wp,6.682_wp,0.669_wp),&
    uff4mof_entry('Al6+3   ',1.22_wp,90.0_wp,4.499_wp,0.505_wp,11.278_wp,1.792_wp,0_wp,1.25_wp,4.06_wp,3.59_wp,1.201_wp),&
    uff4mof_entry('Sc6+3   ',1.44_wp,90.0_wp,3.295_wp,0.019_wp,12.0_wp,2.595_wp,0_wp,0.7_wp,3.395_wp,3.08_wp,1.75_wp),&
    uff4mof_entry('Ti4+2   ',1.38_wp,90.0_wp,3.175_wp,0.017_wp,12.0_wp,2.659_wp,0_wp,0.7_wp,3.47_wp,3.38_wp,1.607_wp),&
    uff4mof_entry('V_4+2   ',1.18_wp,90.0_wp,3.144_wp,0.016_wp,12.0_wp,2.679_wp,0_wp,0.7_wp,3.65_wp,3.41_wp,1.47_wp),&
    uff4mof_entry('V_6+3   ',1.30_wp,90.0_wp,3.144_wp,0.016_wp,12.0_wp,2.679_wp,0_wp,0.7_wp,3.65_wp,3.41_wp,1.47_wp),&
    uff4mof_entry('Cr4+2   ',1.10_wp,90.0_wp,3.023_wp,0.015_wp,12.0_wp,2.463_wp,0_wp,0.7_wp,3.415_wp,3.865_wp,1.402_wp),&
    uff4mof_entry('Cr6f3   ',1.28_wp,90.0_wp,3.023_wp,0.015_wp,12.0_wp,2.463_wp,0_wp,0.7_wp,3.415_wp,3.865_wp,1.402_wp),&
    uff4mof_entry('Mn6+3   ',1.34_wp,90.0_wp,2.961_wp,0.013_wp,12.0_wp,2.43_wp,0_wp,0.7_wp,3.325_wp,4.105_wp,1.533_wp),&
    uff4mof_entry('Mn4+2   ',1.26_wp,90.0_wp,2.961_wp,0.013_wp,12.0_wp,2.43_wp,0_wp,0.7_wp,3.325_wp,4.105_wp,1.533_wp),&
    uff4mof_entry('Fe6+3   ',1.32_wp,90.0_wp,2.912_wp,0.013_wp,12.0_wp,2.43_wp,0_wp,0.7_wp,3.76_wp,4.14_wp,1.393_wp),&
    uff4mof_entry('Fe4+2   ',1.10_wp,90.0_wp,2.912_wp,0.013_wp,12.0_wp,2.43_wp,0_wp,0.7_wp,3.76_wp,4.14_wp,1.393_wp),&
    uff4mof_entry('Co3+2   ',1.24_wp,109.47_wp,2.872_wp,0.014_wp,12.0_wp,1.308_wp,0_wp,0.7_wp,4.105_wp,4.175_wp,1.406_wp),&
    uff4mof_entry('Co4+2   ',1.16_wp,90.0_wp,2.872_wp,0.014_wp,12.0_wp,1.308_wp,0_wp,0.7_wp,4.105_wp,4.175_wp,1.406_wp),&
    uff4mof_entry('Cu4+2   ',1.28_wp,90.0_wp,3.495_wp,0.005_wp,12.0_wp,2.43_wp,0_wp,0.7_wp,4.2_wp,4.22_wp,1.434_wp),&
    uff4mof_entry('Zn4+2   ',1.34_wp,90.0_wp,2.763_wp,0.124_wp,12.0_wp,1.308_wp,0_wp,0.7_wp,5.106_wp,4.285_wp,1.4_wp),&
    uff4mof_entry('Zn3f2   ',1.24_wp,109.47_wp,2.763_wp,0.124_wp,12.0_wp,1.308_wp,0_wp,0.7_wp,5.106_wp,4.285_wp,1.4_wp),&
    uff4mof_entry('Li3f2   ',1.28_wp,109.47_wp,2.451_wp,0.025_wp,12_wp,1.026_wp,0_wp,2_wp,3.006_wp,2.386_wp,1.557_wp),&
    uff4mof_entry('Na3f2   ',1.623_wp,109.47_wp,2.983_wp,0.03_wp,12_wp,1.081_wp,0_wp,1.25_wp,2.843_wp,2.296_wp,2.085_wp),&
    uff4mof_entry('Na4f2   ',1.79_wp,90.0_wp,2.983_wp,0.03_wp,12_wp,1.081_wp,0_wp,1.25_wp,2.843_wp,2.296_wp,2.085_wp),&
    uff4mof_entry('Mg6f3   ',1.525_wp,90.0_wp,3.021_wp,0.111_wp,12_wp,1.787_wp,0_wp,1.25_wp,3.951_wp,3.693_wp,1.5_wp),&
    uff4mof_entry('Al3f2   ',1.28_wp,109.47_wp,4.499_wp,0.505_wp,11.278_wp,1.792_wp,0_wp,1.25_wp,4.06_wp,3.59_wp,1.201_wp),&
    uff4mof_entry('K_3f2   ',2.38_wp,109.47_wp,3.812_wp,0.035_wp,12_wp,1.165_wp,0_wp,0.7_wp,2.421_wp,1.92_wp,2.586_wp),&
    uff4mof_entry('K_4f2   ',2.01_wp,90.0_wp,3.812_wp,0.035_wp,12_wp,1.165_wp,0_wp,0.7_wp,2.421_wp,1.92_wp,2.586_wp),&
    uff4mof_entry('Ca3f2   ',1.705_wp,109.47_wp,3.399_wp,0.238_wp,12_wp,2.141_wp,0_wp,0.7_wp,3.231_wp,2.88_wp,2_wp),&
    uff4mof_entry('V_3f2   ',1.12_wp,109.47_wp,3.144_wp,0.016_wp,12_wp,2.679_wp,0_wp,0.7_wp,3.65_wp,3.41_wp,1.47_wp),&
    uff4mof_entry('Mn1f1   ',1.38_wp,180.0_wp,2.961_wp,0.013_wp,12_wp,2.43_wp,0_wp,0.7_wp,3.325_wp,4.105_wp,1.533_wp),&
    uff4mof_entry('Mn3f2   ',1.18_wp,109.47_wp,2.961_wp,0.013_wp,12_wp,2.43_wp,0_wp,0.7_wp,3.325_wp,4.105_wp,1.533_wp),&
    uff4mof_entry('Mn8f4   ',1.52_wp,109.47_wp,2.961_wp,0.013_wp,12_wp,2.43_wp,0_wp,0.7_wp,3.325_wp,4.105_wp,1.533_wp),&
    uff4mof_entry('Co1f1   ',1.28_wp,180.0_wp,2.872_wp,0.014_wp,12_wp,2.43_wp,0_wp,0.7_wp,4.105_wp,4.175_wp,1.406_wp),&
    uff4mof_entry('Cu1f1   ',1.24_wp,180.0_wp,3.495_wp,0.005_wp,12_wp,1.756_wp,0_wp,0.7_wp,4.2_wp,4.22_wp,1.434_wp),&
    uff4mof_entry('Cu2f2   ',1.11_wp,120.0_wp,3.495_wp,0.005_wp,12_wp,1.756_wp,0_wp,0.7_wp,4.2_wp,4.22_wp,1.434_wp),&
    uff4mof_entry('Cu3f2   ',1.19_wp,109.47_wp,3.495_wp,0.005_wp,12_wp,1.756_wp,0_wp,0.7_wp,4.2_wp,4.22_wp,1.434_wp),&
    uff4mof_entry('Zn1f1   ',1.30_wp,180.0_wp,2.763_wp,0.124_wp,12_wp,1.308_wp,0_wp,0.7_wp,5.106_wp,4.285_wp,1.4_wp),&
    uff4mof_entry('Zn2f2   ',1.30_wp,120.0_wp,2.763_wp,0.124_wp,12_wp,1.308_wp,0_wp,0.7_wp,5.106_wp,4.285_wp,1.4_wp),&
    uff4mof_entry('Ga3f2   ',1.15_wp,109.47_wp,4.383_wp,0.415_wp,11_wp,1.821_wp,0_wp,0.7_wp,3.641_wp,3.16_wp,1.211_wp),&
    uff4mof_entry('Ga6f3   ',1.48_wp,90.0_wp,4.383_wp,0.415_wp,11_wp,1.821_wp,0_wp,0.7_wp,3.641_wp,3.16_wp,1.211_wp),&
    uff4mof_entry('Sr8f4   ',1.82_wp,109.47_wp,3.641_wp,0.235_wp,12_wp,2.449_wp,0_wp,0.2_wp,3.024_wp,2.44_wp,2.415_wp),&
    uff4mof_entry('Y_6f3   ',1.60_wp,90.0_wp,3.345_wp,0.072_wp,12_wp,3.257_wp,0_wp,0.2_wp,3.83_wp,2.81_wp,1.998_wp),&
    uff4mof_entry('Y_8f4   ',1.68_wp,109.47_wp,3.345_wp,0.072_wp,12_wp,3.257_wp,0_wp,0.2_wp,3.83_wp,2.81_wp,1.998_wp),&
    uff4mof_entry('Zr8f4   ',1.68_wp,109.47_wp,3.124_wp,0.069_wp,12_wp,3.667_wp,0_wp,0.2_wp,3.4_wp,3.55_wp,1.758_wp),&
    uff4mof_entry('Nb8f4   ',1.37_wp,109.47_wp,3.165_wp,0.059_wp,12_wp,3.618_wp,0_wp,0.2_wp,3.55_wp,3.38_wp,1.603_wp),&
    uff4mof_entry('Mo3f2   ',1.24_wp,109.47_wp,3.052_wp,0.056_wp,12_wp,3.4_wp,0_wp,0.2_wp,3.465_wp,3.755_wp,1.53_wp),&
    uff4mof_entry('Mo4f2   ',1.40_wp,90.0_wp,3.052_wp,0.056_wp,12_wp,3.4_wp,0_wp,0.2_wp,3.465_wp,3.755_wp,1.53_wp),&
    uff4mof_entry('Mo8f4   ',1.28_wp,109.47_wp,3.052_wp,0.056_wp,12_wp,3.4_wp,0_wp,0.2_wp,3.465_wp,3.755_wp,1.53_wp),&
    uff4mof_entry('Tc4f2   ',1.32_wp,90_wp,2.998_wp,0.048_wp,12_wp,3.4_wp,0_wp,0.2_wp,3.29_wp,3.99_wp,1.5_wp),&
    uff4mof_entry('Ru4f2   ',1.32_wp,90_wp,2.963_wp,0.056_wp,12_wp,3.4_wp,0_wp,0.2_wp,3.575_wp,4.015_wp,1.5_wp),&
    uff4mof_entry('Pd6f3   ',1.19_wp,90_wp,2.899_wp,0.048_wp,12_wp,3.21_wp,0_wp,0.2_wp,4.32_wp,4_wp,1.544_wp),&
    uff4mof_entry('Ag1f1   ',1.22_wp,180_wp,3.148_wp,0.036_wp,12_wp,1.956_wp,0_wp,0.2_wp,4.436_wp,3.134_wp,1.622_wp),&
    uff4mof_entry('Ag2f2   ',1.34_wp,120_wp,3.148_wp,0.036_wp,12_wp,1.956_wp,0_wp,0.2_wp,4.436_wp,3.134_wp,1.622_wp),&
    uff4mof_entry('Ag3f2   ',1.48_wp,109.47_wp,3.148_wp,0.036_wp,12_wp,1.956_wp,0_wp,0.2_wp,4.436_wp,3.134_wp,1.622_wp),&
    uff4mof_entry('Ag4f2   ',1.51_wp,90.0_wp,3.148_wp,0.036_wp,12_wp,1.956_wp,0_wp,0.2_wp,4.436_wp,3.134_wp,1.622_wp),&
    uff4mof_entry('Cd1f1   ',1.40_wp,180.0_wp,2.848_wp,0.228_wp,12_wp,1.65_wp,0_wp,0.2_wp,5.034_wp,3.957_wp,1.6_wp),&
    uff4mof_entry('Cd3f2   ',1.29_wp,109.47_wp,2.848_wp,0.228_wp,12_wp,1.65_wp,0_wp,0.2_wp,5.034_wp,3.957_wp,1.6_wp),&
    uff4mof_entry('Cd4f2   ',1.46_wp,90.0_wp,2.848_wp,0.228_wp,12_wp,1.65_wp,0_wp,0.2_wp,5.034_wp,3.957_wp,1.6_wp),&
    uff4mof_entry('Cd8f4   ',1.64_wp,109.47_wp,2.848_wp,0.228_wp,12_wp,1.65_wp,0_wp,0.2_wp,5.034_wp,3.957_wp,1.6_wp),&
    uff4mof_entry('In3f2   ',1.33_wp,109.47_wp,4.463_wp,0.599_wp,11_wp,2.07_wp,0_wp,0.2_wp,3.506_wp,2.896_wp,1.404_wp),&
    uff4mof_entry('In6f3   ',1.53_wp,90.0_wp,4.463_wp,0.599_wp,11_wp,2.07_wp,0_wp,0.2_wp,3.506_wp,2.896_wp,1.404_wp),&
    uff4mof_entry('In8f4   ',1.53_wp,109.47_wp,4.463_wp,0.599_wp,11_wp,2.07_wp,0_wp,0.2_wp,3.506_wp,2.896_wp,1.404_wp),&
    uff4mof_entry('Ba3f2   ',2.04_wp,109.47_wp,3.703_wp,0.364_wp,12_wp,2.727_wp,0_wp,0.1_wp,2.814_wp,2.396_wp,2.442_wp),&
    uff4mof_entry('La8f4   ',1.66_wp,109.47_wp,3.522_wp,0.017_wp,12_wp,3.3_wp,0_wp,0.1_wp,2.8355_wp,2.7415_wp,2.071_wp),&
    uff4mof_entry('Ce8f4   ',1.76_wp,109.47_wp,3.556_wp,0.013_wp,12_wp,3.3_wp,0_wp,0.1_wp,2.774_wp,2.692_wp,1.925_wp),&
    uff4mof_entry('Pr8f4   ',1.83_wp,109.47_wp,3.606_wp,0.01_wp,12_wp,3.3_wp,0_wp,0.1_wp,2.858_wp,2.564_wp,2.007_wp),&
    uff4mof_entry('Nd8f4   ',1.78_wp,109.47_wp,3.575_wp,0.01_wp,12_wp,3.3_wp,0_wp,0.1_wp,2.8685_wp,2.6205_wp,2.007_wp),&
    uff4mof_entry('Sm8f4   ',1.78_wp,109.47_wp,3.52_wp,0.008_wp,12_wp,3.3_wp,0_wp,0.1_wp,2.9115_wp,2.7195_wp,1.978_wp),&
    uff4mof_entry('Eu6f3   ',1.60_wp,90_wp,3.493_wp,0.008_wp,12_wp,3.3_wp,0_wp,0.1_wp,2.8785_wp,2.7875_wp,2.227_wp),&
    uff4mof_entry('Eu8f4   ',1.74_wp,109.47_wp,3.493_wp,0.008_wp,12_wp,3.3_wp,0_wp,0.1_wp,2.8785_wp,2.7875_wp,2.227_wp),&
    uff4mof_entry('Gd6f3   ',1.55_wp,90_wp,3.368_wp,0.009_wp,12_wp,3.3_wp,0_wp,0.1_wp,3.1665_wp,2.9745_wp,1.968_wp),&
    uff4mof_entry('Gd8f4   ',1.70_wp,109.47_wp,3.368_wp,0.009_wp,12_wp,3.3_wp,0_wp,0.1_wp,3.1665_wp,2.9745_wp,1.968_wp),&
    uff4mof_entry('Tb8f4   ',1.64_wp,109.47_wp,3.451_wp,0.007_wp,12_wp,3.3_wp,0_wp,0.1_wp,3.018_wp,2.834_wp,1.954_wp),&
    uff4mof_entry('Dy6f3   ',1.58_wp,90_wp,3.428_wp,0.007_wp,12_wp,3.3_wp,0_wp,0.1_wp,3.0555_wp,2.8715_wp,1.934_wp),&
    uff4mof_entry('Dy8f4   ',1.70_wp,109.47_wp,3.428_wp,0.007_wp,12_wp,3.3_wp,0_wp,0.1_wp,3.0555_wp,2.8715_wp,1.934_wp),&
    uff4mof_entry('Ho8f4   ',1.70_wp,109.47_wp,3.409_wp,0.007_wp,12_wp,3.416_wp,0_wp,0.1_wp,3.127_wp,2.891_wp,1.925_wp),&
    uff4mof_entry('Er8f4   ',1.64_wp,109.47_wp,3.391_wp,0.007_wp,12_wp,3.3_wp,0_wp,0.1_wp,3.1865_wp,2.9145_wp,1.915_wp),&
    uff4mof_entry('Tm8f4   ',1.67_wp,109.47_wp,3.374_wp,0.006_wp,12_wp,3.3_wp,0_wp,0.1_wp,3.2514_wp,2.9329_wp,2_wp),&
    uff4mof_entry('Yb6f3   ',1.45_wp,90_wp,3.355_wp,0.228_wp,12_wp,2.618_wp,0_wp,0.1_wp,3.2889_wp,2.965_wp,2.158_wp),&
    uff4mof_entry('Yb8f4   ',1.62_wp,109.47_wp,3.355_wp,0.228_wp,12_wp,2.618_wp,0_wp,0.1_wp,3.2889_wp,2.965_wp,2.158_wp),&
    uff4mof_entry('Lu8f4   ',1.66_wp,109.47_wp,3.64_wp,0.041_wp,12_wp,3.271_wp,0_wp,0.1_wp,2.9629_wp,2.4629_wp,1.896_wp),&
    uff4mof_entry('Hf8f4   ',1.46_wp,109.47_wp,3.141_wp,0.072_wp,12_wp,3.921_wp,0_wp,0.1_wp,3.7_wp,3.4_wp,1.759_wp),&
    uff4mof_entry('W_3f2   ',1.16_wp,109.47_wp,3.069_wp,0.067_wp,12_wp,3.7_wp,0_wp,0.1_wp,4.63_wp,3.31_wp,1.538_wp),&
    uff4mof_entry('W_4f2   ',1.345_wp,90_wp,3.069_wp,0.067_wp,12_wp,3.7_wp,0_wp,0.1_wp,4.63_wp,3.31_wp,1.538_wp),&
    uff4mof_entry('W_8f4   ',1.27_wp,109.47_wp,3.069_wp,0.067_wp,12_wp,3.7_wp,0_wp,0.1_wp,4.63_wp,3.31_wp,1.538_wp),&
    uff4mof_entry('Re6f3   ',1.23_wp,90.0_wp,2.954_wp,0.066_wp,12_wp,3.7_wp,0_wp,0.1_wp,3.96_wp,3.92_wp,1.6_wp),&
    uff4mof_entry('Os4f2   ',1.24_wp,90.0_wp,3.12_wp,0.037_wp,12_wp,3.7_wp,0_wp,0.1_wp,5.14_wp,3.63_wp,1.7_wp),&
    uff4mof_entry('Au1f1   ',1.11_wp,180_wp,3.293_wp,0.039_wp,12_wp,2.625_wp,0_wp,0.1_wp,4.894_wp,2.586_wp,1.618_wp),&
    uff4mof_entry('Hg3f2   ',1.248_wp,109.47_wp,2.705_wp,0.385_wp,12_wp,1.75_wp,0_wp,0.1_wp,6.27_wp,4.16_wp,1.6_wp),&
    uff4mof_entry('Pb4f2   ',1.67_wp,90.0_wp,4.297_wp,0.663_wp,12_wp,2.846_wp,0.1_wp,0.1_wp,3.9_wp,3.53_wp,1.444_wp),&
    uff4mof_entry('U_6f3   ',1.65_wp,90_wp,3.395_wp,0.022_wp,12_wp,3.9_wp,0_wp,0_wp,3.341_wp,2.853_wp,1.713_wp),&
    uff4mof_entry('U_8f4   ',1.73_wp,109.47_wp,3.395_wp,0.022_wp,12_wp,3.9_wp,0_wp,0_wp,3.341_wp,2.853_wp,1.713_wp),&
    uff4mof_entry('S_3_f   ',0.854_wp,109.47_wp,4.035_wp,0.274_wp,13.969_wp,2.703_wp,0.484_wp,1.25_wp,6.928_wp,4.486_wp,1.047_wp),&
!>--- CUSTOM entries below: NOT from lammps_interface/UFF4MOF -- added
!>    here because the original table has no octahedral (CN6) Zn type,
!>    which is needed for e.g. tris-chelate Zn(N,N-ligand)3 cage nodes.
!>    theta0=90 (idealized octahedral angle, as for the other 3d-metal
!>    "M6+n" entries e.g. Fe6+2/Cr6+3/Mn6+2). All non-geometric
!>    parameters (x1,D1,zeta,Z1,Vi,Uj,Xi,hard,radius) are IDENTICAL
!>    across every tabulated Zn type, so they are kept unchanged here.
!>    r1 is estimated by applying the same tetrahedral->octahedral
!>    r1 shift observed for Fe (Fe3+2 1.27 -> Fe6+2 1.335, +0.065 Å)
!>    to Zn3+2 (1.193 Å), i.e. NOT a fitted/validated literature value
!>    -- treat as provisional until refit against reference data.
    uff4mof_entry('Zn6+2   ',1.26_wp,90.0_wp,2.763_wp,0.124_wp,12.0_wp,1.308_wp,0_wp,0.7_wp,5.106_wp,4.285_wp,1.4_wp) &
  ]

!&>

contains

!========================================================================================!
!> Extract the leading element symbol from a UFF4MOF label, e.g.
!> "Fe3+2" -> "Fe", "O_3_M" -> "O", "W_6+6" -> "W".
!========================================================================================!
  function uff4mof_element_symbol(label) result(elem)
    character(len=*),intent(in) :: label
    character(len=2) :: elem
    integer :: i,n
    elem = ' '
    n = 0
    do i = 1,len_trim(label)
      if ((label(i:i) >= 'A'.and.label(i:i) <= 'Z').or. &
      &   (label(i:i) >= 'a'.and.label(i:i) <= 'z')) then
        n = n+1
        if (n <= 2) elem(n:n) = label(i:i)
      else
        exit
      end if
    end do
  end function uff4mof_element_symbol

!========================================================================================!
!> Collect all UFF4MOF atom-type labels available for a given element
!> symbol (e.g. "Fe"). Returns nopts=0 and an unallocated array if
!> none are tabulated for that element.
!========================================================================================!
  subroutine uff4mof_options(elementsymbol,labels,nopts)
    character(len=*),intent(in) :: elementsymbol
    character(len=8),allocatable,intent(out) :: labels(:)
    integer,intent(out) :: nopts
    character(len=2) :: query
    integer :: i
    query = adjustl(elementsymbol)
    nopts = 0
    do i = 1,uff4mof_ntypes
      if (trim(uff4mof_element_symbol(uff4mof_table(i)%label)) == trim(query)) then
        nopts = nopts+1
      end if
    end do
    if (nopts == 0) return
    allocate (labels(nopts))
    nopts = 0
    do i = 1,uff4mof_ntypes
      if (trim(uff4mof_element_symbol(uff4mof_table(i)%label)) == trim(query)) then
        nopts = nopts+1
        labels(nopts) = uff4mof_table(i)%label
      end if
    end do
  end subroutine uff4mof_options

!========================================================================================!
!> Coordination number implied by a UFF4MOF label's geometry digit, e.g.
!> "Fe3+2" -> 4 (digit 3 = UFF's tetrahedral code, itself CN4),
!> "Fe4+2" -> 4 (digit 4 = UFF's square-planar code, also CN4),
!> "Fe6+2" -> 6 (digit 6 = UFF's octahedral code, CN6),
!> "Mn8f4" -> 8 (UFF4MOF's extension beyond standard UFF: digit IS the CN).
!> Returns -1 if the label has no geometry digit (e.g. "Na","Cl","Br").
!========================================================================================!
  function uff4mof_label_cn(label) result(cn)
    character(len=*),intent(in) :: label
    integer :: cn
    character(len=2) :: elem
    integer :: i,n,digit
    integer,parameter :: digit2cn(8) = [2,3,4,4,5,6,7,8]
    cn = -1
    elem = uff4mof_element_symbol(label)
    n = len_trim(elem)
    if (len_trim(label) <= n) return !> nothing after the element symbol
    i = n+1
    if (label(i:i) == '_') i = i+1 !> skip the 1-letter-element padding underscore
    if (i > len_trim(label)) return
    digit = -1
    if (label(i:i) >= '1'.and.label(i:i) <= '8') then
      read (label(i:i),'(i1)') digit
    end if
    if (digit >= 1.and.digit <= 8) cn = digit2cn(digit)
  end function uff4mof_label_cn

!========================================================================================!
!> Try to auto-detect which tabulated UFF4MOF type fits a metal atom
!> based only on its (integer) coordination number: if exactly one
!> tabulated type for this element has a matching implied CN, that one
!> is returned with found=.true.; if zero or more than one match (i.e.
!> the geometry is ambiguous from CN alone, e.g. tetrahedral vs.
!> square-planar at CN4), found=.false. and the caller should fall back
!> to asking the user (see [[calculation.nci_metal]] in .toml input).
!========================================================================================!
  subroutine uff4mof_autodetect(elementsymbol,cn,label,found)
    character(len=*),intent(in) :: elementsymbol
    integer,intent(in) :: cn
    character(len=8),intent(out) :: label
    logical,intent(out) :: found
    character(len=8),allocatable :: opts(:)
    integer :: nopts,i,nmatch
    found = .false.
    label = ''
    call uff4mof_options(elementsymbol,opts,nopts)
    if (nopts == 0) return
    nmatch = 0
    do i = 1,nopts
      if (uff4mof_label_cn(opts(i)) == cn) then
        nmatch = nmatch+1
        label = opts(i)
      end if
    end do
    found = (nmatch == 1)
    if (.not.found) label = ''
  end subroutine uff4mof_autodetect

!========================================================================================!
!> All tabulated UFF4MOF labels for an element whose implied CN matches
!> `cn` (i.e. the full tied set uff4mof_autodetect gives up on when
!> nmatch>1). Used to let a geometric tiebreak (comparing actual donor
!> angles against each candidate's theta0) pick among them -- see
!> auto=true handling in search_newnci_metal.f90.
!========================================================================================!
  subroutine uff4mof_candidates_at_cn(elementsymbol,cn,labels,nlabels)
    character(len=*),intent(in) :: elementsymbol
    integer,intent(in) :: cn
    character(len=8),allocatable,intent(out) :: labels(:)
    integer,intent(out) :: nlabels
    character(len=8),allocatable :: opts(:)
    integer :: nopts,i
    nlabels = 0
    call uff4mof_options(elementsymbol,opts,nopts)
    if (nopts == 0) return
    do i = 1,nopts
      if (uff4mof_label_cn(opts(i)) == cn) nlabels = nlabels+1
    end do
    if (nlabels == 0) return
    allocate (labels(nlabels))
    nlabels = 0
    do i = 1,nopts
      if (uff4mof_label_cn(opts(i)) == cn) then
        nlabels = nlabels+1
        labels(nlabels) = opts(i)
      end if
    end do
  end subroutine uff4mof_candidates_at_cn

!========================================================================================!
!> Look up the tabulated bond radius (r1) and ideal coordination angle
!> (theta0) for a resolved UFF4MOF label, e.g. "Zn6+2" -> r1=1.26,
!> theta0=90.0. Used to build metal-ligand distance/angle restraints
!> (see nci_metal_bonds in search_newnci_metal.f90) -- r1 is the UFF
!> "natural bond radius" (Angstrom, uncorrected for bond order/EN), NOT
!> the inflated CN-counting radius crest's own cn_to_bond uses.
!========================================================================================!
  subroutine uff4mof_get_entry(label,r1,theta0,found)
    character(len=*),intent(in) :: label
    real(wp),intent(out) :: r1,theta0
    logical,intent(out) :: found
    integer :: i
    r1 = 0.0_wp
    theta0 = 0.0_wp
    found = .false.
    do i = 1,uff4mof_ntypes
      if (trim(uff4mof_table(i)%label) == trim(adjustl(label))) then
        r1 = uff4mof_table(i)%r1
        theta0 = uff4mof_table(i)%theta0
        found = .true.
        return
      end if
    end do
  end subroutine uff4mof_get_entry

!========================================================================================!
!========================================================================================!
end module uff4mof_data
