!
! Copyright (C) 2004-2024, DWD
! See ./AUTHORS.txt for a list of authors
! See ./LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
!

  TYPE t_CE_boehm
    REAL(wp) :: &
             &  r,        & 
             &  r_e,      & 
             &  v,        & 
             &  m,        & 
             &  alpha,    &
             &  alpha_e,  & 
             &  q_e,      &
             &  k,        & 
             &  DELS,    &
             &  CDI
    REAL(wp) :: SIG = 0.25_wp
    INTEGER  :: j   = 2
  END TYPE t_CE_boehm
