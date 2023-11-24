!-------------------------------------------------------------------------------------------
!
!
!
!
!
!
subroutine THERMO (model, NC, phase_type, calc, guess, T, P, n, Z, lnPhi, DLnPhi_N, DLnPhi_T, DLnPhi_P, H, DH_N, Cp, B, ic)
!
!
!  calc       H        CP        HX        FUG     DFUGT/P    DFUGX      B
!-----------------------------------------------------------------------------
!   0         X
!   1                                       X
!   2         X                             X
!   3         X         X
!   4                                       X         X
!   5                                       X                   X
!   6                                       X         X         X
!   7         X         X                   X         X
!   8         X         X         X         X         X
!   9         X         X         X         X         X         X        X
!  10                                                                    X
!
  implicit none
  
  integer, parameter          :: NCM = 30

  integer                     :: calc, phase_type, ic, guess, model, NC
  real(8)                     :: B, Cp, H, P, T, Z, Z0
  real(8), dimension(NCM)     :: DH_N, DLnPhi_P, DLnPhi_T, lnPhi, n
  real(8), dimension(NCM,NCM) :: DLnPhi_N
  
  integer                     :: common_model
  
  COMMON /MODEL/ common_model
  
  common_model = model
  B = 0; Cp = 0; H = 0
  DH_n(:NC) = 0; lnPhi(:NC) = 0; DLnPhi_N(:NC,:NC) = 0; DLnPhi_P = 0; DLnPhi_T = 0
  if (model == 1) then

!  Call to GCA-EoS
    call GCA_EoS(calc, phase_type, guess, NC, T, P, n, lnPhi, DLnPhi_N, DLnPhi_T, DLnPhi_P, H, Cp, DH_n, Z, B, ic)
  
  elseif (model == 2 .OR. model == 3) then

!  Call to SRK or PR EoS
    call SRKPR_EoS (calc, phase_type, NC, T, P, n, lnPhi, DLnPhi_N, DLnPhi_T, DLnPhi_P, H, Cp, DH_n, Z, B, ic)

  elseif (model == 4 .OR. model == 100) then

    Z0 = Z
    call RKPR_EoS (calc, phase_type, guess, NC, T, P, n, Z0, lnPhi, DLnPhi_N, DLnPhi_T, DLnPhi_P, H, Cp, DH_n, Z, B, ic)
    
  elseif (model == 5) then

    Z0 = Z
    call PCSAFT_EoS (calc, phase_type, guess, NC, T, P, n, Z0, lnPhi, DLnPhi_N, DLnPhi_T, DLnPhi_P, H, Cp, DH_n, Z, B, ic)    

  elseif (model == 6) then

    Z0 = Z
    call SanchezLacombe_EoS (calc, phase_type, guess, NC, T, P, n, Z0, lnPhi, DLnPhi_N, DLnPhi_T, DLnPhi_P, H, Cp, DH_n, Z, B, ic)  

  elseif (model == 7) then

    Z0 = Z
    call MTC_EoS (calc, phase_type, guess, NC, T, P, n, Z0, lnPhi, DLnPhi_N, DLnPhi_T, DLnPhi_P, H, Cp, DH_n, Z, B, ic)  
    
  elseif (model <= 0) then
    
    call gamma_phi (model, calc, phase_type, NC, T, P, n, lnPhi, DLnPhi_N, DLnPhi_T, DLnPhi_P, H, Cp, DH_n, Z, B, IC)
    
  endif

  return
endsubroutine THERMO
!--------------------------------------------------------------------------------
!
!
!
subroutine ReadParameters (inputFile, outputFile, model, NC)

  implicit none
  
  integer, parameter      :: NCM = 30
  
  integer                 :: inputFile, outputFile, model, NC
  integer                 :: NG, NGA, NST
  integer                 :: iCEq = 0, NONHY = 0
  character(20)           :: thermo_name, version

  common /GCA/               NG, NGA, NST
  common /cubicEq/           iCEq, NONHY
  
  if (model == 1) then

!    Call to GCA-EoS
    call GCA_version (version, thermo_name)
    write (outputFile, 1) thermo_name, version
    read (inputFile, *) NG, NGA
    call pargr (NC, NG, NGA, NST, inputFile, outputFile)
    call parmol (NC, NG, NST, inputFile, outputFile)
  
  elseif (model == 2 .OR. model == 3) then

    iCEq = model - 2
    call SRKPR_version (iCEq, version, thermo_name)
    write (outputFile, 1) thermo_name, version    
!    Call to cubic EoS, which only needs Tc, Pc, omega and k
    call indat (iCEq, NC, inputFile, outputFile)

  elseif (model == 4) then

!    RK-PR EOS
    call RKPR_version (version, thermo_name)
    write (outputFile, 1) thermo_name, version
    call readcomp (NC, inputFile, outputFile)

  elseif (model == 5) then
    
!     PC-SAFT EOS
    call PCSAFT_version (2, version, thermo_name)
    write (outputFile, 1) thermo_name, version        
    call readPCSAFT (inputFile, outputFile, NC)

  elseif (model == 6) then
    
!     SL EOS
    call SL_version (version, thermo_name)
    write (outputFile, 1) thermo_name, version        
    call read_SL_parameters (inputFile, outputFile, NC)

  elseif (model == 7) then
    
!     MTC EOS
    call MTC_version (version, thermo_name)
    write (outputFile, 1) thermo_name, version        
    call read_MTC_parameters (inputFile, outputFile, NC)

  elseif (model == 100) then

    call read_hs_mix_parameters (inputFile, outputFile, NC)
    
  elseif (model <= 0) then
    
    call readParametersGammaPhi (model, NC, inputFile, outputFile)
    thermo_name = "A                   "
    version = "A                   "
    
  endif
  write (outputFile, 1) thermo_name, version  
1  format (/, " Thermodynamic model: ", A, /," Package version: ", A, /)
  return

endsubroutine ReadParameters
!--------------------------------------------------------------------------------
!
!
!
!  This subroutine calculates thermodynamic properties using the GCA-EoS model
!
! 
SUBROUTINE GCA_EoS (calc, phase_type, guess, NC, T, P, n, lnPhi, DLnPhi_N, DLnPhi_T, DLnPhi_P, H, Cp, DH_n, Z, B, ic)
!
  implicit none
  integer, parameter           :: NCM = 30
  real(8), parameter           :: R = 82.05D0
!
  integer                      :: calc, guess, n_deriv = 1, NC, phase_type, temp_deriv = 0, ic
  integer                      :: NG, NGA, NST
  
  real(8)                      :: B, Cp, H, HR, P, T, Z
  real(8)                      :: Bvl, Bat, Bas
  real(8), dimension(NCM)      :: DH_N, DLnPhi_P, DLnPhi_T, Hgi, lnPhi, n, x
  real(8), dimension(NCM,NCM)  :: DLnPhi_N

  common /GCA/                    NG, NGA, NST
!
  if ((calc /= 0) .OR. (calc /= 3) .OR. calc /= 10) then
    
!  C�lculo de coeficiente de fugacidad
    if((calc == 5) .OR. (calc == 6) .OR. (calc == 9))n_deriv =  2

    if((calc == 4) .OR. (calc == 6) .OR. (calc == 7) .OR. (calc == 8) .OR. (calc == 9))temp_deriv = 1

    call PARAGC (T, NC, NG, NST, temp_deriv)
    call GCEOS (NC, NG, NST, n_deriv, temp_deriv, T, P, n, lnPhi, DLnPhi_N, DLnPhi_T, DLnPhi_P, Z, guess, phase_type, IC)

    if((calc == 0) .OR. (calc == 2) .OR. (calc == 3) .OR. (calc == 7) .OR. (calc == 8) .OR. (calc == 9))then

!  C�lculo de entalp�as molares parciales
      DH_n(:NC) = -R*T**2*DLnPhi_T(:NC)

    endif

  endif
  if((calc == 0) .OR. (calc == 2) .OR. (calc == 3) .OR. (calc == 7) .OR. (calc == 8) .OR. (calc == 9))then

    call PARAGC (T, NC, NG, NST, 1)
    call NTlpy (NC, NG, NST, T, P, n, HR, Z, guess, phase_type, iC)
    call NTlpy_GI (NC, T, Hgi(:NC))
    H = (HR + dot_product(n(:NC),Hgi(:NC)))/sum(n(:NC))

    if((calc == 8) .OR. (calc == 9))DH_n(:NC) = DH_n(:NC) + Hgi(:NC)
!    if(calc == 3)then
!
!      call Cp_GCA() !algo as�� como una sobrutina para el Cp...
!
!    endif

  endif
  
  if(calc == 9 .OR. calc == 10)then
    
!  Virial coefficient calculation
    call GCAvirial (NC, NG, NST, T, n, x, B, Bvl, Bat, Bas, .FALSE.)
    
  endif

  return
endsubroutine
!--------------------------------------------------------------------------------
!
!
!
!  This subroutine calculates thermodynamic properties using the the SRK or PR models.
!
! 
SUBROUTINE SRKPR_EoS (calc, phase_type, NC, T, P, n, lnPhi, DLnPhi_N, DLnPhi_T, DLnPhi_P, H, Cp, DH_n, Z, B, ic)
!
  implicit none
  integer, parameter           :: NCM = 30
  real(8), parameter           :: R = 82.05D0
!
  integer                      :: calc, ic, n_deriv, NC, phase_type, temp_deriv = 0
  integer                      :: iCEq, NC1, NONHY, NH, NST
  
  real(8)                      :: B, Cp, H, HR = 0, P, ResProp(6), T, Z
  real(8), dimension(NCM)      :: DH_N, DLnPhi_P, DLnPhi_T, Hgi, lnPhi, n
  real(8), dimension(NCM,NCM)  :: DLnPhi_N

  common /cubicEq/                iCEq, NONHY
  COMMON /STYR/                   NC1, NH, n_deriv, temp_deriv, NST
!
  n_deriv = 1
  temp_deriv = 0
  if((calc /= 10))then
    
    if((calc == 5) .OR. (calc == 6) .OR. (calc == 9))n_deriv =  2

    if((calc /= 1) .AND. (calc /= 5))temp_deriv = 1

    call CUBGEN (phase_type, IC, T, P, Z, n, lnPhi, DLnPhi_T, DLnPhi_P, DLnPhi_N, ResProp)
    ResProp = ResProp*R

    if(calc == 8 .OR. calc == 9)then

!  C�lculo de entalp�as molares parciales
      DH_n(:NC) = -R*T**2*DLnPhi_T(:NC)

    endif

  endif
  if((calc == 0) .OR. (calc == 2) .OR. (calc == 3) .OR. (calc == 7) .OR. (calc == 8) .OR. (calc == 9))then
! 
!  Add ideal gas contribution to enthalpy and heat capacity
    call NTlpy_GI (NC, T, Hgi(:NC))
    H = (ResProp(3) + dot_product(n(:NC),Hgi(:NC)))/sum(n(:NC))

    if((calc == 8) .OR. (calc == 9))DH_n(:NC) = DH_n(:NC) + Hgi(:NC)
    if(calc == 3)then

      call Cp_GI(NC, T, n, Hgi, Cp)
      Cp = ResProp(5) + Cp

    endif

  endif
  
  if(calc == 9 .OR. calc == 10)then
    
!  Virial coefficient calculation
    call CubicVirial(NC, n, T, B, h)
    
  endif

  return
endsubroutine
!--------------------------------------------------------------------------------
!
!
!
!  This subroutine calculates thermodynamic properties using the the RK-PR EOS.
!
! 
SUBROUTINE RKPR_EoS (calc, phase_type, guess, NC, T, P, n, Z0, lnPhi, DLnPhi_N, DLnPhi_T, DLnPhi_P, H, Cp, DH_n, Z, B, ic)
!
  implicit none
  integer, parameter           :: NCM = 30
  real(8), parameter           :: R = 82.05D0
!
  integer                      :: calc, guess, ic, n_deriv, NC, phase_type, temp_deriv = 0
  integer                      :: iCEq, NC1, NONHY, NH, NST
  
  real(8)                      :: B, Cp, Cpr = 0, H, HR = 0, P, Sr = 0, T, Ur = 0, V, V0, Z, Z0
  real(8), dimension(NCM)      :: DH_N, DLnPhi_P, DLnPhi_T, Hgi, lnPhi, n
  real(8), dimension(NCM,NCM)  :: DLnPhi_N

  
  Cp = 0; B = 0
  n_deriv = 1
  temp_deriv = 0
  if (calc /= 10)  then
    
    if ((calc == 5) .OR. (calc == 6) .OR. (calc == 9)) n_deriv =  2

    if ((calc /= 1) .AND. (calc /= 5)) temp_deriv = 1

    V0 = sum(n(:NC))*Z0*R*T/P
    call termo (phase_type, n_deriv, temp_deriv, guess, IC, NC, T, P, n, V0, V, lnPhi, DLnPhi_T, DLnPhi_P, DLnPhi_N, Hr, Sr, Cpr)

    if (calc == 8 .OR. calc == 9) then

!  C�lculo de entalp�as molares parciales
      DH_n(:NC) = -R*T**2*DLnPhi_T(:NC)

    endif

  endif
  if((calc == 0) .OR. (calc == 2) .OR. (calc == 3) .OR. (calc == 7) .OR. (calc == 8) .OR. (calc == 9))then
! 
!  Add ideal gas contribution to enthalpy and heat capacity
    call NTlpy_GI (NC, T, Hgi(:NC))
    H = (Hr + dot_product(n(:NC), Hgi(:NC)))/sum(n(:NC))

    if (calc == 8 .OR. calc == 9) DH_n(:NC) = DH_n(:NC) + Hgi(:NC)
    if(calc == 3)then

      call Cp_GI(NC, T, n, Hgi, Cp)
      Cp = cPr + Cp

    endif

  endif
  
!   if(calc == 9 .OR. calc == 10)then
!     
! !  Virial coefficient calculation
! !     call CubicVirial(NC, n, T, B, h)
!     
!   endif
  
  Z = P*V/R/T/sum(n)
  
  return
endsubroutine
!--------------------------------------------------------------------------------
!
!
!
!  This subroutine calculates thermodynamic properties using the the RK-PR EOS.
!
! 
SUBROUTINE PCSAFT_EoS (calc, phase_type, guess, NC, T, P, n, Z0, lnPhi, DLnPhi_N, DLnPhi_T, DLnPhi_P, H, Cp, DH_n, Z, B, ic)
!
  implicit none
  integer, parameter           :: NCM = 30
  real(8), parameter           :: R = 82.05D0
!
  integer                      :: calc, guess, ic, n_deriv, NC, phase_type, temp_deriv = 0
  integer                      :: iCEq, NC1, NONHY, NH, NST
  
  real(8)                      :: B, Cp, Cpr = 0, H, HR = 0, P, Sr = 0, T, Ur = 0, V, V0, Z, Z0
  real(8), dimension(NCM)      :: DH_N, DLnPhi_P, DLnPhi_T, Hgi, lnPhi, n
  real(8), dimension(NCM,NCM)  :: DLnPhi_N

  
  Cp = 0; B = 0
  n_deriv = 1
  temp_deriv = 0
  if (calc /= 10)  then
    
    if ((calc == 5) .OR. (calc == 6) .OR. (calc == 9)) n_deriv =  2

    if ((calc /= 1) .AND. (calc /= 5)) temp_deriv = 1

    V0 = sum(n(:NC))*Z0*R*T/P
    call termo (phase_type, n_deriv, temp_deriv, guess, IC, NC, T, P, n, V0, V, lnPhi, DLnPhi_T, DLnPhi_P, DLnPhi_N, Hr, Sr, Cpr)

    if (calc == 8 .OR. calc == 9) then

!  C�lculo de entalp�as molares parciales
      DH_n(:NC) = -R*T**2*DLnPhi_T(:NC)

    endif

  endif
  if((calc == 0) .OR. (calc == 2) .OR. (calc == 3) .OR. (calc == 7) .OR. (calc == 8) .OR. (calc == 9))then
! 
!  Add ideal gas contribution to enthalpy and heat capacity
    call NTlpy_GI (NC, T, Hgi(:NC))
    H = (Hr + dot_product(n(:NC), Hgi(:NC)))/sum(n(:NC))

    if (calc == 8 .OR. calc == 9) DH_n(:NC) = DH_n(:NC) + Hgi(:NC)
    if(calc == 3)then

      call Cp_GI(NC, T, n, Hgi, Cp)
      Cp = cPr + Cp

    endif

  endif
  
!   if(calc == 9 .OR. calc == 10)then
!     
! !  Virial coefficient calculation
! !     call CubicVirial(NC, n, T, B, h)
!     
!   endif
  
  Z = P*V/R/T/sum(n)
  
  return
endsubroutine
!--------------------------------------------------------------------------------------------------------------------------------------
!
!
!
SUBROUTINE SanchezLacombe_EoS (calc, phase_type, guess, NC, T, P, n, Z0, lnPhi, DLnPhi_N, DLnPhi_T, DLnPhi_P, H, Cp, DH_n, Z, B, ic)
!
  implicit none
  integer, parameter           :: NCM = 30
  real(8), parameter           :: R = 82.05D0
!
  integer                      :: calc, guess, ic, n_deriv, NC, phase_type, temp_deriv = 0
  integer                      :: iCEq, NC1, NONHY, NH, NST
  
  real(8)                      :: B, Cp, Cpr = 0, H, HR = 0, P, Sr = 0, T, Ur = 0, V, V0, Z, Z0
  real(8), dimension(NCM)      :: DH_N, DLnPhi_P, DLnPhi_T, Hgi, lnPhi, n
  real(8), dimension(NCM,NCM)  :: DLnPhi_N

  
  Cp = 0; B = 0
  n_deriv = 1
  temp_deriv = 0
  if (calc /= 10)  then
    
    if ((calc == 5) .OR. (calc == 6) .OR. (calc == 9)) n_deriv =  2

    if ((calc /= 1) .AND. (calc /= 5)) temp_deriv = 1

    V0 = sum(n(:NC))*Z0*R*T/P
    call termo (phase_type, n_deriv, temp_deriv, guess, IC, NC, T, P, n, V0, V, lnPhi, DLnPhi_T, DLnPhi_P, DLnPhi_N, Hr, Sr, Cpr)

    if (calc == 8 .OR. calc == 9) then

!  C�lculo de entalp�as molares parciales
      DH_n(:NC) = -R*T**2*DLnPhi_T(:NC)

    endif

  endif
  if((calc == 0) .OR. (calc == 2) .OR. (calc == 3) .OR. (calc == 7) .OR. (calc == 8) .OR. (calc == 9))then
! 
!  Add ideal gas contribution to enthalpy and heat capacity
    call NTlpy_GI (NC, T, Hgi(:NC))
    H = (Hr + dot_product(n(:NC), Hgi(:NC)))/sum(n(:NC))

    if (calc == 8 .OR. calc == 9) DH_n(:NC) = DH_n(:NC) + Hgi(:NC)
    if(calc == 3)then

      call Cp_GI(NC, T, n, Hgi, Cp)
      Cp = cPr + Cp

    endif

  endif
  
!   if(calc == 9 .OR. calc == 10)then
!     
! !  Virial coefficient calculation
! !     call CubicVirial(NC, n, T, B, h)
!     
!   endif
  
  Z = P*V/R/T/sum(n)
  
  return
endsubroutine
!-----------------------------------------------------------------------------------------------------------------------------
! SUBROUTINE MTC_EoS (calc, phase_type, guess, NC, T, P, n, Z0, lnPhi, DLnPhi_N, DLnPhi_T, DLnPhi_P, H, Cp, DH_n, Z, B, ic)
! !
!   implicit none
!   integer, parameter           :: NCM = 30
!   real(8), parameter           :: R = 82.05D0
! !
!   integer                      :: calc, guess, ic, n_deriv, NC, phase_type, temp_deriv = 0
!   integer                      :: iCEq, NC1, NONHY, NH, NST
!   
!   real(8)                      :: B, Cp, Cpr = 0, H, HR = 0, P, Sr = 0, T, Ur = 0, V, V0, Z, Z0
!   real(8), dimension(NCM)      :: DH_N, DLnPhi_P, DLnPhi_T, Hgi, lnPhi, n
!   real(8), dimension(NCM,NCM)  :: DLnPhi_N
! 
!   
!   Cp = 0; B = 0
!   n_deriv = 1
!   temp_deriv = 0
!   if (calc /= 10)  then
!     
!     if ((calc == 5) .OR. (calc == 6) .OR. (calc == 9)) n_deriv =  2
! 
!     if ((calc /= 1) .AND. (calc /= 5)) temp_deriv = 1
! 
!     V0 = sum(n(:NC))*Z0*R*T/P
!     call termo (phase_type, n_deriv, temp_deriv, guess, IC, NC, T, P, n, V0, V, lnPhi, DLnPhi_T, DLnPhi_P, DLnPhi_N, Hr, Sr, Cpr)
! 
!     if (calc == 8 .OR. calc == 9) then
! 
! !  C�lculo de entalp�as molares parciales
!       DH_n(:NC) = -R*T**2*DLnPhi_T(:NC)
! 
!     endif
! 
!   endif
!   if((calc == 0) .OR. (calc == 2) .OR. (calc == 3) .OR. (calc == 7) .OR. (calc == 8) .OR. (calc == 9))then
! ! 
! !  Add ideal gas contribution to enthalpy and heat capacity
!     call NTlpy_GI (NC, T, Hgi(:NC))
!     H = (Hr + dot_product(n(:NC), Hgi(:NC)))/sum(n(:NC))
! 
!     if (calc == 8 .OR. calc == 9) DH_n(:NC) = DH_n(:NC) + Hgi(:NC)
!     if(calc == 3)then
! 
!       call Cp_GI(NC, T, n, Hgi, Cp)
!       Cp = cPr + Cp
! 
!     endif
! 
!   endif
!   
! !   if(calc == 9 .OR. calc == 10)then
! !     
! ! !  Virial coefficient calculation
! ! !     call CubicVirial(NC, n, T, B, h)
! !     
! !   endif
!   
!   Z = P*V/R/T/sum(n)
!   
!   return
! endsubroutine
!-----------------------------------------------------------------------------------------------------------------------------
!
!
subroutine gamma_phi (model, calc, phase_type, NC, T, P, n, lnPhi, DLnPhi_N, DLnPhi_T, DLnPhi_P, H, Cp, DH_n, Z, B, IC)

  implicit none
  integer, parameter           :: NCM = 30
  real(8), parameter           :: R = 82.05D0

  integer                      :: calc, guess, ic, model, n_deriv, NC, phase_type, p_deriv, temp_deriv = 0
  
  real(8)                      :: B, Cp, Cpr = 0, H, HR = 0, P, Sr = 0, T, Ur = 0, Z
  real(8), dimension(NCM)      :: DH_N, DLnPhi_P, DLnPhi_T, Hgi, lnPhi, n
  real(8), dimension(NCM,NCM)  :: DLnPhi_N
  Cp = 0; B = 0
  n_deriv = 1
  temp_deriv = 0
  p_deriv = 0
  if (calc /= 10)  then
    
    if ((calc == 5) .OR. (calc == 6) .OR. (calc == 9)) n_deriv =  2

    if ((calc /= 1) .AND. (calc /= 5)) then
      temp_deriv = 1
      p_deriv = 1
    endif

    call fugacity_from_gamma_phi (model, NC, phase_type, temp_deriv, p_deriv, n_deriv, T, P, n, lnPhi, DLnPhi_N, DLnPhi_T, DLnPhi_P, Z, IC) 
    


  endif
!   if((calc == 0) .OR. (calc == 2) .OR. (calc == 3) .OR. (calc == 7) .OR. (calc == 8) .OR. (calc == 9))then
! ! 
! !  Add ideal gas contribution to enthalpy and heat capacity
!     call NTlpy_GI (NC, T, Hgi(:NC))
!     H = (Hr + dot_product(n(:NC), Hgi(:NC)))/sum(n(:NC))
! 
!     if (calc == 8 .OR. calc == 9) DH_n(:NC) = DH_n(:NC) + Hgi(:NC)
!     if(calc == 3)then
! 
!       call Cp_GI(NC, T, n, Hgi, Cp)
!       Cp = cPr + Cp
! 
!     endif
! 
!   endif
  
!   if(calc == 9 .OR. calc == 10)then
!     
! !  Virial coefficient calculation
! !     call CubicVirial(NC, n, T, B, h)
!     
!   endif

  
  return
endsubroutine gamma_phi
!--------------------------------------------------------------------------------
!
!  Ideal gas heat capacity at constant pressure.
!
!  Calculation is based on Passut and Danner correlation [7]
!
subroutine Cp_GI (NC, T, n, Cpgi, Cpmix)
!
!  This subroutine calculates ideal gas contribtion to the constant pressure heat
!  capacity of a mixture. Also, it gave back indivitual contribtion.
!
!  Correlation used is based on the Passut & Danner work [1]. Coefficients are red 
!  in input file though there should be a database in the future.
!
  implicit none
  integer, parameter      :: NCM = 30

  integer                 :: NC
  real(8)                 :: Cpmix, T
  real(8), dimension(NCM) :: CpGI, A, B, C, D, E, F, G, M, n, Peneloux
  
!  Unit conversion constants:
  real(8)                 :: BTU_atmCm3 = 1.0412238D4, g_lb = 454.D0, R_K = 1.8D0, TRan
  
  common /GCPROM/            M, Peneloux, A, B, C, D, E, F, G
!
!  Parameters keep original units: BTU/(molLB R)
  TRan = T*R_K  !K -> �R
  Cpgi = 0.D0
  Cpgi(:NC) = B(:NC) + 2.D0*C(:NC)*TRan**2 + 3.D0*D(:NC)*TRan**3 + 4.D0*E(:NC)*TRan**4 + 5.D0*F(:NC)*TRan**5
!
!  Turn into package units: atm*cm3/mol-g/K
  Cpgi(:NC) = Cpgi(:NC)*M(:NC)*BTU_atmCm3/g_lb*R_K  !(BTU/lb �R)*(10412 atm cm3/BTU)*(M lb/mol-lb)*(454 mol-g/mol-lb)*(1.8 �R/K)
  Cpmix = dot_product(n(:NC),Cpgi(:NC))/sum(n(:NC))
!
  return

endsubroutine

!-----------------------------------------------------------
!
!  Optim is a simple optimization subroutine. It doesn't check for positive definitness
!  but only test for decrease in function Q. If Q is increased after step DX, it update is
!  decreased proportionally to alpha.
!  All independent variables must be non-negative. If one becomes negative after an 
!  iteration, it is deactivated. Reactivation is carried out when after a reevaluation of Q, 
!  g(i) < 0.
!
!  NDim is the number of independent variables
!  X are the NDim independient variables
!  FUNC is de name of the subroutine which provides de function to be minimized, gradient and hessian.
!  Q is the value of the function to be minized
!  g is the gradiento of Q
!  H is the hessian of Q
!  DX is the correction step
!  alpha is the correction to the stepsize DX
!  iAct is an integer vector which counts for acivated and deactivated variables
!
subroutine Optim (NDim, iAct, X, func, Q, g, H, Restr, tol)
!
  implicit real(8) (a-h, o-z)
  
!   parameter(tol = 1D-8, maxIt = 100)
  integer, parameter :: maxIt = 1000
    
  real(8) g(NDim), H(NDim,NDim), X(NDim), macheps
    
  real(8) Xnew(NDim), DX(NDim), gNew(NDim), Hnew(NDim,NDim)
  integer iAct(NDim), iPiv(NDim), Restr
  
!  Generate machine epsilon:
  macheps = epsilon(1.D0)  !this should be of the order of 2.2E-16
!  Init variables:  
  ErrQ = 1
  do i = 1, NDim
    g(i) = 1.
    gNew(i) = 1.
    H(i,:NDim) = 0.
    Hnew(i,:NDim) = 0.
  enddo
!
!  Fist objective function evaluation
  call func(NDim, iAct, X, Q, g, H)
!  
  im = 0
!
!  Start external convegence loop
  do in = 1, maxit
!  
!  If... the initial value satisfies g = 0, then it's already 
!  an optimum (it could happen)
    if(maxval(abs(g)) <= tol)then
!    
      return  !*********************
!    
    endif
!  
!  Auxiliary counter:
    im = 1 + im      
    DX = -g
!  Solving linear set of equations:
    call LUDcmp(H, NDim, NDim, iPiv, d)
    call LUBksb(H, NDim, NDim, iPiv, DX)
!  Restart of alfa
    alpha = 1.D0
    Qnew = Q + dabs(Q/10)      
!
!  --------------------------------------------------------------------
    do while(Qnew > Q*(1D0 + 10*macheps))
!      
!  Exit if alpha*|g|_1 < tolerance. 
      if(maxval(abs(alpha*DX)) <= tol)exit
  
!  Calculating new X
      do i = 1, NDim
    
        Xnew(i) = X(i) + alpha*DX(i)
        if(Restr == 1)then
    
!  Negative values not allowed
          if(Xnew(i) < 0.D0)then
    
            Xnew(i) = 2*macheps*X(i)
            iAct(i) = 0
            
          endif
          
        endif
        
      enddo
!
      call func(NDim, iAct, Xnew, Qnew, gnew, Hnew)
!
      if(Qnew < Q*(1D0 + 10*macheps))then

!  Accept step:
        Q = Qnew
        X = Xnew
        g = gnew
        H = Hnew
        exit

      else

!  If it doesn't, cut the step:
        alpha = alpha/3.1623
        if(maxval(abs(alpha*DX)) <= tol)exit        
        im = im + 1

      endif
  
    enddo
!  --------------------------------------------------------------------    
    if(maxval(abs(alpha*DX)) <= tol)exit  
  enddo
!
  return

endsubroutine Optim
!-----------------------------------------------------------
!
!  SAFTQfunction evaluates the function to be minimized in order to obtain the non-bonded fraction
!  of all the associating sites in a phase.
! 
!  La subrutina Qfunction calcula el valor de la funci�n Q, su gra-
!  diente y Hessiano para el c�lculo de segundo orden de la fracci�n
!  no asociada, seg�n lo sugerido por M. L. Michelsen, IECR 2002, 45,
!  8449-8453.
! subroutine SAFTQfunction (NST, iAct, X, Q, g, H)

!
!  Q = Sum(sm(k)�(ln X(k) - X(k) + 1) - 1/2�Sum(Sum(sm(k)�sm(l)�Delta(k,l)�X(l)/V),l=1,NST),k=1,NST)
!
!    = Q1 + Q2
!
!  g(k) = sm(k)/X(k) - sm(k) - Sum(sm(k)*sm(l)*Delta(k,l)*X(l)/V),l=1,NST)
!
!    = sm(k)/X(k) - sm(k) - SUMA
!
!  "H"(k,l) = -(sm(k) + SUMA)/X(k)*d(k,l) - sm(k)*sm(l)*Delta(k,l)/V
!
!  Donde: sm(k)      = es la cantidad de moles de sitios "k"
!         X(k)       = es la fracci�n de sitio "k" no asociada
!         Delta(k,l  = fuerza de asociaci�n entre sitios "k" y "l"
!         NST        = n�mero de sitios totales
!         d(k,l)  = delta de Kronocker, = 1 si "l" = "k", sino = 0.
!         g(k)       = gradiente de Q en la direcci�n "k"
!        "H"(k,l)   = hessiano k,l
!  nota1: tomo a rho = 1/V donde V = vol total. Asume que son iguales porque los moles est�n normalizados (creo).
!
!  Francisco, febrero de 2010.
!
!  ---------------------------------------
!  Siendo que cada vez que se corta un c�lculo, el programa corta casi siempre en esta subrutina.
!  Eso indica que esta es un cuello de botella, por eso, comprim�� el c�lculo de las variables
!  Q, g y H en la menor cantidad de lazos posible. v-1.9.18
!
!-------------------------------------------------------------------
! 
! Function to count blanck spaces at the end of a string.
!
integer function strlen(st)
      
    integer        :: i
    character      :: st*(*)
    
    i = len(st)
    do while (st(i:i) == ' ')
      i = i - 1
    enddo
    strlen = i
    return
end function strlen
!
!-------------------------------------------------------------------
!
subroutine get_critical_properties (NC, Tc, Pc, omega)

  implicit none
  
  integer, parameter                  :: NCM = 30
  
  integer, intent(in)                 :: NC
  
  real(8), dimension(NC), intent(out) :: omega, Pc, Tc
  
  real(8), dimension(NCM)             :: Tc_common, Pc_common, omega_common
  
  common /CRIT/                          Tc_common, Pc_common, omega_common
  
  Tc = Tc_common(:NC)
  Pc = Pc_common(:NC)
  omega = omega_common(:NC)
  
  return
endsubroutine get_critical_properties
!----------------------------------------------------------------------------
!  
!
!  Generic subroutine for calculation of fugacity coefficients at (T,P,n), and
!  its derivatives.
!
!
SUBROUTINE TERMO (MTYP, NDer, NTemp, guess, IC, NC, T, P, rn, V0, &
                  V, PHILOG, DLPHIT, DLPHIP, FUGN, Hr, Sr, Cpr)
      
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  integer, PARAMETER                       :: NCM = 30
  real(8), parameter                       :: RGAS = 82.05D0 !0.08314472d0)
                                           
                                           
  integer, intent(in)                      :: guess, MTyp, NC, NDer, NTemp
  real(8), intent(in)                      :: P, T, V0
  real(8), DIMENSION(NCM), intent(in)      :: rn
                                           
  integer, intent(out)                     :: IC
  real(8), intent(out)                     :: V
  real(8), DIMENSION(NCM), intent(out)     :: PHILOG, DLPHIT, DLPHIP
  real(8), dimension(NCM,NCM), intent(out) :: FUGN
  
  real(8), dimension(NCM)                  :: Arn, ArVn, ArTn, DPDN
  real(8), dimension(NCM,NCM)              :: Arn2
  common /ScndDer/                            dPV, dPDT, dPdn

  TOTN = sum(rn)
  CALL VCALC (MTYP, guess, NC, rn, T, P, V0, V, IC)
  RT = RGAS*T
  Z = P*V/TOTN/RT
  call ArVnder (NDER, NTEMP, NC, rn, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2, ArT, ArTT)
  DPV = -ArV2 - RT*TOTN/V**2
  DPDT = -ArTV + TOTN*RGAS/V
  if (nDer > 0) then
    
    DO 60 I = 1, NC

      PHILOG(I) = Arn(I)/RT - dLOG(Z)
      DPDN(I) = RT/V - ArVn(I)
      DLPHIP(I) = -DPDN(I)/DPV/RT - 1.D0/P
      IF (NTEMP > 0) then

        DLPHIT(I)  = (ArTn(I) - Arn(I)/T)/RT + DPDN(I)*DPDT/DPV/RT + 1.D0/T
      
      endif
60    enddo
62    IF (NDER > 1) then
      DO 63 j = 1, NC
        DO 61 K = j, NC
    
!           Here we use degree-0 functions
          FUGN(j,K) = 1.D0 + TOTN*(Arn2(j,K) + DPDN(j)*DPDN(K)/DPV)/RT
          FUGN(K,j) = FUGN(j,K)

61        enddo
63      enddo
    endif
  endif
  if (nTemp > 0) then
    
    Sr = -ArT
    Ur = Ar + T*Sr
    Hr = Ur + P*V - TOTN*RT
    cVr = -T*ArTT
    cPr = cVr - T*dPdT**2/dPV
    
  endif
  RETURN
endsubroutine TERMO
!---------------------------------------------------------------------------------------
!
!
!
SUBROUTINE ArVnder (NDER, NTD, NC, rn, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2, ArT, ArTT)

  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  
  integer, PARAMETER                       :: NCM = 30
                                           
  integer, intent(in)                      :: NC, NDer, NTD
  real(8), intent(in)                      :: T, V
  real(8), dimension(NCM), intent(in)      :: rn
                                           
  real(8), intent(out)                     :: Ar, ArTV, ArT, ArTT, Arv
  real(8), dimension(NCM), intent(out)     :: Arn, ArTn, ArVn
  real(8), dimension(NCM,NCM), intent(out) :: Arn2
  
  COMMON /MODEL/ n_model
  
!   IF (n_model <= 2) THEN
! !  SRK or PR
!     CALL HelmSRKPR(NDer,NTD,rn,V,T, Ar,ArV,ArTV,ArV2,Arn,ArVn,ArTn,Arn2)
! 
!   ELSE IF (n_model.EQ.3) THEN
  IF (n_model == 4) THEN

    CALL HelmRKPR (NDER, NTD, NC, rn, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2, ArT, ArTT)

   ELSE IF (n_model == 5 ) THEN
! 
    CALL HelmPCSAFT (NDER, NTD, NC, rn(:NC), V, T, Ar, ArV, ArTV, ArV2, Arn(:NC), ArVn(:NC), ArTn(:NC), Arn2(:NC,:NC), ArT, ArTT)
! 
!    ELSE IF (n_model.EQ.6) THEN
!   
!     CALL HelmSPHCT(NDER,NTD,rn,V,T,Ar,ArV,ArTV,ArV2,Arn,ArVn,ArTn,Arn2)
! 
   ELSEif (n_model == 6 ) then  !  
   
        call Helmhotz_SL (NDER, NTD, NC, rn, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2, ArT, ArTT)
!                       (NDER, NTD, n, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2, ArT, ArT2)

   ELSEif (n_model == 7) then  !  
   
        call Helmholtz_MTC (NC, T, V, rn, NDER, NTD, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2, ArT)

!   elseif (n_model == 8) then
!     
! !     call HelmholtzPPCSAFT (nDer, NTD, rn, V, T, Ar, ArV, ArTV, ArV2,Arn, ArVn, ArTn, Arn2)
!     call Helmholtz_CPSAFT (nDer, NTD, rn, V, T, Ar, ArV, ArTV, ArV2,Arn, ArVn, ArTn, Arn2)
!     
  END IF
  return
endsubroutine ArVnder
!-------------------------------------------------------------------------------
!
SUBROUTINE VCALC (ITYP, guess, NC, rn, T, P, V0, V, IC)
!
!  ROUTINE FOR CALCULATION OF VOLUME, GIVEN PRESSURE
!  
!  INPUT:
!  
!  ITYP:        TYPE OF ROOT DESIRED
!  NC:          NO. OF COMPONENTS
!  rn:          FEED MOLES
!  T:           TEMPERATURE
!  P:           PRESSURE
!  
!  OUTPUT:
!  
!  V:           VOLUME
!
!
!  Important note:
!               ZETA, the iteration variable, is *not* the compressibility
!               factor, but b/v. It represents the GREEK letter ZETA.
!
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  
  integer, PARAMETER                  :: NCM = 30
  real(8), parameter                  :: bisecTol = 1, RGAS = 82.05D0, newtol = 1E-10   !0.08314472d0
                                      
  integer, intent(in)                 :: iTyp, guess
  real(8), intent(in)                 :: P, T, V0
  real(8), dimension(NCM), intent(in) :: rn
  
  integer, intent(out)                :: IC
  real(8), intent(out)                :: V
                                      
  real(8)                             :: P_common, T_common
  real(8), dimension(NCM)             :: Arn, ArTN, ArVn, rn_common 
  real(8), dimension(NCM,NCM)         :: Arn2, dBij
  integer                             :: n_deriv, NC_common, t_deriv

  
  common /fromVCALC_real/                rn_common, T_common, P_common
  common /fromVCALC_integer/             n_deriv, NC_common, t_deriv
                                         
  external                               deltapressure
  
  LOGICAL FIRST_RUN
  NDER = 0
  NTemp = 0
  FIRST_RUN = .TRUE.
  TOTN = sum(rn)
  call Bcalc (NC, rn, T, B)
  CPV=B
   S3R = 1.D0/CPV
   DEL = 00D0
   ITER = 0
  alpha = 1      
!
   ZETMIN = 0.D0
   ZETMAX = 1.D0 !- 2D-5*T  !.99D0  This is flexible for low T (V very close to B)
   if (guess == 0) then
    IF (ITYP .GT. 0) THEN
    
      ZETA = .5D0
    
    ELSE
    !    IDEAL GAS ESTIMATE
      ZETA = MIN (.5D0,CPV*P/(TOTN*RGAS*T))
     
    endif
  else
    
    ZETA = B/V0
    IC = ITYP
    
  ENDIF

!  Start volume convergence. Convergence variable is zeta = b/v, "zeta" being the greek 
!  character, not related to the compressibility factor.
!100  CONTINUE
100  del = 1
  iter = 0
  do while (ABS(DEL) .GT. 1D-10)

    V = CPV/ZETA
    ITER = ITER + 1
    call ArVnder (NDER, NTEMP, NC, rn, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2, ArT, ArTT)
    PCALC = TOTN*RGAS*T/V - ArV
    IF (PCALC .GT. P) THEN
    
      ZETMAX =  min(ZETA *1.05D0, 1D0)
      
    ELSE
    
      ZETMIN = ZETA /1.05D0
      
    ENDIF
    AT  = (Ar + V*P) /(T*RGAS) - TOTN*LOG(V)
!    AT is something close to Gr(P,T)
    DER = (ArV2*V**2 + TOTN*RGAS*T)*S3R  ! this is dP/d(rho/B)
    delOld = del
    DEL = -(PCALC-P)/DER*alpha
    if (abs(del) < newtol) then
      
      exit
      
    elseif (abs(del + delOld)/abs(del) < 1D-1) then
      
      del = del/2
      
    endif
    ZETA = ZETA + MAX (MIN (DEL, 0.1D0), -.1D0)
    if (iter >= 30) then      
    
      T_common = T
      P_common = P
      rn_common = rn
      NC_common = NC
      call screening (20, deltapressure, CPV/zetMax, CPV/zetMin, bisecTol, V, DP)
      zeta = CPV/V
      iter = 0
    
    endif
    
      IF (ZETA .GT. ZETMAX .OR. ZETA .LT. ZETMIN) ZETA = .5D0*(ZETMAX+ZETMIN)
  
    IF (ABS(PCALC - P) < 1D-12) exit !GOTO 101
      
!       IF (ABS(DEL) .GT. 1D-10) GOTO 100
  enddo
!       Francisco, 11/08/2014
!   if (abs(DEL) > 1.D-8) goto 100

101  IF (ITYP == 0 ) THEN
!
!    FIRST RUN WAS VAPOUR; RERUN FOR LIQUID
!
    IF (FIRST_RUN) THEN

      VVAP = V
      AVAP = AT
      FIRST_RUN = .FALSE.
      ZETA = 0.5D0
      ZETMAX = 1.D0 !- 0.01*T/500
      alpha = 1

      GOTO 100

    ELSE
    
      IF (AT > AVAP) then
      
        V = VVAP
        IC = -1
        
      else
    
        IC = 1
        
      endif
      
    ENDIF
    
  else
    
    ic = 1 !no error, desired phase or not.
    
  ENDIF
ENDsubroutine VCALC
!-----------------------------------------------------------------------------------
!
!  This general subroutine provides the "co-volume" for specified composition,
!  that will be used by Evalsecond or Vcalc
!
SUBROUTINE Bcalc (NC, x, T, BMIX)

  IMPLICIT real(8) (A-H, O-Z)
  integer, PARAMETER                  :: NCM = 30
  real(8), parameter                  :: RGAS = 82.05D0 !0.08314472d0)
  
  integer, intent(in)                 :: NC
  real(8), intent(in)                 :: T
  real(8), dimension(NCM), intent(in) :: x
  
  real(8), intent(out)                :: Bmix
  
  real(8)                                dBi(NCM), dBij(NCM,NCM), DDB(0:3,NCM)
  DIMENSION                              DD(0:3,NCM),DDT(0:3,NCM),DTT(0:3,NCM),DIA(NCM)
!   
!   real(16)      :: quad_b, quad_T, quad_x(NCo), quad_rho
  
  COMMON /MODEL/   n_model
!   COMMON /MOL/     DC(2),D(2),DT(2),HA(2),HB(2)
!   COMMON /MIXT/    VCPM,CMIX,CVYM,CVYMT,dCVYM(MAXC),dCMIX(MAXC),    &
!                  dCVYMT(MAXC),d2CVYM(MAXC,MAXC),d2CMIX(MAXC,MAXC)
  COMMON /MIXRULE/ NSUB
  COMMON /BMIX/    B
  COMMON /forB/    DDB
!   COMMON /NG/      NGR,NST
  COMMON /rule/    ncomb
  
!   NG=NGR
!   NC=2
!   if (n_model == 5 .OR. n_model == 7) then
!     
!     CALL PARAGC(T,NCO,NG,NST,1)
!     PI=3.1415926536D0
!     XLAM3=0.0d0
!     DO I=1,NCO
!     
!       DGC=D(I)
!       XLAM3=XLAM3+X(I)*DGC**3
!       
!     enddo
!     B=PI/6.D0*XLAM3/1.0D3
!   
!   else if(n_model.EQ.4)then
!     
  if (n_model == 5) then
    
    DD = DDB
    CALL DIAMET (NC, T, DIA, DD, DDT, DTT, NSUB)
!     B = RGAS*(DD(3,1)*X(1)+DD(3,2)*X(2))  !S3
    B = Rgas*sum(dd(3,:NC)*x(:NC))
!   
  elseif (n_model == 10) then
  
    call get_covolume (b)
    
!     
!     CALL Mixture_Param(NSUB,NC,X,T)
!     B=VCPM
!   
!   elseif (n_model == 8)then
!     
!     quad_T = T
!     quad_x = x
!     quad_rho = 1
!     call evaluate_mixture_parameters (NC, quad_T, quad_rho, quad_x)
!     call get_covolume(quad_b)
!     B = quad_b
    
  else
!     if(ncomb < 2)then
      
      call Bnder (NC, x, B, dBi, dBij)  ! Bmix is used in EVALSECOND
      
!     else
!     
!       call Bcubicnder(2,x,B,dBi,dBij)
      
!     end if
  end if
  BMIX=B
  return
endsubroutine Bcalc  
!--------------------------------------------------------------------------------------------------------------------------
!
! This subroutine aims to find a good initial estimate if the conventional has failed. It will be called whenever
! the Newton procedure fails.
! Screening searchs for a root by trial and error, calculating many points in an inverval set by the user. If many
! possible roots are found, the subroutine will work among the closer to the value initially proposed.
!
subroutine screening (n_screen, func, xMin, xMax, tolerance, x, F)
  
  implicit none
  
  integer                               :: i, j, n_candidates = 0
                                        
  real(8)                               :: deltaW = 0, eps, fDscld, fScld, &
&                                           w = 0, w0 = 0, Z = 0
                                        
  integer, dimension(:), allocatable    :: candidate
  real(8), dimension(:), allocatable    :: F_screen, w_screen
                                        
  integer, intent(in)                   :: n_screen
  real(8), intent(in)                   :: xMin, xMax
!   real(8), intent(out)                  :: F, x
  real(8)                 :: F, x, tolerance, absF
  
  eps = dsqrt(epsilon(1D0))
!   w0 = fScld (x0, 1D0, 0D0, 1D0)
  allocate (w_screen(n_screen+1), F_screen(n_screen+1), candidate(n_screen+1))
  n_candidates = 0
  candidate(:n_screen+1) = 0
  
! Set interval

  x = xMin + eps

  call func (x, F_screen(1))
  w_screen(1) = fScld(x, 1D0, xMin, xMax)
  write (*, *) x, F_screen(1)  

  deltaW = (fScld(xMax - eps, 1D0, xMin, xMax) - w_screen(1))/dble(n_screen)
  
! Start screening with logarithmic scaling, within the interval (xMin, xMax).

  do i = 1, n_screen
  
    w_screen(i+1) = w_screen(i) + deltaW
    x = fDscld(w_screen(i+1), 1D0, xMin, xMax)
    call func (x, F_screen(i+1))
    write (*, *) x, F_screen(i+1)
    
! Look for a candidate:

    if (F_screen(i)*F_screen(i+1) < 0D0) then
    
      n_candidates = n_candidates + 1
      candidate(n_candidates) = i
      exit
!       if (n_candidates > 1) then
! 
! !         check if we are moveing away from the initial value originally applied.
!         if (dabs(w0 - w_screen(i+1)) > dabs(w0 - w_screen(candidate(n_candidates-1))) ) then
!           n_candidates = n_candidates - 1
!           exit !we've started to move away from x0
!         endif
!       endif
      
    endif
    
  enddo
  if (n_candidates == 0) then
  
!     No candidates have been found :-(
    return
    
  else
  
!    Screening finished, start Bisection within the interval with candidate closer to x0
    i = candidate(n_candidates)
    do j = 1, n_screen
  
      w = (w_screen(i) + w_screen(i+1))/2
      x = fDscld(w, 1D0, xMin, xMax)
      call func (x, F)
      absF = dabs(F)
      if (absF < tolerance) exit
      if (F*F_screen(i+1) < 0) then
    
        w_screen(i) = w
        F_screen(i) = F
        
      else
    
        w_screen(i+1) = w
        F_screen(i+1) = F          
        
      endif    
      
    enddo
  
  endif
  
  return
endsubroutine screening
!------------------------------------------------------------------------
!
!  difference between specified and calculates pressure
!
subroutine deltapressure (V, deltaP)

  implicit none
  
  integer, parameter          :: NCM = 30
  real(8), parameter          :: R = 82.05D0 !0.08314472d0
                              
  integer                     :: n_deriv, NC, NC_common, t_deriv
  real(8)                     :: Ar, ArT, ArTT, ArTV, Arv, ArV2, Pspec, Pcalc, T, T_common
  real(8), dimension(NCM)     :: Arn, ArTn, ArVn, n, n_common
  real(8), dimension(NCM,NCM) :: Arn2                       
                          
!   real(8), intent(in)       :: V
!   real(8), intent(out)      :: deltaP
  real(8)       :: V
  real(8)       :: deltaP
  
  common /fromVCALC_real/      n_common, T_common, Pspec
  common /fromVCALC_integer/   n_deriv, NC_common, t_deriv

  
  T = T_common 
  n = n_common
  NC = NC_common
  call ArVnder (n_deriv, t_deriv, NC, n,  V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2, ArT, ArTT)
!                (NDER,    NTD,     NC, rn, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2)
  PCALC = sum(n(:NC))*R*T/V - ArV
  deltaP = Pcalc - Pspec
  
  return
endsubroutine deltapressure
!--------------------------------------------------------------------------------------------------------------------------
!
!  Funci�n de escalado de variables: �til para variables con topes m�ximo y m�nimo.
!
!  vMin y vMax son las cotas m�nima y m�xima que puede tomar "x". 
!
!  scl es la variable de escalado, que estira o comprime el dominio.
!
double precision function fScld(x,scl,vMin,vMax)
  implicit none
  real(8) scl, vMax, vMin, x
  
  fScld = scl*dlog((x - vMin)/(vMax - x))

  return
endfunction
!--------------------------------------------------------------------------------------------------------------------------
!
!  Funci�n de desescalado de variables: inversa de la funci�n fScld
!
double precision function fDscld(x,scl,vMin,vMax)
  implicit none
  real(8) ex, scl, vMax, vMin, x

  ex = dexp(x/scl)
  fDscld = (vMin + vMax*ex)/(1.d0 + ex)

  return
endfunction
!--------------------------------------------------------------------------------------------------------------------------
