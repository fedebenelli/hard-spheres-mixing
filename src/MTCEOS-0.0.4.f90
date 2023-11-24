module parametersMTC

  integer, parameter          :: att = 2, fv = 1, NCM = 30, NGM = 60
  real(8), parameter          :: PSI = 1, Rgas = 82.05, z = 10
  
  integer                     :: NG
  integer, dimension(NCM,NGM) :: nu
  
  real(8)                     :: ltot, qtot, rtot, vstr
  real(8), dimension(NGM)     :: Atemp, Btemp, etot, l, nqtot, q, r, u0_R
  real(8), dimension(NGM,NGM) :: e, kij, u_R
  real(8), dimension(NCM,NGM) :: nu_q

endmodule parametersMTC
!----------------------------------------------------------------
subroutine MTC_version (version, thermo_name)

	character(20)  ::thermo_name, version

	version = "0.0.4"
	thermo_name = "MTC-EOS"
	return
endsubroutine MTC_version
!-----------------------------------------------------------------
subroutine Helmholtz_MTC (NC, T, V, n, n_deriv, t_deriv, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2, ArT)

  use parametersMTC
  implicit none
  
  integer, intent(in)                      :: NC, n_deriv, t_deriv
  real(8), intent(in)                      :: n(NCM), T, V
                                           
  real(8), intent(out)                     :: Ar, ArT,  ArTV, ArV, ArV2    !ArTT
  real(8), dimension(NCM), intent(out)     :: Arn, ArVn, ArTn
  real(8), dimension(NCM,NCM), intent(out) :: Arn2
  
  integer                                  :: k, m
  real(8)                                  :: eta, h0, h1, h2(NGM), h3, h4(NGM), h5, h6, nbar_ltot, ntot, r_q_1k, r_q_1m, RT, Vred, z2
  real(8), dimension(2)                    :: Ar_, ArT_, ArTT_, ArTV_, ArV_, ArV2_
  real(8), dimension(NCM,2)                :: Arn_, ArTn_, ArVn_
  real(8), dimension(NGM)                  :: Arnq, r_q
  real(8), dimension(NGM,2)                :: Arnq_
  real(8), dimension(NGM,NGM,2)            :: Arnq2_
  real(8), dimension(NCM,NCM,2)            :: Arn2_

  call MTC_model_parameters (NC, T, n(:NC))
  z2 = z/2
  RT = Rgas*T
  ntot = sum(n(:NC)) 
  Vred = V/vstr !(mol)
  eta = rtot/Vred
  
  h0 = Vred - rtot
  h1 = h0 + qtot
  h2(:NG) = h0 + etot(:NG)
  h3 = 1D0 - eta 
  h4(:NG) = dlog(h1/h2(:NG))
  h5 = dlog(h1/Vred)
  h6 = dlog(h3)
  nbar_ltot = rtot - z2*(rtot - qtot) 
  
  Ar_(fv) = -z2 * h1 * h5  +  h0*h6 + nbar_ltot
  Ar_(att) = PSI*sum( nqtot(:NG)*h4(:NG) )
  Ar = RT*(Ar_(fv) + Ar_(att))
  
  ArV_(fv) =  -z2*h5/vstr + h6/vstr  +  nbar_ltot/V
  ArV_(att) = PSI/vstr /h1 * sum( nqtot(:NG)*(etot(:NG) - qtot)/h2(:NG) )
  ArV = RT*(ArV_(fv) + ArV_(att))
  
  ArV2_(fv) = (-z2*(rtot - qtot)**2 /h1 + rtot*rtot/h0)/V/V
  ArV2_(att) = -PSI/vstr/vstr /h1/h1 * sum( nqtot(:NG)*(etot(:NG) - qtot)*(h1 + h2(:NG))/h2(:NG)**2 )
  ArV2 = RT*(ArV2_(fv) + ArV2_(att))
  if (n_deriv > 0) then
  
    do m = 1, NG
    
      r_q(m) = r(m)/q(m)
      Arnq_(m,att) = PSI*(h4(m) - qtot*(r_q(m) - 1)/h1 + sum( nqtot(:NG)*((r_q(m) - e(:NG,m))/h2(:NG)) ) )
      Arnq_(m,fv) = z2*(r_q(m) - 1)*h5 - r_q(m)*h6 
      Arnq(m) = (Arnq_(m,fv) + Arnq_(m,att))

    enddo
!     Arn(:NC) = matmul(nu_q(:NC,:NG), Arnq_(:NG,fv))
!    Arn(:NC) = matmul(nu_q(:NC,:NG), Arnq_(:NG,att))
    Arn(:NC) = RT*matmul(nu_q(:NC,:NG), Arnq(:NG))
    if (n_deriv > 1) then
    
      do k = 1, NG
        do m = 1, NG
          r_q_1m = (r_q(m) - 1)
          r_q_1k = (r_q(k) - 1)
          Arnq2_(m,k,fv) = -z2*r_q_1m*r_q_1k/h1 + r_q(m)*r_q(k)/h0
          Arnq2_(m,k,att) = -PSI*((r_q_1k + r_q_1m)/h1 - (r_q(m) - e(k,m))/h2(k) - (r_q(k) - e(k,m))/h2(m) &
                            + qtot*r_q_1m*r_q_1k/h1/h1 - sum( nqtot(:NG)*((r_q(m) - e(:NG,m))*(r_q(k) - e(:NG,k))/h2(:NG)**2) ))
        enddo
      enddo
!       Arn2(:NC,:NC) = matmul( nu_q(:NC,:NG), matmul(Arnq2_(:NG,:NG,fv), transpose(nu_q(:NC,:NG)))  )
!      Arn2(:NC,:NC) = matmul( nu_q(:NC,:NG), matmul(Arnq2_(:NG,:NG,att), transpose(nu_q(:NC,:NG)))  )
      Arn2(:NC,:NC) = RT*matmul( nu_q(:NC,:NG), matmul(Arnq2_(:NG,:NG,fv) + Arnq2_(:NG,:NG,att), transpose(nu_q(:NC,:NG)))  )
    
    endif
    
  endif
  if (t_deriv > 0) then
  
    ArT = 0
    ArTV = 0
!     ArTT = 0
    ArTn = 0
    ArVn = 0
  endif
  return
  
endsubroutine Helmholtz_MTC
!----------------------------------------------------------------
!
!

subroutine MTC_model_parameters (NC, T, n)

  use parametersMTC
  implicit none
  
  integer, intent(in) :: NC
  real(8), intent(in) :: n(NC), T
  
  
  integer             :: j, k
  real(8)             :: nbar
  
  qtot = 0
  rtot = 0
  ltot = 0
  do j = 1, NG
    
    nbar = dot_product(n(:NC),nu(:NC,j)) !moles of group "j"
    nqtot(j) = nbar*q(j)                 !moles of surface "j"
    qtot = qtot + nqtot(j)                   !total surface area
    rtot = rtot + nbar*r(j)                  !total volume of segments
    l(j) = z/2 * (r(j) - q(j)) - r(j) + 1
    ltot = ltot + nbar*l(j)                  !"total" bulkiness
    u_R(j,j) = u0_R(j)*(1 + Btemp(j)/T)!self-interaction energy        
    
  enddo
  do j = 1, NG
  
    etot(j) = 0
    do k = 1, NG
    
      if (k /= j) then 
        u_R(k,j) = kij(k,j)*(u_R(k,k) + u_R(j,j))/2
      endif
      e(k,j) = dexp(-u_R(k,j)/T)          !exponential of interaction energy between segments
      etot(j) = etot(j) + nqtot(k)*e(k,j)   !average "total" exponential energy of "j" pondered by its surface.     
      
    enddo
  enddo
  return
  
endsubroutine MTC_model_parameters
!------------------------------------------------------------------
subroutine MTC_covolume (NC, n, b)

  use parametersMTC, only : NG, nu, r, vstr
  implicit none
  
  integer, intent(in) :: NC
  real(8), intent(in) :: n(NC)
  real(8), intent(out) :: b
  
  
  integer             :: j, nbar
  
  b = vstr *dot_product(n(:NC),matmul(nu(:NC,:NG),(r(:NG)))) 

  return
  
endsubroutine MTC_covolume
!----------------------------------------------------------------
!
!
subroutine read_MTC_parameters (inputFile, outputFile, NC)

  use parametersMTC

  implicit none
  
  character, parameter    :: TAB = achar(9)
  
  integer, intent(in)     :: inputFile, NC, outputFile
  
  integer                 :: i, j, k, nown_ij
                      
  real(8), dimension(NCM) :: HHA, HHB, HHC, HHD, HHE, HHF, HHG, M, omega, Pc, Peneloux, Tc
  character(10)           :: common_groupName(NGM), compoundName(NCM), groupName(NGM)
  common /NAME/              compoundName, common_groupName
  common /GCPROM/            M, Peneloux, HHA, HHB, HHC, HHD, HHE, HHF, HHG
  common /CRIT/              TC, PC, omega
  
  do i = 1, NC
    read (inputFile, *) compoundName(i)
    read (inputFile, *) Tc(i), Pc(i), omega(i)
    read (inputFile, *) M(i), Peneloux(i), HHA(i), HHB(i), HHC(i), HHD(i), HHE(i), HHF(i), HHG(i)
  enddo
    
!   read (inputFile, *) 

  read (inputFile, *) NG

  do i = 1, NC
    
    read (inputFile, *) nu(i,:NG)

  enddo
  read (inputFile, *) vstr
  do j = 1, NG
  
    read (inputFile, *) groupName(j),  u0_R(j), Atemp(j), Btemp(j), r(j), q(j)
    kij(:NG,i) = 1
    nu_q(:NC,j) = nu(:NC,j)*q(j)
    
  enddo
  common_groupName(:NG) = groupName(:NG)
  write (outputFile, *) "Compound", TAB, "Tc(K)", TAB, "Pc(atm)", TAB, "omega", (TAB, groupName(j), j = 1, NG)
  do i = 1, NC
    write (outputFile, '(2A, F7.2, A, G10.4, A, F6.4, <NG>(A, I))')compoundName(i), TAB, Tc(i), TAB, Pc(i), TAB, omega(i), (TAB, nu(i,j), j = 1, NG)
  enddo
  read (inputFile, *) nown_ij
  do k = 1, nown_ij
    
    read (inputFile, *) i, j, kij(i,j)
    kij(j,i) = kij(i,j)
    
  enddo
  
  return  

endsubroutine read_MTC_parameters
