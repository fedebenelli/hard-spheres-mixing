	SUBROUTINE PHUB (NC, NMIX, P, TEMP, Y, BB, BF, FUG, ZN)
!  ************************************					   
!  Calculation of fugacity coefficients for pure substances and	   
!  multicomponent mixtures. When specified, the program accounts for    
!  dimerisation in the vapor phase.						 
!  NMIX = 0: only pure component calculation.				   
!  NDIM = 0: no dimerisation are included.					
!							CORRECTED 3/6 - 1982  /TJ   
!
!							Se formateó el código (goto -> then, identación)
!							Ahora se devuelve ln(phi) en lugar de Phi
!							Se incluyerons las derivadas del ln(phi)
!							agosto/2020 /FAS
!  ************************************

	IMPLICIT REAL*8 (A-H, O-Z)
	
	integer, parameter         :: NCM = 30
	real(8), parameter         :: Rg = 82.05D0
	
	integer, intent(in)        :: NC, NMIX
	real(8), intent(in)        :: P, TEMP, Y(NCM), BB(NCM,NCM), BF(NCM,NCM)
	real(8), intent(out)       :: FUG(NCM), ZN
	                           
	integer                    :: I, J, K, MDIM(NCM,NCM), NDIM,  pivot_vector(NCM)
	
	REAL(8)                    :: AS, SDZ, DFDZ(NCM,NCM+1), ResZ
	REAL(8), dimension(NCM)    :: DNDZ, SF, SKJ, SSJ, SQ, Z
	real(8), dimension(NCM,NCM):: EQK, SFF
     
	integer                    :: common_MDIM(NCM,NCM), common_NDIM
! 	real(8), dimension(NCM)    :: common_DMU, common_eta, common_Pc, common_Rd, common_Tc, common_Zc
! 	COMMON /VIRDAT/               common_PC, common_RD, common_DMU, common_ETA, common_TC,ZC
	COMMON /MPHU/                 common_MDIM, common_NDIM
	
	MDIM = common_MDIM
	NDIM = common_NDIM
	RESZ = 1e-8
	if (NMIX > 0) then
	
!	  Calculation for mixture:							   
	  if (NDIM > 0) then
	  
!	    Calculation for mixture with dimerisation:	
	    NC1 = NC+1	
!	    Calculation of true fugacity koefficients:				   
	    DO I = 1, NC	
	    
	      SF(I) = DEXP(BF(I,I)*P/(RG*TEMP))					     
	      AS=DABS(BF(I,I))								  
	      SSJ(I)=AS/BF(I,I)*(AS**(1.D0/3))
	      
	    enddo
	    DO  I=1,NC									
	      DO  J=1,NC	
	      
	        SFF(I,J)=DEXP(((SSJ(I)+SSJ(J))**3)*P/(8.D0*RG*TEMP))	!(39)	  
              SFF(J,I)=SFF(I,J)
              
	      enddo
	    enddo
	    
!	    Calculation of  equlibrium constants:					  
	    DO I = 1, NC									
	      DO J = I, NC
	      
	        EQK(I,J) = 0.D0								     
	        ND=MDIM(I,J)									
	        IF (ND /= 0) EQK(I,J) = (BF(I,J) - BB(I,J)) / (RG*TEMP) * (SF(I)*SF(J)*P/SFF(I,J))	!(37+38+39)							   
	        EQK(J,I) = EQK(I,J)								 
	        
	      enddo
	    enddo
	    
!	    Calculation of true molefractions:					     
!	    Initial guess:									 
	    DO I = 1, NC
	    
	      Z(I)=Y(I)
	      IF (Y(I) == 0.D0) Z(I) = 1.D-24
	      
	    enddo
!	    Calculation of mass - balances and derivatives for the jacobian:     
	    NP=0
	    SDZ = 1
	    do while (SDZ > RESZ)
	    
100         SKIJ=0.D0									   
	      DO I=1,NC									
	        DO J=I,NC									
	        
	          SKIJ = SKIJ + EQK(I,J)*Z(I)*Z(J)
	          
	        enddo
	      enddo
	      ZN=1.D0/(1.D0 + SKIJ)							    
	      DO I=1,NC
	      
	        SKJ(I)=0.D0									 
	        SSJ(I)=0.D0									 
	        DO J=1,NC
	        
	          SKJ(I)=SKJ(I)+EQK(I,J)*Z(I)*Z(J)					    
	          SSJ(I)=SSJ(I)+EQK(I,J)*Z(J)	
	          
	        enddo
	        
	      enddo
            
	      DO I=1,NC
	      
	        SQ(I)=Z(I)+EQK(I,I)*Z(I)**2+SKJ(I)					  
	        DFDZ(I,NC+1)=-(ZN*SQ(I)-Y(I))						 
	        DNDZ(I)=-(EQK(I,I)*Z(I)+SSJ(I))/((1.D0+SKIJ)**2)			
	        DFDZ(I,I)=DNDZ(I)*SQ(I)+ZN*(1.D0+3.D0*EQK(I,I)*Z(I)+SSJ(I))	 
	        
	      enddo
	      DO I = 1, NC									
	        DO K = 1, NC	
	        
	          IF (K /= I) DFDZ(I,K)=DNDZ(K)*SQ(I)+ZN*EQK(I,K)*Z(I)
	          
	        enddo
	      enddo
	      call LUDcmp (dFdz(:NC,:NC), NC, NC, pivot_vector, d)
	      call LUBksb (dFdz(:NC,:NC), NC, NC, pivot_vector, dFdz(:NC,NC1))
	      SDZ = 0.D0									    
	      DO I=1, NC	
	      
	        DZ = DFDZ(I,NC1)								    
	        Q = DABS(DZ/(0.5D0*Z(I)))							 
	        IF (Q > 1.) DZ=DZ/Q							    
	        SDZ = SDZ + DABS(DZ)								  
	        Z(I) = Z(I) + DZ									
	        
	      enddo
	      NP=NP+1									     
	      IF (NP > 15) exit 

	    enddo
	    do i = 1, NC
		
	      FUG(i) = dlog(Z(i)/Y(i)) + (BF(i,i)*P/(RG*TEMP)) !(41)
! 	      if (n_deriv > 0) then
! 	        Derivatives of z wrt y, for number of moles derivatives of ln(phi). 
!	        They were deprecated in original version, but they should do no harm
! 	        for j = 1, NC
! 		    
! 	          dZdY(j,i) = -ZN*SQ(j)
! 		    
! 	        enddo 
! 	        dZdY(i,i) = dZdY(i,i) - 1
! 	        call LUBksb (dFdz(:NC,:NC), NC, NC, pivot_vector, dZdY(:NC,NC1))
! 		  
! 	      endif
		
	    enddo
	    
	  else
	  
!	    Calculation for mixture without dimerisation				 
	    B=0.D0										
	    DO I=1,NC									
	    
	      SSJ(I)=0.D0									 
	      DO J=1,NC									
	      
	        SSJ(I)=SSJ(I)+Y(J)*BB(I,J)						    
	        B = B + Y(I)*Y(J)*BB(I,J)	!mixture virial coefficient						   
	        
	      enddo
	    enddo
	    do i = 1, NC
		
	      FUG(i) = (P*(2.D0*SSJ(i)-B)/(RG*TEMP))
		
	    enddo

	  endif
	
	else
	
!	  Pure component calculation							 
	  if (NDIM > 0) then
	  
!	    Dimerisation included:							     
	    DO I=1,NC
	    
	      ND = MDIM(I,I)									
! 	      IF (ND == 0) GO TO 90
	      if (ND > 0) then
	      
	        SF(I) = DEXP(BF(I,I)*P/(RG*TEMP))					     
	        EQK(I,I) = (BF(I,I)-BB(I,I))/(RG*TEMP)*SF(I)*P			    
	        Z(I) = -(0.5D0/EQK(I,I))+(0.25D0/(EQK(I,I)**2)+1.D0/EQK(I,I))**0.5D0
	        
	        FUG(I) = dlog(Z(I)) + BF(I,I)*P/(RG*TEMP)
            
            else
            
!	        Without dimerisation:
	        FUG(I) = BB(I,I)*P/(RG*TEMP)
	        
	      endif
	      
	    enddo

	  else

!	    Pure component calculation, - without dimerisation for ALL COMPOUNDS			 
	    do i = 1, NC
	      
	      FUG(i) = BB(i,i)*P/(RG*TEMP)
		
	    enddo

	  endif

	endif
1000  RETURN										
	ENDSUBROUTINE PHUB
