	SUBROUTINE SVIP (NC, NMIX, TEMP, BB, BF)	
!	***************************
!	Calculation of pure component and cross virial coefficients		 
!	for multicomponent systems at temperature temp	
!	from Hayden and O'Connell
!	IEC PROC.DES.DEV.14(3)209(1975)
!
!	
!	
!	 NC = NUMBER OF COMPONENT
!	 BF = BFREE  
!	 BB = BTOTAL 
!	 Common block /VIRDAT/ contains
!		RD = MEAN RADIUS OF GYRATION IN A
!		DMU = DIPOLE MOMENT IN DEBYE
!		ETA(1) AND ETA(2) = ASSOCIATION PARAMETERS (PURE COMPONENTS)
!		ETA(I,J) = SOLVATION PARAMETER (CROSS INTERACTION)
	
!	
!	 For given common values of virdat the subroutine will
!	 return values of Bfree and Btotal
!	 ***************************
!	
	implicit real(8) (a-h, o-z) !Original version was single precision
	
	integer, parameter :: NCM = 30	
	
	DIMENSION BF(NCM,NCM),BB(NCM,NCM),W(NCM,NCM),EPSI(NCM,NCM),SIGM3(NCM,NCM),RDMU(NCM,NCM),RDMM(NCM,NCM),A(NCM,NCM), &
&	          DELH(NCM,NCM), D(NCM,NCM),B0(NCM,NCM)									
	COMMON /VIRDAT/       RD(NCM), DMU(NCM), ETA(NCM,NCM)
	COMMON /CRIT/         TC(NCM), PC(NCM), omega(NCM)

	
!	The equation numbers refer to the article by Hayden and O'Connell	 
!
!	Calculation of component parameters	
!	EQ'S 15,30,17,24,25,23,21,22,10		 
!	
	DO 101 I = 1, NC											 		  
	W(I,I)=0.006*RD(I)+0.02087*RD(I)**2-0.00136*RD(I)**3		  
	EPSI(I,I)=TC(I)*(0.748+0.91*W(I,I)-0.4*ETA(I,I)/(2.+20.*W(I,I)))  
	SIGM3(I,I)=(2.44-W(I,I))**3*(TC(I)/PC(I))				 
	IF(DMU(I)-1.45)102,102,103 
103	PN=16.+400.*W(I,I)
	C=2.882-1.882*W(I,I)/(0.03+W(I,I))
! 	write(*,*) C
	XI=DMU(I)**4/(C*EPSI(I,I)*(SIGM3(I,I)**2)*TC(I)*5.723E-8)
! 	write(*, '( 5G)') xi, epsi(i,i), sigm3(i,i), Tc(i), DMU(i)
	
	PPN=PN/(PN-6.)
	EPSI(I,I)=EPSI(I,I)*(1.-XI*PPN+PPN*(PPN+1.)*(XI**2)/2.)		  
	SIGM3(I,I)=SIGM3(I,I)*(1.+3.*XI/(PN-6.))	
102 	RDMU(I,I)=(DMU(I)**2)*7243.8/(EPSI(I,I)*SIGM3(I,I))			
	IF (RDMU(I,I)-0.04) 801,802,802
801	RDMM(I,I)=RDMU(I,I)
	GO TO 850
802	IF(RDMU(I,I)-0.25)803,804,804
803	RDMM(I,I)=0.00
	GO TO 850	
804	RDMM(I,I)=RDMU(I,I)-0.25
850	CONTINUE
!	
!	Last parameters (pure component), eq's 7,8,9,29
!
	B0(I,I)=1.2618*SIGM3(I,I)
	A(I,I)=-0.3-0.05*RDMU(I,I)							 
	DELH(I,I)=1.99+0.2*RDMU(I,I)**2						  
	IF (ETA(I,I)-4.0)806,806,807						  
806	D(I,I)=650./(EPSI(I,I)+300.)						  
	GO TO 805										
807	D(I,I)=42800./(EPSI(I,I)+22400.)						 
805	CONTINUE										 
!	 											  
!	Calculation of virial coefficients (pure components), eq's 14  	 											  
!		13,26,6,29									  
!	 											  
	TSTR=EPSI(I,I)/TEMP-1.6*W(I,I)						
	BFN=0.94-1.47*TSTR-0.85*TSTR**2+1.015*TSTR**3					!(13)/b0
	BFP=(0.75-3.*TSTR+2.1*TSTR**2+2.1*TSTR**3)*RDMM(I,I)				!(26)/b0 (Bfnp - Bf)/b0 
	BF(I,I)=(BFN-BFP)*B0(I,I)							  
	BB(I,I)=BF(I,I)+B0(I,I)*A(I,I)*EXP(DELH(I,I)*EPSI(I,I)/TEMP)	
	IF (ETA(I,I))101,101,809							
809	BCHEM=B0(I,I)*EXP(ETA(I,I)*(D(I,I)-4.27))*(1.-EXP(1500.*ETA(I,I)/TEMP)) !   /21 !******								
! 	write(*,*) bchem
	BB(I,I)=BB(I,I)+BCHEM
101	CONTINUE
	IF (NMIX.EQ.0) GO TO 991							
	IF (NC.GT.1) GO TO 400							  
	GO TO 991										
400	NS=NC-1										  
	DO 990 I=1,NS									  
	K = I + 1										
	DO 990 J=K, NC
	
!		Parameters for mixture calculation					  
!		nonpolar-nonpolar, eq's 32,3NCM,NCM4					 
!	 											  
401	EPSI(I,J)=0.7* SQRT(EPSI(I,I)*EPSI(J,J))+0.60/(1./EPSI(I,I)+1./EPSI(J,J))										  
	SIGM3(I,J)= SQRT(SIGM3(I,I)*SIGM3(J,J))					
	W(I,J)=0.5*(W(I,I)+W(J,J))							 
	IF(DMU(I)*DMU(J))500,501,500						  
!	 
!		Polar-nonpolar, eq's 38,24,36,37
!	 
501	IF(DMU(I)+DMU(J)-2.)500,500,19						
19	XI38=( DMU(I)**2 * (EPSI(J,J)**2*SIGM3(J,J))**(1./3.) * SIGM3(J,J) + &
&	      DMU(J)**2*(EPSI(I,I)**2* SIGM3(I,I))**(1./3.)*SIGM3(I,I) )/(EPSI(I,J)*SIGM3(I,J)**2)									 
	PN=16. + 400.*W(I,J)								
	EPSI(I,J)=EPSI(I,J)*(1.+XI38*PN/(PN-6.))				  
	SIGM3(I,J)=SIGM3(I,J)*(1.-3.*XI38/(PN-6.))				
!	 											  
!		POLAR-POLAR, EQ'S 35,27							 
!	 											  
500	RDMU(I,J)=7243.8*DMU(I)*DMU(J)/(EPSI(I,J)*SIGM3(I,J))		 
	IF(RDMU(I,J)-0.04)14,15,15							 
14	RDMM(I,J)=RDMU(I,J)								  
	GO TO 600										
15	IF(RDMU(I,J)-0.25)16,17,17							 
16	RDMM(I,J)=0.									
	GO TO 600										
17	RDMM(I,J)=RDMU(I,J)-0.25							
600	CONTINUE	
!
!		Last parameters, eq's 7,8,9,29						
!	 											  
!	 											  
	B0(I,J)=1.2618*SIGM3(I,J)							  
	A(I,J)=-0.3-0.05*RDMU(I,J)							 
	DELH(I,J)=1.99+0.2*RDMU(I,J)**2						  
	IF(ETA(I,J)-4.)604,604,605							 
604	D(I,J)=650./(EPSI(I,J)+300.)						  
	GO TO 609										
605	D(I,J)=42800./(EPSI(I,J)+22400.)						 
609	CONTINUE										 
!	
!	Calculation of virial coefficients, eq's 14,13,26,6,29		
!	  
	TSTR=EPSI(I,J)/TEMP-1.6*W(I,J)						
	BFN=0.94-1.47*TSTR-0.85*TSTR**2+1.015*TSTR**3				
	BFP=(0.75-3.*TSTR+2.1*TSTR**2+2.1*TSTR**3)*RDMM(I,J)		  
	BF(I,J)=(BFN-BFP)*B0(I,J)							  
	BF(J,I)=BF(I,J)									
	BB(I,J)=BF(I,J)+B0(I,J)*A(I,J)* EXP(DELH(I,J)*EPSI(I,J)/TEMP)	  
	BB(J,I)=BB(I,J)									
	IF(ETA(I,J))651,651,653							 
653	BCHEM=B0(I,J)*EXP(ETA(I,J)*(D(I,J)-4.27))*(1.-EXP(1500.*ETA(I,J)/TEMP))										 
	BB(I,J)=BB(I,J)+BCHEM								
	BB(J,I)=BB(I,J)									
651	CONTINUE										 
990  	CONTINUE										 
991 	RETURN										
	END											 
! *******************************************************************
