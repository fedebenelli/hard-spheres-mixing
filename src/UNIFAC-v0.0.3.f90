!
!  *********************************************************************
	SUBROUTINE PARINL (NC, NG, MODEL, inputFile, outputFile)				  
!  *********************************************************************

	IMPLICIT REAL*8(A-H,O-Z)
	integer, parameter :: NCM = 30, NGMax = 30
	
	integer, intent(in) :: inputFile, NC, MODEL, outputFile
	integer, intent(inout) :: NG
	
	COMMON /UNIF/ QT(NGMax,NGMax),TAU(NGMax,NGMax),S(NGMax,NGMax),F(NGMax),Q(NGMax),R(NGMax),P(NGMax,NGMax)
	
	DIMENSION RT(NGMax,NGmax),A(32,32),NGM(NGMax),MAINSG(57)			 
	DIMENSION MS(NCM,NGMax,2),NY(NCM,2*NGMax), JH(57), IH(2*NGMax)			   
	REAL*4 RR(57), QQ(57)							    
	REAL*4    A1(32),A2(32),A3(32),A4(32),A5(32),A6(32),A7(32),A8(32),&
     &A9(32),A10(32),A11(32),A12(32),A13(32),A14(32),A15(32),A16(32),A17&
     &(32),A18(32),A19(32),A20(32),A21(32),A22(32),A23(32),A24(32),A25(3&
     &2),A26(32),A27(32),A28(32),A29(32),A30(32),A31(32),A32(32)        
      DATA MAINSG/4*1,4*2,2*3,3*4,5,6,7,8,9,2*10,11,12,2*13,2*14,4*15,3*&
     &16,3*17,2*18,19,20,2*21,22,3*23,24,25,26,3*27,28,29,30,31,32/     
      DATA RR/.9011,.6744,.4469,.2195,1.3454,1.1167,.8886,1.1173,.5313,.&
     &3652,1.2663,1.0396,.8121,1.,3.2499,3.2491,.92,.8952,1.6724,1.4457,&
     &.998,3.168,1.3013,1.528,1.9031,1.6764,1.145,.9183,.6908,.9183,1.46&
     &54,1.238,1.006,2.2564,2.0606,1.8016,2.87,2.6401,3.39,1.1562,1.8701&
     &,1.6434,1.06,2.0086,1.7818,1.5544,1.4199,2.4088,4.0013,2.9993,2.83&
     &32,2.667,3.3092,2.4317,3.0856,4.0358,2.8266/   	
      DATA QQ/.848,.54,.228,0.,1.176,.867,.676,.988,.4,.12,.968,.66,.348&
     &,1.2,3.128,3.124,1.4,.68,1.488,1.18,.948,2.484,1.224,1.532,1.728,1&
     &.42,1.088,.78,.468,1.1,1.264,.952,.724,1.988,1.684,1.448,2.41,2.18&
     &4,2.91,.844,1.724,1.416,.816,1.868,1.56,1.248,1.104,2.248,3.568,2.&
     &113,1.833,1.553,2.86,2.192,2.736,3.2,2.472/			     
	DATA A1/0.,292.3,156.5,104.4,328.2,-136.7,-131.9,342.4,-159.8,66.5&
     &6,146.1,14.78,1744.,-320.1,1571.,73.8,27.9,21.23,89.97,-59.06,29.0&
     &8,175.8,94.34,193.6,108.5,81.49,-128.8,147.3,-11.91,14.91,67.84,36&
     &.42/										  
	DATA A2/74.54,0.,-94.78,-269.7,470.7,-135.7,-135.7,220.6,1.,306.1,&
     &517.,1.,-48.52,485.6,76.44,-24.36,-52.71,-185.1,-293.7,1.,34.78,1.&
     &,375.4,5*1.,176.7,132.1,42.73,60.82/					
	DATA A3/-114.8,340.7,0.,-146.8,-9.21,-223.,-252.,372.8,-473.2,-78.&
     &31,-75.3,-10.44,75.49,114.8,52.13,4.68,1.,288.5,-4.7,777.8,56.41,-&
     &218.9,113.6,7.18,247.3,-50.71,-255.3,1.,-80.48,-17.78,59.16,29.77/
	DATA A4/-115.7,4102.,167.,0.,1.27,-162.6,-273.6,203.7,-470.4,-73.8&
     &7,223.2,-184.9,147.3,-170.,65.69,122.9,1.,33.61,134.7,-47.13,-53.2&
     &9,-15.41,-97.05,-127.1,453.4,-30.28,-124.6,3*1.,26.59,55.97	/
	DATA A5/644.6,724.4,703.9,4000.,0.,-281.1,-268.8,-122.4,-63.15,216&
     &.,-431.3,444.7,118.4,180.6,137.1,455.1,669.2,418.4,713.5,1989.,201&
     &1.,529.,483.8,332.6,-289.3,-99.56,-319.2,837.9,4*1.		  /
	DATA A6/329.6,1731.,511.5,136.6,937.3,2*0.,247.,-547.,401.7,643.4,&
     &-94.64,728.7,-76.64,-218.1,351.5,-186.1,-465.7,-260.3,3*1.,264.7,9&
     &*1./										  
	DATA A7/310.7,1731.,577.3,906.8,991.3,2*0.,104.9,-547.2,-127.6,231&
     &.4,732.3,349.1,-152.8,-218.1,351.5,-401.6,-465.7,512.2,3*1.,264.7,&
     &9*1./										 
	DATA A8/1300.,896.,859.4,5695.,28.73,-61.29,5.89,0.,-595.9,634.8,6&
     &23.7,211.6,652.3,385.9,212.8,770.,740.4,793.2,1205.,390.7,63.48,-2&
     &39.8,13.32,439.9,-424.3,1.,203.,1153.,-311.,-262.6,1.11,1.	 /
	DATA A9/2255.,1.,1649.,292.6,-195.5,-153.2,-153.2,344.5,0.,-568.,3&
     &*1.,-337.3,4*1.,1616.,2*1.,-860.3,1.,-230.4,523.,1.,-222.7,5*1.  /
	DATA A10/472.6,343.7,593.7,916.7,67.07,-47.41,353.8,-171.8,-825.7,&
     &0.,128.,48.93,-101.3,58.84,52.38,483.9,550.6,342.2,550.,190.5,-349&
     &.2,857.7,377.,211.6,82.77,2*1.,417.4,4*1.	    /		  
	DATA A11/158.1,-214.7,362.3,1218.,1409.,-344.1,-338.6,-349.9,1.,-3&
     &7.36,0.,-311.6,1051.,1090.,1.,-47.51,16*1./			     
	DATA A12/383.,1.,31.14,715.6,-140.3,299.3,-241.8,66.95,1.,120.3,17&
     &24.,0.,-115.7,-46.13,2*1.,808.8,203.1,70.14,5*1.,-75.23,1.,-201.9,&
     &123.2,1.,-281.9,2*1./							   
	DATA A13/139.4,1647.,461.8,339.1,-104.,244.4,-57.98,-465.7,1.,1247&
     &.,.75,1919.,0.,1417.,1402.,337.1,437.7,370.4,438.1,1349.,1.,681.4,&
     &152.4,1.,-1707.,2*1.,639.7,4*1.	    /				
	DATA A14/972.4,-577.5,6.,5688.,195.6,19.57,487.1,-6.32,-898.3,258.&
     &70,-245.8,57.7,-117.6,0.,461.3,1.,-132.9,176.5,129.5,-246.3,2.41,3&
     &*1.,29.86,7*1.				     /			    
	DATA A15/662.1,289.3,32.14,213.1,262.5,1970.,1970.,64.42,1.,5.202,&
     &2*1.,-96.62,-235.7,0.,225.4,-197.7,-20.93,113.9,3*1.,-94.49,9*1. /
	DATA A16/42.14,99.61,-18.81,-114.1,62.05,-166.4,-166.4,315.9,1.,10&
     &00.,751.8,1.,19.77,1.,301.1,0.,-21.35,-157.1,11.8,13*1.	 /   
	DATA A17/-243.9,337.1,2*1.,272.2,128.6,507.8,370.7,1.,-301.,1.,-34&
     &7.9,1670.,108.9,137.8,110.5,0.,1.,17.97,13*1.	   /	     
	DATA A18/7.5,4583.,-231.9,-12.14,-61.57,2*1544.,356.8,1.,12.01,1.,&
     &-249.3,48.15,-209.7,-154.3,249.2,1.,0.,51.9,1.,-15.62,-216.3,4*1.,&
     &-114.7,5*1. /								     
	DATA A19/-5.55,5831.,3000.,-141.3,-41.75,224.6,-207.,502.9,4894.,-&
     &10.88,1.,61.59,43.83,54.57,47.67,62.42,56.33,-30.1,0.,-255.4,-54.8&
     &6,8455.,-34.68,514.6,8*1. /						   
	DATA A20/924.8,1.,-878.1,-107.3,-597.1,2*1.,-97.27,1.,902.6,2*1.,8&
     &74.3,629.,4*1.,475.8,0.,-465.2,1.,794.4,1.,-241.7,1.,-906.5,5*1. /
	DATA A21/696.8,405.9,29.13,1208.,-189.3,2*1.,198.3,1.,430.6,3*1.,-&
     &149.2,3*1.,70.04,492.,346.2,0.,5*1.,-169.7,5*1.	  /	    
	DATA A22/902.2,1.,1.64,689.6,-348.2,1.,1.,-109.8,-851.6,1010.,2*1.&
     &,942.2,4*1.,-75.5,1302.,2*1.,0.,1.,175.8,164.4,1.,-944.9,5*1./    
	DATA A23/556.7,425.7,-1.77,3629.,-30.7,150.8,150.8,1538.6,1.,400.,&
     &2*1.,446.3,1.,95.18,3*1.,490.9,-154.5,2*1.,0.,1.,481.3,7*1.	/
	DATA A24/575.7,1.,-11.19,-175.6,-159.,2*1.,32.92,-16.13,-328.6,8*1&
     &.,534.7,2*1.,179.9,1.,0.,-246.,7*1.	     /			 
	DATA A25/527.5,1.,358.9,337.7,536.6,2*1.,-269.2,-538.6,211.6,1.,-2&
     &78.2,572.7,343.1,5*1.,124.8,1.,125.3,139.8,963.,0.,7*1.    /	
	DATA A26/269.2,1.,363.5,1023.,53.37,20*1.,0.,6*1.	    /	
	DATA A27/-300.,1.,-578.2,-390.7,183.3,2*1.,-873.6,-637.3,2*1.,-208&
     &.4,5*1.,18.98,1.,-387.7,134.3,924.5,4*1.,0.,5*1.     /		
	DATA A28/-63.6,3*1.,-44.44,2*1.,1429.,1.,148.,1.,-13.91,-2.16,14*1&
     &.,0.,4*1./									  
	DATA A29/928.3,500.7,364.2,4*1.,-364.2,20*1.,0.,3*1.    /	   
	DATA A30/331.,115.4,-58.1,4*1.,-117.4,3*1.,173.8,17*1.,0.,2*1.   /
	DATA A31/561.4,784.4,21.97,238.,3*1.,18.41,22*1.,0.,1.	 /    
	DATA A32/956.5,265.4,84.16,132.2,27*1.,0./   

	CHARACTER(10) :: compoundName(NCM)
	COMMON /NAME/    compoundName
	
	iout = 1
	DO 5 I=1,32		
	
	  A(I,1)=A1(I)	
	  A(I,2)=A2(I)	
	  A(I,3)=A3(I)	
	  A(I,4)=A4(I)	
	  A(I,5)=A5(I)	
	  A(I,6)=A6(I)	
	  A(I,7)=A7(I)	
	  A(I,8)=A8(I)	
	  A(I,9)=A9(I)	
	  A(I,10)=A10(I)	
	  A(I,11)=A11(I)	
	  A(I,12)=A12(I)	
	  A(I,13)=A13(I)	
	  A(I,14)=A14(I)	
	  A(I,15)=A15(I)	
	  A(I,16)=A16(I)	
	  A(I,17)=A17(I)	
	  A(I,18)=A18(I)	
	  A(I,19)=A19(I)	
	  A(I,20)=A20(I)	
	  A(I,21)=A21(I)	
	  A(I,22)=A22(I)	
	  A(I,23)=A23(I)	
	  A(I,24)=A24(I)	
	  A(I,25)=A25(I)	
	  A(I,26)=A26(I)	
	  A(I,27)=A27(I)	
	  A(I,28)=A28(I)	
	  A(I,29)=A29(I)	
	  A(I,30)=A30(I)	
	  A(I,31)=A31(I)
	  A(I,32)=A32(I)
	  
5	enddo  		
	if (IOUT == 0) IOUT=6		
	NK=NC					
	DO 15 I = 1, NGMax
	  DO 150 J = 1, NGMax			
	    P(I,J)=0.D0			
	    QT(I,J)=0.D0			
	    RT(I,J)=0.D0			
!    15 RT(I,J)=0.D0			
150	  enddo
15	enddo
	if (MODEL == 1) then ! /= 1) GOTO 19
	  
!	  UNIQUAC model 
	  NG = NK										 
	  DO 16 I = 1, NK									

	    read (inputFile, *) RT(I,I), QT(I,I), (P(I,J), J = 1, NK)			     

16	  enddo

	endif
19	CONTINUE									    
! 	if (MODEL == 1) GOTO 21							  
	if (model /=1) then
	  
! 	  UNIFAC model

	  read (inputFile, *) NG	  
	  read (inputFile, *) IOWNRQ, IOWNP	!# of parameters from user	
! 	  if (IOWNRQ == 0) GOTO 10							 
	  if (IOWNRQ > 0) then	    
	    DO 6 I=1,IOWNRQ
	    
	      read (inputFile, *) K, RR(K), QQ(K)
		
6	    enddo
	  endif

!10	  if (IOWNP == 0) GOTO 14							  
	  if (IOWNP > 0) then	    
	    DO 11 I = 1, IOWNP 
	    
	      read (inputFile, *) J, K, P(J,K), P(k,j)							  
	
11	    enddo
14	  endif
! 	  DO 48 I=1,NK									
! 	    DO 48 J=1,10									
! 	      DO 48 K=1,2									 
! 48	        MS(I,J,K)=0									 	  
! 	  DO 49 I=1,57									
! 49	    JH(I)=0	
	  MS(:NK,:NGMax,:2) = 0	! ATENCIÓN!!! :10 -> :NG?	   
	  JH(:57) = 0		! ATENCIÓN!!!
	
!    
!	  Read group contribution per compound
	  DO 50 I = 1, NK							
	  
	    read (inputFile, *) (MS(I,J,1), MS(I,J,2), J = 1, NGMax)	!10				  
	  
50	  enddo
	  IC = 1										  
	  DO 71 I = 1, NK									
	    DO 70 J = 1, NGMax !10									
	      if (MS(I,J,1)  == 0) GOTO 71						    
		  
	        IH(IC) = MS(I,J,1)  
		  
	        if (IC == 1) GOTO 69
	        
	          if (IH(IC) == IH(IC-1)) GOTO 70						
	          
	          if (IH(IC) > IH(IC-1)) GOTO 69						
	          
	          if (IC > 2) GOTO 55							     
	          
	            IHH=IH(1)									   
	            IH(1)=IH(2)									 
	            IH(2)=IHH	
	          
	            GOTO 69 
	            
55	          I1 = IC - 1									     
	          DO 65 I2 = 1, I1								     
	          
	            if (IH(IC) > IH(I2)) GOTO 65	
			  
	            if (IH(IC) == IH(I2)) GOTO 70						  
			  
	              I4 = IC - I2									    
	              DO 61 I3 = 1, I4	
	              
61	                IH(IC+1-I3) = IH(IC-I3)							   
                    
	              IH(I2)=MS(I,J,1)	
			  
65	          CONTINUE

69	        IC=IC+1									     
	        if (IC > 20) then
	          WRITE(*,607)	  
	          write (outputFile, 607)	
	        endif
607	        FORMAT (" ** WARNING: number of sub groups must not exceed 20 **")
		    
70 	    continue									    

71 	  continue
	  IC = IC - 1									     
	  DO 73 I = 1, IC									
	    JH(IH(I))=I								 
73	  enddo
! 	  DO 72 I=1,10									
! 	    DO 72 J=1,20									
! 72	      NY(I,J) = 0									   
	  NY(:NC,:NG) = 0
	  DO 75 I=1,NK									
	    DO 74 J=1, NGMax !10									
	    
	      if (MS(I,J,1) == 0) GOTO 75	
		  
	      N1 = MS(I,J,1)									
	      N2 = MS(I,J,2)									
		
	      if (N1 == 0) GOTO 75							     
		  
	      N3 = JH(N1)									   
	
74	    NY(I,N3) = N2									 
75	  CONTINUE									    
	  I = 0										   
	  NGMGL = 0									     
	  DO 80 K = 1, IC
	  
	    NSG = IH(K)									   
	    NGMNY = MAINSG(NSG)								 
	    if (NGMNY /= NGMGL) I = I + 1							
	    NGM(I) = NGMNY									
	    NGMGL = NGMNY									 
	    DO 801 J = 1, NK									
	    
	      RT(I,J) = RT(I,J) + NY(J,K)*RR(NSG)					     
	      QT(I,J) = QT(I,J) + NY(J,K)*QQ(NSG)					     
		
801	    enddo
80	  enddo
	  NG=I										  
	  WRITE (outputFile, '(//," Sub groups: ", 6X, <2*NG>I3)') (IH(K), K = 1, IC)						   
	  WRITE (outputFile, '(" Main groups:", 6X, <2*NG>I3)') (MAINSG(IH(K)), K = 1, IC)					 
	  WRITE (outputFile, '(/, " Component")')	
  								
	  DO 90 I = 1, NK
	  
	    WRITE (outputFile, '(6X, I2, X, A10, <2*NG>I3)') I, compoundName(i), (NY(I,K), K = 1, IC)					     
	    
90	  enddo
	  WRITE (outputFile, 699)									
	  DO 20 I = 1, NG									
	    DO 201 J = 1, NG									
	      NI = NGM(I)									   
	      NJ = NGM(J)									   
	      if (DABS(P(I,J)) > 0.D0) GO TO 201					   
	      
	        P(I,J) = A(NI,NJ)
	  
201	    enddo
20      enddo									    
	  WRITE (outputFile, '(" Group R- ans Q-values",/)')									
	  DO 95 K = 1, IC									
	    NN=IH(K)									    
	    WRITE (outputFile, 613) NN, RR(NN), QQ(NN)						 
95	  enddo
	  WRITE (outputFile, 699)	    									    

21	endif !CONTINUE									    
	WRITE (outputFile, 604)									
	DO 25 I = 1, NG
	
	  WRITE (outputFile, 603) (P(I,J), J = 1, NG)						  
	  
25	enddo
	WRITE (outputFile, 699)									
	if (MODEL == 0) WRITE (outputFile, 605)						   
	if (MODEL == 1) WRITE (outputFile, 627)						   
		    
	DO 30 I=1,NK									
	  Q(I)=0.D0									   
	  R(I)=0.D0									   
	  DO 301 K=1,NG									
	    Q(I)=Q(I)+QT(K,I)								 
	    R(I)=R(I)+RT(K,I)								 
301	  enddo
30	enddo
	DO 40 I=1,NK									
	
	  WRITE (outputFile, 606) I, R(I), Q(I)							
	  
40	enddo								    
501	FORMAT(20I3)									
502	FORMAT(8F10.2)								    
503	FORMAT(I3,2F10.2)								 
504	FORMAT(2I3,F10.2)								 
603	FORMAT(1X,10F12.3)								
604	FORMAT('  Interaction parameters',/)					
605	FORMAT(' UNIFAC molecular R and Q',/)				     
606	FORMAT(I5,2F15.4)								 
		  
613	FORMAT(1X,I3,2F10.4)							    
627	FORMAT(' Specified UNIQUAC R and Q',/)				    
699	FORMAT(//)									  
	RETURN										
	ENDsubroutine PARINL
!
!  ********************************************************************
!
!	Cálculo de parámetros independientes de la composición.
!	Útil sólo cuando T = cte en un extractor...
!
	SUBROUTINE PARAM (NC, NG, T)						
!  ********************************************************************
!
	IMPLICIT REAL*8(A-H,O-Z)
    integer, parameter :: NCM = 30, NGM = 30							
	COMMON/UNIF/QT(NGM,NGM),TAU(NGM,NGM),S(NGM,NGM),F(NGM),Q(NGM),R(NGM),P(NGM,NGM)
	DO 30 I=1,NG									
	  DO 301 J=1,NG									
          TAU(I,J)=DEXP(-P(I,J)/T)							
301	  enddo
30	enddo
	DO 50 I=1,NC									
	DO 50 K=1,NG									
	S(K,I)=0.D0									 
	DO 50 M=1,NG									
   50 S(K,I)=S(K,I)+QT(M,I)*TAU(M,K)						
	DO 60 I=1,NC									
	F(I)=1.D0									   
	DO 60 J=1,NG									
   60 F(I)=F(I)+QT(J,I)*DLOG(S(J,I))						
	RETURN										
	ENDsubroutine PARAM
!
!   *******************************************************************
	SUBROUTINE UNIQUACFAC (NDIF, NACT, NC, NG, T, X, ACT, DACT, TACT)		   
!   *******************************************************************
!
	IMPLICIT REAL*8(A-H,O-Z)
	integer, parameter :: NCM = 30, NGM = 30
	COMMON /UNIF/ QT(NGM,NGM),TAU(NGM,NGM),S(NGM,NGM),F(NGM),Q(NGM),R(NGM),P(NGM,NGM)								    
	DIMENSION     X(NCM),GAM(NCM),ACT(NCM),DACT(NCM,NCM),THETA(NCM),PHI(NCM),RI(NCM),QI(NCM),QIL(NCM),RIL(NCM),QID(NCM),ETAL(NCM),TACT(NCM),U(NGM,NCM),V(NGM,NCM)									    
	DIMENSION     DETA(NGM),DS(NGM,NCM),ETA(NGM),TETAR(NGM),H3(NGM,NCM)	   
	
	THETS = 0.D0									  
	PHS = 0.D0
	call param (NC, NG, T)
	DO 10 I = 1,NC									
	
	  THETA(I) = X(I)*Q(I)								
	  PHI(I) = R(I)*X(I)								  
	  THETS = THETS + THETA(I)							    
	  PHS = PHS + PHI(I)
	  
10	enddo

	DO 20 I=1,NC									
	
	  THETA(I) = THETA(I)/THETS							 
	  PHI(I) = PHI(I)/PHS								 
	  RI(I) = R(I)/PHS								    
	  RIL(I) = DLOG(RI(I))								
	  QI(I) = Q(I)/THETS								  
	  QID(I) = 1.D0 - RI(I)/QI(I)							   
	  QIL(I) = DLOG(QI(I))
	  
20	enddo
	DO 30 I=1,NC									
	
	  XX = F(I) + Q(I)*(1.D0 - QIL(I)) - RI(I) + RIL(I)				     
	  XX = XX - 5*Q(I)*(QID(I) + RIL(I) - QIL(I))					
	  GAM(I) = XX
	  
30	enddo
	DO 40 I=1,NG									
	TETAR(I)=0.									 
	ETA(I)=0.									   
	DO 45 J=1,NC									
	ETA(I)=ETA(I)+S(I,J)*X(J)						     
   45 TETAR(I)=TETAR(I)+QT(I,J)*X(J)						
   40 ETAL(I)=DLOG(ETA(I))							    
	DO 50 I=1,NC									
	  DO 60 J=1,NG									
	    U(J,I)=S(J,I)/ETA(J)							    
	    V(J,I)=U(J,I)*TETAR(J)							  
   60     GAM(I)=GAM(I)-V(J,I)-QT(J,I)*ETAL(J)					
	  ACT(I)=DEXP(GAM(I))							     
   50   if (NACT == 1) ACT(I)=ACT(I)*X(I)					    
	if (NDIF == 0) GOTO 1000							 
! 	IF (NDIF == 2) GO TO 81							 
	DO 70 I=1,NC									
	DO 70 J=I,NC									
	XX=Q(I)*QI(J)*(1.-5.*QID(I)*QID(J))+(1.-RI(I))*(1.-RI(J))	   
	DO 75 K=1,NG									
   75 XX=XX+U(K,I)*(V(K,J)-QT(K,J))-U(K,J)*QT(K,I)			    
	DACT(I,J)=XX									
	DACT(J,I)=XX									
	if (NACT == 1) GOTO 70							   
	DACT(I,J)=DACT(I,J)*ACT(I)						    
	IF (J == I) GO TO 70							    
	DACT(J,I)=DACT(J,I)*ACT(J)						    
   70 CONTINUE									    
	if (NACT == 0) GOTO 81							   
	DO 80 I=1,NC									
	DO 80 J=1,NC									
	DACT(I,J)=ACT(I)*(DACT(I,J)-1.D0)					   
	if (J == I) DACT(I,J)=DACT(I,J)+DEXP(GAM(I))			     
   80 CONTINUE									    
   81 CONTINUE									    
	if (NDIF == 1) GOTO 1000							 
	DO 150 K=1,NG								     
	DETA(K)=0.D0									
	DO 150 I=1,NC								     
	DS(K,I)=0.D0									
	DO 151 M=1,NG								     
	if (QT(M,I).EQ.0.D0) GOTO 151						  
	DS(K,I)=DS(K,I)-QT(M,I)*DLOG(TAU(M,K))*TAU(M,K)/T		     
  151 CONTINUE									    
  150 DETA(K)=DETA(K)+DS(K,I)*X(I)						  
	DO 152 I=1,NC								     
	TACT(I)=0.D0									
	DO 153 K=1,NG								     
	H3(K,I)=(-S(K,I)*DETA(K)/ETA(K)+DS(K,I))/ETA(K)			 
	HH=H3(K,I)*(TETAR(K)-QT(K,I)*ETA(K)/S(K,I))			     
  153 TACT(I)=TACT(I)-HH								
  152 TACT(I)=TACT(I)*ACT(I)							  
 1000 CONTINUE									    
	RETURN										
	ENDsubroutine
