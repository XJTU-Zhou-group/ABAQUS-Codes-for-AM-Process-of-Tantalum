      SUBROUTINE UMATHT(U,DUDT,DUDG,FLUX,DFDT,DFDG,
     1 STATEV,TEMP,DTEMP,DTEMDX,TIME,DTIME,PREDEF,DPRED,
     2 CMNAME,NTGRD,NSTATV,PROPS,NPROPS,COORDS,PNEWDT,
     3 NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'

C
      CHARACTER*80 CMNAME
      DIMENSION DUDG(NTGRD),FLUX(NTGRD),DFDT(NTGRD),
     1 DFDG(NTGRD,NTGRD),STATEV(NSTATV),DTEMDX(NTGRD),
     2 TIME(2),PREDEF(1),DPRED(1),PROPS(NPROPS),COORDS(3),
     3 TABLE_K(2,6),TABLE_C(2,7),TABLE_KP(2,6),Mat_ID(9999999,8)

      common /Mat_ID/Mat_ID
      REAL*8 TempatTdT,TMEP,DTEMP,TABLE_K,TABLE_KP,TABLE_C,
     1 asoltemp, aliqtemp, alatht, TempatT, deltu, DTIME, DTDTIME,
     2 Mat_ID


C     Solidus Temperature
      asoltemp = PROPS(1)
C     Liquidus Temperature
      aliqtemp = PROPS(2)
C     Latent Heat of Fusion
      alatht = PROPS(3)


C     Conductivity of Solid
      DATA TABLE_K/57.1969,139.554, 57.2936,315.805, 65.3774,2571.24,
     &66.3789,3160.37, 66.4,3253, 66.9686,3900.6/
C     Conductivity of Powder
      DATA TABLE_KP/8.916114315,139.554, 8.931188353,315.805,
     &10.19132806,2571.24, 10.34744646,3160.37, 66.4,3253, 66.96,3900.6/
C     Specific Heat
      DATA TABLE_C/138693500,322, 158794000,1454, 178894500,2265,
     &204020100,2792, 232914600,3153, 250502500,3281, 384924600,4182/ 
     
      INC = 0
      INC1 = 0
      INC2 = 0

c                    
      TempatTdT = temp+dtemp
      TempatT = temp 
C      DTDTIME = DTEMP/DTIME

      Mat_ID(NOEL,NPT) = statev(1)     
      if (TempatTdT .GT. aliqtemp .AND. DTEMP .LT. 0.0) then
            Mat_ID(NOEL,NPT) = 1
      endif 

      if (TempatTdT .GT. 1300.0) then
            statev(2) = 1
      endif

      statev(1) = Mat_ID(NOEL, NPT)

C     Conductivity of Powder
      IF (statev(1) .EQ. 0) THEN
        IF (TempatTdT .LE. TABLE_KP(2,1)) THEN
          COND = TABLE_KP(1,1)
          DCOND = 0.0d0
        ELSEIF (TempatTdT .GE. TABLE_KP(2,6)) THEN
          COND = TABLE_KP(1,6)
          DCOND = 0.d0
        ELSEIF (TempatTdT .GT. TABLE_KP(2,1) .AND. TempatTdT .LT. TABLE_KP(2,6)) THEN
          DO K1 = 1,5
            TL1 = TABLE_KP(2, K1+1)
            IF (TempatTdT .LT. TL1 .AND. INC .EQ. 0) THEN
              TL0 = TABLE_KP(2,K1)
              DT = TL1-TL0
              C0 = TABLE_KP(1, K1)
              C1 = TABLE_KP(1, K1+1)
              DC = C1-C0
              DCOND = DC/DT
              COND = DCOND*(TempatTdT-TL0)+C0
              INC = 1
            ENDIF
          END DO
        END IF
C     Conductivity of Solid
      ELSEIF (statev(1) .EQ. 1) THEN
        IF (TempatTdT .LE. TABLE_K(2,1)) THEN
          COND = TABLE_K(1,1)
          DCOND = 0.d0
        ELSEIF (TempatTdT .GE. TABLE_K(2,6)) THEN
          COND = TABLE_K(1,6)
          DCOND = 0.d0
        ELSEIF (TempatTdT .GT. TABLE_K(2,1) .AND. TempatTdT .LT. TABLE_K(2,6)) THEN
          DO K1=1,5
            TL1 = TABLE_K(2, K1+1)
            IF (TempatTdT .LT. TL1 .AND. INC1 .EQ. 0) THEN
              TL0 = TABLE_K(2,K1)
              DT = TL1-TL0
              C0 = TABLE_K(1,K1)
              C1 = TABLE_K(1,K1+1)
              DC = C1-C0
              DCOND = DC/DT
              COND = DCOND*(TempatTdT-TL0)+C0
              INC1 = 1
            ENDIF
          END DO
        END IF
      END IF
C     Specific Heat
      IF (TempatTdT .LE. TABLE_C(2,1)) THEN
        SPECHT = TABLE_C(1,1)
      ELSEIF (TempatTdT .GE. TABLE_C(2,7)) THEN
        SPECHT = TABLE_C(1,7)
      ELSEIF (TempatTdT .GT. TABLE_C(2,1) .AND. TempatTdT .LT. TABLE_C(2,7)) THEN
        DO K1=1,6
          TL1 = TABLE_C(2, K1+1)
          IF (TempatTdT .LT. TL1 .AND. INC2 .EQ. 0) THEN
            TL0 = TABLE_C(2, K1)
            DT = TL1-TL0
            C2 = TABLE_C(1, K1)
            C3 = TABLE_C(1, K1+1)
            DCC = C3-C2
            SPECHT = (DCC/DT)*(TempatTdT-TL0)+C2
            INC2 = 1
          ENDIF
        END DO
      END IF

      DUDT = SPECHT
      deltu = DUDT*DTEMP
      
C     
      ulatn1 = 0.0d0
      ulatn2 = 0.0d0
      ulatnp = 0.0d0
      slope = 0.0d0
      frac = 0.25d0
C
c                    
c     account for latent heat effects
c                                        
      if (TempatT  .gt. asoltemp .and. TempatT .lt. aliqtemp) then
         ulatn1 = (TempatT-asoltemp)*alatht/(aliqtemp-asoltemp)
      else if (TempatT .gt. aliqtemp) then
         ulatn1 = alatht
      end if
c                    
      if (TempatTdT .gt. asoltemp .and. TempatTdT .lt. aliqtemp) then
         ulatn2 = (TempatTdT-asoltemp)*alatht/(aliqtemp-asoltemp)
         slope = alatht/(aliqtemp-asoltemp)
      else if (TempatTdT .gt. aliqtemp) then
         ulatn2 = alatht
         slope = 0.0d0
      end if
c                    
      if (ulatn2 .ne. ulatn1) then
         deltu = deltu+ulatn2-ulatn1
         dudt = dudt+slope
         if (slope .eq. 0.d0) then
            tempp = TempatTdT-frac*dtemp
            if (tempp .gt. asoltemp .and. tempp .lt. aliqtemp) then
               ulatnp = (tempp-asoltemp)*alatht/(aliqtemp-asoltemp)
               slope = alatht/(aliqtemp-asoltemp)
            else if (tempp .gt. aliqtemp) then
               ulatnp = alatht
               slope=0.0d0
            end if
c                          
            if (ulatnp .ne. ulatn2) then
               dudt = dudt+slope
            end if
         end if
      end if
c      

      
      U = U+deltu

C     Heat FLUX
      DO I=1, NTGRD
        FLUX(I) = -COND*DTEMDX(I)
        DFDG(I, I) = -COND
        DFDT(I) = -DCOND*DTEMDX(I)
      END DO

      RETURN
      END


      SUBROUTINE DFLUX(FLUX,SOL,JSTEP,JINC,TIME,NOEL,NPT,COORDS,JLTYP,
     1 TEMP,PRESS,SNAME)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION COORDS(3),FLUX(2),TIME(2)
      CHARACTER*80 SNAME


      real*8 t,pi,phi,eta,P,Q,v,S,x,y,z,x0,y0,z0,Hs,Iz,a,b,c,ff
      t=time(2)
      pi=3.1415926
      eta = 0.8
c      P = 10000.0
      P = 600000.0
      v = 1000
C     defining position parameters
      x=coords(1)
      y=coords(2)
      z=coords(3)
      Q = P*eta
      a = 0.25
      b = 0.2
      cf = 0.25
      cr = 0.25
      ff = 1
      fr = 1
      ox = 1

      IF (t<10.002000 .AND. t>10.000000) THEN
        x0 = 1.000000 + v*(t-10.000000)
        y0 =1.000000
        z0 =2.05
      ELSE IF (t<10.004000.AND. t>10.002000) THEN
        x0 = 3.000000 - v*(t-10.002000)
        y0 =1.200000
        z0 =2.05
      ELSE IF (t<10.006000 .AND. t>10.004000) THEN
        x0 = 1.000000 + v*(t-10.004000)
        y0 =1.400000
        z0 =2.05
      ELSE IF (t<10.008000.AND. t>10.006000) THEN
        x0 = 3.000000 - v*(t-10.006000)
        y0 =1.600000
        z0 =2.05
      ELSE IF (t<10.010000 .AND. t>10.008000) THEN
        x0 = 1.000000 + v*(t-10.008000)
        y0 =1.800000
        z0 =2.05
      ELSE IF (t<10.012000.AND. t>10.010000) THEN
        x0 = 3.000000 - v*(t-10.010000)
        y0 =2.000000
        z0 =2.05
      ELSE IF (t<10.014000 .AND. t>10.012000) THEN
        x0 = 1.000000 + v*(t-10.012000)
        y0 =2.200000
        z0 =2.05
      ELSE IF (t<10.016000.AND. t>10.014000) THEN
        x0 = 3.000000 - v*(t-10.014000)
        y0 =2.400000
        z0 =2.05
      ELSE IF (t<10.018000 .AND. t>10.016000) THEN
        x0 = 1.000000 + v*(t-10.016000)
        y0 =2.600000
        z0 =2.05
      ELSE IF (t<10.020000.AND. t>10.018000) THEN
        x0 = 3.000000 - v*(t-10.018000)
        y0 =2.800000
        z0 =2.05
      ELSE IF (t<20.022000 .AND. t>20.020000) THEN
        x0 = 1.000000 + v*(t-20.020000)
        y0 =1.000000
        z0 =2.10
      ELSE IF (t<20.024000.AND. t>20.022000) THEN
        x0 = 3.000000 - v*(t-20.022000)
        y0 =1.200000
        z0 =2.10
      ELSE IF (t<20.026000 .AND. t>20.024000) THEN
        x0 = 1.000000 + v*(t-20.024000)
        y0 =1.400000
        z0 =2.10
      ELSE IF (t<20.028000.AND. t>20.026000) THEN
        x0 = 3.000000 - v*(t-20.026000)
        y0 =1.600000
        z0 =2.10
      ELSE IF (t<20.030000 .AND. t>20.028000) THEN
        x0 = 1.000000 + v*(t-20.028000)
        y0 =1.800000
        z0 =2.10
      ELSE IF (t<20.032000.AND. t>20.030000) THEN
        x0 = 3.000000 - v*(t-20.030000)
        y0 =2.000000
        z0 =2.10
      ELSE IF (t<20.034000 .AND. t>20.032000) THEN
        x0 = 1.000000 + v*(t-20.032000)
        y0 =2.200000
        z0 =2.10
      ELSE IF (t<20.036000.AND. t>20.034000) THEN
        x0 = 3.000000 - v*(t-20.034000)
        y0 =2.400000
        z0 =2.10
      ELSE IF (t<20.038000 .AND. t>20.036000) THEN
        x0 = 1.000000 + v*(t-20.036000)
        y0 =2.600000
        z0 =2.10
      ELSE IF (t<20.040000.AND. t>20.038000) THEN
        x0 = 3.000000 - v*(t-20.038000)
        y0 =2.800000
        z0 =2.10
      ELSE
        ox =0
      END IF

      heatf = 6.0*sqrt(3.0)*ff*Q/(a*b*cf*pi*sqrt(pi))
      shapef = exp(-3.0*(x-x0)**2/cf**2-3.0*(y-y0)**2/a**2-3.0*(z-z0)**2/b**2)
      heatr = 6.0*sqrt(3.0)*fr*Q/(a*b*cr*pi*sqrt(pi))
      shaper = exp(-3.0*(x-x0)**2/cr**2-3.0*(y-y0)**2/a**2-3.0*(z-z0)**2/b**2)
      IF (x >= x0) THEN
         flux(1) = heatf*shapef*ox
      ELSE IF (x < x0) THEN
         flux(1) = heatr*shaper*ox
      END IF
         FLUX(2) = 0
      RETURN
      END
