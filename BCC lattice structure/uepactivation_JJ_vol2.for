      subroutine uepactivationvol(
     * lFlags,
     * epaName,
     * noel,
     * nElemNodes,
     * iElemNodes,
     * mcrd,
     * coordNodes,
     * uNodes,
     * kstep,
     * kinc,
     * time,
     * dtime,
     * temp,
     * npredef,
     * predef,
     * nsvars,
     * svars,
     * sol,
     * solinc,
     * volFract,
     * nVolumeAddEvents,
     * volFractAdded,
     * csiAdded)     
C
      include 'aba_param.inc'
C
      dimension   
     * lFlags(*),
     * iElemNodes(nElemNodes),
     * coordNodes(mcrd,nElemNodes),
     * uNodes(mcrd,nElemNodes),
     * time(2),
     * svars(2,nsvars),
     * temp(2,nElemNodes),
     * predef(2,npredef,nElemNodes),
     * sol(nElemNodes),
     * solinc(nElemNodes),
     * volFract(*),
     * volFractAdded(*),
     * csiAdded(3,*)

      character*80 epaName

c      user coding to define material volume fraction added.
      t = time(2)
C      IF (t .GT. 0.8) THEN
C      volFractAdded(1)=1.0
c      ELSE IF (t .LT. 0.7 .AND. t .GT. 0.5) THEN
c      volFractAdded(1)=1.0
c      ELSE
c      volFractAdded(1)=0.0      
C      END IF

C      x1 = coordNodes(3,1)
C      x2 = coordNodes(1,1)
C      x3 = coordNodes(2,1)
      tp = temp(1,nElemNodes)
      dt = temp(2,nElemNodes)
      volFractAdded(1)=0.0
      IF (tp .GT. 3250 .AND.  tp .LT. 3300 .AND. dt .LT. 0.0) THEN
C      IF (tp .GT. 3250 .AND. dt .LT. 0.0) THEN
C        IF (x3 .EQ. 2.05) THEN
           volFractAdded(1)=1.0
C        ELSE IF (x3 .EQ. 2.1) THEN
C           volFractAdded(1)=1.0      
C        END IF
      END IF
      return
      end  