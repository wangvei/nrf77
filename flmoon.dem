	PROGRAM D1R1
C	Driver for routine FLMOON
	PARAMETER(TZONE=5.0)
	CHARACTER PHASE(4)*15,TIMSTR(2)*3
	DATA PHASE/'new moon','first quarter',
     *			'full moon','last quarter'/
	DATA TIMSTR/' AM',' PM'/
	WRITE(*,*) 'Date of the next few phases of the moon'
	WRITE(*,*) 'Enter today''s date (e.g. 1,31,1982)'
	TIMZON=-TZONE/24.0
	READ(*,*) IM,ID,IY
C	Approximate number of full moons since January 1900
	N=12.37*(IY-1900+(IM-0.5)/12.0)
	NPH=2
	J1=JULDAY(IM,ID,IY)
	CALL FLMOON(N,NPH,J2,FRAC)
	N=N+(J1-J2)/28.0
	WRITE(*,'(/1X,T6,A,T19,A,T32,A)') 'Date','Time(EST)','Phase'
	DO 11 I=1,20
		CALL FLMOON(N,NPH,J2,FRAC)
		IFRAC=NINT(24.*(FRAC+TIMZON))
		IF (IFRAC.LT.0) THEN
			J2=J2-1
			IFRAC=IFRAC+24
		ENDIF
		IF (IFRAC.GE.12) THEN
			J2=J2+1
			IFRAC=IFRAC-12
		ELSE
			IFRAC=IFRAC+12
		ENDIF
		IF (IFRAC.GT.12) THEN
			IFRAC=IFRAC-12
			ISTR=2
		ELSE
			ISTR=1
		ENDIF
		CALL CALDAT(J2,IM,ID,IY)
		WRITE(*,'(1X,2I3,I5,T20,I2,A,5X,A)') IM,ID,IY,
     *			IFRAC,TIMSTR(ISTR),PHASE(NPH+1)
		IF (NPH.EQ.3) THEN
			NPH=0
			N=N+1
		ELSE
			NPH=NPH+1
		ENDIF
11	CONTINUE
	END
