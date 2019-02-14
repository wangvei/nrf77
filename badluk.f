      PROGRAM BADLUK
      DATA IYBEG,IYEND /1900,2000/
      TIMZON=-5./24
      WRITE (*,'(1X,A,I5,A,I5)') 'Full moons on Friday the 13th from',
     *    IYBEG,' to',IYEND
      DO 12 IYYY=IYBEG,IYEND
        DO 11 IM=1,12
          JDAY=JULDAY(IM,13,IYYY)
          IDWK=MOD(JDAY+1,7)
          IF(IDWK.EQ.5) THEN
            N=12.37*(IYYY-1900+(IM-0.5)/12.)
            ICON=0
1           CALL FLMOON(N,2,JD,FRAC)
            IFRAC=NINT(24.*(FRAC+TIMZON))
            IF(IFRAC.LT.0)THEN
              JD=JD-1
              IFRAC=IFRAC+24
            ENDIF
            IF(IFRAC.GT.12)THEN
              JD=JD+1
              IFRAC=IFRAC-12
            ELSE
              IFRAC=IFRAC+12
            ENDIF
            IF(JD.EQ.JDAY)THEN
              WRITE (*,'(/1X,I2,A,I2,A,I4)') IM,'/',13,'/',IYYY
              WRITE (*,'(1X,A,I2,A)') 'Full moon ',IFRAC,
     *            ' hrs after midnight (EST).'
              GOTO 2
            ELSE
              IC=ISIGN(1,JDAY-JD)
              IF(IC.EQ.-ICON) GOTO 2
              ICON=IC
              N=N+IC
            ENDIF
            GOTO 1
2           CONTINUE
          ENDIF
11      CONTINUE
12    CONTINUE
      END
