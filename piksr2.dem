	PROGRAM D8R2
C	Driver for routine PIKSR2
	DIMENSION A(100),B(100)
	OPEN(5,FILE='TARRAY.DAT',STATUS='OLD')
	READ(5,*) (A(I),I=1,100)
	CLOSE(5)
C	Generate B-array
	DO 11 I=1,100
		B(I)=I
11	CONTINUE
C	Sort A and mix B
	CALL PIKSR2(100,A,B)
	WRITE(*,*) 'After sorting A and mixing B, array A is:'
	DO 12 I=1,10
		WRITE(*,'(1X,10F6.2)') (A(10*(I-1)+J), J=1,10)
12	CONTINUE
	WRITE(*,*) '...and array B is:'
	DO 13 I=1,10
		WRITE(*,'(1X,10F6.2)') (B(10*(I-1)+J), J=1,10)
13	CONTINUE
	WRITE(*,*) 'press RETURN to continue...'
	READ(*,*)
C	Sort B and mix A
	CALL PIKSR2(100,B,A)
	WRITE(*,*) 'After sorting B and mixing A, array A is:'
	DO 14 I=1,10
		WRITE(*,'(1X,10F6.2)') (A(10*(I-1)+J), J=1,10)
14	CONTINUE
	WRITE(*,*) '...and array B is:'
	DO 15 I=1,10
		WRITE(*,'(1X,10F6.2)') (B(10*(I-1)+J), J=1,10)
15	CONTINUE
	END
