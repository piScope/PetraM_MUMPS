C
C  This file is read destirbuted Z matrix
C
      PROGRAM MUMPS_TEST
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'zmumps_struc.h'
      TYPE (ZMUMPS_STRUC) mumps_par
      LOGICAL LOAD_RHS, FILE_EXIST1, FILE_EXIST2      
      INTEGER IERR, I
      INTEGER(8) I8
      CHARACTER(LEN=13) FILENAME
      CHARACTER(LEN=16) ROW_T, COL_T, REAL_T, IMAG_T
      CHARACTER(LEN=200) LINE, L
      INTEGER(8) NNZ_LOCAL, NNZ_TOTAL, COUNT, ASIZE, ASIZE_LOCAL
      INTEGER(8) D1, D2, ROW, COL, STAT
      REAL RPART, IPART
      
      CALL MPI_INIT(IERR)
C Define a communicator for the package.
      mumps_par%COMM = MPI_COMM_WORLD
C  Initialize an instance of the package
C  for L U factorization (sym = 0, with working host)
      mumps_par%JOB = -1
      mumps_par%SYM = 0
      mumps_par%PAR = 1
      CALL ZMUMPS(mumps_par)
      IF (mumps_par%INFOG(1).LT.0) THEN
       WRITE(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",
     &            "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1), 
     &            "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2) 
       GOTO 500
      END IF
      
C     Read file (filename is like matrxi.000000)
      WRITE(FILENAME, "(A7, I6.6)") "matrix.", mumps_par%MYID

C     Check if input file exist
      INQUIRE(FILE=FILENAME, EXIST=FILE_EXIST1)
      INQUIRE(FILE= "rhs_0.000000", EXIST=FILE_EXIST2)
      
      IF (FILE_EXIST1) THEN
         write(*,*) "reading:", FILENAME
      ENDIF
      IF (FILE_EXIST2) THEN
         LOAD_RHS = .TRUE.
      ELSE
         LOAD_RHS = .FALSE.
      ENDIF
      

C     Find Number of NNZ      
      OPEN (1, file = FILENAME)
      DO
         READ (1,*, END=10)
         NNZ_LOCAL = NNZ_LOCAL + 1
      END DO
 10   CLOSE (1)
      mumps_par%NNZ_loc = NNZ_LOCAL
      
C     Total NNZ (not used as input)      
      CALL MPI_ALLREDUCE(NNZ_LOCAL, NNZ_TOTAL, 1, MPI_LONG, MPI_SUM,
     &                MPI_COMM_WORLD, ierr)
      
C     Distributed matrix, parmetis      
      mumps_par%ICNTL(18) = 3
C     ParMETIS (28 = 2, 29 = 2)
      mumps_par%ICNTL(28) = 2
      mumps_par%ICNTL(29) = 2
C     OUT-OF-CORE 22 = 1      
      mumps_par%ICNTL(22) = 1
C     Central RHS 20 = 0
      mumps_par%ICNTL(20) = 0

      ALLOCATE( mumps_par%IRN_loc( mumps_par%NNZ_loc ) )
      ALLOCATE( mumps_par%JCN_loc ( mumps_par%NNZ_loc ) )
      ALLOCATE( mumps_par%A_loc( mumps_par%NNZ_loc ) )

      COUNT = 1
      ASIZE = 0
      IF (NNZ_LOCAL .GT. 0) THEN      
       OPEN (1, file = FILENAME)
       DO
         READ (1,"(A)", END=20) LINE
         D1 = 1
         D2 = 1
         DO 
            IF (LINE(D2:D2) .EQ. " ")  THEN
               ROW_T = LINE(D1:D2-1)
               EXIT
            END IF
            D2 = D2 + 1
         END DO

         D1 = D2+1
         D2 = D2+1
         DO 
            IF (LINE(D2:D2) .EQ. " ")  THEN
               COL_T = LINE(D1:D2-1)
               EXIT
            END IF
            D2 = D2 + 1
         END DO

         D1 = D2+1
         D2 = D2+1         
         DO 
            IF (LINE(D2:D2) .EQ. " ")  THEN
               REAL_T = LINE(D1:D2-1)
               EXIT
            END IF
            D2 = D2 + 1
         END DO
         D2 = D2+1                  
         IMAG_T = LINE(D2:)
         
         READ(ROW_T, *, iostat=STAT)  ROW
         READ(COL_T, *, iostat=STAT)  COL         
         READ(REAL_T, *, iostat=STAT) RPART
         READ(IMAG_T, *, iostat=STAT) IPART         
         COL = COL + 1
         ROW = ROW + 1
         
         IF (ASIZE_LOCAL .LT. ROW) THEN
            ASIZE_LOCAL = ROW
         END IF
         IF (ASIZE_LOCAL .LT. COL) THEN
            ASIZE_LOCAL = COL
         END IF
         
C         write(*,*) ROW, COL, RPART, IPART
         
         mumps_par%IRN_loc(COUNT) = ROW
         mumps_par%JCN_loc(COUNT) = COL
         mumps_par%A_loc(COUNT) = cmplx(RPART, IPART)


C         write(*,*) LINE, ROW, COL, mumps_par%A_loc(COUNT)
C         IF (COUNT .EQ. 30) THEN
C            EXIT
C         ENDIF
         COUNT = COUNT + 1      
       END DO
      END IF
 20   CLOSE (1)
      
C     SIZE of MATRIX
      CALL MPI_ALLREDUCE(ASIZE_LOCAL, ASIZE, 1, MPI_LONG, MPI_MAX,
     &     MPI_COMM_WORLD, ierr)
      IF ( mumps_par%MYID .eq. 0 ) THEN            
          write(*,*) "Done reading matrix, SIZE = ", ASIZE 
          write(*,*) "total nnz:", NNZ_TOTAL
      ENDIF
      
      mumps_par%N = ASIZE

      IF ( mumps_par%MYID .eq. 0 ) THEN
         ALLOCATE( mumps_par%RHS ( mumps_par%N  ) )

         IF (LOAD_RHS) THEN         
          OPEN (1, file = "rhs_0.000000") 
          DO
            READ (1,"(A)", END=30) LINE
            D1 = 1
            D2 = 1
            DO 
               IF (LINE(D2:D2) .EQ. " ")  THEN
                  COL_T = LINE(D1:D2-1)
                  EXIT
               END IF
               D2 = D2 + 1
            END DO
            D1 = D2+1
            D2 = D2+1
            DO 
               IF (LINE(D2:D2) .EQ. " ")  THEN
                  REAL_T = LINE(D1:D2-1)
                  EXIT
               END IF
               D2 = D2 + 1
            END DO
            D2 = D2+1
            IMAG_T = LINE(D2:)
            
            READ(COL_T, *, iostat=STAT)  COL                        
            READ(REAL_T, *, iostat=STAT) RPART
            READ(IMAG_T, *, iostat=STAT) IPART
            COL = COL + 1
            mumps_par%RHS(COL) = cmplx(RPART, IPART)
C            IF ((RPART .NE. 0) .OR. (IPART .NE. 0)) THEN
C               write(*,*)  COL, mumps_par%RHS(COL)
C            ENDIF
          END DO
         END IF
 30      CLOSE (1)
      END IF
      
C  Call package for solution
      mumps_par%JOB = 6
      CALL ZMUMPS(mumps_par)
      IF (mumps_par%INFOG(1).LT.0) THEN
       WRITE(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",
     &            "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1), 
     &            "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2) 
       GOTO 500
      END IF
C  Solution has been assembled on the host
      IF ( mumps_par%MYID .eq. 0 ) THEN
c     WRITE( 6, * ) ' Solution is ',(mumps_par%RHS(I),I=1,mumps_par%N)
        WRITE(6, *) ' Solution (1) ',mumps_par%RHS(1)
        WRITE(6, *) ' Solution (2) ',mumps_par%RHS(2)
        WRITE(6, *) ' Solution (3) ',mumps_par%RHS(3)
        WRITE(6, *) ' Solution (-2) ',mumps_par%RHS(ASIZE-1)
        WRITE(6, *) ' Solution (-1) ',mumps_par%RHS(ASIZE)
      END IF
      
C  Deallocate user data
      IF ( mumps_par%MYID .eq. 0 )THEN
         DEALLOCATE( mumps_par%RHS)
      ENDIF

      DEALLOCATE( mumps_par%IRN_loc )      
      DEALLOCATE( mumps_par%JCN_loc )
      DEALLOCATE( mumps_par%A_loc   )


C  Destroy the instance (deallocate internal data structures)
      mumps_par%JOB = -2
      CALL ZMUMPS(mumps_par)
      IF (mumps_par%INFOG(1).LT.0) THEN
       WRITE(6,'(A,A,I6,A,I9)') " ERROR RETURN: ",
     &            "  mumps_par%INFOG(1)= ", mumps_par%INFOG(1), 
     &            "  mumps_par%INFOG(2)= ", mumps_par%INFOG(2) 
       GOTO 500
      END IF
 500  CALL MPI_FINALIZE(IERR)
      STOP
      END
