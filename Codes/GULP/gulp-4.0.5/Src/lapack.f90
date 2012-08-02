  SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
  use datatypes
  integer(i4)        INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
  integer(i4)        IPIV( * )
  real(dp)           A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DGETRF computes an LU factorization of a general M-by-N matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 3 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the M-by-N matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
!                has been completed, but the factor U is exactly
!                singular, and division by zero will occur if it is used
!                to solve a system of equations.
!
!  =====================================================================
!
!     .. Parameters ..
  real(dp)           ONE
  PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
  integer(i4)        I, IINFO, J, JB, NB
!     ..
!     .. External Subroutines ..
  EXTERNAL           DGEMM, DGETF2, DLASWP, DTRSM, XERBLA
!     ..
!     .. External Functions ..
  integer(i4)        ILAENV
  EXTERNAL           ILAENV
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
  INFO = 0
  IF( M.LT.0 ) THEN
     INFO = -1
  ELSE IF( N.LT.0 ) THEN
     INFO = -2
  ELSE IF( LDA.LT.MAX( 1_i4, M ) ) THEN
     INFO = -4
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DGETRF', -INFO )
     RETURN
  END IF
!
!     Quick return if possible
!
  IF( M.EQ.0 .OR. N.EQ.0 ) &
     RETURN
!
!     Determine the block size for this environment.
!
  NB = ILAENV( 1_i4, 'DGETRF', ' ', M, N, -1_i4, -1_i4 )
  IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
!
!        Use unblocked code.
!
     CALL DGETF2( M, N, A, LDA, IPIV, INFO )
  ELSE
!
!        Use blocked code.
!
     DO 20 J = 1, MIN( M, N ), NB
        JB = MIN( MIN( M, N )-J+1_i4, NB )
!
!           Factor diagonal and subdiagonal blocks and test for exact
!           singularity.
!
        CALL DGETF2( M-J+1_i4, JB, A( J, J ), LDA, IPIV( J ), IINFO )
!
!           Adjust INFO and the pivot indices.
!
        IF( INFO.EQ.0 .AND. IINFO.GT.0 ) &
           INFO = IINFO + J - 1
        DO 10 I = J, MIN( M, J+JB-1 )
           IPIV( I ) = J - 1 + IPIV( I )
10       CONTINUE
!
!           Apply interchanges to columns 1:J-1.
!
        CALL DLASWP( J-1_i4, A, LDA, J, J+JB-1_i4, IPIV, 1_i4 )
!
        IF( J+JB.LE.N ) THEN
!
!              Apply interchanges to columns J+JB:N.
!
           CALL DLASWP( N-J-JB+1_i4, A( 1_i4, J+JB ), LDA, J, J+JB-1_i4, &
                        IPIV, 1_i4 )
!
!              Compute block row of U.
!
           CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB, &
                       N-J-JB+1_i4, ONE, A( J, J ), LDA, A( J, J+JB ), &
                       LDA )
           IF( J+JB.LE.M ) THEN
!
!                 Update trailing submatrix.
!
              CALL DGEMM( 'No transpose', 'No transpose', M-J-JB+1_i4, &
                          N-J-JB+1_i4, JB, -ONE, A( J+JB, J ), LDA, &
                          A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ), &
                          LDA )
           END IF
        END IF
20    CONTINUE
  END IF
  RETURN
!
!     End of DGETRF
!
  END
  SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
  use datatypes
  integer(i4)        INFO, LDA, LWORK, N
!     ..
!     .. Array Arguments ..
  integer(i4)        IPIV( * )
  real(dp)           A( LDA, * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DGETRI computes the inverse of a matrix using the LU factorization
!  computed by DGETRF.
!
!  This method inverts U and then computes inv(A) by solving the system
!  inv(A)*L = inv(U) for inv(A).
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the factors L and U from the factorization
!          A = P*L*U as computed by DGETRF.
!          On exit, if INFO = 0, the inverse of the original matrix A.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  IPIV    (input) INTEGER array, dimension (N)
!          The pivot indices from DGETRF; for 1<=i<=N, row i of the
!          matrix was interchanged with row IPIV(i).
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
!          On exit, if INFO=0, then WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.  LWORK >= max(1,N).
!          For optimal performance LWORK >= N*NB, where NB is
!          the optimal blocksize returned by ILAENV.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
!                singular and its inverse could not be computed.
!
!  =====================================================================
!
!     .. Parameters ..
  real(dp)           ZERO, ONE
  PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
  LOGICAL            LQUERY
  integer(i4)        I, IWS, J, JB, JJ, JP, LDWORK, LWKOPT, NB, &
                     NBMIN, NN
!     ..
!     .. External Functions ..
  integer(i4)        ILAENV
  EXTERNAL           ILAENV
!     ..
!     .. External Subroutines ..
  EXTERNAL           DGEMM, DGEMV, DSWAP, DTRSM, DTRTRI, XERBLA
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
  INFO = 0
  NB = ILAENV( 1_i4, 'DGETRI', ' ', N, -1_i4, -1_i4, -1_i4 )
  LWKOPT = N*NB
  WORK( 1 ) = LWKOPT
  LQUERY = ( LWORK.EQ.-1 )
  IF( N.LT.0 ) THEN
     INFO = -1
  ELSE IF( LDA.LT.MAX( 1_i4, N ) ) THEN
     INFO = -3
  ELSE IF( LWORK.LT.MAX( 1_i4, N ) .AND. .NOT.LQUERY ) THEN
     INFO = -6
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DGETRI', -INFO )
     RETURN
  ELSE IF( LQUERY ) THEN
     RETURN
  END IF
!
!     Quick return if possible
!
  IF( N.EQ.0 ) &
     RETURN
!
!     Form inv(U).  If INFO > 0 from DTRTRI, then U is singular,
!     and the inverse is not computed.
!
  CALL DTRTRI( 'Upper', 'Non-unit', N, A, LDA, INFO )
  IF( INFO.GT.0 ) &
     RETURN
!
  NBMIN = 2
  LDWORK = N
  IF( NB.GT.1 .AND. NB.LT.N ) THEN
     IWS = MAX( LDWORK*NB, 1 )
     IF( LWORK.LT.IWS ) THEN
        NB = LWORK / LDWORK
        NBMIN = MAX( 2_i4, ILAENV( 2_i4, 'DGETRI', ' ', N, -1_i4, -1_i4, -1_i4 ) )
     END IF
  ELSE
     IWS = N
  END IF
!
!     Solve the equation inv(A)*L = inv(U) for inv(A).
!
  IF( NB.LT.NBMIN .OR. NB.GE.N ) THEN
!
!        Use unblocked code.
!
     DO 20 J = N, 1, -1
!
!           Copy current column of L to WORK and replace with zeros.
!
        DO 10 I = J + 1, N
           WORK( I ) = A( I, J )
           A( I, J ) = ZERO
10       CONTINUE
!
!           Compute current column of inv(A).
!
        IF( J.LT.N ) &
           CALL DGEMV( 'No transpose', N, N-J, -ONE, A( 1, J+1 ), &
                       LDA, WORK( J+1 ), 1_i4, ONE, A( 1, J ), 1_i4 )
20    CONTINUE
  ELSE
!
!        Use blocked code.
!
     NN = ( ( N-1 ) / NB )*NB + 1
     DO 50 J = NN, 1, -NB
        JB = MIN( NB, N-J+1 )
!
!           Copy current block column of L to WORK and replace with
!           zeros.
!
        DO 40 JJ = J, J + JB - 1
           DO 30 I = JJ + 1, N
              WORK( I+( JJ-J )*LDWORK ) = A( I, JJ )
              A( I, JJ ) = ZERO
30          CONTINUE
40       CONTINUE
!
!           Compute current block column of inv(A).
!
        IF( J+JB.LE.N ) &
           CALL DGEMM( 'No transpose', 'No transpose', N, JB, &
                       N-J-JB+1_i4, -ONE, A( 1, J+JB ), LDA, &
                       WORK( J+JB ), LDWORK, ONE, A( 1, J ), LDA )
        CALL DTRSM( 'Right', 'Lower', 'No transpose', 'Unit', N, JB, &
                    ONE, WORK( J ), LDWORK, A( 1, J ), LDA )
50    CONTINUE
  END IF
!
!     Apply column interchanges.
!
  DO 60 J = N - 1, 1, -1
     JP = IPIV( J )
     IF( JP.NE.J ) &
        CALL DSWAP( N, A( 1, J ), 1_i4, A( 1, JP ), 1_i4 )
60 CONTINUE
!
  WORK( 1 ) = IWS
  RETURN
!
!     End of DGETRI
!
  END
  SUBROUTINE DLASWP( N, A, LDA, K1, K2, IPIV, INCX )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
  use datatypes
  integer(i4)        INCX, K1, K2, LDA, N
!     ..
!     .. Array Arguments ..
  integer(i4)        IPIV( * )
  real(dp)           A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DLASWP performs a series of row interchanges on the matrix A.
!  One row interchange is initiated for each of rows K1 through K2 of A.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the matrix of column dimension N to which the row
!          interchanges will be applied.
!          On exit, the permuted matrix.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!
!  K1      (input) INTEGER
!          The first element of IPIV for which a row interchange will
!          be done.
!
!  K2      (input) INTEGER
!          The last element of IPIV for which a row interchange will
!          be done.
!
!  IPIV    (input) INTEGER array, dimension (M*abs(INCX))
!          The vector of pivot indices.  Only the elements in positions
!          K1 through K2 of IPIV are accessed.
!          IPIV(K) = L implies rows K and L are to be interchanged.
!
!  INCX    (input) INTEGER
!          The increment between successive values of IPIV.  If IPIV
!          is negative, the pivots are applied in reverse order.
!
!  Further Details
!  ===============
!
!  Modified by
!   R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA
!
! =====================================================================
!
!     .. Local Scalars ..
  integer(i4)        I, I1, I2L, INC, IP, IX, IX0, J, K, N32
  real(dp)           TEMP
!     ..
!     .. Executable Statements ..
!
!     Interchange row I with row IPIV(I) for each of rows K1 through K2.
!
  IF( INCX.GT.0 ) THEN
     IX0 = K1
     I1 = K1
     I2L = K2
     INC = 1
  ELSE IF( INCX.LT.0 ) THEN
     IX0 = 1 + ( 1-K2 )*INCX
     I1 = K2
     I2L = K1
     INC = -1
  ELSE
     RETURN
  END IF
!
  N32 = ( N / 32 )*32
  IF( N32.NE.0 ) THEN
     DO 30 J = 1, N32, 32
        IX = IX0
        DO 20 I = I1, I2L, INC
           IP = IPIV( IX )
           IF( IP.NE.I ) THEN
              DO 10 K = J, J + 31
                 TEMP = A( I, K )
                 A( I, K ) = A( IP, K )
                 A( IP, K ) = TEMP
10             CONTINUE
           END IF
           IX = IX + INCX
20       CONTINUE
30    CONTINUE
  END IF
  IF( N32.NE.N ) THEN
     N32 = N32 + 1
     IX = IX0
     DO 50 I = I1, I2L, INC
        IP = IPIV( IX )
        IF( IP.NE.I ) THEN
           DO 40 K = N32, N
              TEMP = A( I, K )
              A( I, K ) = A( IP, K )
              A( IP, K ) = TEMP
40          CONTINUE
        END IF
        IX = IX + INCX
50    CONTINUE
  END IF
!
  RETURN
!
!     End of DLASWP
!
  END
  FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
  use datatypes
  CHARACTER(len=*)   NAME, OPTS
  integer(i4)        ISPEC, N1, N2, N3, N4
  integer(i4)        ILAENV
!     ..
!
!  Purpose
!  =======
!
!  ILAENV is called from the LAPACK routines to choose problem-dependent
!  parameters for the local environment.  See ISPEC for a description of
!  the parameters.
!
!  This version provides a set of parameters which should give good,
!  but not optimal, performance on many of the currently available
!  computers.  Users are encouraged to modify this subroutine to set
!  the tuning parameters for their particular machine using the option
!  and problem size information in the arguments.
!
!  This routine will not function correctly if it is converted to all
!  lower case.  Converting it to all upper case is allowed.
!
!  Arguments
!  =========
!
!  ISPEC   (input) INTEGER
!          Specifies the parameter to be returned as the value of
!          ILAENV.
!          = 1: the optimal blocksize; if this value is 1, an unblocked
!               algorithm will give the best performance.
!          = 2: the minimum block size for which the block routine
!               should be used; if the usable block size is less than
!               this value, an unblocked routine should be used.
!          = 3: the crossover point (in a block routine, for N less
!               than this value, an unblocked routine should be used)
!          = 4: the number of shifts, used in the nonsymmetric
!               eigenvalue routines
!          = 5: the minimum column dimension for blocking to be used;
!               rectangular blocks must have dimension at least k by m,
!               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
!          = 6: the crossover point for the SVD (when reducing an m by n
!               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
!               this value, a QR factorization is used first to reduce
!               the matrix to a triangular form.)
!          = 7: the number of processors
!          = 8: the crossover point for the multishift QR and QZ methods
!               for nonsymmetric eigenvalue problems.
!          = 9: maximum size of the subproblems at the bottom of the
!               computation tree in the divide-and-conquer algorithm
!               (used by xGELSD and xGESDD)
!          =10: ieee NaN arithmetic can be trusted not to trap
!          =11: infinity arithmetic can be trusted not to trap
!
!  NAME    (input) CHARACTER*(*)
!          The name of the calling subroutine, in either upper case or
!          lower case.
!
!  OPTS    (input) CHARACTER*(*)
!          The character options to the subroutine NAME, concatenated
!          into a single character string.  For example, UPLO = 'U',
!          TRANS = 'T', and DIAG = 'N' for a triangular routine would
!          be specified as OPTS = 'UTN'.
!
!  N1      (input) INTEGER
!  N2      (input) INTEGER
!  N3      (input) INTEGER
!  N4      (input) INTEGER
!          Problem dimensions for the subroutine NAME; these may not all
!          be required.
!
! (ILAENV) (output) INTEGER
!          >= 0: the value of the parameter specified by ISPEC
!          < 0:  if ILAENV = -k, the k-th argument had an illegal value.
!
!  Further Details
!  ===============
!
!  The following conventions have been used when calling ILAENV from the
!  LAPACK routines:
!  1)  OPTS is a concatenation of all of the character options to
!      subroutine NAME, in the same order that they appear in the
!      argument list for NAME, even if they are not used in determining
!      the value of the parameter specified by ISPEC.
!  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
!      that they appear in the argument list for NAME.  N1 is used
!      first, N2 second, and so on, and unused problem dimensions are
!      passed a value of -1.
!  3)  The parameter value returned by ILAENV is checked for validity in
!      the calling subroutine.  For example, ILAENV is used to retrieve
!      the optimal blocksize for STRTRI as follows:
!
!      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
!      IF( NB.LE.1 ) NB = MAX( 1, N )
!
!  =====================================================================
!
!     .. Local Scalars ..
  LOGICAL            CNAME, SNAME
  CHARACTER(len=1)   C1
  CHARACTER(len=2)   C2, C4
  CHARACTER(len=3)   C3
  CHARACTER(len=6)   SUBNAM
  integer(i4)        I, IC, IZ, NB, NBMIN, NX
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
!     ..
!     .. External Functions ..
  integer(i4)        IEEECK
  EXTERNAL           IEEECK
!     ..
!     .. Executable Statements ..
!
  GO TO ( 100, 100, 100, 400, 500, 600, 700, 800, 900, 1000, &
          1100 ) ISPEC
!
!     Invalid value for ISPEC
!
  ILAENV = -1
  RETURN
!
100 CONTINUE
!
!     Convert NAME to upper case if the first character is lower case.
!
  ILAENV = 1
  SUBNAM = NAME
  IC = ICHAR( SUBNAM( 1:1 ) )
  IZ = ICHAR( 'Z' )
  IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
!
!        ASCII character set
!
     IF( IC.GE.97 .AND. IC.LE.122 ) THEN
        SUBNAM( 1:1 ) = CHAR( IC-32 )
        DO 10 I = 2, 6
           IC = ICHAR( SUBNAM( I:I ) )
           IF( IC.GE.97 .AND. IC.LE.122 ) &
              SUBNAM( I:I ) = CHAR( IC-32 )
10       CONTINUE
     END IF
!
  ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
!
!        EBCDIC character set
!
     IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR. &
         ( IC.GE.145 .AND. IC.LE.153 ) .OR. &
         ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
        SUBNAM( 1:1 ) = CHAR( IC+64 )
        DO 20 I = 2, 6
           IC = ICHAR( SUBNAM( I:I ) )
           IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR. &
               ( IC.GE.145 .AND. IC.LE.153 ) .OR. &
               ( IC.GE.162 .AND. IC.LE.169 ) ) &
              SUBNAM( I:I ) = CHAR( IC+64 )
20       CONTINUE
     END IF
!
  ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
!
!        Prime machines:  ASCII+128
!
     IF( IC.GE.225 .AND. IC.LE.250 ) THEN
        SUBNAM( 1:1 ) = CHAR( IC-32 )
        DO 30 I = 2, 6
           IC = ICHAR( SUBNAM( I:I ) )
           IF( IC.GE.225 .AND. IC.LE.250 ) &
              SUBNAM( I:I ) = CHAR( IC-32 )
30       CONTINUE
     END IF
  END IF
!
  C1 = SUBNAM( 1:1 )
  SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
  CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
  IF( .NOT.( CNAME .OR. SNAME ) ) &
     RETURN
  C2 = SUBNAM( 2:3 )
  C3 = SUBNAM( 4:6 )
  C4 = C3( 2:3 )
!
  GO TO ( 110, 200, 300 ) ISPEC
!
110 CONTINUE
!
!     ISPEC = 1:  block size
!
!     In these examples, separate code is provided for setting NB for
!     real and complex.  We assume that NB will take the same value in
!     single or double precision.
!
  NB = 1
!
  IF( C2.EQ.'GE' ) THEN
     IF( C3.EQ.'TRF' ) THEN
        IF( SNAME ) THEN
           NB = 64
        ELSE
           NB = 64
        END IF
     ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. &
              C3.EQ.'QLF' ) THEN
        IF( SNAME ) THEN
           NB = 32
        ELSE
           NB = 32
        END IF
     ELSE IF( C3.EQ.'HRD' ) THEN
        IF( SNAME ) THEN
           NB = 32
        ELSE
           NB = 32
        END IF
     ELSE IF( C3.EQ.'BRD' ) THEN
        IF( SNAME ) THEN
           NB = 32
        ELSE
           NB = 32
        END IF
     ELSE IF( C3.EQ.'TRI' ) THEN
        IF( SNAME ) THEN
           NB = 64
        ELSE
           NB = 64
        END IF
     END IF
  ELSE IF( C2.EQ.'PO' ) THEN
     IF( C3.EQ.'TRF' ) THEN
        IF( SNAME ) THEN
           NB = 64
        ELSE
           NB = 64
        END IF
     END IF
  ELSE IF( C2.EQ.'SY' ) THEN
     IF( C3.EQ.'TRF' ) THEN
        IF( SNAME ) THEN
           NB = 64
        ELSE
           NB = 64
        END IF
     ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
        NB = 32
     ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
        NB = 64
     END IF
  ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
     IF( C3.EQ.'TRF' ) THEN
        NB = 64
     ELSE IF( C3.EQ.'TRD' ) THEN
        NB = 32
     ELSE IF( C3.EQ.'GST' ) THEN
        NB = 64
     END IF
  ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
     IF( C3( 1:1 ).EQ.'G' ) THEN
        IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
            C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
            C4.EQ.'BR' ) THEN
           NB = 32
        END IF
     ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
        IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
            C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
            C4.EQ.'BR' ) THEN
           NB = 32
        END IF
     END IF
  ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
     IF( C3( 1:1 ).EQ.'G' ) THEN
        IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
            C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
            C4.EQ.'BR' ) THEN
           NB = 32
        END IF
     ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
        IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
            C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
            C4.EQ.'BR' ) THEN
           NB = 32
        END IF
     END IF
  ELSE IF( C2.EQ.'GB' ) THEN
     IF( C3.EQ.'TRF' ) THEN
        IF( SNAME ) THEN
           IF( N4.LE.64 ) THEN
              NB = 1
           ELSE
              NB = 32
           END IF
        ELSE
           IF( N4.LE.64 ) THEN
              NB = 1
           ELSE
              NB = 32
           END IF
        END IF
     END IF
  ELSE IF( C2.EQ.'PB' ) THEN
     IF( C3.EQ.'TRF' ) THEN
        IF( SNAME ) THEN
           IF( N2.LE.64 ) THEN
              NB = 1
           ELSE
              NB = 32
           END IF
        ELSE
           IF( N2.LE.64 ) THEN
              NB = 1
           ELSE
              NB = 32
           END IF
        END IF
     END IF
  ELSE IF( C2.EQ.'TR' ) THEN
     IF( C3.EQ.'TRI' ) THEN
        IF( SNAME ) THEN
           NB = 64
        ELSE
           NB = 64
        END IF
     END IF
  ELSE IF( C2.EQ.'LA' ) THEN
     IF( C3.EQ.'UUM' ) THEN
        IF( SNAME ) THEN
           NB = 64
        ELSE
           NB = 64
        END IF
     END IF
  ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
     IF( C3.EQ.'EBZ' ) THEN
        NB = 1
     END IF
  END IF
  ILAENV = NB
  RETURN
!
200 CONTINUE
!
!     ISPEC = 2:  minimum block size
!
  NBMIN = 2
  IF( C2.EQ.'GE' ) THEN
     IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. &
         C3.EQ.'QLF' ) THEN
        IF( SNAME ) THEN
           NBMIN = 2
        ELSE
           NBMIN = 2
        END IF
     ELSE IF( C3.EQ.'HRD' ) THEN
        IF( SNAME ) THEN
           NBMIN = 2
        ELSE
           NBMIN = 2
        END IF
     ELSE IF( C3.EQ.'BRD' ) THEN
        IF( SNAME ) THEN
           NBMIN = 2
        ELSE
           NBMIN = 2
        END IF
     ELSE IF( C3.EQ.'TRI' ) THEN
        IF( SNAME ) THEN
           NBMIN = 2
        ELSE
           NBMIN = 2
        END IF
     END IF
  ELSE IF( C2.EQ.'SY' ) THEN
     IF( C3.EQ.'TRF' ) THEN
        IF( SNAME ) THEN
           NBMIN = 8
        ELSE
           NBMIN = 8
        END IF
     ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
        NBMIN = 2
     END IF
  ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
     IF( C3.EQ.'TRD' ) THEN
        NBMIN = 2
     END IF
  ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
     IF( C3( 1:1 ).EQ.'G' ) THEN
        IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
            C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
            C4.EQ.'BR' ) THEN
           NBMIN = 2
        END IF
     ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
        IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
            C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
            C4.EQ.'BR' ) THEN
           NBMIN = 2
        END IF
     END IF
  ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
     IF( C3( 1:1 ).EQ.'G' ) THEN
        IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
            C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
            C4.EQ.'BR' ) THEN
           NBMIN = 2
        END IF
     ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
        IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
            C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
            C4.EQ.'BR' ) THEN
           NBMIN = 2
        END IF
     END IF
  END IF
  ILAENV = NBMIN
  RETURN
!
300 CONTINUE
!
!     ISPEC = 3:  crossover point
!
  NX = 0
  IF( C2.EQ.'GE' ) THEN
     IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. &
         C3.EQ.'QLF' ) THEN
        IF( SNAME ) THEN
           NX = 128
        ELSE
           NX = 128
        END IF
     ELSE IF( C3.EQ.'HRD' ) THEN
        IF( SNAME ) THEN
           NX = 128
        ELSE
           NX = 128
        END IF
     ELSE IF( C3.EQ.'BRD' ) THEN
        IF( SNAME ) THEN
           NX = 128
        ELSE
           NX = 128
        END IF
     END IF
  ELSE IF( C2.EQ.'SY' ) THEN
     IF( SNAME .AND. C3.EQ.'TRD' ) THEN
        NX = 32
     END IF
  ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
     IF( C3.EQ.'TRD' ) THEN
        NX = 32
     END IF
  ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
     IF( C3( 1:1 ).EQ.'G' ) THEN
        IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
            C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
            C4.EQ.'BR' ) THEN
           NX = 128
        END IF
     END IF
  ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
     IF( C3( 1:1 ).EQ.'G' ) THEN
        IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
            C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
            C4.EQ.'BR' ) THEN
           NX = 128
        END IF
     END IF
  END IF
  ILAENV = NX
  RETURN
!
400 CONTINUE
!
!     ISPEC = 4:  number of shifts (used by xHSEQR)
!
  ILAENV = 6
  RETURN
!
500 CONTINUE
!
!     ISPEC = 5:  minimum column dimension (not used)
!
  ILAENV = 2
  RETURN
!
600 CONTINUE 
!
!     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
!
  ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
  RETURN
!
700 CONTINUE
!
!     ISPEC = 7:  number of processors (not used)
!
  ILAENV = 1
  RETURN
!
800 CONTINUE
!
!     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
!
  ILAENV = 50
  RETURN
!
900 CONTINUE
!
!     ISPEC = 9:  maximum size of the subproblems at the bottom of the
!                 computation tree in the divide-and-conquer algorithm
!                 (used by xGELSD and xGESDD)
!
  ILAENV = 25
  RETURN
!
1000 CONTINUE
!
!     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap
!
!     ILAENV = 0
  ILAENV = 1
  IF( ILAENV.EQ.1 ) THEN
     ILAENV = IEEECK( 0_i4, 0.0_dp, 1.0_dp ) 
  END IF
  RETURN
!
1100 CONTINUE
!
!     ISPEC = 11: infinity arithmetic can be trusted not to trap
!
!     ILAENV = 0
  ILAENV = 1
  IF( ILAENV.EQ.1 ) THEN
     ILAENV = IEEECK( 1_i4, 0.0_dp, 1.0_dp ) 
  END IF
  RETURN
!
!     End of ILAENV
!
  END
  SUBROUTINE DGETF2( M, N, A, LDA, IPIV, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1992
!
!     .. Scalar Arguments ..
  use datatypes
  integer(i4)        INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
  integer(i4)        IPIV( * )
  real(dp)           A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DGETF2 computes an LU factorization of a general m-by-n matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 2 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the m by n matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.
!
!  =====================================================================
!
!     .. Parameters ..
  real(dp)           ONE, ZERO
  PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
  integer(i4)        J, JP
!     ..
!     .. External Functions ..
  integer(i4)        IDAMAX
  EXTERNAL           IDAMAX
!     ..
!     .. External Subroutines ..
  EXTERNAL           DGER, DSCAL, DSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
  INFO = 0
  IF( M.LT.0 ) THEN
     INFO = -1
  ELSE IF( N.LT.0 ) THEN
     INFO = -2
  ELSE IF( LDA.LT.MAX( 1_i4, M ) ) THEN
     INFO = -4
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DGETF2', -INFO )
     RETURN
  END IF
!
!     Quick return if possible
!
  IF( M.EQ.0 .OR. N.EQ.0 ) &
     RETURN
!
  DO 10 J = 1, MIN( M, N )
!
!        Find pivot and test for singularity.
!
     JP = J - 1 + IDAMAX( M-J+1_i4, A( J, J ), 1_i4 )
     IPIV( J ) = JP
     IF( A( JP, J ).NE.ZERO ) THEN
!
!           Apply the interchange to columns 1:N.
!
        IF( JP.NE.J ) &
           CALL DSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
!
!           Compute elements J+1:M of J-th column.
!
        IF( J.LT.M ) &
           CALL DSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1_i4 )
!
     ELSE IF( INFO.EQ.0 ) THEN
!
        INFO = J
     END IF
!
     IF( J.LT.MIN( M, N ) ) THEN
!
!           Update trailing submatrix.
!
        CALL DGER( M-J, N-J, -ONE, A( J+1, J ), 1_i4, A( J, J+1 ), LDA, &
                   A( J+1, J+1 ), LDA )
     END IF
10 CONTINUE
  RETURN
!
!     End of DGETF2
!
  END
  SUBROUTINE DTRTRI( UPLO, DIAG, N, A, LDA, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
  use datatypes
  CHARACTER          DIAG, UPLO
  integer(i4)        INFO, LDA, N
!     ..
!     .. Array Arguments ..
  real(dp)           A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DTRTRI computes the inverse of a real upper or lower triangular
!  matrix A.
!
!  This is the Level 3 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  A is upper triangular;
!          = 'L':  A is lower triangular.
!
!  DIAG    (input) CHARACTER*1
!          = 'N':  A is non-unit triangular;
!          = 'U':  A is unit triangular.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the triangular matrix A.  If UPLO = 'U', the
!          leading N-by-N upper triangular part of the array A contains
!          the upper triangular matrix, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading N-by-N lower triangular part of the array A contains
!          the lower triangular matrix, and the strictly upper
!          triangular part of A is not referenced.  If DIAG = 'U', the
!          diagonal elements of A are also not referenced and are
!          assumed to be 1.
!          On exit, the (triangular) inverse of the original matrix, in
!          the same storage format.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = i, A(i,i) is exactly zero.  The triangular
!               matrix is singular and its inverse can not be computed.
!
!  =====================================================================
!
!     .. Parameters ..
  real(dp)           ONE, ZERO
  PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
  LOGICAL            NOUNIT, UPPER
  integer(i4)        J, JB, NB, NN
!     ..
!     .. External Functions ..
  LOGICAL            LSAME
  integer(i4)        ILAENV
  EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
  EXTERNAL           DTRMM, DTRSM, DTRTI2, XERBLA
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
  INFO = 0
  UPPER = LSAME( UPLO, 'U' )
  NOUNIT = LSAME( DIAG, 'N' )
  IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
     INFO = -1
  ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
     INFO = -2
  ELSE IF( N.LT.0 ) THEN
     INFO = -3
  ELSE IF( LDA.LT.MAX( 1_i4, N ) ) THEN
     INFO = -5
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DTRTRI', -INFO )
     RETURN
  END IF
!
!     Quick return if possible
!
  IF( N.EQ.0 ) &
     RETURN
!
!     Check for singularity if non-unit.
!
  IF( NOUNIT ) THEN
     DO 10 INFO = 1, N
        IF( A( INFO, INFO ).EQ.ZERO ) &
           RETURN
10    CONTINUE
     INFO = 0
  END IF
!
!     Determine the block size for this environment.
!
  NB = ILAENV( 1_i4, 'DTRTRI', UPLO // DIAG, N, -1_i4, -1_i4, -1_i4 )
  IF( NB.LE.1 .OR. NB.GE.N ) THEN
!
!        Use unblocked code
!
     CALL DTRTI2( UPLO, DIAG, N, A, LDA, INFO )
  ELSE
!
!        Use blocked code
!
     IF( UPPER ) THEN
!
!           Compute inverse of upper triangular matrix
!
        DO 20 J = 1, N, NB
           JB = MIN( NB, N-J+1 )
!
!              Compute rows 1:j-1 of current block column
!
           CALL DTRMM( 'Left', 'Upper', 'No transpose', DIAG, J-1_i4, &
                       JB, ONE, A, LDA, A( 1, J ), LDA )
           CALL DTRSM( 'Right', 'Upper', 'No transpose', DIAG, J-1_i4, &
                       JB, -ONE, A( J, J ), LDA, A( 1, J ), LDA )
!
!              Compute inverse of current diagonal block
!
           CALL DTRTI2( 'Upper', DIAG, JB, A( J, J ), LDA, INFO )
20       CONTINUE
     ELSE
!
!           Compute inverse of lower triangular matrix
!
        NN = ( ( N-1 ) / NB )*NB + 1
        DO 30 J = NN, 1, -NB
           JB = MIN( NB, N-J+1 )
           IF( J+JB.LE.N ) THEN
!
!                 Compute rows j+jb:n of current block column
!
              CALL DTRMM( 'Left', 'Lower', 'No transpose', DIAG, &
                          N-J-JB+1_i4, JB, ONE, A( J+JB, J+JB ), LDA, &
                          A( J+JB, J ), LDA )
              CALL DTRSM( 'Right', 'Lower', 'No transpose', DIAG, &
                          N-J-JB+1_i4, JB, -ONE, A( J, J ), LDA, &
                          A( J+JB, J ), LDA )
           END IF
!
!              Compute inverse of current diagonal block
!
           CALL DTRTI2( 'Lower', DIAG, JB, A( J, J ), LDA, INFO )
30       CONTINUE
     END IF
  END IF
!
  RETURN
!
!     End of DTRTRI
!
  END
  SUBROUTINE DTRTI2( UPLO, DIAG, N, A, LDA, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
  use datatypes
  CHARACTER          DIAG, UPLO
  integer(i4)        INFO, LDA, N
!     ..
!     .. Array Arguments ..
  real(dp)           A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DTRTI2 computes the inverse of a real upper or lower triangular
!  matrix.
!
!  This is the Level 2 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the matrix A is upper or lower triangular.
!          = 'U':  Upper triangular
!          = 'L':  Lower triangular
!
!  DIAG    (input) CHARACTER*1
!          Specifies whether or not the matrix A is unit triangular.
!          = 'N':  Non-unit triangular
!          = 'U':  Unit triangular
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the triangular matrix A.  If UPLO = 'U', the
!          leading n by n upper triangular part of the array A contains
!          the upper triangular matrix, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading n by n lower triangular part of the array A contains
!          the lower triangular matrix, and the strictly upper
!          triangular part of A is not referenced.  If DIAG = 'U', the
!          diagonal elements of A are also not referenced and are
!          assumed to be 1.
!
!          On exit, the (triangular) inverse of the original matrix, in
!          the same storage format.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
  real(dp)           ONE
  PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
  LOGICAL            NOUNIT, UPPER
  integer(i4)        J
  real(dp)           AJJ
!     ..
!     .. External Functions ..
  LOGICAL            LSAME
  EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
  EXTERNAL           DSCAL, DTRMV, XERBLA
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
  INFO = 0
  UPPER = LSAME( UPLO, 'U' )
  NOUNIT = LSAME( DIAG, 'N' )
  IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
     INFO = -1
  ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
     INFO = -2
  ELSE IF( N.LT.0 ) THEN
     INFO = -3
  ELSE IF( LDA.LT.MAX( 1_i4, N ) ) THEN
     INFO = -5
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DTRTI2', -INFO )
     RETURN
  END IF
!
  IF( UPPER ) THEN
!
!        Compute inverse of upper triangular matrix.
!
     DO 10 J = 1, N
        IF( NOUNIT ) THEN
           A( J, J ) = ONE / A( J, J )
           AJJ = -A( J, J )
        ELSE
           AJJ = -ONE
        END IF
!
!           Compute elements 1:j-1 of j-th column.
!
        CALL DTRMV( 'Upper', 'No transpose', DIAG, J-1_i4, A, LDA, &
                    A( 1, J ), 1_i4 )
        CALL DSCAL( J-1_i4, AJJ, A( 1, J ), 1_i4 )
10    CONTINUE
  ELSE
!
!        Compute inverse of lower triangular matrix.
!
     DO 20 J = N, 1, -1
        IF( NOUNIT ) THEN
           A( J, J ) = ONE / A( J, J )
           AJJ = -A( J, J )
        ELSE
           AJJ = -ONE
        END IF
        IF( J.LT.N ) THEN
!
!              Compute elements j+1:n of j-th column.
!
           CALL DTRMV( 'Lower', 'No transpose', DIAG, N-J, &
                       A( J+1, J+1 ), LDA, A( J+1, J ), 1_i4 )
           CALL DSCAL( N-J, AJJ, A( J+1, J ), 1_i4 )
        END IF
20    CONTINUE
  END IF
!
  RETURN
!
!     End of DTRTI2
!
  END
  FUNCTION IEEECK( ISPEC, ZERO, ONE )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1998
!
!     .. Scalar Arguments ..
  use datatypes
  integer(i4)        ISPEC
  integer(i4)        IEEECK
  real(dp)           ONE, ZERO
!     ..
!
!  Purpose
!  =======
!
!  IEEECK is called from the ILAENV to verify that Infinity and
!  possibly NaN arithmetic is safe (i.e. will not trap).
!
!  Arguments
!  =========
!
!  ISPEC   (input) INTEGER
!          Specifies whether to test just for inifinity arithmetic
!          or whether to test for infinity and NaN arithmetic.
!          = 0: Verify infinity arithmetic only.
!          = 1: Verify infinity and NaN arithmetic.
!
!  ZERO    (input) REAL
!          Must contain the value 0.0
!          This is passed to prevent the compiler from optimizing
!          away this code.
!
!  ONE     (input) REAL
!          Must contain the value 1.0
!          This is passed to prevent the compiler from optimizing
!          away this code.
!
!  RETURN VALUE:  INTEGER
!          = 0:  Arithmetic failed to produce the correct answers
!          = 1:  Arithmetic produced the correct answers
!
!     .. Local Scalars ..
  real(dp)           NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF, &
                     NEGZRO, NEWZRO, POSINF
!     ..
!     .. Executable Statements ..
  IEEECK = 1
!
  POSINF = ONE / ZERO
  IF( POSINF.LE.ONE ) THEN
     IEEECK = 0
     RETURN
  END IF
!
  NEGINF = -ONE / ZERO
  IF( NEGINF.GE.ZERO ) THEN
     IEEECK = 0
     RETURN
  END IF
!
  NEGZRO = ONE / ( NEGINF+ONE )
  IF( NEGZRO.NE.ZERO ) THEN
     IEEECK = 0
     RETURN
  END IF
!
  NEGINF = ONE / NEGZRO
  IF( NEGINF.GE.ZERO ) THEN
     IEEECK = 0
     RETURN
  END IF
!
  NEWZRO = NEGZRO + ZERO
  IF( NEWZRO.NE.ZERO ) THEN
     IEEECK = 0
     RETURN
  END IF
!
  POSINF = ONE / NEWZRO
  IF( POSINF.LE.ONE ) THEN
     IEEECK = 0
     RETURN
  END IF
!
  NEGINF = NEGINF*POSINF
  IF( NEGINF.GE.ZERO ) THEN
     IEEECK = 0
     RETURN
  END IF
!
  POSINF = POSINF*POSINF
  IF( POSINF.LE.ONE ) THEN
     IEEECK = 0
     RETURN
  END IF
!
!
!
!
!     Return if we were only asked to check infinity arithmetic
!
  IF( ISPEC.EQ.0 ) &
     RETURN
!
  NAN1 = POSINF + NEGINF
!
  NAN2 = POSINF / NEGINF
!
  NAN3 = POSINF / POSINF
!
  NAN4 = POSINF*ZERO
!
  NAN5 = NEGINF*NEGZRO
!
  NAN6 = NAN5*0.0
!
  IF( NAN1.EQ.NAN1 ) THEN
     IEEECK = 0
     RETURN
  END IF
!
  IF( NAN2.EQ.NAN2 ) THEN
     IEEECK = 0
     RETURN
  END IF
!
  IF( NAN3.EQ.NAN3 ) THEN
     IEEECK = 0
     RETURN
  END IF
!
  IF( NAN4.EQ.NAN4 ) THEN
     IEEECK = 0
     RETURN
  END IF
!
  IF( NAN5.EQ.NAN5 ) THEN
     IEEECK = 0
     RETURN
  END IF
!
  IF( NAN6.EQ.NAN6 ) THEN
     IEEECK = 0
     RETURN
  END IF
!
  RETURN
  END
!
  SUBROUTINE DSPTRF( UPLO, N, AP, IPIV, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
  use datatypes
  CHARACTER          UPLO
  integer(i4)        INFO, N
!     ..
!     .. Array Arguments ..
  integer(i4)        IPIV( * )
  real(dp)           AP( * )
!     ..
!
!  Purpose
!  =======
!
!  DSPTRF computes the factorization of a real symmetric matrix A stored
!  in packed format using the Bunch-Kaufman diagonal pivoting method:
!
!     A = U*D*U**T  or  A = L*D*L**T
!
!  where U (or L) is a product of permutation and unit upper (lower)
!  triangular matrices, and D is symmetric and block diagonal with
!  1-by-1 and 2-by-2 diagonal blocks.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  AP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2)
!          On entry, the upper or lower triangle of the symmetric matrix
!          A, packed columnwise in a linear array.  The j-th column of A
!          is stored in the array AP as follows:
!          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
!
!          On exit, the block diagonal matrix D and the multipliers used
!          to obtain the factor U or L, stored as a packed triangular
!
!          matrix overwriting A (see below for further details).
!
!  IPIV    (output) INTEGER array, dimension (N)
!          Details of the interchanges and the block structure of D.
!          If IPIV(k) > 0, then rows and columns k and IPIV(k) were
!          interchanged and D(k,k) is a 1-by-1 diagonal block.
!          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
!          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
!          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
!          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
!          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = i, D(i,i) is exactly zero.  The factorization
!               has been completed, but the block diagonal matrix D is
!               exactly singular, and division by zero will occur if it
!               is used to solve a system of equations.
!
!  Further Details
!  ===============
!
!  5-96 - Based on modifications by J. Lewis, Boeing Computer Services
!         Company
!
!  If UPLO = 'U', then A = U*D*U', where
!     U = P(n)*U(n)* ... *P(k)U(k)* ...,
!  i.e., U is a product of terms P(k)*U(k), where k decreases from n to
!  1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
!  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
!  defined by IPIV(k), and U(k) is a unit upper triangular matrix, such
!  that if the diagonal block D(k) is of order s (s = 1 or 2), then
!
!             (   I    v    0   )   k-s
!     U(k) =  (   0    I    0   )   s
!             (   0    0    I   )   n-k
!                k-s   s   n-k
!
!  If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k).
!  If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k),
!  and A(k,k), and v overwrites A(1:k-2,k-1:k).
!
!  If UPLO = 'L', then A = L*D*L', where
!     L = P(1)*L(1)* ... *P(k)*L(k)* ...,
!  i.e., L is a product of terms P(k)*L(k), where k increases from 1 to
!  n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
!  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
!  defined by IPIV(k), and L(k) is a unit lower triangular matrix, such
!  that if the diagonal block D(k) is of order s (s = 1 or 2), then
!
!             (   I    0     0   )  k-1
!     L(k) =  (   0    I     0   )  s
!             (   0    v     I   )  n-k-s+1
!                k-1   s  n-k-s+1
!
!  If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k).
!  If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k),
!  and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1).
!
!  =====================================================================
!
!     .. Parameters ..
  real(dp)           ZERO, ONE
  PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
  real(dp)           EIGHT, SEVTEN
  PARAMETER          ( EIGHT = 8.0D+0, SEVTEN = 17.0D+0 )
!     ..
!     .. Local Scalars ..
  LOGICAL            UPPER
  integer(i4)        I, IMAX, J, JMAX, K, KC, KK, KNC, KP, KPC, &
                     KSTEP, KX, NPP
  real(dp)           ABSAKK, ALPHA, COLMAX, D11, D12, D21, D22, R1, &
                     ROWMAX, T, WK, WKM1, WKP1
!     ..
!     .. External Functions ..
  LOGICAL            LSAME
  integer(i4)        IDAMAX
  EXTERNAL           LSAME, IDAMAX
!     ..
!     .. External Subroutines ..
  EXTERNAL           DSCAL, DSPR, DSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          ABS, MAX, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
  INFO = 0
  UPPER = LSAME( UPLO, 'U' )
  IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
     INFO = -1
  ELSE IF( N.LT.0 ) THEN
     INFO = -2
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DSPTRF', -INFO )
     RETURN
  END IF
!
!     Initialize ALPHA for use in choosing pivot block size.
!
  ALPHA = ( ONE+SQRT( SEVTEN ) ) / EIGHT
!
  IF( UPPER ) THEN
!
!        Factorize A as U*D*U' using the upper triangle of A
!
!        K is the main loop index, decreasing from N to 1 in steps of
!        1 or 2
!
     K = N
     KC = ( N-1 )*N / 2 + 1
10    CONTINUE
     KNC = KC
!
!        If K < 1, exit from loop
!
     IF( K.LT.1 ) &
        GO TO 110
     KSTEP = 1
!
!        Determine rows and columns to be interchanged and whether
!        a 1-by-1 or 2-by-2 pivot block will be used
!
     ABSAKK = ABS( AP( KC+K-1 ) )
!
!        IMAX is the row-index of the largest off-diagonal element in
!        column K, and COLMAX is its absolute value
!
     IF( K.GT.1 ) THEN
        IMAX = IDAMAX( K-1_i4, AP( KC ), 1_i4 )
        COLMAX = ABS( AP( KC+IMAX-1 ) )
     ELSE
        COLMAX = ZERO
     END IF
!
     IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN
!
!           Column K is zero: set INFO and continue
!
        IF( INFO.EQ.0 ) &
           INFO = K
        KP = K
     ELSE
        IF( ABSAKK.GE.ALPHA*COLMAX ) THEN
!
!              no interchange, use 1-by-1 pivot block
!
           KP = K
        ELSE
!
!              JMAX is the column-index of the largest off-diagonal
!              element in row IMAX, and ROWMAX is its absolute value
!
           ROWMAX = ZERO
           JMAX = IMAX
           KX = IMAX*( IMAX+1 ) / 2 + IMAX
           DO 20 J = IMAX + 1, K
              IF( ABS( AP( KX ) ).GT.ROWMAX ) THEN
                 ROWMAX = ABS( AP( KX ) )
                 JMAX = J
              END IF
              KX = KX + J
20          CONTINUE
           KPC = ( IMAX-1 )*IMAX / 2 + 1
           IF( IMAX.GT.1 ) THEN
              JMAX = IDAMAX( IMAX-1_i4, AP( KPC ), 1_i4 )
              ROWMAX = MAX( ROWMAX, ABS( AP( KPC+JMAX-1 ) ) )
           END IF
!
           IF( ABSAKK.GE.ALPHA*COLMAX*( COLMAX / ROWMAX ) ) THEN
!
!                 no interchange, use 1-by-1 pivot block
!
              KP = K
           ELSE IF( ABS( AP( KPC+IMAX-1 ) ).GE.ALPHA*ROWMAX ) THEN
!
!                 interchange rows and columns K and IMAX, use 1-by-1
!                 pivot block
!
              KP = IMAX
           ELSE
!
!                 interchange rows and columns K-1 and IMAX, use 2-by-2
!                 pivot block
!
              KP = IMAX
              KSTEP = 2
           END IF
        END IF
!
        KK = K - KSTEP + 1
        IF( KSTEP.EQ.2 ) &
           KNC = KNC - K + 1
        IF( KP.NE.KK ) THEN
!
!              Interchange rows and columns KK and KP in the leading
!              submatrix A(1:k,1:k)
!
           CALL DSWAP( KP-1_i4, AP( KNC ), 1_i4, AP( KPC ), 1_i4 )
           KX = KPC + KP - 1
           DO 30 J = KP + 1, KK - 1
              KX = KX + J - 1
              T = AP( KNC+J-1 )
              AP( KNC+J-1 ) = AP( KX )
              AP( KX ) = T
30          CONTINUE
           T = AP( KNC+KK-1 )
           AP( KNC+KK-1 ) = AP( KPC+KP-1 )
           AP( KPC+KP-1 ) = T
           IF( KSTEP.EQ.2 ) THEN
              T = AP( KC+K-2 )
              AP( KC+K-2 ) = AP( KC+KP-1 )
              AP( KC+KP-1 ) = T
           END IF
        END IF
!
!           Update the leading submatrix
!
        IF( KSTEP.EQ.1 ) THEN
!
!              1-by-1 pivot block D(k): column k now holds
!
!              W(k) = U(k)*D(k)
!
!              where U(k) is the k-th column of U
!
!              Perform a rank-1 update of A(1:k-1,1:k-1) as
!
!              A := A - U(k)*D(k)*U(k)' = A - W(k)*1/D(k)*W(k)'
!
           R1 = ONE / AP( KC+K-1 )
           CALL DSPR( UPLO, K-1_i4, -R1, AP( KC ), 1_i4, AP )
!
!              Store U(k) in column k
!
           CALL DSCAL( K-1_i4, R1, AP( KC ), 1_i4 )
        ELSE
!
!              2-by-2 pivot block D(k): columns k and k-1 now hold
!
!              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)
!
!              where U(k) and U(k-1) are the k-th and (k-1)-th columns
!              of U
!
!              Perform a rank-2 update of A(1:k-2,1:k-2) as
!
!              A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )'
!                 = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )'
!
           IF( K.GT.2 ) THEN
!
              D12 = AP( K-1+( K-1 )*K / 2 )
              D22 = AP( K-1+( K-2 )*( K-1 ) / 2 ) / D12
              D11 = AP( K+( K-1 )*K / 2 ) / D12
              T = ONE / ( D11*D22-ONE )
              D12 = T / D12
!
              DO 50 J = K - 2, 1, -1
                 WKM1 = D12*( D11*AP( J+( K-2 )*( K-1 ) / 2 )- &
                        AP( J+( K-1 )*K / 2 ) )
                 WK = D12*( D22*AP( J+( K-1 )*K / 2 )- &
                      AP( J+( K-2 )*( K-1 ) / 2 ) )
                 DO 40 I = J, 1, -1
                    AP( I+( J-1 )*J / 2 ) = AP( I+( J-1 )*J / 2 ) - &
                       AP( I+( K-1 )*K / 2 )*WK - &
                       AP( I+( K-2 )*( K-1 ) / 2 )*WKM1
40                CONTINUE
                 AP( J+( K-1 )*K / 2 ) = WK
                 AP( J+( K-2 )*( K-1 ) / 2 ) = WKM1
50             CONTINUE
!
           END IF
!
        END IF
     END IF
!
!        Store details of the interchanges in IPIV
!
     IF( KSTEP.EQ.1 ) THEN
        IPIV( K ) = KP
     ELSE
        IPIV( K ) = -KP
        IPIV( K-1 ) = -KP
     END IF
!
!        Decrease K and return to the start of the main loop
!
     K = K - KSTEP
     KC = KNC - K
     GO TO 10
!
  ELSE
!
!        Factorize A as L*D*L' using the lower triangle of A
!
!        K is the main loop index, increasing from 1 to N in steps of
!        1 or 2
!
     K = 1
     KC = 1
     NPP = N*( N+1 ) / 2
60    CONTINUE
     KNC = KC
!
!        If K > N, exit from loop
!
     IF( K.GT.N ) &
        GO TO 110
     KSTEP = 1
!
!        Determine rows and columns to be interchanged and whether
!        a 1-by-1 or 2-by-2 pivot block will be used
!
     ABSAKK = ABS( AP( KC ) )
!
!        IMAX is the row-index of the largest off-diagonal element in
!        column K, and COLMAX is its absolute value
!
     IF( K.LT.N ) THEN
        IMAX = K + IDAMAX( N-K, AP( KC+1 ), 1_i4 )
        COLMAX = ABS( AP( KC+IMAX-K ) )
     ELSE
        COLMAX = ZERO
     END IF
!
     IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN
!
!           Column K is zero: set INFO and continue
!
        IF( INFO.EQ.0 ) &
           INFO = K
        KP = K
     ELSE
        IF( ABSAKK.GE.ALPHA*COLMAX ) THEN
!
!              no interchange, use 1-by-1 pivot block
!
           KP = K
        ELSE
!
!              JMAX is the column-index of the largest off-diagonal
!              element in row IMAX, and ROWMAX is its absolute value
!
           ROWMAX = ZERO
           KX = KC + IMAX - K
           DO 70 J = K, IMAX - 1
              IF( ABS( AP( KX ) ).GT.ROWMAX ) THEN
                 ROWMAX = ABS( AP( KX ) )
                 JMAX = J
              END IF
              KX = KX + N - J
70          CONTINUE
           KPC = NPP - ( N-IMAX+1 )*( N-IMAX+2 ) / 2 + 1
           IF( IMAX.LT.N ) THEN
              JMAX = IMAX + IDAMAX( N-IMAX, AP( KPC+1 ), 1_i4 )
              ROWMAX = MAX( ROWMAX, ABS( AP( KPC+JMAX-IMAX ) ) )
           END IF
!
           IF( ABSAKK.GE.ALPHA*COLMAX*( COLMAX / ROWMAX ) ) THEN
!
!                 no interchange, use 1-by-1 pivot block
!
              KP = K
           ELSE IF( ABS( AP( KPC ) ).GE.ALPHA*ROWMAX ) THEN
!
!                 interchange rows and columns K and IMAX, use 1-by-1
!                 pivot block
!
              KP = IMAX
           ELSE
!
!                 interchange rows and columns K+1 and IMAX, use 2-by-2
!                 pivot block
!
              KP = IMAX
              KSTEP = 2
           END IF
        END IF
!
        KK = K + KSTEP - 1
        IF( KSTEP.EQ.2 ) &
           KNC = KNC + N - K + 1
        IF( KP.NE.KK ) THEN
!
!              Interchange rows and columns KK and KP in the trailing
!              submatrix A(k:n,k:n)
!
           IF( KP.LT.N ) &
              CALL DSWAP( N-KP, AP( KNC+KP-KK+1 ), 1_i4, AP( KPC+1 ), &
                          1_i4 )
           KX = KNC + KP - KK
           DO 80 J = KK + 1, KP - 1
              KX = KX + N - J + 1
              T = AP( KNC+J-KK )
              AP( KNC+J-KK ) = AP( KX )
              AP( KX ) = T
80          CONTINUE
           T = AP( KNC )
           AP( KNC ) = AP( KPC )
           AP( KPC ) = T
           IF( KSTEP.EQ.2 ) THEN
              T = AP( KC+1 )
              AP( KC+1 ) = AP( KC+KP-K )
              AP( KC+KP-K ) = T
           END IF
        END IF
!
!           Update the trailing submatrix
!
        IF( KSTEP.EQ.1 ) THEN
!
!              1-by-1 pivot block D(k): column k now holds
!
!              W(k) = L(k)*D(k)
!
!              where L(k) is the k-th column of L
!
           IF( K.LT.N ) THEN
!
!                 Perform a rank-1 update of A(k+1:n,k+1:n) as
!
!                 A := A - L(k)*D(k)*L(k)' = A - W(k)*(1/D(k))*W(k)'
!
              R1 = ONE / AP( KC )
              CALL DSPR( UPLO, N-K, -R1, AP( KC+1 ), 1_i4, &
                         AP( KC+N-K+1 ) )
!
!                 Store L(k) in column K
!
              CALL DSCAL( N-K, R1, AP( KC+1 ), 1_i4 )
           END IF
        ELSE
!
!              2-by-2 pivot block D(k): columns K and K+1 now hold
!
!              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)
!
!              where L(k) and L(k+1) are the k-th and (k+1)-th columns
!              of L
!
           IF( K.LT.N-1 ) THEN
!
!                 Perform a rank-2 update of A(k+2:n,k+2:n) as
!
!                 A := A - ( L(k) L(k+1) )*D(k)*( L(k) L(k+1) )'
!                    = A - ( W(k) W(k+1) )*inv(D(k))*( W(k) W(k+1) )'
!
              D21 = AP( K+1+( K-1 )*( 2*N-K ) / 2 )
              D11 = AP( K+1+K*( 2*N-K-1 ) / 2 ) / D21
              D22 = AP( K+( K-1 )*( 2*N-K ) / 2 ) / D21
              T = ONE / ( D11*D22-ONE )
              D21 = T / D21
!
              DO 100 J = K + 2, N
                 WK = D21*( D11*AP( J+( K-1 )*( 2*N-K ) / 2 )- &
                      AP( J+K*( 2*N-K-1 ) / 2 ) )
                 WKP1 = D21*( D22*AP( J+K*( 2*N-K-1 ) / 2 )- &
                        AP( J+( K-1 )*( 2*N-K ) / 2 ) )
!
                 DO 90 I = J, N
                    AP( I+( J-1 )*( 2*N-J ) / 2 ) = AP( I+( J-1 )* &
                       ( 2*N-J ) / 2 ) - AP( I+( K-1 )*( 2*N-K ) / &
                       2 )*WK - AP( I+K*( 2*N-K-1 ) / 2 )*WKP1
90                CONTINUE
!
                 AP( J+( K-1 )*( 2*N-K ) / 2 ) = WK
                 AP( J+K*( 2*N-K-1 ) / 2 ) = WKP1
!
100             CONTINUE
           END IF
        END IF
     END IF
!
!        Store details of the interchanges in IPIV
!
     IF( KSTEP.EQ.1 ) THEN
        IPIV( K ) = KP
     ELSE
        IPIV( K ) = -KP
        IPIV( K+1 ) = -KP
     END IF
!
!        Increase K and return to the start of the main loop
!
     K = K + KSTEP
     KC = KNC + N - K + 2
     GO TO 60
!
  END IF
!
110 CONTINUE
  RETURN
!
!     End of DSPTRF
!
  END
!
  SUBROUTINE DSPTRI( UPLO, N, AP, IPIV, WORK, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
  use datatypes
  CHARACTER          UPLO
  integer(i4)        INFO, N
!     ..
!     .. Array Arguments ..
  integer(i4)        IPIV( * )
  real(dp)           AP( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DSPTRI computes the inverse of a real symmetric indefinite matrix
!  A in packed storage using the factorization A = U*D*U**T or
!  A = L*D*L**T computed by DSPTRF.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the details of the factorization are stored
!          as an upper or lower triangular matrix.
!          = 'U':  Upper triangular, form is A = U*D*U**T;
!          = 'L':  Lower triangular, form is A = L*D*L**T.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  AP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2)
!          On entry, the block diagonal matrix D and the multipliers
!          used to obtain the factor U or L as computed by DSPTRF,
!          stored as a packed triangular matrix.
!
!          On exit, if INFO = 0, the (symmetric) inverse of the original
!          matrix, stored as a packed triangular matrix. The j-th column
!          of inv(A) is stored in the array AP as follows:
!          if UPLO = 'U', AP(i + (j-1)*j/2) = inv(A)(i,j) for 1<=i<=j;
!          if UPLO = 'L',
!             AP(i + (j-1)*(2n-j)/2) = inv(A)(i,j) for j<=i<=n.
!
!  IPIV    (input) INTEGER array, dimension (N)
!          Details of the interchanges and the block structure of D
!          as determined by DSPTRF.
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = i, D(i,i) = 0; the matrix is singular and its
!               inverse could not be computed.
!
!  =====================================================================
!
!     .. Parameters ..
  real(dp)           ONE, ZERO
  PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
  LOGICAL            UPPER
  integer(i4)        J, K, KC, KCNEXT, KP, KPC, KSTEP, KX, NPP
  real(dp)           AK, AKKP1, AKP1, D, T, TEMP
!     ..
!     .. External Functions ..
  LOGICAL            LSAME
  real(dp)           DDOT
  EXTERNAL           LSAME, DDOT
!     ..
!     .. External Subroutines ..
  EXTERNAL           DCOPY, DSPMV, DSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          ABS
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
  INFO = 0
  UPPER = LSAME( UPLO, 'U' )
  IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
     INFO = -1
  ELSE IF( N.LT.0 ) THEN
     INFO = -2
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DSPTRI', -INFO )
     RETURN
  END IF
!
!     Quick return if possible
!
  IF( N.EQ.0 ) &
     RETURN
!
!     Check that the diagonal matrix D is nonsingular.
!
  IF( UPPER ) THEN
!
!        Upper triangular storage: examine D from bottom to top
!
     KP = N*( N+1 ) / 2
     DO 10 INFO = N, 1, -1
        IF( IPIV( INFO ).GT.0 .AND. AP( KP ).EQ.ZERO ) &
           RETURN
        KP = KP - INFO
10    CONTINUE
  ELSE
!
!        Lower triangular storage: examine D from top to bottom.
!
     KP = 1
     DO 20 INFO = 1, N
        IF( IPIV( INFO ).GT.0 .AND. AP( KP ).EQ.ZERO ) &
           RETURN
        KP = KP + N - INFO + 1
20    CONTINUE
  END IF
  INFO = 0
!
  IF( UPPER ) THEN
!
!        Compute inv(A) from the factorization A = U*D*U'.
!
!        K is the main loop index, increasing from 1 to N in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
     K = 1
     KC = 1
30    CONTINUE
!
!        If K > N, exit from loop.
!
     IF( K.GT.N ) &
        GO TO 50
!
     KCNEXT = KC + K
     IF( IPIV( K ).GT.0 ) THEN
!
!           1 x 1 diagonal block
!
!           Invert the diagonal block.
!
        AP( KC+K-1 ) = ONE / AP( KC+K-1 )
!
!           Compute column K of the inverse.
!
        IF( K.GT.1 ) THEN
           CALL DCOPY( K-1_i4, AP( KC ), 1_i4, WORK, 1_i4 )
           CALL DSPMV( UPLO, K-1_i4, -ONE, AP, WORK, 1_i4, ZERO, AP( KC ), &
                       1_i4 )
           AP( KC+K-1 ) = AP( KC+K-1 ) - &
                          DDOT( K-1_i4, WORK, 1_i4, AP( KC ), 1_i4 )
        END IF
        KSTEP = 1
     ELSE
!
!           2 x 2 diagonal block
!
!           Invert the diagonal block.
!
        T = ABS( AP( KCNEXT+K-1 ) )
        AK = AP( KC+K-1 ) / T
        AKP1 = AP( KCNEXT+K ) / T
        AKKP1 = AP( KCNEXT+K-1 ) / T
        D = T*( AK*AKP1-ONE )
        AP( KC+K-1 ) = AKP1 / D
        AP( KCNEXT+K ) = AK / D
        AP( KCNEXT+K-1 ) = -AKKP1 / D
!
!           Compute columns K and K+1 of the inverse.
!
        IF( K.GT.1 ) THEN
           CALL DCOPY( K-1_i4, AP( KC ), 1_i4, WORK, 1_i4 )
           CALL DSPMV( UPLO, K-1_i4, -ONE, AP, WORK, 1_i4, ZERO, AP( KC ), &
                       1_i4 )
           AP( KC+K-1 ) = AP( KC+K-1 ) - &
                          DDOT( K-1_i4, WORK, 1_i4, AP( KC ), 1_i4 )
           AP( KCNEXT+K-1 ) = AP( KCNEXT+K-1 ) - &
                              DDOT( K-1_i4, AP( KC ), 1_i4, AP( KCNEXT ), &
                              1_i4 )
           CALL DCOPY( K-1_i4, AP( KCNEXT ), 1_i4, WORK, 1_i4 )
           CALL DSPMV( UPLO, K-1_i4, -ONE, AP, WORK, 1_i4, ZERO, &
                       AP( KCNEXT ), 1_i4 )
           AP( KCNEXT+K ) = AP( KCNEXT+K ) - &
                            DDOT( K-1_i4, WORK, 1_i4, AP( KCNEXT ), 1_i4 )
        END IF
        KSTEP = 2
        KCNEXT = KCNEXT + K + 1
     END IF
!
     KP = ABS( IPIV( K ) )
     IF( KP.NE.K ) THEN
!
!           Interchange rows and columns K and KP in the leading
!           submatrix A(1:k+1,1:k+1)
!
        KPC = ( KP-1 )*KP / 2 + 1
        CALL DSWAP( KP-1_i4, AP( KC ), 1_i4, AP( KPC ), 1_i4 )
        KX = KPC + KP - 1
        DO 40 J = KP + 1, K - 1
           KX = KX + J - 1
           TEMP = AP( KC+J-1 )
           AP( KC+J-1 ) = AP( KX )
           AP( KX ) = TEMP
40       CONTINUE
        TEMP = AP( KC+K-1 )
        AP( KC+K-1 ) = AP( KPC+KP-1 )
        AP( KPC+KP-1 ) = TEMP
        IF( KSTEP.EQ.2 ) THEN
           TEMP = AP( KC+K+K-1 )
           AP( KC+K+K-1 ) = AP( KC+K+KP-1 )
           AP( KC+K+KP-1 ) = TEMP
        END IF
     END IF
!
     K = K + KSTEP
     KC = KCNEXT
     GO TO 30
50    CONTINUE
!
  ELSE
!
!        Compute inv(A) from the factorization A = L*D*L'.
!
!        K is the main loop index, increasing from 1 to N in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
     NPP = N*( N+1 ) / 2
     K = N
     KC = NPP
60    CONTINUE
!
!        If K < 1, exit from loop.
!
     IF( K.LT.1 ) &
        GO TO 80
!
     KCNEXT = KC - ( N-K+2 )
     IF( IPIV( K ).GT.0 ) THEN
!
!           1 x 1 diagonal block
!
!           Invert the diagonal block.
!
        AP( KC ) = ONE / AP( KC )
!
!           Compute column K of the inverse.
!
        IF( K.LT.N ) THEN
           CALL DCOPY( N-K, AP( KC+1 ), 1_i4, WORK, 1_i4 )
           CALL DSPMV( UPLO, N-K, -ONE, AP( KC+N-K+1 ), WORK, 1_i4, &
                       ZERO, AP( KC+1 ), 1_i4 )
           AP( KC ) = AP( KC ) - DDOT( N-K, WORK, 1_i4, AP( KC+1 ), 1_i4 )
        END IF
        KSTEP = 1
     ELSE
!
!           2 x 2 diagonal block
!
!           Invert the diagonal block.
!
        T = ABS( AP( KCNEXT+1 ) )
        AK = AP( KCNEXT ) / T
        AKP1 = AP( KC ) / T
        AKKP1 = AP( KCNEXT+1 ) / T
        D = T*( AK*AKP1-ONE )
        AP( KCNEXT ) = AKP1 / D
        AP( KC ) = AK / D
        AP( KCNEXT+1 ) = -AKKP1 / D
!
!           Compute columns K-1 and K of the inverse.
!
        IF( K.LT.N ) THEN
           CALL DCOPY( N-K, AP( KC+1 ), 1_i4, WORK, 1_i4 )
           CALL DSPMV( UPLO, N-K, -ONE, AP( KC+( N-K+1 ) ), WORK, 1_i4, &
                       ZERO, AP( KC+1 ), 1_i4 )
           AP( KC ) = AP( KC ) - DDOT( N-K, WORK, 1_i4, AP( KC+1 ), 1_i4 )
           AP( KCNEXT+1 ) = AP( KCNEXT+1 ) - &
                            DDOT( N-K, AP( KC+1 ), 1_i4, &
                            AP( KCNEXT+2 ), 1_i4 )
           CALL DCOPY( N-K, AP( KCNEXT+2 ), 1_i4, WORK, 1_i4 )
           CALL DSPMV( UPLO, N-K, -ONE, AP( KC+( N-K+1 ) ), WORK, 1_i4, &
                       ZERO, AP( KCNEXT+2 ), 1_i4 )
           AP( KCNEXT ) = AP( KCNEXT ) - &
                          DDOT( N-K, WORK, 1_i4, AP( KCNEXT+2 ), 1_i4 )
        END IF
        KSTEP = 2
        KCNEXT = KCNEXT - ( N-K+3 )
     END IF
!
     KP = ABS( IPIV( K ) )
     IF( KP.NE.K ) THEN
!
!           Interchange rows and columns K and KP in the trailing
!           submatrix A(k-1:n,k-1:n)
!
        KPC = NPP - ( N-KP+1 )*( N-KP+2 ) / 2 + 1
        IF( KP.LT.N ) &
           CALL DSWAP( N-KP, AP( KC+KP-K+1 ), 1_i4, AP( KPC+1 ), 1_i4 )
        KX = KC + KP - K
        DO 70 J = K + 1, KP - 1
           KX = KX + N - J + 1
           TEMP = AP( KC+J-K )
           AP( KC+J-K ) = AP( KX )
           AP( KX ) = TEMP
70       CONTINUE
        TEMP = AP( KC )
        AP( KC ) = AP( KPC )
        AP( KPC ) = TEMP
        IF( KSTEP.EQ.2 ) THEN
           TEMP = AP( KC-N+K-1 )
           AP( KC-N+K-1 ) = AP( KC-N+KP-1 )
           AP( KC-N+KP-1 ) = TEMP
        END IF
     END IF
!
     K = K - KSTEP
     KC = KCNEXT
     GO TO 60
80    CONTINUE
  END IF
!
  RETURN
!
!     End of DSPTRI
!
  END
!
  SUBROUTINE DSPEV( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, INFO )
!
!  -- LAPACK driver routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
  use datatypes
  CHARACTER          JOBZ, UPLO
  integer(i4)        INFO, LDZ, N
!     ..
!     .. Array Arguments ..
  real(dp)           AP( * ), W( * ), WORK( * ), Z( LDZ, * )
!     ..
!
!  Purpose
!  =======
!
!  DSPEV computes all the eigenvalues and, optionally, eigenvectors of a
!  real symmetric matrix A in packed storage.
!
!  Arguments
!  =========
!
!  JOBZ    (input) CHARACTER*1
!          = 'N':  Compute eigenvalues only;
!          = 'V':  Compute eigenvalues and eigenvectors.
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  AP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2)
!          On entry, the upper or lower triangle of the symmetric matrix
!          A, packed columnwise in a linear array.  The j-th column of A
!          is stored in the array AP as follows:
!          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
!
!          On exit, AP is overwritten by values generated during the
!          reduction to tridiagonal form.  If UPLO = 'U', the diagonal
!          and first superdiagonal of the tridiagonal matrix T overwrite
!          the corresponding elements of A, and if UPLO = 'L', the
!          diagonal and first subdiagonal of T overwrite the
!          corresponding elements of A.
!
!  W       (output) DOUBLE PRECISION array, dimension (N)
!          If INFO = 0, the eigenvalues in ascending order.
!
!  Z       (output) DOUBLE PRECISION array, dimension (LDZ, N)
!          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal
!          eigenvectors of the matrix A, with the i-th column of Z
!          holding the eigenvector associated with W(i).
!          If JOBZ = 'N', then Z is not referenced.
!
!  LDZ     (input) INTEGER
!          The leading dimension of the array Z.  LDZ >= 1, and if
!          JOBZ = 'V', LDZ >= max(1,N).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (3*N)
!
!  INFO    (output) INTEGER
!          = 0:  successful exit.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          > 0:  if INFO = i, the algorithm failed to converge; i
!                off-diagonal elements of an intermediate tridiagonal
!                form did not converge to zero.
!
!  =====================================================================
!
!     .. Parameters ..
  real(dp)           ZERO, ONE
  PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
  LOGICAL            WANTZ
  integer(i4)        IINFO, IMAX, INDE, INDTAU, INDWRK, ISCALE
  real(dp)           ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, &
                     SMLNUM
!     ..
!     .. External Functions ..
  LOGICAL            LSAME
  real(dp)           DLAMCH, DLANSP
  EXTERNAL           LSAME, DLAMCH, DLANSP
!     ..
!     .. External Subroutines ..
  EXTERNAL           DOPGTR, DSCAL, DSPTRD, DSTEQR, DSTERF, XERBLA
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
  WANTZ = LSAME( JOBZ, 'V' )
!
  INFO = 0
  IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
     INFO = -1
  ELSE IF( .NOT.( LSAME( UPLO, 'U' ) .OR. LSAME( UPLO, 'L' ) ) ) &
            THEN
     INFO = -2
  ELSE IF( N.LT.0 ) THEN
     INFO = -3
  ELSE IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) THEN
     INFO = -7
  END IF
!
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DSPEV ', -INFO )
     RETURN
  END IF
!
!     Quick return if possible
!
  IF( N.EQ.0 ) &
     RETURN
!
  IF( N.EQ.1 ) THEN
     W( 1 ) = AP( 1 )
     IF( WANTZ ) &
        Z( 1, 1 ) = ONE
     RETURN
  END IF
!
!     Get machine constants.
!
  SAFMIN = DLAMCH( 'Safe minimum' )
  EPS = DLAMCH( 'Precision' )
  SMLNUM = SAFMIN / EPS
  BIGNUM = ONE / SMLNUM
  RMIN = SQRT( SMLNUM )
  RMAX = SQRT( BIGNUM )
!
!     Scale matrix to allowable range, if necessary.
!
  ANRM = DLANSP( 'M', UPLO, N, AP, WORK )
  ISCALE = 0
  IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
     ISCALE = 1
     SIGMA = RMIN / ANRM
  ELSE IF( ANRM.GT.RMAX ) THEN
     ISCALE = 1
     SIGMA = RMAX / ANRM
  END IF
  IF( ISCALE.EQ.1 ) THEN
     CALL DSCAL( ( N*( N+1_i4 ) ) / 2_i4, SIGMA, AP, 1_i4 )
  END IF
!
!     Call DSPTRD to reduce symmetric packed matrix to tridiagonal form.
!
  INDE = 1
  INDTAU = INDE + N
  CALL DSPTRD( UPLO, N, AP, W, WORK( INDE ), WORK( INDTAU ), IINFO )
!
!     For eigenvalues only, call DSTERF.  For eigenvectors, first call
!     DOPGTR to generate the orthogonal matrix, then call DSTEQR.
!
  IF( .NOT.WANTZ ) THEN
     CALL DSTERF( N, W, WORK( INDE ), INFO )
  ELSE
     INDWRK = INDTAU + N
     CALL DOPGTR( UPLO, N, AP, WORK( INDTAU ), Z, LDZ, &
                  WORK( INDWRK ), IINFO )
     CALL DSTEQR( JOBZ, N, W, WORK( INDE ), Z, LDZ, WORK( INDTAU ), &
                  INFO )
  END IF
!
!     If matrix was scaled, then rescale eigenvalues appropriately.
!
  IF( ISCALE.EQ.1 ) THEN
     IF( INFO.EQ.0 ) THEN
        IMAX = N
     ELSE
        IMAX = INFO - 1
     END IF
     CALL DSCAL( IMAX, ONE / SIGMA, W, 1_i4 )
  END IF
!
  RETURN
!
!     End of DSPEV
!
  END
!
  FUNCTION DLANSP( NORM, UPLO, N, AP, WORK )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
  use datatypes
  CHARACTER          NORM, UPLO
  integer(i4)        N
!     ..
!     .. Array Arguments ..
  real(dp)           AP( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DLANSP  returns the value of the one norm,  or the Frobenius norm, or
!  the  infinity norm,  or the  element of  largest absolute value  of a
!  real symmetric matrix A,  supplied in packed form.
!
!  Description
!  ===========
!
!  DLANSP returns the value
!
!     DLANSP = ( max(abs(A(i,j))), NORM = 'M' or 'm'
!              (
!              ( norm1(A),         NORM = '1', 'O' or 'o'
!              (
!              ( normI(A),         NORM = 'I' or 'i'
!              (
!              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
!
!  where  norm1  denotes the  one norm of a matrix (maximum column sum),
!  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
!  normF  denotes the  Frobenius norm of a matrix (square root of sum of
!  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.
!
!  Arguments
!  =========
!
!  NORM    (input) CHARACTER*1
!          Specifies the value to be returned in DLANSP as described
!          above.
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the upper or lower triangular part of the
!          symmetric matrix A is supplied.
!          = 'U':  Upper triangular part of A is supplied
!          = 'L':  Lower triangular part of A is supplied
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.  When N = 0, DLANSP is
!          set to zero.
!
!  AP      (input) DOUBLE PRECISION array, dimension (N*(N+1)/2)
!          The upper or lower triangle of the symmetric matrix A, packed
!          columnwise in a linear array.  The j-th column of A is stored
!          in the array AP as follows:
!          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK),
!          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,
!          WORK is not referenced.
!
! =====================================================================
!
!     .. Parameters ..
  real(dp)           ONE, ZERO
  PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
  integer(i4)        I, J, K
  real(dp)           ABSA, SCALE, SUM, VALUE
!     ..
!     .. External Subroutines ..
  EXTERNAL           DLASSQ
!     ..
!     .. External Functions ..
  LOGICAL            LSAME
  EXTERNAL           LSAME
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          ABS, MAX, SQRT
  real(dp)           DLANSP
!     ..
!     .. Executable Statements ..
!
  IF( N.EQ.0 ) THEN
     VALUE = ZERO
  ELSE IF( LSAME( NORM, 'M' ) ) THEN
!
!        Find max(abs(A(i,j))).
!
     VALUE = ZERO
     IF( LSAME( UPLO, 'U' ) ) THEN
        K = 1
        DO 20 J = 1, N
           DO 10 I = K, K + J - 1
              VALUE = MAX( VALUE, ABS( AP( I ) ) )
10          CONTINUE
           K = K + J
20       CONTINUE
     ELSE
        K = 1
        DO 40 J = 1, N
           DO 30 I = K, K + N - J
              VALUE = MAX( VALUE, ABS( AP( I ) ) )
30          CONTINUE
           K = K + N - J + 1
40       CONTINUE
     END IF
  ELSE IF( ( LSAME( NORM, 'I' ) ) .OR. ( LSAME( NORM, 'O' ) ) .OR. &
           ( NORM.EQ.'1' ) ) THEN
!
!        Find normI(A) ( = norm1(A), since A is symmetric).
!
     VALUE = ZERO
     K = 1
     IF( LSAME( UPLO, 'U' ) ) THEN
        DO 60 J = 1, N
           SUM = ZERO
           DO 50 I = 1, J - 1
              ABSA = ABS( AP( K ) )
              SUM = SUM + ABSA
              WORK( I ) = WORK( I ) + ABSA
              K = K + 1
50          CONTINUE
           WORK( J ) = SUM + ABS( AP( K ) )
           K = K + 1
60       CONTINUE
        DO 70 I = 1, N
           VALUE = MAX( VALUE, WORK( I ) )
70       CONTINUE
     ELSE
        DO 80 I = 1, N
           WORK( I ) = ZERO
80       CONTINUE
        DO 100 J = 1, N
           SUM = WORK( J ) + ABS( AP( K ) )
           K = K + 1
           DO 90 I = J + 1, N
              ABSA = ABS( AP( K ) )
              SUM = SUM + ABSA
              WORK( I ) = WORK( I ) + ABSA
              K = K + 1
90          CONTINUE
           VALUE = MAX( VALUE, SUM )
100       CONTINUE
     END IF
  ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
!
!        Find normF(A).
!
     SCALE = ZERO
     SUM = ONE
     K = 2
     IF( LSAME( UPLO, 'U' ) ) THEN
        DO 110 J = 2, N
           CALL DLASSQ( J-1_i4, AP( K ), 1_i4, SCALE, SUM )
           K = K + J
110       CONTINUE
     ELSE
        DO 120 J = 1, N - 1
           CALL DLASSQ( N-J, AP( K ), 1_i4, SCALE, SUM )
           K = K + N - J + 1
120       CONTINUE
     END IF
     SUM = 2*SUM
     K = 1
     DO 130 I = 1, N
        IF( AP( K ).NE.ZERO ) THEN
           ABSA = ABS( AP( K ) )
           IF( SCALE.LT.ABSA ) THEN
              SUM = ONE + SUM*( SCALE / ABSA )**2
              SCALE = ABSA
           ELSE
              SUM = SUM + ( ABSA / SCALE )**2
           END IF
        END IF
        IF( LSAME( UPLO, 'U' ) ) THEN
           K = K + I + 1
        ELSE
           K = K + N - I + 1
        END IF
130    CONTINUE
     VALUE = SCALE*SQRT( SUM )
  END IF
!
  DLANSP = VALUE
  RETURN
!
!     End of DLANSP
!
  END
  SUBROUTINE DOPGTR( UPLO, N, AP, TAU, Q, LDQ, WORK, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
  use datatypes
  CHARACTER          UPLO
  integer(i4)        INFO, LDQ, N
!     ..
!     .. Array Arguments ..
  real(dp)           AP( * ), Q( LDQ, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DOPGTR generates a real orthogonal matrix Q which is defined as the
!  product of n-1 elementary reflectors H(i) of order n, as returned by
!  DSPTRD using packed storage:
!
!  if UPLO = 'U', Q = H(n-1) . . . H(2) H(1),
!
!  if UPLO = 'L', Q = H(1) H(2) . . . H(n-1).
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U': Upper triangular packed storage used in previous
!                 call to DSPTRD;
!          = 'L': Lower triangular packed storage used in previous
!                 call to DSPTRD.
!
!  N       (input) INTEGER
!          The order of the matrix Q. N >= 0.
!
!  AP      (input) DOUBLE PRECISION array, dimension (N*(N+1)/2)
!          The vectors which define the elementary reflectors, as
!          returned by DSPTRD.
!
!  TAU     (input) DOUBLE PRECISION array, dimension (N-1)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DSPTRD.
!
!  Q       (output) DOUBLE PRECISION array, dimension (LDQ,N)
!          The N-by-N orthogonal matrix Q.
!
!  LDQ     (input) INTEGER
!          The leading dimension of the array Q. LDQ >= max(1,N).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (N-1)
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
  real(dp)           ZERO, ONE
  PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
  LOGICAL            UPPER
  integer(i4)        I, IINFO, IJ, J
!     ..
!     .. External Functions ..
  LOGICAL            LSAME
  EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
  EXTERNAL           DORG2L, DORG2R, XERBLA
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
  INFO = 0
  UPPER = LSAME( UPLO, 'U' )
  IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
     INFO = -1
  ELSE IF( N.LT.0 ) THEN
     INFO = -2
  ELSE IF( LDQ.LT.MAX( 1_i4, N ) ) THEN
     INFO = -6
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DOPGTR', -INFO )
     RETURN
  END IF
!
!     Quick return if possible
!
  IF( N.EQ.0 ) &
     RETURN
!
  IF( UPPER ) THEN
!
!        Q was determined by a call to DSPTRD with UPLO = 'U'
!
!        Unpack the vectors which define the elementary reflectors and
!        set the last row and column of Q equal to those of the unit
!        matrix
!
     IJ = 2
     DO 20 J = 1, N - 1
        DO 10 I = 1, J - 1
           Q( I, J ) = AP( IJ )
           IJ = IJ + 1
10       CONTINUE
        IJ = IJ + 2
        Q( N, J ) = ZERO
20    CONTINUE
     DO 30 I = 1, N - 1
        Q( I, N ) = ZERO
30    CONTINUE
     Q( N, N ) = ONE
!
!        Generate Q(1:n-1,1:n-1)
!
     CALL DORG2L( N-1_i4, N-1_i4, N-1_i4, Q, LDQ, TAU, WORK, IINFO )
!
  ELSE
!
!        Q was determined by a call to DSPTRD with UPLO = 'L'.
!
!        Unpack the vectors which define the elementary reflectors and
!        set the first row and column of Q equal to those of the unit
!        matrix
!
     Q( 1, 1 ) = ONE
     DO 40 I = 2, N
        Q( I, 1 ) = ZERO
40    CONTINUE
     IJ = 3
     DO 60 J = 2, N
        Q( 1, J ) = ZERO
        DO 50 I = J + 1, N
           Q( I, J ) = AP( IJ )
           IJ = IJ + 1
50       CONTINUE
        IJ = IJ + 2
60    CONTINUE
     IF( N.GT.1 ) THEN
!
!           Generate Q(2:n,2:n)
!
        CALL DORG2R( N-1_i4, N-1_i4, N-1_i4, Q( 2, 2 ), LDQ, TAU, WORK, &
                     IINFO )
     END IF
  END IF
  RETURN
!
!     End of DOPGTR
!
  END
  SUBROUTINE DSPTRD( UPLO, N, AP, D, E, TAU, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
  use datatypes
  CHARACTER          UPLO
  integer(i4)        INFO, N
!     ..
!     .. Array Arguments ..
  real(dp)           AP( * ), D( * ), E( * ), TAU( * )
!     ..
!
!  Purpose
!  =======
!
!  DSPTRD reduces a real symmetric matrix A stored in packed form to
!  symmetric tridiagonal form T by an orthogonal similarity
!  transformation: Q**T * A * Q = T.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  AP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2)
!          On entry, the upper or lower triangle of the symmetric matrix
!          A, packed columnwise in a linear array.  The j-th column of A
!          is stored in the array AP as follows:
!          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
!          On exit, if UPLO = 'U', the diagonal and first superdiagonal
!          of A are overwritten by the corresponding elements of the
!          tridiagonal matrix T, and the elements above the first
!          superdiagonal, with the array TAU, represent the orthogonal
!          matrix Q as a product of elementary reflectors; if UPLO
!          = 'L', the diagonal and first subdiagonal of A are over-
!          written by the corresponding elements of the tridiagonal
!          matrix T, and the elements below the first subdiagonal, with
!          the array TAU, represent the orthogonal matrix Q as a product
!          of elementary reflectors. See Further Details.
!
!  D       (output) DOUBLE PRECISION array, dimension (N)
!          The diagonal elements of the tridiagonal matrix T:
!          D(i) = A(i,i).
!
!  E       (output) DOUBLE PRECISION array, dimension (N-1)
!          The off-diagonal elements of the tridiagonal matrix T:
!          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
!
!  TAU     (output) DOUBLE PRECISION array, dimension (N-1)
!          The scalar factors of the elementary reflectors (see Further
!          Details).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  Further Details
!  ===============
!
!  If UPLO = 'U', the matrix Q is represented as a product of elementary
!  reflectors
!
!     Q = H(n-1) . . . H(2) H(1).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in AP,
!  overwriting A(1:i-1,i+1), and tau is stored in TAU(i).
!
!  If UPLO = 'L', the matrix Q is represented as a product of elementary
!  reflectors
!
!     Q = H(1) H(2) . . . H(n-1).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in AP,
!  overwriting A(i+2:n,i), and tau is stored in TAU(i).
!
!  =====================================================================
!
!     .. Parameters ..
  real(dp)           ONE, ZERO, HALF
  PARAMETER          ( ONE = 1.0D0, ZERO = 0.0D0, &
                     HALF = 1.0D0 / 2.0D0 )
!     ..
!     .. Local Scalars ..
  LOGICAL            UPPER
  integer(i4)        I, I1, I1I1, II
  real(dp)           ALPHA, TAUI
!     ..
!     .. External Subroutines ..
  EXTERNAL           DAXPY, DLARFG, DSPMV, DSPR2, XERBLA
!     ..
!     .. External Functions ..
  LOGICAL            LSAME
  real(dp)           DDOT
  EXTERNAL           LSAME, DDOT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
  INFO = 0
  UPPER = LSAME( UPLO, 'U' )
  IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
     INFO = -1
  ELSE IF( N.LT.0 ) THEN
     INFO = -2
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DSPTRD', -INFO )
     RETURN
  END IF
!
!     Quick return if possible
!
  IF( N.LE.0 ) &
     RETURN
!
  IF( UPPER ) THEN
!
!        Reduce the upper triangle of A.
!        I1 is the index in AP of A(1,I+1).
!
     I1 = N*( N-1 ) / 2 + 1
     DO 10 I = N - 1, 1, -1
!
!           Generate elementary reflector H(i) = I - tau * v * v'
!           to annihilate A(1:i-1,i+1)
!
        CALL DLARFG( I, AP( I1+I-1 ), AP( I1 ), 1_i4, TAUI )
        E( I ) = AP( I1+I-1 )
!
        IF( TAUI.NE.ZERO ) THEN
!
!              Apply H(i) from both sides to A(1:i,1:i)
!
           AP( I1+I-1 ) = ONE
!
!              Compute  y := tau * A * v  storing y in TAU(1:i)
!
           CALL DSPMV( UPLO, I, TAUI, AP, AP( I1 ), 1_i4, ZERO, TAU, &
                       1_i4 )
!
!              Compute  w := y - 1/2 * tau * (y'*v) * v
!
           ALPHA = -HALF*TAUI*DDOT( I, TAU, 1_i4, AP( I1 ), 1_i4 )
           CALL DAXPY( I, ALPHA, AP( I1 ), 1_i4, TAU, 1_i4 )
!
!              Apply the transformation as a rank-2 update:
!                 A := A - v * w' - w * v'
!
           CALL DSPR2( UPLO, I, -ONE, AP( I1 ), 1_i4, TAU, 1_i4, AP )
!
           AP( I1+I-1 ) = E( I )
        END IF
        D( I+1 ) = AP( I1+I )
        TAU( I ) = TAUI
        I1 = I1 - I
10    CONTINUE
     D( 1 ) = AP( 1 )
  ELSE
!
!        Reduce the lower triangle of A. II is the index in AP of
!        A(i,i) and I1I1 is the index of A(i+1,i+1).
!
     II = 1
     DO 20 I = 1, N - 1
        I1I1 = II + N - I + 1
!
!           Generate elementary reflector H(i) = I - tau * v * v'
!           to annihilate A(i+2:n,i)
!
        CALL DLARFG( N-I, AP( II+1 ), AP( II+2 ), 1_i4, TAUI )
        E( I ) = AP( II+1 )
!
        IF( TAUI.NE.ZERO ) THEN
!
!              Apply H(i) from both sides to A(i+1:n,i+1:n)
!
           AP( II+1 ) = ONE
!
!              Compute  y := tau * A * v  storing y in TAU(i:n-1)
!
           CALL DSPMV( UPLO, N-I, TAUI, AP( I1I1 ), AP( II+1 ), 1_i4, &
                       ZERO, TAU( I ), 1_i4 )
!
!              Compute  w := y - 1/2 * tau * (y'*v) * v
!
           ALPHA = -HALF*TAUI*DDOT( N-I, TAU( I ), 1_i4, AP( II+1 ), &
                   1_i4 )
           CALL DAXPY( N-I, ALPHA, AP( II+1 ), 1_i4, TAU( I ), 1_i4 )
!
!              Apply the transformation as a rank-2 update:
!                 A := A - v * w' - w * v'
!
           CALL DSPR2( UPLO, N-I, -ONE, AP( II+1 ), 1_i4, TAU( I ), 1_i4, &
                       AP( I1I1 ) )
!
           AP( II+1 ) = E( I )
        END IF
        D( I ) = AP( II )
        TAU( I ) = TAUI
        II = I1I1
20    CONTINUE
     D( N ) = AP( II )
  END IF
!
  RETURN
!
!     End of DSPTRD
!
  END
  SUBROUTINE DSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
  use datatypes
  CHARACTER          COMPZ
  integer(i4)        INFO, LDZ, N
!     ..
!     .. Array Arguments ..
  real(dp)           D( * ), E( * ), WORK( * ), Z( LDZ, * )
!     ..
!
!  Purpose
!  =======
!
!  DSTEQR computes all eigenvalues and, optionally, eigenvectors of a
!  symmetric tridiagonal matrix using the implicit QL or QR method.
!  The eigenvectors of a full or band symmetric matrix can also be found
!  if DSYTRD or DSPTRD or DSBTRD has been used to reduce this matrix to
!  tridiagonal form.
!
!  Arguments
!  =========
!
!  COMPZ   (input) CHARACTER*1
!          = 'N':  Compute eigenvalues only.
!          = 'V':  Compute eigenvalues and eigenvectors of the original
!                  symmetric matrix.  On entry, Z must contain the
!                  orthogonal matrix used to reduce the original matrix
!                  to tridiagonal form.
!          = 'I':  Compute eigenvalues and eigenvectors of the
!                  tridiagonal matrix.  Z is initialized to the identity
!                  matrix.
!
!  N       (input) INTEGER
!          The order of the matrix.  N >= 0.
!
!  D       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, the diagonal elements of the tridiagonal matrix.
!          On exit, if INFO = 0, the eigenvalues in ascending order.
!
!  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
!          On entry, the (n-1) subdiagonal elements of the tridiagonal
!          matrix.
!          On exit, E has been destroyed.
!
!  Z       (input/output) DOUBLE PRECISION array, dimension (LDZ, N)
!          On entry, if  COMPZ = 'V', then Z contains the orthogonal
!          matrix used in the reduction to tridiagonal form.
!          On exit, if INFO = 0, then if  COMPZ = 'V', Z contains the
!          orthonormal eigenvectors of the original symmetric matrix,
!          and if COMPZ = 'I', Z contains the orthonormal eigenvectors
!          of the symmetric tridiagonal matrix.
!          If COMPZ = 'N', then Z is not referenced.
!
!  LDZ     (input) INTEGER
!          The leading dimension of the array Z.  LDZ >= 1, and if
!          eigenvectors are desired, then  LDZ >= max(1,N).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (max(1,2*N-2))
!          If COMPZ = 'N', then WORK is not referenced.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  the algorithm has failed to find all the eigenvalues in
!                a total of 30*N iterations; if INFO = i, then i
!                elements of E have not converged to zero; on exit, D
!                and E contain the elements of a symmetric tridiagonal
!                matrix which is orthogonally similar to the original
!                matrix.
!
!  =====================================================================
!
!     .. Parameters ..
  real(dp)           ZERO, ONE, TWO, THREE
  PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, &
                     THREE = 3.0D0 )
  integer(i4)        MAXIT
  PARAMETER          ( MAXIT = 30 )
!     ..
!     .. Local Scalars ..
  integer(i4)        I, ICOMPZ, II, ISCALE, J, JTOT, K, L, L1, LEND, &
                     LENDM1, LENDP1, LENDSV, LM1, LSV, M, MM, MM1, &
                     NM1, NMAXIT
  real(dp)           ANORM, B, C, EPS, EPS2, F, G, P, R, RT1, RT2, &
                     S, SAFMAX, SAFMIN, SSFMAX, SSFMIN, TST
!     ..
!     .. External Functions ..
  LOGICAL            LSAME
  real(dp)           DLAMCH, DLANST, DLAPY2
  EXTERNAL           LSAME, DLAMCH, DLANST, DLAPY2
!     ..
!     .. External Subroutines ..
  EXTERNAL           DLAE2, DLAEV2, DLARTG, DLASCL, DLASET, DLASR, &
                     DLASRT, DSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          ABS, MAX, SIGN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
  INFO = 0
!
  IF( LSAME( COMPZ, 'N' ) ) THEN
     ICOMPZ = 0
  ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
     ICOMPZ = 1
  ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
     ICOMPZ = 2
  ELSE
     ICOMPZ = -1
  END IF
  IF( ICOMPZ.LT.0 ) THEN
     INFO = -1
  ELSE IF( N.LT.0 ) THEN
     INFO = -2
  ELSE IF( ( LDZ.LT.1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1_i4, &
           N ) ) ) THEN
     INFO = -6
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DSTEQR', -INFO )
     RETURN
  END IF
!
!     Quick return if possible
!
  IF( N.EQ.0 ) &
     RETURN
!
  IF( N.EQ.1 ) THEN
     IF( ICOMPZ.EQ.2 ) &
        Z( 1, 1 ) = ONE
     RETURN
  END IF
!
!     Determine the unit roundoff and over/underflow thresholds.
!
  EPS = DLAMCH( 'E' )
  EPS2 = EPS**2
  SAFMIN = DLAMCH( 'S' )
  SAFMAX = ONE / SAFMIN
  SSFMAX = SQRT( SAFMAX ) / THREE
  SSFMIN = SQRT( SAFMIN ) / EPS2
!
!     Compute the eigenvalues and eigenvectors of the tridiagonal
!     matrix.
!
  IF( ICOMPZ.EQ.2 ) &
     CALL DLASET( 'Full', N, N, ZERO, ONE, Z, LDZ )
!
  NMAXIT = N*MAXIT
  JTOT = 0
!
!     Determine where the matrix splits and choose QL or QR iteration
!     for each block, according to whether top or bottom diagonal
!     element is smaller.
!
  L1 = 1
  NM1 = N - 1
!
10 CONTINUE
  IF( L1.GT.N ) &
     GO TO 160
  IF( L1.GT.1 ) &
     E( L1-1 ) = ZERO
  IF( L1.LE.NM1 ) THEN
     DO 20 M = L1, NM1
        TST = ABS( E( M ) )
        IF( TST.EQ.ZERO ) &
           GO TO 30
        IF( TST.LE.( SQRT( ABS( D( M ) ) )*SQRT( ABS( D( M+ &
            1 ) ) ) )*EPS ) THEN
           E( M ) = ZERO
           GO TO 30
        END IF
20    CONTINUE
  END IF
  M = N
!
30 CONTINUE
  L = L1
  LSV = L
  LEND = M
  LENDSV = LEND
  L1 = M + 1
  IF( LEND.EQ.L ) &
     GO TO 10
!
!     Scale submatrix in rows and columns L to LEND
!
  ANORM = DLANST( 'I', LEND-L+1_i4, D( L ), E( L ) )
  ISCALE = 0
  IF( ANORM.EQ.ZERO ) &
     GO TO 10
  IF( ANORM.GT.SSFMAX ) THEN
     ISCALE = 1
     CALL DLASCL( 'G', 0_i4, 0_i4, ANORM, SSFMAX, LEND-L+1_i4, 1_i4, D( L ), N, &
                  INFO )
     CALL DLASCL( 'G', 0_i4, 0_i4, ANORM, SSFMAX, LEND-L, 1_i4, E( L ), N, &
                  INFO )
  ELSE IF( ANORM.LT.SSFMIN ) THEN
     ISCALE = 2
     CALL DLASCL( 'G', 0_i4, 0_i4, ANORM, SSFMIN, LEND-L+1_i4, 1_i4, D( L ), N, &
                  INFO )
     CALL DLASCL( 'G', 0_i4, 0_i4, ANORM, SSFMIN, LEND-L, 1_i4, E( L ), N, &
                  INFO )
  END IF
!
!     Choose between QL and QR iteration
!
  IF( ABS( D( LEND ) ).LT.ABS( D( L ) ) ) THEN
     LEND = LSV
     L = LENDSV
  END IF
!
  IF( LEND.GT.L ) THEN
!
!        QL Iteration
!
!        Look for small subdiagonal element.
!
40    CONTINUE
     IF( L.NE.LEND ) THEN
        LENDM1 = LEND - 1
        DO 50 M = L, LENDM1
           TST = ABS( E( M ) )**2
           IF( TST.LE.( EPS2*ABS( D( M ) ) )*ABS( D( M+1 ) )+ &
               SAFMIN )GO TO 60
50       CONTINUE
     END IF
!
     M = LEND
!
60    CONTINUE
     IF( M.LT.LEND ) &
        E( M ) = ZERO
     P = D( L )
     IF( M.EQ.L ) &
        GO TO 80
!
!        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
!        to compute its eigensystem.
!
     IF( M.EQ.L+1 ) THEN
        IF( ICOMPZ.GT.0 ) THEN
           CALL DLAEV2( D( L ), E( L ), D( L+1 ), RT1, RT2, C, S )
           WORK( L ) = C
           WORK( N-1+L ) = S
           CALL DLASR( 'R', 'V', 'B', N, 2_i4, WORK( L ), &
                       WORK( N-1+L ), Z( 1, L ), LDZ )
        ELSE
           CALL DLAE2( D( L ), E( L ), D( L+1 ), RT1, RT2 )
        END IF
        D( L ) = RT1
        D( L+1 ) = RT2
        E( L ) = ZERO
        L = L + 2
        IF( L.LE.LEND ) &
           GO TO 40
        GO TO 140
     END IF
!
     IF( JTOT.EQ.NMAXIT ) &
        GO TO 140
     JTOT = JTOT + 1
!
!        Form shift.
!
     G = ( D( L+1 )-P ) / ( TWO*E( L ) )
     R = DLAPY2( G, ONE )
     G = D( M ) - P + ( E( L ) / ( G+SIGN( R, G ) ) )
!
     S = ONE
     C = ONE
     P = ZERO
!
!        Inner loop
!
     MM1 = M - 1
     DO 70 I = MM1, L, -1
        F = S*E( I )
        B = C*E( I )
        CALL DLARTG( G, F, C, S, R )
        IF( I.NE.M-1 ) &
           E( I+1 ) = R
        G = D( I+1 ) - P
        R = ( D( I )-G )*S + TWO*C*B
        P = S*R
        D( I+1 ) = G + P
        G = C*R - B
!
!           If eigenvectors are desired, then save rotations.
!
        IF( ICOMPZ.GT.0 ) THEN
           WORK( I ) = C
           WORK( N-1+I ) = -S
        END IF
!
70    CONTINUE
!
!        If eigenvectors are desired, then apply saved rotations.
!
     IF( ICOMPZ.GT.0 ) THEN
        MM = M - L + 1
        CALL DLASR( 'R', 'V', 'B', N, MM, WORK( L ), WORK( N-1+L ), &
                    Z( 1, L ), LDZ )
     END IF
!
     D( L ) = D( L ) - P
     E( L ) = G
     GO TO 40
!
!        Eigenvalue found.
!
80    CONTINUE
     D( L ) = P
!
     L = L + 1
     IF( L.LE.LEND ) &
        GO TO 40
     GO TO 140
!
  ELSE
!
!        QR Iteration
!
!        Look for small superdiagonal element.
!
90    CONTINUE
     IF( L.NE.LEND ) THEN
        LENDP1 = LEND + 1
        DO 100 M = L, LENDP1, -1
           TST = ABS( E( M-1 ) )**2
           IF( TST.LE.( EPS2*ABS( D( M ) ) )*ABS( D( M-1 ) )+ &
               SAFMIN )GO TO 110
100       CONTINUE
     END IF
!
     M = LEND
!
110    CONTINUE
     IF( M.GT.LEND ) &
        E( M-1 ) = ZERO
     P = D( L )
     IF( M.EQ.L ) &
        GO TO 130
!
!        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
!        to compute its eigensystem.
!
     IF( M.EQ.L-1 ) THEN
        IF( ICOMPZ.GT.0 ) THEN
           CALL DLAEV2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2, C, S )
           WORK( M ) = C
           WORK( N-1+M ) = S
           CALL DLASR( 'R', 'V', 'F', N, 2_i4, WORK( M ), &
                       WORK( N-1+M ), Z( 1, L-1 ), LDZ )
        ELSE
           CALL DLAE2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2 )
        END IF
        D( L-1 ) = RT1
        D( L ) = RT2
        E( L-1 ) = ZERO
        L = L - 2
        IF( L.GE.LEND ) &
           GO TO 90
        GO TO 140
     END IF
!
     IF( JTOT.EQ.NMAXIT ) &
        GO TO 140
     JTOT = JTOT + 1
!
!        Form shift.
!
     G = ( D( L-1 )-P ) / ( TWO*E( L-1 ) )
     R = DLAPY2( G, ONE )
     G = D( M ) - P + ( E( L-1 ) / ( G+SIGN( R, G ) ) )
!
     S = ONE
     C = ONE
     P = ZERO
!
!        Inner loop
!
     LM1 = L - 1
     DO 120 I = M, LM1
        F = S*E( I )
        B = C*E( I )
        CALL DLARTG( G, F, C, S, R )
        IF( I.NE.M ) &
           E( I-1 ) = R
        G = D( I ) - P
        R = ( D( I+1 )-G )*S + TWO*C*B
        P = S*R
        D( I ) = G + P
        G = C*R - B
!
!           If eigenvectors are desired, then save rotations.
!
        IF( ICOMPZ.GT.0 ) THEN
           WORK( I ) = C
           WORK( N-1+I ) = S
        END IF
!
120    CONTINUE
!
!        If eigenvectors are desired, then apply saved rotations.
!
     IF( ICOMPZ.GT.0 ) THEN
        MM = L - M + 1
        CALL DLASR( 'R', 'V', 'F', N, MM, WORK( M ), WORK( N-1+M ), &
                    Z( 1, M ), LDZ )
     END IF
!
     D( L ) = D( L ) - P
     E( LM1 ) = G
     GO TO 90
!
!        Eigenvalue found.
!
130    CONTINUE
     D( L ) = P
!
     L = L - 1
     IF( L.GE.LEND ) &
        GO TO 90
     GO TO 140
!
  END IF
!
!     Undo scaling if necessary
!
140 CONTINUE
  IF( ISCALE.EQ.1 ) THEN
     CALL DLASCL( 'G', 0_i4, 0_i4, SSFMAX, ANORM, LENDSV-LSV+1_i4, 1_i4, &
                  D( LSV ), N, INFO )
     CALL DLASCL( 'G', 0_i4, 0_i4, SSFMAX, ANORM, LENDSV-LSV, 1_i4, E( LSV ), &
                  N, INFO )
  ELSE IF( ISCALE.EQ.2 ) THEN
     CALL DLASCL( 'G', 0_i4, 0_i4, SSFMIN, ANORM, LENDSV-LSV+1_i4, 1_i4, &
                  D( LSV ), N, INFO )
     CALL DLASCL( 'G', 0_i4, 0_i4, SSFMIN, ANORM, LENDSV-LSV, 1_i4, E( LSV ), &
                  N, INFO )
  END IF
!
!     Check for no convergence to an eigenvalue after a total
!     of N*MAXIT iterations.
!
  IF( JTOT.LT.NMAXIT ) &
     GO TO 10
  DO 150 I = 1, N - 1
     IF( E( I ).NE.ZERO ) &
        INFO = INFO + 1
150 CONTINUE
  GO TO 190
!
!     Order eigenvalues and eigenvectors.
!
160 CONTINUE
  IF( ICOMPZ.EQ.0 ) THEN
!
!        Use Quick Sort
!
     CALL DLASRT( 'I', N, D, INFO )
!
  ELSE
!
!        Use Selection Sort to minimize swaps of eigenvectors
!
     DO 180 II = 2, N
        I = II - 1
        K = I
        P = D( I )
        DO 170 J = II, N
           IF( D( J ).LT.P ) THEN
              K = J
              P = D( J )
           END IF
170       CONTINUE
        IF( K.NE.I ) THEN
           D( K ) = D( I )
           D( I ) = P
           CALL DSWAP( N, Z( 1, I ), 1_i4, Z( 1, K ), 1_i4 )
        END IF
180    CONTINUE
  END IF
!
190 CONTINUE
  RETURN
!
!     End of DSTEQR
!
  END
  SUBROUTINE DSTERF( N, D, E, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
  use datatypes
  integer(i4)        INFO, N
!     ..
!     .. Array Arguments ..
  real(dp)           D( * ), E( * )
!     ..
!
!  Purpose
!  =======
!
!  DSTERF computes all eigenvalues of a symmetric tridiagonal matrix
!  using the Pal-Walker-Kahan variant of the QL or QR algorithm.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the matrix.  N >= 0.
!
!  D       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, the n diagonal elements of the tridiagonal matrix.
!          On exit, if INFO = 0, the eigenvalues in ascending order.
!
!  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
!          On entry, the (n-1) subdiagonal elements of the tridiagonal
!          matrix.
!          On exit, E has been destroyed.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  the algorithm failed to find all of the eigenvalues in
!                a total of 30*N iterations; if INFO = i, then i
!                elements of E have not converged to zero.
!
!  =====================================================================
!
!     .. Parameters ..
  real(dp)           ZERO, ONE, TWO, THREE
  PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, &
                     THREE = 3.0D0 )
  integer(i4)        MAXIT
  PARAMETER          ( MAXIT = 30 )
!     ..
!     .. Local Scalars ..
  integer(i4)        I, ISCALE, JTOT, L, L1, LEND, LENDSV, LSV, M, &
                     NMAXIT
  real(dp)           ALPHA, ANORM, BB, C, EPS, EPS2, GAMMA, OLDC, &
                     OLDGAM, P, R, RT1, RT2, RTE, S, SAFMAX, SAFMIN, &
                     SIGMA, SSFMAX, SSFMIN
!     ..
!     .. External Functions ..
  real(dp)           DLAMCH, DLANST, DLAPY2
  EXTERNAL           DLAMCH, DLANST, DLAPY2
!     ..
!     .. External Subroutines ..
  EXTERNAL           DLAE2, DLASCL, DLASRT, XERBLA
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          ABS, SIGN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
  INFO = 0
!
!     Quick return if possible
!
  IF( N.LT.0 ) THEN
     INFO = -1
     CALL XERBLA( 'DSTERF', -INFO )
     RETURN
  END IF
  IF( N.LE.1 ) &
     RETURN
!
!     Determine the unit roundoff for this environment.
!
  EPS = DLAMCH( 'E' )
  EPS2 = EPS**2
  SAFMIN = DLAMCH( 'S' )
  SAFMAX = ONE / SAFMIN
  SSFMAX = SQRT( SAFMAX ) / THREE
  SSFMIN = SQRT( SAFMIN ) / EPS2
!
!     Compute the eigenvalues of the tridiagonal matrix.
!
  NMAXIT = N*MAXIT
  SIGMA = ZERO
  JTOT = 0
!
!     Determine where the matrix splits and choose QL or QR iteration
!     for each block, according to whether top or bottom diagonal
!     element is smaller.
!
  L1 = 1
!
10 CONTINUE
  IF( L1.GT.N ) &
     GO TO 170
  IF( L1.GT.1 ) &
     E( L1-1 ) = ZERO
  DO 20 M = L1, N - 1
     IF( ABS( E( M ) ).LE.( SQRT( ABS( D( M ) ) )*SQRT( ABS( D( M+ &
         1 ) ) ) )*EPS ) THEN
        E( M ) = ZERO
        GO TO 30
     END IF
20 CONTINUE
  M = N
!
30 CONTINUE
  L = L1
  LSV = L
  LEND = M
  LENDSV = LEND
  L1 = M + 1
  IF( LEND.EQ.L ) &
     GO TO 10
!
!     Scale submatrix in rows and columns L to LEND
!
  ANORM = DLANST( 'I', LEND-L+1_i4, D( L ), E( L ) )
  ISCALE = 0
  IF( ANORM.GT.SSFMAX ) THEN
     ISCALE = 1
     CALL DLASCL( 'G', 0_i4, 0_i4, ANORM, SSFMAX, LEND-L+1_i4, 1_i4, D( L ), N, &
                  INFO )
     CALL DLASCL( 'G', 0_i4, 0_i4, ANORM, SSFMAX, LEND-L, 1_i4, E( L ), N, &
                  INFO )
  ELSE IF( ANORM.LT.SSFMIN ) THEN
     ISCALE = 2
     CALL DLASCL( 'G', 0_i4, 0_i4, ANORM, SSFMIN, LEND-L+1_i4, 1_i4, D( L ), N, &
                  INFO )
     CALL DLASCL( 'G', 0_i4, 0_i4, ANORM, SSFMIN, LEND-L, 1_i4, E( L ), N, &
                  INFO )
  END IF
!
  DO 40 I = L, LEND - 1
     E( I ) = E( I )**2
40 CONTINUE
!
!     Choose between QL and QR iteration
!
  IF( ABS( D( LEND ) ).LT.ABS( D( L ) ) ) THEN
     LEND = LSV
     L = LENDSV
  END IF
!
  IF( LEND.GE.L ) THEN
!
!        QL Iteration
!
!        Look for small subdiagonal element.
!
50    CONTINUE
     IF( L.NE.LEND ) THEN
        DO 60 M = L, LEND - 1
           IF( ABS( E( M ) ).LE.EPS2*ABS( D( M )*D( M+1 ) ) ) &
              GO TO 70
60       CONTINUE
     END IF
     M = LEND
!
70    CONTINUE
     IF( M.LT.LEND ) &
        E( M ) = ZERO
     P = D( L )
     IF( M.EQ.L ) &
        GO TO 90
!
!        If remaining matrix is 2 by 2, use DLAE2 to compute its
!        eigenvalues.
!
     IF( M.EQ.L+1 ) THEN
        RTE = SQRT( E( L ) )
        CALL DLAE2( D( L ), RTE, D( L+1 ), RT1, RT2 )
        D( L ) = RT1
        D( L+1 ) = RT2
        E( L ) = ZERO
        L = L + 2
        IF( L.LE.LEND ) &
           GO TO 50
        GO TO 150
     END IF
!
     IF( JTOT.EQ.NMAXIT ) &
        GO TO 150
     JTOT = JTOT + 1
!
!        Form shift.
!
     RTE = SQRT( E( L ) )
     SIGMA = ( D( L+1 )-P ) / ( TWO*RTE )
     R = DLAPY2( SIGMA, ONE )
     SIGMA = P - ( RTE / ( SIGMA+SIGN( R, SIGMA ) ) )
!
     C = ONE
     S = ZERO
     GAMMA = D( M ) - SIGMA
     P = GAMMA*GAMMA
!
!        Inner loop
!
     DO 80 I = M - 1, L, -1
        BB = E( I )
        R = P + BB
        IF( I.NE.M-1 ) &
           E( I+1 ) = S*R
        OLDC = C
        C = P / R
        S = BB / R
        OLDGAM = GAMMA
        ALPHA = D( I )
        GAMMA = C*( ALPHA-SIGMA ) - S*OLDGAM
        D( I+1 ) = OLDGAM + ( ALPHA-GAMMA )
        IF( C.NE.ZERO ) THEN
           P = ( GAMMA*GAMMA ) / C
        ELSE
           P = OLDC*BB
        END IF
80    CONTINUE
!
     E( L ) = S*P
     D( L ) = SIGMA + GAMMA
     GO TO 50
!
!        Eigenvalue found.
!
90    CONTINUE
     D( L ) = P
!
     L = L + 1
     IF( L.LE.LEND ) &
        GO TO 50
     GO TO 150
!
  ELSE
!
!        QR Iteration
!
!        Look for small superdiagonal element.
!
100    CONTINUE
     DO 110 M = L, LEND + 1, -1
        IF( ABS( E( M-1 ) ).LE.EPS2*ABS( D( M )*D( M-1 ) ) ) &
           GO TO 120
110    CONTINUE
     M = LEND
!
120    CONTINUE
     IF( M.GT.LEND ) &
        E( M-1 ) = ZERO
     P = D( L )
     IF( M.EQ.L ) &
        GO TO 140
!
!        If remaining matrix is 2 by 2, use DLAE2 to compute its
!        eigenvalues.
!
     IF( M.EQ.L-1 ) THEN
        RTE = SQRT( E( L-1 ) )
        CALL DLAE2( D( L ), RTE, D( L-1 ), RT1, RT2 )
        D( L ) = RT1
        D( L-1 ) = RT2
        E( L-1 ) = ZERO
        L = L - 2
        IF( L.GE.LEND ) &
           GO TO 100
        GO TO 150
     END IF
!
     IF( JTOT.EQ.NMAXIT ) &
        GO TO 150
     JTOT = JTOT + 1
!
!        Form shift.
!
     RTE = SQRT( E( L-1 ) )
     SIGMA = ( D( L-1 )-P ) / ( TWO*RTE )
     R = DLAPY2( SIGMA, ONE )
     SIGMA = P - ( RTE / ( SIGMA+SIGN( R, SIGMA ) ) )
!
     C = ONE
     S = ZERO
     GAMMA = D( M ) - SIGMA
     P = GAMMA*GAMMA
!
!        Inner loop
!
     DO 130 I = M, L - 1
        BB = E( I )
        R = P + BB
        IF( I.NE.M ) &
           E( I-1 ) = S*R
        OLDC = C
        C = P / R
        S = BB / R
        OLDGAM = GAMMA
        ALPHA = D( I+1 )
        GAMMA = C*( ALPHA-SIGMA ) - S*OLDGAM
        D( I ) = OLDGAM + ( ALPHA-GAMMA )
        IF( C.NE.ZERO ) THEN
           P = ( GAMMA*GAMMA ) / C
        ELSE
           P = OLDC*BB
        END IF
130    CONTINUE
!
     E( L-1 ) = S*P
     D( L ) = SIGMA + GAMMA
     GO TO 100
!
!        Eigenvalue found.
!
140    CONTINUE
     D( L ) = P
!
     L = L - 1
     IF( L.GE.LEND ) &
        GO TO 100
     GO TO 150
!
  END IF
!
!     Undo scaling if necessary
!
150 CONTINUE
  IF( ISCALE.EQ.1 ) &
     CALL DLASCL( 'G', 0_i4, 0_i4, SSFMAX, ANORM, LENDSV-LSV+1_i4, 1_i4, &
                  D( LSV ), N, INFO )
  IF( ISCALE.EQ.2 ) &
     CALL DLASCL( 'G', 0_i4, 0_i4, SSFMIN, ANORM, LENDSV-LSV+1_i4, 1_i4, &
                  D( LSV ), N, INFO )
!
!     Check for no convergence to an eigenvalue after a total
!     of N*MAXIT iterations.
!
  IF( JTOT.LT.NMAXIT ) &
     GO TO 10
  DO 160 I = 1, N - 1
     IF( E( I ).NE.ZERO ) &
        INFO = INFO + 1
160 CONTINUE
  GO TO 180
!
!     Sort eigenvalues in increasing order.
!
170 CONTINUE
  CALL DLASRT( 'I', N, D, INFO )
!
180 CONTINUE
  RETURN
!
!     End of DSTERF
!
  END
  SUBROUTINE DLAE2( A, B, C, RT1, RT2 )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
  use datatypes
  real(dp)           A, B, C, RT1, RT2
!     ..
!
!  Purpose
!  =======
!
!  DLAE2  computes the eigenvalues of a 2-by-2 symmetric matrix
!     [  A   B  ]
!     [  B   C  ].
!  On return, RT1 is the eigenvalue of larger absolute value, and RT2
!  is the eigenvalue of smaller absolute value.
!
!  Arguments
!  =========
!
!  A       (input) DOUBLE PRECISION
!          The (1,1) element of the 2-by-2 matrix.
!
!  B       (input) DOUBLE PRECISION
!          The (1,2) and (2,1) elements of the 2-by-2 matrix.
!
!  C       (input) DOUBLE PRECISION
!          The (2,2) element of the 2-by-2 matrix.
!
!  RT1     (output) DOUBLE PRECISION
!          The eigenvalue of larger absolute value.
!
!  RT2     (output) DOUBLE PRECISION
!          The eigenvalue of smaller absolute value.
!
!  Further Details
!  ===============
!
!  RT1 is accurate to a few ulps barring over/underflow.
!
!  RT2 may be inaccurate if there is massive cancellation in the
!  determinant A*C-B*B; higher precision or correctly rounded or
!  correctly truncated arithmetic would be needed to compute RT2
!  accurately in all cases.
!
!  Overflow is possible only if RT1 is within a factor of 5 of overflow.
!  Underflow is harmless if the input data is 0 or exceeds
!     underflow_threshold / macheps.
!
! =====================================================================
!
!     .. Parameters ..
  real(dp)           ONE
  PARAMETER          ( ONE = 1.0D0 )
  real(dp)           TWO
  PARAMETER          ( TWO = 2.0D0 )
  real(dp)           ZERO
  PARAMETER          ( ZERO = 0.0D0 )
  real(dp)           HALF
  PARAMETER          ( HALF = 0.5D0 )
!     ..
!     .. Local Scalars ..
  real(dp)           AB, ACMN, ACMX, ADF, DF, RT, SM, TB
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          ABS, SQRT
!     ..
!     .. Executable Statements ..
!
!     Compute the eigenvalues
!
  SM = A + C
  DF = A - C
  ADF = ABS( DF )
  TB = B + B
  AB = ABS( TB )
  IF( ABS( A ).GT.ABS( C ) ) THEN
     ACMX = A
     ACMN = C
  ELSE
     ACMX = C
     ACMN = A
  END IF
  IF( ADF.GT.AB ) THEN
     RT = ADF*SQRT( ONE+( AB / ADF )**2 )
  ELSE IF( ADF.LT.AB ) THEN
     RT = AB*SQRT( ONE+( ADF / AB )**2 )
  ELSE
!
!        Includes case AB=ADF=0
!
     RT = AB*SQRT( TWO )
  END IF
  IF( SM.LT.ZERO ) THEN
     RT1 = HALF*( SM-RT )
!
!        Order of execution important.
!        To get fully accurate smaller eigenvalue,
!        next line needs to be executed in higher precision.
!
     RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
  ELSE IF( SM.GT.ZERO ) THEN
     RT1 = HALF*( SM+RT )
!
!        Order of execution important.
!        To get fully accurate smaller eigenvalue,
!        next line needs to be executed in higher precision.
!
     RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
  ELSE
!
!        Includes case RT1 = RT2 = 0
!
     RT1 = HALF*RT
     RT2 = -HALF*RT
  END IF
  RETURN
!
!     End of DLAE2
!
  END
  SUBROUTINE DLAEV2( A, B, C, RT1, RT2, CS1, SN1 )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
  use datatypes
  real(dp)           A, B, C, CS1, RT1, RT2, SN1
!     ..
!
!  Purpose
!  =======
!
!  DLAEV2 computes the eigendecomposition of a 2-by-2 symmetric matrix
!     [  A   B  ]
!     [  B   C  ].
!  On return, RT1 is the eigenvalue of larger absolute value, RT2 is the
!  eigenvalue of smaller absolute value, and (CS1,SN1) is the unit right
!  eigenvector for RT1, giving the decomposition
!
!     [ CS1  SN1 ] [  A   B  ] [ CS1 -SN1 ]  =  [ RT1  0  ]
!     [-SN1  CS1 ] [  B   C  ] [ SN1  CS1 ]     [  0  RT2 ].
!
!  Arguments
!  =========
!
!  A       (input) DOUBLE PRECISION
!          The (1,1) element of the 2-by-2 matrix.
!
!  B       (input) DOUBLE PRECISION
!          The (1,2) element and the conjugate of the (2,1) element of
!          the 2-by-2 matrix.
!
!  C       (input) DOUBLE PRECISION
!          The (2,2) element of the 2-by-2 matrix.
!
!  RT1     (output) DOUBLE PRECISION
!          The eigenvalue of larger absolute value.
!
!  RT2     (output) DOUBLE PRECISION
!          The eigenvalue of smaller absolute value.
!
!  CS1     (output) DOUBLE PRECISION
!  SN1     (output) DOUBLE PRECISION
!          The vector (CS1, SN1) is a unit right eigenvector for RT1.
!
!  Further Details
!  ===============
!
!  RT1 is accurate to a few ulps barring over/underflow.
!
!  RT2 may be inaccurate if there is massive cancellation in the
!  determinant A*C-B*B; higher precision or correctly rounded or
!  correctly truncated arithmetic would be needed to compute RT2
!  accurately in all cases.
!
!  CS1 and SN1 are accurate to a few ulps barring over/underflow.
!
!  Overflow is possible only if RT1 is within a factor of 5 of overflow.
!  Underflow is harmless if the input data is 0 or exceeds
!     underflow_threshold / macheps.
!
! =====================================================================
!
!     .. Parameters ..
  real(dp)           ONE
  PARAMETER          ( ONE = 1.0D0 )
  real(dp)           TWO
  PARAMETER          ( TWO = 2.0D0 )
  real(dp)           ZERO
  PARAMETER          ( ZERO = 0.0D0 )
  real(dp)           HALF
  PARAMETER          ( HALF = 0.5D0 )
!     ..
!     .. Local Scalars ..
  integer(i4)        SGN1, SGN2
  real(dp)           AB, ACMN, ACMX, ACS, ADF, CS, CT, DF, RT, SM, &
                     TB, TN
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          ABS, SQRT
!     ..
!     .. Executable Statements ..
!
!     Compute the eigenvalues
!
  SM = A + C
  DF = A - C
  ADF = ABS( DF )
  TB = B + B
  AB = ABS( TB )
  IF( ABS( A ).GT.ABS( C ) ) THEN
     ACMX = A
     ACMN = C
  ELSE
     ACMX = C
     ACMN = A
  END IF
  IF( ADF.GT.AB ) THEN
     RT = ADF*SQRT( ONE+( AB / ADF )**2 )
  ELSE IF( ADF.LT.AB ) THEN
     RT = AB*SQRT( ONE+( ADF / AB )**2 )
  ELSE
!
!        Includes case AB=ADF=0
!
     RT = AB*SQRT( TWO )
  END IF
  IF( SM.LT.ZERO ) THEN
     RT1 = HALF*( SM-RT )
     SGN1 = -1
!
!        Order of execution important.
!        To get fully accurate smaller eigenvalue,
!        next line needs to be executed in higher precision.
!
     RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
  ELSE IF( SM.GT.ZERO ) THEN
     RT1 = HALF*( SM+RT )
     SGN1 = 1
!
!        Order of execution important.
!        To get fully accurate smaller eigenvalue,
!        next line needs to be executed in higher precision.
!
     RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
  ELSE
!
!        Includes case RT1 = RT2 = 0
!
     RT1 = HALF*RT
     RT2 = -HALF*RT
     SGN1 = 1
  END IF
!
!     Compute the eigenvector
!
  IF( DF.GE.ZERO ) THEN
     CS = DF + RT
     SGN2 = 1
  ELSE
     CS = DF - RT
     SGN2 = -1
  END IF
  ACS = ABS( CS )
  IF( ACS.GT.AB ) THEN
     CT = -TB / CS
     SN1 = ONE / SQRT( ONE+CT*CT )
     CS1 = CT*SN1
  ELSE
     IF( AB.EQ.ZERO ) THEN
        CS1 = ONE
        SN1 = ZERO
     ELSE
        TN = -CS / TB
        CS1 = ONE / SQRT( ONE+TN*TN )
        SN1 = TN*CS1
     END IF
  END IF
  IF( SGN1.EQ.SGN2 ) THEN
     TN = CS1
     CS1 = -SN1
     SN1 = TN
  END IF
  RETURN
!
!     End of DLAEV2
!
  END
  FUNCTION DLANST( NORM, N, D, E )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
  use datatypes
  CHARACTER          NORM
  integer(i4)        N
!     ..
!     .. Array Arguments ..
  real(dp)           D( * ), E( * )
!     ..
!
!  Purpose
!  =======
!
!  DLANST  returns the value of the one norm,  or the Frobenius norm, or
!  the  infinity norm,  or the  element of  largest absolute value  of a
!  real symmetric tridiagonal matrix A.
!
!  Description
!  ===========
!
!  DLANST returns the value
!
!     DLANST = ( max(abs(A(i,j))), NORM = 'M' or 'm'
!              (
!              ( norm1(A),         NORM = '1', 'O' or 'o'
!              (
!              ( normI(A),         NORM = 'I' or 'i'
!              (
!              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
!
!  where  norm1  denotes the  one norm of a matrix (maximum column sum),
!  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
!  normF  denotes the  Frobenius norm of a matrix (square root of sum of
!  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.
!
!  Arguments
!  =========
!
!  NORM    (input) CHARACTER*1
!          Specifies the value to be returned in DLANST as described
!          above.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.  When N = 0, DLANST is
!          set to zero.
!
!  D       (input) DOUBLE PRECISION array, dimension (N)
!          The diagonal elements of A.
!
!  E       (input) DOUBLE PRECISION array, dimension (N-1)
!          The (n-1) sub-diagonal or super-diagonal elements of A.
!
!  =====================================================================
!
!     .. Parameters ..
  real(dp)           ONE, ZERO
  PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
  integer(i4)        I
  real(dp)           ANORM, SCALE, SUM
!     ..
!     .. External Functions ..
  LOGICAL            LSAME
  EXTERNAL           LSAME
  real(dp)           DLANST
!     ..
!     .. External Subroutines ..
  EXTERNAL           DLASSQ
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          ABS, MAX, SQRT
!     ..
!     .. Executable Statements ..
!
  IF( N.LE.0 ) THEN
     ANORM = ZERO
  ELSE IF( LSAME( NORM, 'M' ) ) THEN
!
!        Find max(abs(A(i,j))).
!
     ANORM = ABS( D( N ) )
     DO 10 I = 1, N - 1
        ANORM = MAX( ANORM, ABS( D( I ) ) )
        ANORM = MAX( ANORM, ABS( E( I ) ) )
10    CONTINUE
  ELSE IF( LSAME( NORM, 'O' ) .OR. NORM.EQ.'1' .OR. &
           LSAME( NORM, 'I' ) ) THEN
!
!        Find norm1(A).
!
     IF( N.EQ.1 ) THEN
        ANORM = ABS( D( 1 ) )
     ELSE
        ANORM = MAX( ABS( D( 1 ) )+ABS( E( 1 ) ), &
                ABS( E( N-1 ) )+ABS( D( N ) ) )
        DO 20 I = 2, N - 1
           ANORM = MAX( ANORM, ABS( D( I ) )+ABS( E( I ) )+ &
                   ABS( E( I-1 ) ) )
20       CONTINUE
     END IF
  ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
!
!        Find normF(A).
!
     SCALE = ZERO
     SUM = ONE
     IF( N.GT.1 ) THEN
        CALL DLASSQ( N-1_i4, E, 1_i4, SCALE, SUM )
        SUM = 2*SUM
     END IF
     CALL DLASSQ( N, D, 1_i4, SCALE, SUM )
     ANORM = SCALE*SQRT( SUM )
  END IF
!
  DLANST = ANORM
  RETURN
!
!     End of DLANST
!
  END
  FUNCTION DLAPY2( X, Y )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
  use datatypes
  real(dp)           X, Y
!     ..
!
!  Purpose
!  =======
!
!  DLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
!  overflow.
!
!  Arguments
!  =========
!
!  X       (input) DOUBLE PRECISION
!  Y       (input) DOUBLE PRECISION
!          X and Y specify the values x and y.
!
!  =====================================================================
!
!     .. Parameters ..
  real(dp)           ZERO
  PARAMETER          ( ZERO = 0.0D0 )
  real(dp)           ONE
  PARAMETER          ( ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
  real(dp)           W, XABS, YABS, Z
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          ABS, MAX, MIN, SQRT
  real(dp)           DLAPY2
!     ..
!     .. Executable Statements ..
!
  XABS = ABS( X )
  YABS = ABS( Y )
  W = MAX( XABS, YABS )
  Z = MIN( XABS, YABS )
  IF( Z.EQ.ZERO ) THEN
     DLAPY2 = W
  ELSE
     DLAPY2 = W*SQRT( ONE+( Z / W )**2 )
  END IF
  RETURN
!
!     End of DLAPY2
!
  END
  SUBROUTINE DLARTG( F, G, CS, SN, R )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
  use datatypes
  real(dp)           CS, F, G, R, SN
!     ..
!
!  Purpose
!  =======
!
!  DLARTG generate a plane rotation so that
!
!     [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.
!     [ -SN  CS  ]     [ G ]     [ 0 ]
!
!  This is a slower, more accurate version of the BLAS1 routine DROTG,
!  with the following other differences:
!     F and G are unchanged on return.
!     If G=0, then CS=1 and SN=0.
!     If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any
!        floating point operations (saves work in DBDSQR when
!        there are zeros on the diagonal).
!
!  If F exceeds G in magnitude, CS will be positive.
!
!  Arguments
!  =========
!
!  F       (input) DOUBLE PRECISION
!          The first component of vector to be rotated.
!
!  G       (input) DOUBLE PRECISION
!          The second component of vector to be rotated.
!
!  CS      (output) DOUBLE PRECISION
!          The cosine of the rotation.
!
!  SN      (output) DOUBLE PRECISION
!          The sine of the rotation.
!
!  R       (output) DOUBLE PRECISION
!          The nonzero component of the rotated vector.
!
!  =====================================================================
!
!     .. Parameters ..
  real(dp)           ZERO
  PARAMETER          ( ZERO = 0.0D0 )
  real(dp)           ONE
  PARAMETER          ( ONE = 1.0D0 )
  real(dp)           TWO
  PARAMETER          ( TWO = 2.0D0 )
!     ..
!     .. Local Scalars ..
  LOGICAL            FIRST
  integer(i4)        COUNT, I
  real(dp)           EPS, F1, G1, SAFMIN, SAFMN2, SAFMX2, SCALE
!     ..
!     .. External Functions ..
  real(dp)           DLAMCH
  EXTERNAL           DLAMCH
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          ABS, INT, LOG, MAX, SQRT
!     ..
!     .. Save statement ..
  SAVE               FIRST, SAFMX2, SAFMIN, SAFMN2
!     ..
!     .. Data statements ..
  DATA               FIRST / .TRUE. /
!     ..
!     .. Executable Statements ..
!
  IF( FIRST ) THEN
     FIRST = .FALSE.
     SAFMIN = DLAMCH( 'S' )
     EPS = DLAMCH( 'E' )
     SAFMN2 = DLAMCH( 'B' )**INT( LOG( SAFMIN / EPS ) / &
              LOG( DLAMCH( 'B' ) ) / TWO )
     SAFMX2 = ONE / SAFMN2
  END IF
  IF( G.EQ.ZERO ) THEN
     CS = ONE
     SN = ZERO
     R = F
  ELSE IF( F.EQ.ZERO ) THEN
     CS = ZERO
     SN = ONE
     R = G
  ELSE
     F1 = F
     G1 = G
     SCALE = MAX( ABS( F1 ), ABS( G1 ) )
     IF( SCALE.GE.SAFMX2 ) THEN
        COUNT = 0
10       CONTINUE
        COUNT = COUNT + 1
        F1 = F1*SAFMN2
        G1 = G1*SAFMN2
        SCALE = MAX( ABS( F1 ), ABS( G1 ) )
        IF( SCALE.GE.SAFMX2 ) &
           GO TO 10
        R = SQRT( F1**2+G1**2 )
        CS = F1 / R
        SN = G1 / R
        DO 20 I = 1, COUNT
           R = R*SAFMX2
20       CONTINUE
     ELSE IF( SCALE.LE.SAFMN2 ) THEN
        COUNT = 0
30       CONTINUE
        COUNT = COUNT + 1
        F1 = F1*SAFMX2
        G1 = G1*SAFMX2
        SCALE = MAX( ABS( F1 ), ABS( G1 ) )
        IF( SCALE.LE.SAFMN2 ) &
           GO TO 30
        R = SQRT( F1**2+G1**2 )
        CS = F1 / R
        SN = G1 / R
        DO 40 I = 1, COUNT
           R = R*SAFMN2
40       CONTINUE
     ELSE
        R = SQRT( F1**2+G1**2 )
        CS = F1 / R
        SN = G1 / R
     END IF
     IF( ABS( F ).GT.ABS( G ) .AND. CS.LT.ZERO ) THEN
        CS = -CS
        SN = -SN
        R = -R
     END IF
  END IF
  RETURN
!
!     End of DLARTG
!
  END
  SUBROUTINE DLARFG( N, ALPHA, X, INCX, TAU )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
  use datatypes
  integer(i4)        INCX, N
  real(dp)           ALPHA, TAU
!     ..
!     .. Array Arguments ..
  real(dp)           X( * )
!     ..
!
!  Purpose
!  =======
!
!  DLARFG generates a real elementary reflector H of order n, such
!  that
!
!        H * ( alpha ) = ( beta ),   H' * H = I.
!            (   x   )   (   0  )
!
!  where alpha and beta are scalars, and x is an (n-1)-element real
!  vector. H is represented in the form
!
!        H = I - tau * ( 1 ) * ( 1 v' ) ,
!                      ( v )
!
!  where tau is a real scalar and v is a real (n-1)-element
!  vector.
!
!  If the elements of x are all zero, then tau = 0 and H is taken to be
!  the unit matrix.
!
!  Otherwise  1 <= tau <= 2.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the elementary reflector.
!
!  ALPHA   (input/output) DOUBLE PRECISION
!          On entry, the value alpha.
!          On exit, it is overwritten with the value beta.
!
!  X       (input/output) DOUBLE PRECISION array, dimension
!                         (1+(N-2)*abs(INCX))
!          On entry, the vector x.
!          On exit, it is overwritten with the vector v.
!
!  INCX    (input) INTEGER
!          The increment between elements of X. INCX > 0.
!
!  TAU     (output) DOUBLE PRECISION
!          The value tau.
!
!  =====================================================================
!
!     .. Parameters ..
  real(dp)           ONE, ZERO
  PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
  integer(i4)        J, KNT
  real(dp)           BETA, RSAFMN, SAFMIN, XNORM
!     ..
!     .. External Functions ..
  real(dp)           DLAMCH, DLAPY2, DNRM2
  EXTERNAL           DLAMCH, DLAPY2, DNRM2
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          ABS, SIGN
!     ..
!     .. External Subroutines ..
  EXTERNAL           DSCAL
!     ..
!     .. Executable Statements ..
!
  IF( N.LE.1 ) THEN
     TAU = ZERO
     RETURN
  END IF
!
  XNORM = DNRM2( N-1_i4, X, INCX )
!
  IF( XNORM.EQ.ZERO ) THEN
!
!        H  =  I
!
     TAU = ZERO
  ELSE
!
!        general case
!
     BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
     SAFMIN = DLAMCH( 'S' ) / DLAMCH( 'E' )
     IF( ABS( BETA ).LT.SAFMIN ) THEN
!
!           XNORM, BETA may be inaccurate; scale X and recompute them
!
        RSAFMN = ONE / SAFMIN
        KNT = 0
10       CONTINUE
        KNT = KNT + 1
        CALL DSCAL( N-1_i4, RSAFMN, X, INCX )
        BETA = BETA*RSAFMN
        ALPHA = ALPHA*RSAFMN
        IF( ABS( BETA ).LT.SAFMIN ) &
           GO TO 10
!
!           New BETA is at most 1, at least SAFMIN
!
        XNORM = DNRM2( N-1_i4, X, INCX )
        BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
        TAU = ( BETA-ALPHA ) / BETA
        CALL DSCAL( N-1_i4, ONE / ( ALPHA-BETA ), X, INCX )
!
!           If ALPHA is subnormal, it may lose relative accuracy
!
        ALPHA = BETA
        DO 20 J = 1, KNT
           ALPHA = ALPHA*SAFMIN
20       CONTINUE
     ELSE
        TAU = ( BETA-ALPHA ) / BETA
        CALL DSCAL( N-1_i4, ONE / ( ALPHA-BETA ), X, INCX )
        ALPHA = BETA
     END IF
  END IF
!
  RETURN
!
!     End of DLARFG
!
  END
  SUBROUTINE DLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
  use datatypes
  CHARACTER          TYPE
  integer(i4)        INFO, KL, KU, LDA, M, N
  real(dp)           CFROM, CTO
!     ..
!     .. Array Arguments ..
  real(dp)           A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DLASCL multiplies the M by N real matrix A by the real scalar
!  CTO/CFROM.  This is done without over/underflow as long as the final
!  result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
!  A may be full, upper triangular, lower triangular, upper Hessenberg,
!  or banded.
!
!  Arguments
!  =========
!
!  TYPE    (input) CHARACTER*1
!          TYPE indices the storage type of the input matrix.
!          = 'G':  A is a full matrix.
!          = 'L':  A is a lower triangular matrix.
!          = 'U':  A is an upper triangular matrix.
!          = 'H':  A is an upper Hessenberg matrix.
!          = 'B':  A is a symmetric band matrix with lower bandwidth KL
!                  and upper bandwidth KU and with the only the lower
!                  half stored.
!          = 'Q':  A is a symmetric band matrix with lower bandwidth KL
!                  and upper bandwidth KU and with the only the upper
!                  half stored.
!          = 'Z':  A is a band matrix with lower bandwidth KL and upper
!                  bandwidth KU.
!
!  KL      (input) INTEGER
!          The lower bandwidth of A.  Referenced only if TYPE = 'B',
!          'Q' or 'Z'.
!
!  KU      (input) INTEGER
!          The upper bandwidth of A.  Referenced only if TYPE = 'B',
!          'Q' or 'Z'.
!
!  CFROM   (input) DOUBLE PRECISION
!  CTO     (input) DOUBLE PRECISION
!          The matrix A is multiplied by CTO/CFROM. A(I,J) is computed
!          without over/underflow if the final result CTO*A(I,J)/CFROM
!          can be represented without over/underflow.  CFROM must be
!          nonzero.
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,M)
!          The matrix to be multiplied by CTO/CFROM.  See TYPE for the
!          storage type.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  INFO    (output) INTEGER
!          0  - successful exit
!          <0 - if INFO = -i, the i-th argument had an illegal value.
!
!  =====================================================================
!
!     .. Parameters ..
  real(dp)           ZERO, ONE
  PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
  LOGICAL            DONE
  integer(i4)        I, ITYPE, J, K1, K2, K3, K4
  real(dp)           BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM
!     ..
!     .. External Functions ..
  LOGICAL            LSAME
  real(dp)           DLAMCH
  EXTERNAL           LSAME, DLAMCH
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          ABS, MAX, MIN
!     ..
!     .. External Subroutines ..
  EXTERNAL           XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
  INFO = 0
!
  IF( LSAME( TYPE, 'G' ) ) THEN
     ITYPE = 0
  ELSE IF( LSAME( TYPE, 'L' ) ) THEN
     ITYPE = 1
  ELSE IF( LSAME( TYPE, 'U' ) ) THEN
     ITYPE = 2
  ELSE IF( LSAME( TYPE, 'H' ) ) THEN
     ITYPE = 3
  ELSE IF( LSAME( TYPE, 'B' ) ) THEN
     ITYPE = 4
  ELSE IF( LSAME( TYPE, 'Q' ) ) THEN
     ITYPE = 5
  ELSE IF( LSAME( TYPE, 'Z' ) ) THEN
     ITYPE = 6
  ELSE
     ITYPE = -1
  END IF
!
  IF( ITYPE.EQ.-1 ) THEN
     INFO = -1
  ELSE IF( CFROM.EQ.ZERO ) THEN
     INFO = -4
  ELSE IF( M.LT.0 ) THEN
     INFO = -6
  ELSE IF( N.LT.0 .OR. ( ITYPE.EQ.4 .AND. N.NE.M ) .OR. &
           ( ITYPE.EQ.5 .AND. N.NE.M ) ) THEN
     INFO = -7
  ELSE IF( ITYPE.LE.3 .AND. LDA.LT.MAX( 1, M ) ) THEN
     INFO = -9
  ELSE IF( ITYPE.GE.4 ) THEN
     IF( KL.LT.0 .OR. KL.GT.MAX( M-1, 0 ) ) THEN
        INFO = -2
     ELSE IF( KU.LT.0 .OR. KU.GT.MAX( N-1, 0 ) .OR. &
              ( ( ITYPE.EQ.4 .OR. ITYPE.EQ.5 ) .AND. KL.NE.KU ) ) &
               THEN
        INFO = -3
     ELSE IF( ( ITYPE.EQ.4 .AND. LDA.LT.KL+1 ) .OR. &
              ( ITYPE.EQ.5 .AND. LDA.LT.KU+1 ) .OR. &
              ( ITYPE.EQ.6 .AND. LDA.LT.2*KL+KU+1 ) ) THEN
        INFO = -9
     END IF
  END IF
!
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DLASCL', -INFO )
     RETURN
  END IF
!
!     Quick return if possible
!
  IF( N.EQ.0 .OR. M.EQ.0 ) &
     RETURN
!
!     Get machine parameters
!
  SMLNUM = DLAMCH( 'S' )
  BIGNUM = ONE / SMLNUM
!
  CFROMC = CFROM
  CTOC = CTO
!
10 CONTINUE
  CFROM1 = CFROMC*SMLNUM
  CTO1 = CTOC / BIGNUM
  IF( ABS( CFROM1 ).GT.ABS( CTOC ) .AND. CTOC.NE.ZERO ) THEN
     MUL = SMLNUM
     DONE = .FALSE.
     CFROMC = CFROM1
  ELSE IF( ABS( CTO1 ).GT.ABS( CFROMC ) ) THEN
     MUL = BIGNUM
     DONE = .FALSE.
     CTOC = CTO1
  ELSE
     MUL = CTOC / CFROMC
     DONE = .TRUE.
  END IF
!
  IF( ITYPE.EQ.0 ) THEN
!
!        Full matrix
!
     DO 30 J = 1, N
        DO 20 I = 1, M
           A( I, J ) = A( I, J )*MUL
20       CONTINUE
30    CONTINUE
!
  ELSE IF( ITYPE.EQ.1 ) THEN
!
!        Lower triangular matrix
!
     DO 50 J = 1, N
        DO 40 I = J, M
           A( I, J ) = A( I, J )*MUL
40       CONTINUE
50    CONTINUE
!
  ELSE IF( ITYPE.EQ.2 ) THEN
!
!        Upper triangular matrix
!
     DO 70 J = 1, N
        DO 60 I = 1, MIN( J, M )
           A( I, J ) = A( I, J )*MUL
60       CONTINUE
70    CONTINUE
!
  ELSE IF( ITYPE.EQ.3 ) THEN
!
!        Upper Hessenberg matrix
!
     DO 90 J = 1, N
        DO 80 I = 1, MIN( J+1, M )
           A( I, J ) = A( I, J )*MUL
80       CONTINUE
90    CONTINUE
!
  ELSE IF( ITYPE.EQ.4 ) THEN
!
!        Lower half of a symmetric band matrix
!
     K3 = KL + 1
     K4 = N + 1
     DO 110 J = 1, N
        DO 100 I = 1, MIN( K3, K4-J )
           A( I, J ) = A( I, J )*MUL
100       CONTINUE
110    CONTINUE
!
  ELSE IF( ITYPE.EQ.5 ) THEN
!
!        Upper half of a symmetric band matrix
!
     K1 = KU + 2
     K3 = KU + 1
     DO 130 J = 1, N
        DO 120 I = MAX( K1-J, 1 ), K3
           A( I, J ) = A( I, J )*MUL
120       CONTINUE
130    CONTINUE
!
  ELSE IF( ITYPE.EQ.6 ) THEN
!
!        Band matrix
!
     K1 = KL + KU + 2
     K2 = KL + 1
     K3 = 2*KL + KU + 1
     K4 = KL + KU + 1 + M
     DO 150 J = 1, N
        DO 140 I = MAX( K1-J, K2 ), MIN( K3, K4-J )
           A( I, J ) = A( I, J )*MUL
140       CONTINUE
150    CONTINUE
!
  END IF
!
  IF( .NOT.DONE ) &
     GO TO 10
!
  RETURN
!
!     End of DLASCL
!
  END
  SUBROUTINE DLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
  use datatypes
  CHARACTER          UPLO
  integer(i4)        LDA, M, N
  real(dp)           ALPHA, BETA
!     ..
!     .. Array Arguments ..
  real(dp)           A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DLASET initializes an m-by-n matrix A to BETA on the diagonal and
!  ALPHA on the offdiagonals.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies the part of the matrix A to be set.
!          = 'U':      Upper triangular part is set; the strictly lower
!                      triangular part of A is not changed.
!          = 'L':      Lower triangular part is set; the strictly upper
!                      triangular part of A is not changed.
!          Otherwise:  All of the matrix A is set.
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  ALPHA   (input) DOUBLE PRECISION
!          The constant to which the offdiagonal elements are to be set.
!
!  BETA    (input) DOUBLE PRECISION
!          The constant to which the diagonal elements are to be set.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On exit, the leading m-by-n submatrix of A is set as follows:
!
!          if UPLO = 'U', A(i,j) = ALPHA, 1<=i<=j-1, 1<=j<=n,
!          if UPLO = 'L', A(i,j) = ALPHA, j+1<=i<=m, 1<=j<=n,
!          otherwise,     A(i,j) = ALPHA, 1<=i<=m, 1<=j<=n, i.ne.j,
!
!          and, for all UPLO, A(i,i) = BETA, 1<=i<=min(m,n).
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
! =====================================================================
!
!     .. Local Scalars ..
  integer(i4)        I, J
!     ..
!     .. External Functions ..
  LOGICAL            LSAME
  EXTERNAL           LSAME
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          MIN
!     ..
!     .. Executable Statements ..
!
  IF( LSAME( UPLO, 'U' ) ) THEN
!
!        Set the strictly upper triangular or trapezoidal part of the
!        array to ALPHA.
!
     DO 20 J = 2, N
        DO 10 I = 1, MIN( J-1, M )
           A( I, J ) = ALPHA
10       CONTINUE
20    CONTINUE
!
  ELSE IF( LSAME( UPLO, 'L' ) ) THEN
!
!        Set the strictly lower triangular or trapezoidal part of the
!        array to ALPHA.
!
     DO 40 J = 1, MIN( M, N )
        DO 30 I = J + 1, M
           A( I, J ) = ALPHA
30       CONTINUE
40    CONTINUE
!
  ELSE
!
!        Set the leading m-by-n submatrix to ALPHA.
!
     DO 60 J = 1, N
        DO 50 I = 1, M
           A( I, J ) = ALPHA
50       CONTINUE
60    CONTINUE
  END IF
!
!     Set the first min(M,N) diagonal elements to BETA.
!
  DO 70 I = 1, MIN( M, N )
     A( I, I ) = BETA
70 CONTINUE
!
  RETURN
!
!     End of DLASET
!
  END
  SUBROUTINE DLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
  use datatypes
  CHARACTER          DIRECT, PIVOT, SIDE
  integer(i4)        LDA, M, N
!     ..
!     .. Array Arguments ..
  real(dp)           A( LDA, * ), C( * ), S( * )
!     ..
!
!  Purpose
!  =======
!
!  DLASR   performs the transformation
!
!     A := P*A,   when SIDE = 'L' or 'l'  (  Left-hand side )
!
!     A := A*P',  when SIDE = 'R' or 'r'  ( Right-hand side )
!
!  where A is an m by n real matrix and P is an orthogonal matrix,
!  consisting of a sequence of plane rotations determined by the
!  parameters PIVOT and DIRECT as follows ( z = m when SIDE = 'L' or 'l'
!  and z = n when SIDE = 'R' or 'r' ):
!
!  When  DIRECT = 'F' or 'f'  ( Forward sequence ) then
!
!     P = P( z - 1 )*...*P( 2 )*P( 1 ),
!
!  and when DIRECT = 'B' or 'b'  ( Backward sequence ) then
!
!     P = P( 1 )*P( 2 )*...*P( z - 1 ),
!
!  where  P( k ) is a plane rotation matrix for the following planes:
!
!     when  PIVOT = 'V' or 'v'  ( Variable pivot ),
!        the plane ( k, k + 1 )
!
!     when  PIVOT = 'T' or 't'  ( Top pivot ),
!        the plane ( 1, k + 1 )
!
!     when  PIVOT = 'B' or 'b'  ( Bottom pivot ),
!        the plane ( k, z )
!
!  c( k ) and s( k )  must contain the  cosine and sine that define the
!  matrix  P( k ).  The two by two plane rotation part of the matrix
!  P( k ), R( k ), is assumed to be of the form
!
!     R( k ) = (  c( k )  s( k ) ).
!              ( -s( k )  c( k ) )
!
!  This version vectorises across rows of the array A when SIDE = 'L'.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          Specifies whether the plane rotation matrix P is applied to
!          A on the left or the right.
!          = 'L':  Left, compute A := P*A
!          = 'R':  Right, compute A:= A*P'
!
!  DIRECT  (input) CHARACTER*1
!          Specifies whether P is a forward or backward sequence of
!          plane rotations.
!          = 'F':  Forward, P = P( z - 1 )*...*P( 2 )*P( 1 )
!          = 'B':  Backward, P = P( 1 )*P( 2 )*...*P( z - 1 )
!
!  PIVOT   (input) CHARACTER*1
!          Specifies the plane for which P(k) is a plane rotation
!          matrix.
!          = 'V':  Variable pivot, the plane (k,k+1)
!          = 'T':  Top pivot, the plane (1,k+1)
!          = 'B':  Bottom pivot, the plane (k,z)
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  If m <= 1, an immediate
!          return is effected.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  If n <= 1, an
!          immediate return is effected.
!
!  C, S    (input) DOUBLE PRECISION arrays, dimension
!                  (M-1) if SIDE = 'L'
!                  (N-1) if SIDE = 'R'
!          c(k) and s(k) contain the cosine and sine that define the
!          matrix P(k).  The two by two plane rotation part of the
!          matrix P(k), R(k), is assumed to be of the form
!          R( k ) = (  c( k )  s( k ) ).
!                   ( -s( k )  c( k ) )
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          The m by n matrix A.  On exit, A is overwritten by P*A if
!          SIDE = 'R' or by A*P' if SIDE = 'L'.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  =====================================================================
!
!     .. Parameters ..
  real(dp)           ONE, ZERO
  PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
  integer(i4)        I, INFO, J
  real(dp)           CTEMP, STEMP, TEMP
!     ..
!     .. External Functions ..
  LOGICAL            LSAME
  EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
  EXTERNAL           XERBLA
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
  INFO = 0
  IF( .NOT.( LSAME( SIDE, 'L' ) .OR. LSAME( SIDE, 'R' ) ) ) THEN
     INFO = 1
  ELSE IF( .NOT.( LSAME( PIVOT, 'V' ) .OR. LSAME( PIVOT, &
           'T' ) .OR. LSAME( PIVOT, 'B' ) ) ) THEN
     INFO = 2
  ELSE IF( .NOT.( LSAME( DIRECT, 'F' ) .OR. LSAME( DIRECT, 'B' ) ) ) &
            THEN
     INFO = 3
  ELSE IF( M.LT.0 ) THEN
     INFO = 4
  ELSE IF( N.LT.0 ) THEN
     INFO = 5
  ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
     INFO = 9
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DLASR ', INFO )
     RETURN
  END IF
!
!     Quick return if possible
!
  IF( ( M.EQ.0 ) .OR. ( N.EQ.0 ) ) &
     RETURN
  IF( LSAME( SIDE, 'L' ) ) THEN
!
!        Form  P * A
!
     IF( LSAME( PIVOT, 'V' ) ) THEN
        IF( LSAME( DIRECT, 'F' ) ) THEN
           DO 20 J = 1, M - 1
              CTEMP = C( J )
              STEMP = S( J )
              IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                 DO 10 I = 1, N
                    TEMP = A( J+1, I )
                    A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                    A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
10                CONTINUE
              END IF
20          CONTINUE
        ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
           DO 40 J = M - 1, 1, -1
              CTEMP = C( J )
              STEMP = S( J )
              IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                 DO 30 I = 1, N
                    TEMP = A( J+1, I )
                    A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                    A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
30                CONTINUE
              END IF
40          CONTINUE
        END IF
     ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
        IF( LSAME( DIRECT, 'F' ) ) THEN
           DO 60 J = 2, M
              CTEMP = C( J-1 )
              STEMP = S( J-1 )
              IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                 DO 50 I = 1, N
                    TEMP = A( J, I )
                    A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                    A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
50                CONTINUE
              END IF
60          CONTINUE
        ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
           DO 80 J = M, 2, -1
              CTEMP = C( J-1 )
              STEMP = S( J-1 )
              IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                 DO 70 I = 1, N
                    TEMP = A( J, I )
                    A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                    A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
70                CONTINUE
              END IF
80          CONTINUE
        END IF
     ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
        IF( LSAME( DIRECT, 'F' ) ) THEN
           DO 100 J = 1, M - 1
              CTEMP = C( J )
              STEMP = S( J )
              IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                 DO 90 I = 1, N
                    TEMP = A( J, I )
                    A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                    A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
90                CONTINUE
              END IF
100          CONTINUE
        ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
           DO 120 J = M - 1, 1, -1
              CTEMP = C( J )
              STEMP = S( J )
              IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                 DO 110 I = 1, N
                    TEMP = A( J, I )
                    A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                    A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
110                CONTINUE
              END IF
120          CONTINUE
        END IF
     END IF
  ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!        Form A * P'
!
     IF( LSAME( PIVOT, 'V' ) ) THEN
        IF( LSAME( DIRECT, 'F' ) ) THEN
           DO 140 J = 1, N - 1
              CTEMP = C( J )
              STEMP = S( J )
              IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                 DO 130 I = 1, M
                    TEMP = A( I, J+1 )
                    A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                    A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
130                CONTINUE
              END IF
140          CONTINUE
        ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
           DO 160 J = N - 1, 1, -1
              CTEMP = C( J )
              STEMP = S( J )
              IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                 DO 150 I = 1, M
                    TEMP = A( I, J+1 )
                    A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                    A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
150                CONTINUE
              END IF
160          CONTINUE
        END IF
     ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
        IF( LSAME( DIRECT, 'F' ) ) THEN
           DO 180 J = 2, N
              CTEMP = C( J-1 )
              STEMP = S( J-1 )
              IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                 DO 170 I = 1, M
                    TEMP = A( I, J )
                    A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                    A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
170                CONTINUE
              END IF
180          CONTINUE
        ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
           DO 200 J = N, 2, -1
              CTEMP = C( J-1 )
              STEMP = S( J-1 )
              IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                 DO 190 I = 1, M
                    TEMP = A( I, J )
                    A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                    A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
190                CONTINUE
              END IF
200          CONTINUE
        END IF
     ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
        IF( LSAME( DIRECT, 'F' ) ) THEN
           DO 220 J = 1, N - 1
              CTEMP = C( J )
              STEMP = S( J )
              IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                 DO 210 I = 1, M
                    TEMP = A( I, J )
                    A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                    A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
210                CONTINUE
              END IF
220          CONTINUE
        ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
           DO 240 J = N - 1, 1, -1
              CTEMP = C( J )
              STEMP = S( J )
              IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                 DO 230 I = 1, M
                    TEMP = A( I, J )
                    A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                    A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
230                CONTINUE
              END IF
240          CONTINUE
        END IF
     END IF
  END IF
!
  RETURN
!
!     End of DLASR
!
  END
  SUBROUTINE DLASRT( ID, N, D, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
  use datatypes
  CHARACTER          ID
  integer(i4)        INFO, N
!     ..
!     .. Array Arguments ..
  real(dp)           D( * )
!     ..
!
!  Purpose
!  =======
!
!  Sort the numbers in D in increasing order (if ID = 'I') or
!  in decreasing order (if ID = 'D' ).
!
!  Use Quick Sort, reverting to Insertion sort on arrays of
!  size <= 20. Dimension of STACK limits N to about 2**32.
!
!  Arguments
!  =========
!
!  ID      (input) CHARACTER*1
!          = 'I': sort D in increasing order;
!          = 'D': sort D in decreasing order.
!
!  N       (input) INTEGER
!          The length of the array D.
!
!  D       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, the array to be sorted.
!          On exit, D has been sorted into increasing order
!          (D(1) <= ... <= D(N) ) or into decreasing order
!          (D(1) >= ... >= D(N) ), depending on ID.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
  integer(i4)        SELECT
  PARAMETER          ( SELECT = 20 )
!     ..
!     .. Local Scalars ..
  integer(i4)        DIR, ENDD, I, J, START, STKPNT
  real(dp)           D1, D2, D3, DMNMX, TMP
!     ..
!     .. Local Arrays ..
  integer(i4)        STACK( 2, 32 )
!     ..
!     .. External Functions ..
  LOGICAL            LSAME
  EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
  EXTERNAL           XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input paramters.
!
  INFO = 0
  DIR = -1
  IF( LSAME( ID, 'D' ) ) THEN
     DIR = 0
  ELSE IF( LSAME( ID, 'I' ) ) THEN
     DIR = 1
  END IF
  IF( DIR.EQ.-1 ) THEN
     INFO = -1
  ELSE IF( N.LT.0 ) THEN
     INFO = -2
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DLASRT', -INFO )
     RETURN
  END IF
!
!     Quick return if possible
!
  IF( N.LE.1 ) &
     RETURN
!
  STKPNT = 1
  STACK( 1, 1 ) = 1
  STACK( 2, 1 ) = N
10 CONTINUE
  START = STACK( 1, STKPNT )
  ENDD = STACK( 2, STKPNT )
  STKPNT = STKPNT - 1
  IF( ENDD-START.LE.SELECT .AND. ENDD-START.GT.0 ) THEN
!
!        Do Insertion sort on D( START:ENDD )
!
     IF( DIR.EQ.0 ) THEN
!
!           Sort into decreasing order
!
        DO 30 I = START + 1, ENDD
           DO 20 J = I, START + 1, -1
              IF( D( J ).GT.D( J-1 ) ) THEN
                 DMNMX = D( J )
                 D( J ) = D( J-1 )
                 D( J-1 ) = DMNMX
              ELSE
                 GO TO 30
              END IF
20          CONTINUE
30       CONTINUE
!
     ELSE
!
!           Sort into increasing order
!
        DO 50 I = START + 1, ENDD
           DO 40 J = I, START + 1, -1
              IF( D( J ).LT.D( J-1 ) ) THEN
                 DMNMX = D( J )
                 D( J ) = D( J-1 )
                 D( J-1 ) = DMNMX
              ELSE
                 GO TO 50
              END IF
40          CONTINUE
50       CONTINUE
!
     END IF
!
  ELSE IF( ENDD-START.GT.SELECT ) THEN
!
!        Partition D( START:ENDD ) and stack parts, largest one first
!
!        Choose partition entry as median of 3
!
     D1 = D( START )
     D2 = D( ENDD )
     I = ( START+ENDD ) / 2
     D3 = D( I )
     IF( D1.LT.D2 ) THEN
        IF( D3.LT.D1 ) THEN
           DMNMX = D1
        ELSE IF( D3.LT.D2 ) THEN
           DMNMX = D3
        ELSE
           DMNMX = D2
        END IF
     ELSE
        IF( D3.LT.D2 ) THEN
           DMNMX = D2
        ELSE IF( D3.LT.D1 ) THEN
           DMNMX = D3
        ELSE
           DMNMX = D1
        END IF
     END IF
!
     IF( DIR.EQ.0 ) THEN
!
!           Sort into decreasing order
!
        I = START - 1
        J = ENDD + 1
60       CONTINUE
70       CONTINUE
        J = J - 1
        IF( D( J ).LT.DMNMX ) &
           GO TO 70
80       CONTINUE
        I = I + 1
        IF( D( I ).GT.DMNMX ) &
           GO TO 80
        IF( I.LT.J ) THEN
           TMP = D( I )
           D( I ) = D( J )
           D( J ) = TMP
           GO TO 60
        END IF
        IF( J-START.GT.ENDD-J-1 ) THEN
           STKPNT = STKPNT + 1
           STACK( 1, STKPNT ) = START
           STACK( 2, STKPNT ) = J
           STKPNT = STKPNT + 1
           STACK( 1, STKPNT ) = J + 1
           STACK( 2, STKPNT ) = ENDD
        ELSE
           STKPNT = STKPNT + 1
           STACK( 1, STKPNT ) = J + 1
           STACK( 2, STKPNT ) = ENDD
           STKPNT = STKPNT + 1
           STACK( 1, STKPNT ) = START
           STACK( 2, STKPNT ) = J
        END IF
     ELSE
!
!           Sort into increasing order
!
        I = START - 1
        J = ENDD + 1
90       CONTINUE
100       CONTINUE
        J = J - 1
        IF( D( J ).GT.DMNMX ) &
           GO TO 100
110       CONTINUE
        I = I + 1
        IF( D( I ).LT.DMNMX ) &
           GO TO 110
        IF( I.LT.J ) THEN
           TMP = D( I )
           D( I ) = D( J )
           D( J ) = TMP
           GO TO 90
        END IF
        IF( J-START.GT.ENDD-J-1 ) THEN
           STKPNT = STKPNT + 1
           STACK( 1, STKPNT ) = START
           STACK( 2, STKPNT ) = J
           STKPNT = STKPNT + 1
           STACK( 1, STKPNT ) = J + 1
           STACK( 2, STKPNT ) = ENDD
        ELSE
           STKPNT = STKPNT + 1
           STACK( 1, STKPNT ) = J + 1
           STACK( 2, STKPNT ) = ENDD
           STKPNT = STKPNT + 1
           STACK( 1, STKPNT ) = START
           STACK( 2, STKPNT ) = J
        END IF
     END IF
  END IF
  IF( STKPNT.GT.0 ) &
     GO TO 10
  RETURN
!
!     End of DLASRT
!
  END
  SUBROUTINE DLASSQ( N, X, INCX, SCALE, SUMSQ )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
  use datatypes
  integer(i4)        INCX, N
  real(dp)           SCALE, SUMSQ
!     ..
!     .. Array Arguments ..
  real(dp)           X( * )
!     ..
!
!  Purpose
!  =======
!
!  DLASSQ  returns the values  scl  and  smsq  such that
!
!     ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
!
!  where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is
!  assumed to be non-negative and  scl  returns the value
!
!     scl = max( scale, abs( x( i ) ) ).
!
!  scale and sumsq must be supplied in SCALE and SUMSQ and
!  scl and smsq are overwritten on SCALE and SUMSQ respectively.
!
!  The routine makes only one pass through the vector x.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of elements to be used from the vector X.
!
!  X       (input) DOUBLE PRECISION array, dimension (N)
!          The vector for which a scaled sum of squares is computed.
!             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.
!
!  INCX    (input) INTEGER
!          The increment between successive values of the vector X.
!          INCX > 0.
!
!  SCALE   (input/output) DOUBLE PRECISION
!          On entry, the value  scale  in the equation above.
!          On exit, SCALE is overwritten with  scl , the scaling factor
!          for the sum of squares.
!
!  SUMSQ   (input/output) DOUBLE PRECISION
!          On entry, the value  sumsq  in the equation above.
!          On exit, SUMSQ is overwritten with  smsq , the basic sum of
!          squares from which  scl  has been factored out.
!
! =====================================================================
!
!     .. Parameters ..
  real(dp)           ZERO
  PARAMETER          ( ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
  integer(i4)        IX
  real(dp)           ABSXI
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          ABS
!     ..
!     .. Executable Statements ..
!
  IF( N.GT.0 ) THEN
     DO 10 IX = 1, 1 + ( N-1 )*INCX, INCX
        IF( X( IX ).NE.ZERO ) THEN
           ABSXI = ABS( X( IX ) )
           IF( SCALE.LT.ABSXI ) THEN
              SUMSQ = 1 + SUMSQ*( SCALE / ABSXI )**2
              SCALE = ABSXI
           ELSE
              SUMSQ = SUMSQ + ( ABSXI / SCALE )**2
           END IF
        END IF
10    CONTINUE
  END IF
  RETURN
!
!     End of DLASSQ
!
  END
  SUBROUTINE DORG2L( M, N, K, A, LDA, TAU, WORK, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
  use datatypes
  integer(i4)        INFO, K, LDA, M, N
!     ..
!     .. Array Arguments ..
  real(dp)           A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DORG2L generates an m by n real matrix Q with orthonormal columns,
!  which is defined as the last n columns of a product of k elementary
!  reflectors of order m
!
!        Q  =  H(k) . . . H(2) H(1)
!
!  as returned by DGEQLF.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix Q. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix Q. M >= N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines the
!          matrix Q. N >= K >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the (n-k+i)-th column must contain the vector which
!          defines the elementary reflector H(i), for i = 1,2,...,k, as
!          returned by DGEQLF in the last k columns of its array
!          argument A.
!          On exit, the m by n matrix Q.
!
!  LDA     (input) INTEGER
!          The first dimension of the array A. LDA >= max(1,M).
!
!  TAU     (input) DOUBLE PRECISION array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGEQLF.
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument has an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
  real(dp)           ONE, ZERO
  PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
  integer(i4)        I, II, J, L
!     ..
!     .. External Subroutines ..
  EXTERNAL           DLARF, DSCAL, XERBLA
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
  INFO = 0
  IF( M.LT.0 ) THEN
     INFO = -1
  ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
     INFO = -2
  ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
     INFO = -3
  ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
     INFO = -5
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DORG2L', -INFO )
     RETURN
  END IF
!
!     Quick return if possible
!
  IF( N.LE.0 ) &
     RETURN
!
!     Initialise columns 1:n-k to columns of the unit matrix
!
  DO 20 J = 1, N - K
     DO 10 L = 1, M
        A( L, J ) = ZERO
10    CONTINUE
     A( M-N+J, J ) = ONE
20 CONTINUE
!
  DO 40 I = 1, K
     II = N - K + I
!
!        Apply H(i) to A(1:m-k+i,1:n-k+i) from the left
!
     A( M-N+II, II ) = ONE
     CALL DLARF( 'Left', M-N+II, II-1_i4, A( 1, II ), 1_i4, TAU( I ), A, &
                 LDA, WORK )
     CALL DSCAL( M-N+II-1_i4, -TAU( I ), A( 1, II ), 1_i4 )
     A( M-N+II, II ) = ONE - TAU( I )
!
!        Set A(m-k+i+1:m,n-k+i) to zero
!
     DO 30 L = M - N + II + 1, M
        A( L, II ) = ZERO
30    CONTINUE
40 CONTINUE
  RETURN
!
!     End of DORG2L
!
  END
  SUBROUTINE DORG2R( M, N, K, A, LDA, TAU, WORK, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
  use datatypes
  integer(i4)        INFO, K, LDA, M, N
!     ..
!     .. Array Arguments ..
  real(dp)           A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DORG2R generates an m by n real matrix Q with orthonormal columns,
!  which is defined as the first n columns of a product of k elementary
!  reflectors of order m
!
!        Q  =  H(1) H(2) . . . H(k)
!
!  as returned by DGEQRF.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix Q. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix Q. M >= N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines the
!          matrix Q. N >= K >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the i-th column must contain the vector which
!          defines the elementary reflector H(i), for i = 1,2,...,k, as
!          returned by DGEQRF in the first k columns of its array
!          argument A.
!          On exit, the m-by-n matrix Q.
!
!  LDA     (input) INTEGER
!          The first dimension of the array A. LDA >= max(1,M).
!
!  TAU     (input) DOUBLE PRECISION array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGEQRF.
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument has an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
  real(dp)           ONE, ZERO
  PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
  integer(i4)        I, J, L
!     ..
!     .. External Subroutines ..
  EXTERNAL           DLARF, DSCAL, XERBLA
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
  INFO = 0
  IF( M.LT.0 ) THEN
     INFO = -1
  ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
     INFO = -2
  ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
     INFO = -3
  ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
     INFO = -5
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DORG2R', -INFO )
     RETURN
  END IF
!
!     Quick return if possible
!
  IF( N.LE.0 ) &
     RETURN
!
!     Initialise columns k+1:n to columns of the unit matrix
!
  DO 20 J = K + 1, N
     DO 10 L = 1, M
        A( L, J ) = ZERO
10    CONTINUE
     A( J, J ) = ONE
20 CONTINUE
!
  DO 40 I = K, 1, -1
!
!        Apply H(i) to A(i:m,i:n) from the left
!
     IF( I.LT.N ) THEN
        A( I, I ) = ONE
        CALL DLARF( 'Left', M-I+1_i4, N-I, A( I, I ), 1_i4, TAU( I ), &
                    A( I, I+1 ), LDA, WORK )
     END IF
     IF( I.LT.M ) &
        CALL DSCAL( M-I, -TAU( I ), A( I+1, I ), 1_i4 )
     A( I, I ) = ONE - TAU( I )
!
!        Set A(1:i-1,i) to zero
!
     DO 30 L = 1, I - 1
        A( L, I ) = ZERO
30    CONTINUE
40 CONTINUE
  RETURN
!
!     End of DORG2R
!
  END
  SUBROUTINE DLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
  use datatypes
  CHARACTER          SIDE
  integer(i4)        INCV, LDC, M, N
  real(dp)           TAU
!     ..
!     .. Array Arguments ..
  real(dp)           C( LDC, * ), V( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DLARF applies a real elementary reflector H to a real m by n matrix
!  C, from either the left or the right. H is represented in the form
!
!        H = I - tau * v * v'
!
!  where tau is a real scalar and v is a real vector.
!
!  If tau = 0, then H is taken to be the unit matrix.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': form  H * C
!          = 'R': form  C * H
!
!  M       (input) INTEGER
!          The number of rows of the matrix C.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C.
!
!  V       (input) DOUBLE PRECISION array, dimension
!                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
!                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
!          The vector v in the representation of H. V is not used if
!          TAU = 0.
!
!  INCV    (input) INTEGER
!          The increment between elements of v. INCV <> 0.
!
!  TAU     (input) DOUBLE PRECISION
!          The value tau in the representation of H.
!
!  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
!          On entry, the m by n matrix C.
!          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
!          or C * H if SIDE = 'R'.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension
!                         (N) if SIDE = 'L'
!                      or (M) if SIDE = 'R'
!
!  =====================================================================
!
!     .. Parameters ..
  real(dp)           ONE, ZERO
  PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. External Subroutines ..
  EXTERNAL           DGEMV, DGER
!     ..
!     .. External Functions ..
  LOGICAL            LSAME
  EXTERNAL           LSAME
!     ..
!     .. Executable Statements ..
!
  IF( LSAME( SIDE, 'L' ) ) THEN
!
!        Form  H * C
!
     IF( TAU.NE.ZERO ) THEN
!
!           w := C' * v
!
        CALL DGEMV( 'Transpose', M, N, ONE, C, LDC, V, INCV, ZERO, &
                    WORK, 1_i4 )
!
!           C := C - v * w'
!
        CALL DGER( M, N, -TAU, V, INCV, WORK, 1_i4, C, LDC )
     END IF
  ELSE
!
!        Form  C * H
!
     IF( TAU.NE.ZERO ) THEN
!
!           w := C * v
!
        CALL DGEMV( 'No transpose', M, N, ONE, C, LDC, V, INCV, &
                    ZERO, WORK, 1_i4 )
!
!           C := C - w * v'
!
        CALL DGER( M, N, -TAU, WORK, 1_i4, V, INCV, C, LDC )
     END IF
  END IF
  RETURN
!
!     End of DLARF
!
  END
  SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
!
!  -- LAPACK driver routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
  use datatypes
  CHARACTER          JOBZ, UPLO
  integer(i4)        INFO, LDA, LWORK, N
!     ..
!     .. Array Arguments ..
  real(dp)           A( LDA, * ), W( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DSYEV computes all eigenvalues and, optionally, eigenvectors of a
!  real symmetric matrix A.
!
!  Arguments
!  =========
!
!  JOBZ    (input) CHARACTER*1
!          = 'N':  Compute eigenvalues only;
!          = 'V':  Compute eigenvalues and eigenvectors.
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the
!          leading N-by-N upper triangular part of A contains the
!          upper triangular part of the matrix A.  If UPLO = 'L',
!          the leading N-by-N lower triangular part of A contains
!          the lower triangular part of the matrix A.
!          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
!          orthonormal eigenvectors of the matrix A.
!          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
!          or the upper triangle (if UPLO='U') of A, including the
!          diagonal, is destroyed.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  W       (output) DOUBLE PRECISION array, dimension (N)
!          If INFO = 0, the eigenvalues in ascending order.
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The length of the array WORK.  LWORK >= max(1,3*N-1).
!          For optimal efficiency, LWORK >= (NB+2)*N,
!          where NB is the blocksize for DSYTRD returned by ILAENV.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, the algorithm failed to converge; i
!                off-diagonal elements of an intermediate tridiagonal
!                form did not converge to zero.
!
!  =====================================================================
!
!     .. Parameters ..
  real(dp)           ZERO, ONE
  PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
  LOGICAL            LOWER, LQUERY, WANTZ
  integer(i4)        IINFO, IMAX, INDE, INDTAU, INDWRK, ISCALE, &
                     LLWORK, LOPT, LWKOPT, NB
  real(dp)           ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, &
                     SMLNUM
!     ..
!     .. External Functions ..
  LOGICAL            LSAME
  integer(i4)        ILAENV
  real(dp)           DLAMCH, DLANSY
  EXTERNAL           LSAME, ILAENV, DLAMCH, DLANSY
!     ..
!     .. External Subroutines ..
  EXTERNAL           DLASCL, DORGTR, DSCAL, DSTEQR, DSTERF, DSYTRD, &
                     XERBLA
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          MAX, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
  WANTZ = LSAME( JOBZ, 'V' )
  LOWER = LSAME( UPLO, 'L' )
  LQUERY = ( LWORK.EQ.-1 )
!
  INFO = 0
  IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
     INFO = -1
  ELSE IF( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) THEN
     INFO = -2
  ELSE IF( N.LT.0 ) THEN
     INFO = -3
  ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
     INFO = -5
  ELSE IF( LWORK.LT.MAX( 1, 3*N-1 ) .AND. .NOT.LQUERY ) THEN
     INFO = -8
  END IF
!
  IF( INFO.EQ.0 ) THEN
     NB = ILAENV( 1_i4, 'DSYTRD', UPLO, N, -1_i4, -1_i4, -1_i4 )
     LWKOPT = MAX( 1_i4, ( NB+2_i4 )*N )
     WORK( 1 ) = LWKOPT
  END IF
!
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DSYEV ', -INFO )
     RETURN
  ELSE IF( LQUERY ) THEN
     RETURN
  END IF
!
!     Quick return if possible
!
  IF( N.EQ.0 ) THEN
     WORK( 1 ) = 1
     RETURN
  END IF
!
  IF( N.EQ.1 ) THEN
     W( 1 ) = A( 1, 1 )
     WORK( 1 ) = 3
     IF( WANTZ ) &
        A( 1, 1 ) = ONE
     RETURN
  END IF
!
!     Get machine constants.
!
  SAFMIN = DLAMCH( 'Safe minimum' )
  EPS = DLAMCH( 'Precision' )
  SMLNUM = SAFMIN / EPS
  BIGNUM = ONE / SMLNUM
  RMIN = SQRT( SMLNUM )
  RMAX = SQRT( BIGNUM )
!
!     Scale matrix to allowable range, if necessary.
!
  ANRM = DLANSY( 'M', UPLO, N, A, LDA, WORK )
  ISCALE = 0
  IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
     ISCALE = 1
     SIGMA = RMIN / ANRM
  ELSE IF( ANRM.GT.RMAX ) THEN
     ISCALE = 1
     SIGMA = RMAX / ANRM
  END IF
  IF( ISCALE.EQ.1 ) &
     CALL DLASCL( UPLO, 0_i4, 0_i4, ONE, SIGMA, N, N, A, LDA, INFO )
!
!     Call DSYTRD to reduce symmetric matrix to tridiagonal form.
!
  INDE = 1
  INDTAU = INDE + N
  INDWRK = INDTAU + N
  LLWORK = LWORK - INDWRK + 1
  CALL DSYTRD( UPLO, N, A, LDA, W, WORK( INDE ), WORK( INDTAU ), &
               WORK( INDWRK ), LLWORK, IINFO )
  LOPT = 2*N + NINT(WORK( INDWRK ))
!
!     For eigenvalues only, call DSTERF.  For eigenvectors, first call
!     DORGTR to generate the orthogonal matrix, then call DSTEQR.
!
  IF( .NOT.WANTZ ) THEN
     CALL DSTERF( N, W, WORK( INDE ), INFO )
  ELSE
     CALL DORGTR( UPLO, N, A, LDA, WORK( INDTAU ), WORK( INDWRK ), &
                  LLWORK, IINFO )
     CALL DSTEQR( JOBZ, N, W, WORK( INDE ), A, LDA, WORK( INDTAU ), &
                  INFO )
  END IF
!
!     If matrix was scaled, then rescale eigenvalues appropriately.
!
  IF( ISCALE.EQ.1 ) THEN
     IF( INFO.EQ.0 ) THEN
        IMAX = N
     ELSE
        IMAX = INFO - 1
     END IF
     CALL DSCAL( IMAX, ONE / SIGMA, W, 1_i4 )
  END IF
!
!     Set WORK(1) to optimal workspace size.
!
  WORK( 1 ) = LWKOPT
!
  RETURN
!
!     End of DSYEV
!
  END
  SUBROUTINE DLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV, &
                     T, LDT, C, LDC, WORK, LDWORK )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
  use datatypes
  CHARACTER          DIRECT, SIDE, STOREV, TRANS
  integer(i4)        K, LDC, LDT, LDV, LDWORK, M, N
!     ..
!     .. Array Arguments ..
  real(dp)           C( LDC, * ), T( LDT, * ), V( LDV, * ), &
                     WORK( LDWORK, * )
!     ..
!
!  Purpose
!  =======
!
!  DLARFB applies a real block reflector H or its transpose H' to a
!  real m by n matrix C, from either the left or the right.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': apply H or H' from the Left
!          = 'R': apply H or H' from the Right
!
!  TRANS   (input) CHARACTER*1
!          = 'N': apply H (No transpose)
!          = 'T': apply H' (Transpose)
!
!  DIRECT  (input) CHARACTER*1
!          Indicates how H is formed from a product of elementary
!          reflectors
!          = 'F': H = H(1) H(2) . . . H(k) (Forward)
!          = 'B': H = H(k) . . . H(2) H(1) (Backward)
!
!  STOREV  (input) CHARACTER*1
!          Indicates how the vectors which define the elementary
!          reflectors are stored:
!          = 'C': Columnwise
!          = 'R': Rowwise
!
!  M       (input) INTEGER
!          The number of rows of the matrix C.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C.
!
!  K       (input) INTEGER
!          The order of the matrix T (= the number of elementary
!          reflectors whose product defines the block reflector).
!
!  V       (input) DOUBLE PRECISION array, dimension
!                                (LDV,K) if STOREV = 'C'
!                                (LDV,M) if STOREV = 'R' and SIDE = 'L'
!                                (LDV,N) if STOREV = 'R' and SIDE = 'R'
!          The matrix V. See further details.
!
!  LDV     (input) INTEGER
!          The leading dimension of the array V.
!          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);
!          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);
!          if STOREV = 'R', LDV >= K.
!
!  T       (input) DOUBLE PRECISION array, dimension (LDT,K)
!          The triangular k by k matrix T in the representation of the
!          block reflector.
!
!  LDT     (input) INTEGER
!          The leading dimension of the array T. LDT >= K.
!
!  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
!          On entry, the m by n matrix C.
!          On exit, C is overwritten by H*C or H'*C or C*H or C*H'.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDA >= max(1,M).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (LDWORK,K)
!
!  LDWORK  (input) INTEGER
!          The leading dimension of the array WORK.
!          If SIDE = 'L', LDWORK >= max(1,N);
!          if SIDE = 'R', LDWORK >= max(1,M).
!
!  =====================================================================
!
!     .. Parameters ..
  real(dp)           ONE
  PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
  CHARACTER          TRANST
  integer(i4)        I, J
!     ..
!     .. External Functions ..
  LOGICAL            LSAME
  EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
  EXTERNAL           DCOPY, DGEMM, DTRMM
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
  IF( M.LE.0 .OR. N.LE.0 ) &
     RETURN
!
  IF( LSAME( TRANS, 'N' ) ) THEN
     TRANST = 'T'
  ELSE
     TRANST = 'N'
  END IF
!
  IF( LSAME( STOREV, 'C' ) ) THEN
!
     IF( LSAME( DIRECT, 'F' ) ) THEN
!
!           Let  V =  ( V1 )    (first K rows)
!                     ( V2 )
!           where  V1  is unit lower triangular.
!
        IF( LSAME( SIDE, 'L' ) ) THEN
!
!              Form  H * C  or  H' * C  where  C = ( C1 )
!                                                  ( C2 )
!
!              W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)
!
!              W := C1'
!
           DO 10 J = 1, K
              CALL DCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1_i4 )
10          CONTINUE
!
!              W := W * V1
!
           CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N, &
                       K, ONE, V, LDV, WORK, LDWORK )
           IF( M.GT.K ) THEN
!
!                 W := W + C2'*V2
!
              CALL DGEMM( 'Transpose', 'No transpose', N, K, M-K, &
                          ONE, C( K+1, 1 ), LDC, V( K+1, 1 ), LDV, &
                          ONE, WORK, LDWORK )
           END IF
!
!              W := W * T'  or  W * T
!
           CALL DTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K, &
                       ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V * W'
!
           IF( M.GT.K ) THEN
!
!                 C2 := C2 - V2 * W'
!
              CALL DGEMM( 'No transpose', 'Transpose', M-K, N, K, &
                          -ONE, V( K+1, 1 ), LDV, WORK, LDWORK, ONE, &
                          C( K+1, 1 ), LDC )
           END IF
!
!              W := W * V1'
!
           CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K, &
                       ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W'
!
           DO 30 J = 1, K
              DO 20 I = 1, N
                 C( J, I ) = C( J, I ) - WORK( I, J )
20             CONTINUE
30          CONTINUE
!
        ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H'  where  C = ( C1  C2 )
!
!              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
!
!              W := C1
!
           DO 40 J = 1, K
              CALL DCOPY( M, C( 1, J ), 1_i4, WORK( 1, J ), 1_i4 )
40          CONTINUE
!
!              W := W * V1
!
           CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M, &
                       K, ONE, V, LDV, WORK, LDWORK )
           IF( N.GT.K ) THEN
!
!                 W := W + C2 * V2
!
              CALL DGEMM( 'No transpose', 'No transpose', M, K, N-K, &
                          ONE, C( 1, K+1 ), LDC, V( K+1, 1 ), LDV, &
                          ONE, WORK, LDWORK )
           END IF
!
!              W := W * T  or  W * T'
!
           CALL DTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K, &
                       ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V'
!
           IF( N.GT.K ) THEN
!
!                 C2 := C2 - W * V2'
!
              CALL DGEMM( 'No transpose', 'Transpose', M, N-K, K, &
                          -ONE, WORK, LDWORK, V( K+1, 1 ), LDV, ONE, &
                          C( 1, K+1 ), LDC )
           END IF
!
!              W := W * V1'
!
           CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', M, K, &
                       ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W
!
           DO 60 J = 1, K
              DO 50 I = 1, M
                 C( I, J ) = C( I, J ) - WORK( I, J )
50             CONTINUE
60          CONTINUE
        END IF
!
     ELSE
!
!           Let  V =  ( V1 )
!                     ( V2 )    (last K rows)
!           where  V2  is unit upper triangular.
!
        IF( LSAME( SIDE, 'L' ) ) THEN
!
!              Form  H * C  or  H' * C  where  C = ( C1 )
!                                                  ( C2 )
!
!              W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)
!
!              W := C2'
!
           DO 70 J = 1, K
              CALL DCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1_i4 )
70          CONTINUE
!
!              W := W * V2
!
           CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N, &
                       K, ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
           IF( M.GT.K ) THEN
!
!                 W := W + C1'*V1
!
              CALL DGEMM( 'Transpose', 'No transpose', N, K, M-K, &
                          ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
           END IF
!
!              W := W * T'  or  W * T
!
           CALL DTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K, &
                       ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V * W'
!
           IF( M.GT.K ) THEN
!
!                 C1 := C1 - V1 * W'
!
              CALL DGEMM( 'No transpose', 'Transpose', M-K, N, K, &
                          -ONE, V, LDV, WORK, LDWORK, ONE, C, LDC )
           END IF
!
!              W := W * V2'
!
           CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', N, K, &
                       ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
!
!              C2 := C2 - W'
!
           DO 90 J = 1, K
              DO 80 I = 1, N
                 C( M-K+J, I ) = C( M-K+J, I ) - WORK( I, J )
80             CONTINUE
90          CONTINUE
!
        ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H'  where  C = ( C1  C2 )
!
!              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
!
!              W := C2
!
           DO 100 J = 1, K
              CALL DCOPY( M, C( 1, N-K+J ), 1_i4, WORK( 1, J ), 1_i4 )
100          CONTINUE
!
!              W := W * V2
!
           CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M, &
                       K, ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
           IF( N.GT.K ) THEN
!
!                 W := W + C1 * V1
!
              CALL DGEMM( 'No transpose', 'No transpose', M, K, N-K, &
                          ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
           END IF
!
!              W := W * T  or  W * T'
!
           CALL DTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K, &
                       ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V'
!
           IF( N.GT.K ) THEN
!
!                 C1 := C1 - W * V1'
!
              CALL DGEMM( 'No transpose', 'Transpose', M, N-K, K, &
                          -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
           END IF
!
!              W := W * V2'
!
           CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', M, K, &
                       ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
!
!              C2 := C2 - W
!
           DO 120 J = 1, K
              DO 110 I = 1, M
                 C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
110             CONTINUE
120          CONTINUE
        END IF
     END IF
!
  ELSE IF( LSAME( STOREV, 'R' ) ) THEN
!
     IF( LSAME( DIRECT, 'F' ) ) THEN
!
!           Let  V =  ( V1  V2 )    (V1: first K columns)
!           where  V1  is unit upper triangular.
!
        IF( LSAME( SIDE, 'L' ) ) THEN
!
!              Form  H * C  or  H' * C  where  C = ( C1 )
!                                                  ( C2 )
!
!              W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)
!
!              W := C1'
!
           DO 130 J = 1, K
              CALL DCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1_i4 )
130          CONTINUE
!
!              W := W * V1'
!
           CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', N, K, &
                       ONE, V, LDV, WORK, LDWORK )
           IF( M.GT.K ) THEN
!
!                 W := W + C2'*V2'
!
              CALL DGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE, &
                          C( K+1, 1 ), LDC, V( 1, K+1 ), LDV, ONE, &
                          WORK, LDWORK )
           END IF
!
!              W := W * T'  or  W * T
!
           CALL DTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K, &
                       ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V' * W'
!
           IF( M.GT.K ) THEN
!
!                 C2 := C2 - V2' * W'
!
              CALL DGEMM( 'Transpose', 'Transpose', M-K, N, K, -ONE, &
                          V( 1, K+1 ), LDV, WORK, LDWORK, ONE, &
                          C( K+1, 1 ), LDC )
           END IF
!
!              W := W * V1
!
           CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N, &
                       K, ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W'
!
           DO 150 J = 1, K
              DO 140 I = 1, N
                 C( J, I ) = C( J, I ) - WORK( I, J )
140             CONTINUE
150          CONTINUE
!
        ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H'  where  C = ( C1  C2 )
!
!              W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)
!
!              W := C1
!
           DO 160 J = 1, K
              CALL DCOPY( M, C( 1, J ), 1_i4, WORK( 1, J ), 1_i4 )
160          CONTINUE
!
!              W := W * V1'
!
           CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', M, K, &
                       ONE, V, LDV, WORK, LDWORK )
           IF( N.GT.K ) THEN
!
!                 W := W + C2 * V2'
!
              CALL DGEMM( 'No transpose', 'Transpose', M, K, N-K, &
                          ONE, C( 1, K+1 ), LDC, V( 1, K+1 ), LDV, &
                          ONE, WORK, LDWORK )
           END IF
!
!              W := W * T  or  W * T'
!
           CALL DTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K, &
                       ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V
!
           IF( N.GT.K ) THEN
!
!                 C2 := C2 - W * V2
!
              CALL DGEMM( 'No transpose', 'No transpose', M, N-K, K, &
                          -ONE, WORK, LDWORK, V( 1, K+1 ), LDV, ONE, &
                          C( 1, K+1 ), LDC )
           END IF
!
!              W := W * V1
!
           CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M, &
                       K, ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W
!
           DO 180 J = 1, K
              DO 170 I = 1, M
                 C( I, J ) = C( I, J ) - WORK( I, J )
170             CONTINUE
180          CONTINUE
!
        END IF
!
     ELSE
!
!           Let  V =  ( V1  V2 )    (V2: last K columns)
!           where  V2  is unit lower triangular.
!
        IF( LSAME( SIDE, 'L' ) ) THEN
!
!              Form  H * C  or  H' * C  where  C = ( C1 )
!                                                  ( C2 )
!
!              W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)
!
!              W := C2'
!
           DO 190 J = 1, K
              CALL DCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1_i4 )
190          CONTINUE
!
!              W := W * V2'
!
           CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K, &
                       ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
           IF( M.GT.K ) THEN
!
!                 W := W + C1'*V1'
!
              CALL DGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE, &
                          C, LDC, V, LDV, ONE, WORK, LDWORK )
           END IF
!
!              W := W * T'  or  W * T
!
           CALL DTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K, &
                       ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V' * W'
!
           IF( M.GT.K ) THEN
!
!                 C1 := C1 - V1' * W'
!
              CALL DGEMM( 'Transpose', 'Transpose', M-K, N, K, -ONE, &
                          V, LDV, WORK, LDWORK, ONE, C, LDC )
           END IF
!
!              W := W * V2
!
           CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N, &
                       K, ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
!
!              C2 := C2 - W'
!
           DO 210 J = 1, K
              DO 200 I = 1, N
                 C( M-K+J, I ) = C( M-K+J, I ) - WORK( I, J )
200             CONTINUE
210          CONTINUE
!
        ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H'  where  C = ( C1  C2 )
!
!              W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)
!
!              W := C2
!
           DO 220 J = 1, K
              CALL DCOPY( M, C( 1, N-K+J ), 1_i4, WORK( 1, J ), 1_i4 )
220          CONTINUE
!
!              W := W * V2'
!
           CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', M, K, &
                       ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
           IF( N.GT.K ) THEN
!
!                 W := W + C1 * V1'
!
              CALL DGEMM( 'No transpose', 'Transpose', M, K, N-K, &
                          ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
           END IF
!
!              W := W * T  or  W * T'
!
           CALL DTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K, &
                       ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V
!
           IF( N.GT.K ) THEN
!
!                 C1 := C1 - W * V1
!
              CALL DGEMM( 'No transpose', 'No transpose', M, N-K, K, &
                          -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
           END IF
!
!              W := W * V2
!
           CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M, &
                       K, ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
!
!              C1 := C1 - W
!
           DO 240 J = 1, K
              DO 230 I = 1, M
                 C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
230             CONTINUE
240          CONTINUE
!
        END IF
!
     END IF
  END IF
!
  RETURN
!
!     End of DLARFB
!
  END
  SUBROUTINE DLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
  use datatypes
  CHARACTER          DIRECT, STOREV
  integer(i4)        K, LDT, LDV, N
!     ..
!     .. Array Arguments ..
  real(dp)           T( LDT, * ), TAU( * ), V( LDV, * )
!     ..
!
!  Purpose
!  =======
!
!  DLARFT forms the triangular factor T of a real block reflector H
!  of order n, which is defined as a product of k elementary reflectors.
!
!  If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
!
!  If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
!
!  If STOREV = 'C', the vector which defines the elementary reflector
!  H(i) is stored in the i-th column of the array V, and
!
!     H  =  I - V * T * V'
!
!  If STOREV = 'R', the vector which defines the elementary reflector
!  H(i) is stored in the i-th row of the array V, and
!
!     H  =  I - V' * T * V
!
!  Arguments
!  =========
!
!  DIRECT  (input) CHARACTER*1
!          Specifies the order in which the elementary reflectors are
!          multiplied to form the block reflector:
!          = 'F': H = H(1) H(2) . . . H(k) (Forward)
!          = 'B': H = H(k) . . . H(2) H(1) (Backward)
!
!  STOREV  (input) CHARACTER*1
!          Specifies how the vectors which define the elementary
!          reflectors are stored (see also Further Details):
!          = 'C': columnwise
!          = 'R': rowwise
!
!  N       (input) INTEGER
!          The order of the block reflector H. N >= 0.
!
!  K       (input) INTEGER
!          The order of the triangular factor T (= the number of
!          elementary reflectors). K >= 1.
!
!  V       (input/output) DOUBLE PRECISION array, dimension
!                               (LDV,K) if STOREV = 'C'
!                               (LDV,N) if STOREV = 'R'
!          The matrix V. See further details.
!
!  LDV     (input) INTEGER
!          The leading dimension of the array V.
!          If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K.
!
!  TAU     (input) DOUBLE PRECISION array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i).
!
!  T       (output) DOUBLE PRECISION array, dimension (LDT,K)
!          The k by k triangular factor T of the block reflector.
!          If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is
!          lower triangular. The rest of the array is not used.
!
!  LDT     (input) INTEGER
!          The leading dimension of the array T. LDT >= K.
!
!  Further Details
!  ===============
!
!  The shape of the matrix V and the storage of the vectors which define
!  the H(i) is best illustrated by the following example with n = 5 and
!  k = 3. The elements equal to 1 are not stored; the corresponding
!  array elements are modified but restored on exit. The rest of the
!  array is not used.
!
!  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
!
!               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
!                   ( v1  1    )                     (     1 v2 v2 v2 )
!                   ( v1 v2  1 )                     (        1 v3 v3 )
!                   ( v1 v2 v3 )
!                   ( v1 v2 v3 )
!
!  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
!
!               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
!                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
!                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
!                   (     1 v3 )
!                   (        1 )
!
!  =====================================================================
!
!     .. Parameters ..
  real(dp)           ONE, ZERO
  PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
  integer(i4)        I, J
  real(dp)           VII
!     ..
!     .. External Subroutines ..
  EXTERNAL           DGEMV, DTRMV
!     ..
!     .. External Functions ..
  LOGICAL            LSAME
  EXTERNAL           LSAME
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
  IF( N.EQ.0 ) &
     RETURN
!
  IF( LSAME( DIRECT, 'F' ) ) THEN
     DO 20 I = 1, K
        IF( TAU( I ).EQ.ZERO ) THEN
!
!              H(i)  =  I
!
           DO 10 J = 1, I
              T( J, I ) = ZERO
10          CONTINUE
        ELSE
!
!              general case
!
           VII = V( I, I )
           V( I, I ) = ONE
           IF( LSAME( STOREV, 'C' ) ) THEN
!
!                 T(1:i-1,i) := - tau(i) * V(i:n,1:i-1)' * V(i:n,i)
!
              CALL DGEMV( 'Transpose', N-I+1_i4, I-1_i4, -TAU( I ), &
                          V( I, 1 ), LDV, V( I, I ), 1_i4, ZERO, &
                          T( 1, I ), 1_i4 )
           ELSE
!
!                 T(1:i-1,i) := - tau(i) * V(1:i-1,i:n) * V(i,i:n)'
!
              CALL DGEMV( 'No transpose', I-1_i4, N-I+1_i4, -TAU( I ), &
                          V( 1, I ), LDV, V( I, I ), LDV, ZERO, &
                          T( 1, I ), 1_i4 )
           END IF
           V( I, I ) = VII
!
!              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)
!
           CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I-1_i4, T, &
                       LDT, T( 1, I ), 1_i4 )
           T( I, I ) = TAU( I )
        END IF
20    CONTINUE
  ELSE
     DO 40 I = K, 1, -1
        IF( TAU( I ).EQ.ZERO ) THEN
!
!              H(i)  =  I
!
           DO 30 J = I, K
              T( J, I ) = ZERO
30          CONTINUE
        ELSE
!
!              general case
!
           IF( I.LT.K ) THEN
              IF( LSAME( STOREV, 'C' ) ) THEN
                 VII = V( N-K+I, I )
                 V( N-K+I, I ) = ONE
!
!                    T(i+1:k,i) :=
!                            - tau(i) * V(1:n-k+i,i+1:k)' * V(1:n-k+i,i)
!
                 CALL DGEMV( 'Transpose', N-K+I, K-I, -TAU( I ), &
                             V( 1, I+1 ), LDV, V( 1, I ), 1_i4, ZERO, &
                             T( I+1, I ), 1_i4 )
                 V( N-K+I, I ) = VII
              ELSE
                 VII = V( I, N-K+I )
                 V( I, N-K+I ) = ONE
!
!                    T(i+1:k,i) :=
!                            - tau(i) * V(i+1:k,1:n-k+i) * V(i,1:n-k+i)'
!
                 CALL DGEMV( 'No transpose', K-I, N-K+I, -TAU( I ), &
                             V( I+1, 1 ), LDV, V( I, 1 ), LDV, ZERO, &
                             T( I+1, I ), 1_i4 )
                 V( I, N-K+I ) = VII
              END IF
!
!                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)
!
              CALL DTRMV( 'Lower', 'No transpose', 'Non-unit', K-I, &
                          T( I+1, I+1 ), LDT, T( I+1, I ), 1_i4 )
           END IF
           T( I, I ) = TAU( I )
        END IF
40    CONTINUE
  END IF
  RETURN
!
!     End of DLARFT
!
  END
  SUBROUTINE DLATRD( UPLO, N, NB, A, LDA, E, TAU, W, LDW )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
  use datatypes
  CHARACTER          UPLO
  integer(i4)        LDA, LDW, N, NB
!     ..
!     .. Array Arguments ..
  real(dp)           A( LDA, * ), E( * ), TAU( * ), W( LDW, * )
!     ..
!
!  Purpose
!  =======
!
!  DLATRD reduces NB rows and columns of a real symmetric matrix A to
!  symmetric tridiagonal form by an orthogonal similarity
!  transformation Q' * A * Q, and returns the matrices V and W which are
!  needed to apply the transformation to the unreduced part of A.
!
!  If UPLO = 'U', DLATRD reduces the last NB rows and columns of a
!  matrix, of which the upper triangle is supplied;
!  if UPLO = 'L', DLATRD reduces the first NB rows and columns of a
!  matrix, of which the lower triangle is supplied.
!
!  This is an auxiliary routine called by DSYTRD.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER
!          Specifies whether the upper or lower triangular part of the
!          symmetric matrix A is stored:
!          = 'U': Upper triangular
!          = 'L': Lower triangular
!
!  N       (input) INTEGER
!          The order of the matrix A.
!
!  NB      (input) INTEGER
!          The number of rows and columns to be reduced.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!          n-by-n upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading n-by-n lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!          On exit:
!          if UPLO = 'U', the last NB columns have been reduced to
!            tridiagonal form, with the diagonal elements overwriting
!            the diagonal elements of A; the elements above the diagonal
!            with the array TAU, represent the orthogonal matrix Q as a
!            product of elementary reflectors;
!          if UPLO = 'L', the first NB columns have been reduced to
!            tridiagonal form, with the diagonal elements overwriting
!            the diagonal elements of A; the elements below the diagonal
!            with the array TAU, represent the  orthogonal matrix Q as a
!            product of elementary reflectors.
!          See Further Details.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= (1,N).
!
!  E       (output) DOUBLE PRECISION array, dimension (N-1)
!          If UPLO = 'U', E(n-nb:n-1) contains the superdiagonal
!          elements of the last NB columns of the reduced matrix;
!          if UPLO = 'L', E(1:nb) contains the subdiagonal elements of
!          the first NB columns of the reduced matrix.
!
!  TAU     (output) DOUBLE PRECISION array, dimension (N-1)
!          The scalar factors of the elementary reflectors, stored in
!          TAU(n-nb:n-1) if UPLO = 'U', and in TAU(1:nb) if UPLO = 'L'.
!          See Further Details.
!
!  W       (output) DOUBLE PRECISION array, dimension (LDW,NB)
!          The n-by-nb matrix W required to update the unreduced part
!          of A.
!
!  LDW     (input) INTEGER
!          The leading dimension of the array W. LDW >= max(1,N).
!
!  Further Details
!  ===============
!
!  If UPLO = 'U', the matrix Q is represented as a product of elementary
!  reflectors
!
!     Q = H(n) H(n-1) . . . H(n-nb+1).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(i:n) = 0 and v(i-1) = 1; v(1:i-1) is stored on exit in A(1:i-1,i),
!  and tau in TAU(i-1).
!
!  If UPLO = 'L', the matrix Q is represented as a product of elementary
!  reflectors
!
!     Q = H(1) H(2) . . . H(nb).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(1:i) = 0 and v(i+1) = 1; v(i+1:n) is stored on exit in A(i+1:n,i),
!  and tau in TAU(i).
!
!  The elements of the vectors v together form the n-by-nb matrix V
!  which is needed, with W, to apply the transformation to the unreduced
!  part of the matrix, using a symmetric rank-2k update of the form:
!  A := A - V*W' - W*V'.
!
!  The contents of A on exit are illustrated by the following examples
!  with n = 5 and nb = 2:
!
!  if UPLO = 'U':                       if UPLO = 'L':
!
!    (  a   a   a   v4  v5 )              (  d                  )
!    (      a   a   v4  v5 )              (  1   d              )
!    (          a   1   v5 )              (  v1  1   a          )
!    (              d   1  )              (  v1  v2  a   a      )
!    (                  d  )              (  v1  v2  a   a   a  )
!
!  where d denotes a diagonal element of the reduced matrix, a denotes
!  an element of the original matrix that is unchanged, and vi denotes
!  an element of the vector defining H(i).
!
!  =====================================================================
!
!     .. Parameters ..
  real(dp)           ZERO, ONE, HALF
  PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, HALF = 0.5D+0 )
!     ..
!     .. Local Scalars ..
  integer(i4)        I, IW
  real(dp)           ALPHA
!     ..
!     .. External Subroutines ..
  EXTERNAL           DAXPY, DGEMV, DLARFG, DSCAL, DSYMV
!     ..
!     .. External Functions ..
  LOGICAL            LSAME
  real(dp)           DDOT
  EXTERNAL           LSAME, DDOT
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          MIN
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
  IF( N.LE.0 ) &
     RETURN
!
  IF( LSAME( UPLO, 'U' ) ) THEN
!
!        Reduce last NB columns of upper triangle
!
     DO 10 I = N, N - NB + 1, -1
        IW = I - N + NB
        IF( I.LT.N ) THEN
!
!              Update A(1:i,i)
!
           CALL DGEMV( 'No transpose', I, N-I, -ONE, A( 1, I+1 ), &
                       LDA, W( I, IW+1 ), LDW, ONE, A( 1, I ), 1_i4 )
           CALL DGEMV( 'No transpose', I, N-I, -ONE, W( 1, IW+1 ), &
                       LDW, A( I, I+1 ), LDA, ONE, A( 1, I ), 1_i4 )
        END IF
        IF( I.GT.1 ) THEN
!
!              Generate elementary reflector H(i) to annihilate
!              A(1:i-2,i)
!
           CALL DLARFG( I-1_i4, A( I-1, I ), A( 1, I ), 1_i4, TAU( I-1 ) )
           E( I-1 ) = A( I-1, I )
           A( I-1, I ) = ONE
!
!              Compute W(1:i-1,i)
!
           CALL DSYMV( 'Upper', I-1_i4, ONE, A, LDA, A( 1, I ), 1_i4, &
                       ZERO, W( 1, IW ), 1_i4 )
           IF( I.LT.N ) THEN
              CALL DGEMV( 'Transpose', I-1_i4, N-I, ONE, W( 1, IW+1 ), &
                          LDW, A( 1, I ), 1_i4, ZERO, W( I+1, IW ), 1_i4 )
              CALL DGEMV( 'No transpose', I-1_i4, N-I, -ONE, &
                          A( 1, I+1 ), LDA, W( I+1, IW ), 1_i4, ONE, &
                          W( 1, IW ), 1_i4 )
              CALL DGEMV( 'Transpose', I-1_i4, N-I, ONE, A( 1, I+1 ), &
                          LDA, A( 1, I ), 1_i4, ZERO, W( I+1, IW ), 1_i4 )
              CALL DGEMV( 'No transpose', I-1_i4, N-I, -ONE, &
                          W( 1, IW+1 ), LDW, W( I+1, IW ), 1_i4, ONE, &
                          W( 1, IW ), 1_i4 )
           END IF
           CALL DSCAL( I-1_i4, TAU( I-1 ), W( 1, IW ), 1_i4 )
           ALPHA = -HALF*TAU( I-1 )*DDOT( I-1_i4, W( 1, IW ), 1_i4, &
                   A( 1, I ), 1_i4 )
           CALL DAXPY( I-1_i4, ALPHA, A( 1, I ), 1_i4, W( 1, IW ), 1_i4 )
        END IF
!
10    CONTINUE
  ELSE
!
!        Reduce first NB columns of lower triangle
!
     DO 20 I = 1, NB
!
!           Update A(i:n,i)
!
        CALL DGEMV( 'No transpose', N-I+1_i4, I-1_i4, -ONE, A( I, 1 ), &
                    LDA, W( I, 1 ), LDW, ONE, A( I, I ), 1_i4 )
        CALL DGEMV( 'No transpose', N-I+1_i4, I-1_i4, -ONE, W( I, 1 ), &
                    LDW, A( I, 1 ), LDA, ONE, A( I, I ), 1_i4 )
        IF( I.LT.N ) THEN
!
!              Generate elementary reflector H(i) to annihilate
!              A(i+2:n,i)
!
           CALL DLARFG( N-I, A( I+1, I ), A( MIN( I+2, N ), I ), 1_i4, &
                        TAU( I ) )
           E( I ) = A( I+1, I )
           A( I+1, I ) = ONE
!
!              Compute W(i+1:n,i)
!
           CALL DSYMV( 'Lower', N-I, ONE, A( I+1, I+1 ), LDA, &
                       A( I+1, I ), 1_i4, ZERO, W( I+1, I ), 1_i4 )
           CALL DGEMV( 'Transpose', N-I, I-1_i4, ONE, W( I+1, 1 ), LDW, &
                       A( I+1, I ), 1_i4, ZERO, W( 1, I ), 1_i4 )
           CALL DGEMV( 'No transpose', N-I, I-1_i4, -ONE, A( I+1, 1 ), &
                       LDA, W( 1, I ), 1_i4, ONE, W( I+1, I ), 1_i4 )
           CALL DGEMV( 'Transpose', N-I, I-1_i4, ONE, A( I+1, 1 ), LDA, &
                       A( I+1, I ), 1_i4, ZERO, W( 1, I ), 1_i4 )
           CALL DGEMV( 'No transpose', N-I, I-1_i4, -ONE, W( I+1, 1 ), &
                       LDW, W( 1, I ), 1_i4, ONE, W( I+1, I ), 1_i4 )
           CALL DSCAL( N-I, TAU( I ), W( I+1, I ), 1_i4 )
           ALPHA = -HALF*TAU( I )*DDOT( N-I, W( I+1, I ), 1_i4, &
                   A( I+1, I ), 1_i4 )
           CALL DAXPY( N-I, ALPHA, A( I+1, I ), 1_i4, W( I+1, I ), 1_i4 )
        END IF
!
20    CONTINUE
  END IF
!
  RETURN
!
!     End of DLATRD
!
  END
  SUBROUTINE DORGQL( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
  use datatypes
  integer(i4)        INFO, K, LDA, LWORK, M, N
!     ..
!     .. Array Arguments ..
  real(dp)           A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DORGQL generates an M-by-N real matrix Q with orthonormal columns,
!  which is defined as the last N columns of a product of K elementary
!  reflectors of order M
!
!        Q  =  H(k) . . . H(2) H(1)
!
!  as returned by DGEQLF.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix Q. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix Q. M >= N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines the
!          matrix Q. N >= K >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the (n-k+i)-th column must contain the vector which
!          defines the elementary reflector H(i), for i = 1,2,...,k, as
!          returned by DGEQLF in the last k columns of its array
!          argument A.
!          On exit, the M-by-N matrix Q.
!
!  LDA     (input) INTEGER
!          The first dimension of the array A. LDA >= max(1,M).
!
!  TAU     (input) DOUBLE PRECISION array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGEQLF.
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK. LWORK >= max(1,N).
!          For optimum performance LWORK >= N*NB, where NB is the
!          optimal blocksize.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument has an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
  real(dp)           ZERO
  PARAMETER          ( ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
  LOGICAL            LQUERY
  integer(i4)        I, IB, IINFO, IWS, J, KK, L, LDWORK, LWKOPT, &
                     NB, NBMIN, NX
!     ..
!     .. External Subroutines ..
  EXTERNAL           DLARFB, DLARFT, DORG2L, XERBLA
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          MAX, MIN
!     ..
!     .. External Functions ..
  integer(i4)        ILAENV
  EXTERNAL           ILAENV
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
  INFO = 0
  NB = ILAENV( 1_i4, 'DORGQL', ' ', M, N, K, -1_i4 )
  LWKOPT = MAX( 1_i4, N )*NB
  WORK( 1 ) = LWKOPT
  LQUERY = ( LWORK.EQ.-1 )
  IF( M.LT.0 ) THEN
     INFO = -1
  ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
     INFO = -2
  ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
     INFO = -3
  ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
     INFO = -5
  ELSE IF( LWORK.LT.MAX( 1_i4, N ) .AND. .NOT.LQUERY ) THEN
     INFO = -8
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DORGQL', -INFO )
     RETURN
  ELSE IF( LQUERY ) THEN
     RETURN
  END IF
!
!     Quick return if possible
!
  IF( N.LE.0 ) THEN
     WORK( 1 ) = 1
     RETURN
  END IF
!
  NBMIN = 2
  NX = 0
  IWS = N
  IF( NB.GT.1 .AND. NB.LT.K ) THEN
!
!        Determine when to cross over from blocked to unblocked code.
!
     NX = MAX( 0_i4, ILAENV( 3_i4, 'DORGQL', ' ', M, N, K, -1_i4 ) )
     IF( NX.LT.K ) THEN
!
!           Determine if workspace is large enough for blocked code.
!
        LDWORK = N
        IWS = LDWORK*NB
        IF( LWORK.LT.IWS ) THEN
!
!              Not enough workspace to use optimal NB:  reduce NB and
!              determine the minimum value of NB.
!
           NB = LWORK / LDWORK
           NBMIN = MAX( 2_i4, ILAENV( 2_i4, 'DORGQL', ' ', M, N, K, -1_i4 ) )
        END IF
     END IF
  END IF
!
  IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
!
!        Use blocked code after the first block.
!        The last kk columns are handled by the block method.
!
     KK = MIN( K, ( ( K-NX+NB-1 ) / NB )*NB )
!
!        Set A(m-kk+1:m,1:n-kk) to zero.
!
     DO 20 J = 1, N - KK
        DO 10 I = M - KK + 1, M
           A( I, J ) = ZERO
10       CONTINUE
20    CONTINUE
  ELSE
     KK = 0
  END IF
!
!     Use unblocked code for the first or only block.
!
  CALL DORG2L( M-KK, N-KK, K-KK, A, LDA, TAU, WORK, IINFO )
!
  IF( KK.GT.0 ) THEN
!
!        Use blocked code
!
     DO 50 I = K - KK + 1, K, NB
        IB = MIN( NB, K-I+1 )
        IF( N-K+I.GT.1 ) THEN
!
!              Form the triangular factor of the block reflector
!              H = H(i+ib-1) . . . H(i+1) H(i)
!
           CALL DLARFT( 'Backward', 'Columnwise', M-K+I+IB-1_i4, IB, &
                        A( 1, N-K+I ), LDA, TAU( I ), WORK, LDWORK )
!
!              Apply H to A(1:m-k+i+ib-1,1:n-k+i-1) from the left
!
           CALL DLARFB( 'Left', 'No transpose', 'Backward', &
                        'Columnwise', M-K+I+IB-1_i4, N-K+I-1_i4, IB, &
                        A( 1, N-K+I ), LDA, WORK, LDWORK, A, LDA, &
                        WORK( IB+1 ), LDWORK )
        END IF
!
!           Apply H to rows 1:m-k+i+ib-1 of current block
!
        CALL DORG2L( M-K+I+IB-1_i4, IB, IB, A( 1, N-K+I ), LDA, &
                     TAU( I ), WORK, IINFO )
!
!           Set rows m-k+i+ib:m of current block to zero
!
        DO 40 J = N - K + I, N - K + I + IB - 1
           DO 30 L = M - K + I + IB, M
              A( L, J ) = ZERO
30          CONTINUE
40       CONTINUE
50    CONTINUE
  END IF
!
  WORK( 1 ) = IWS
  RETURN
!
!     End of DORGQL
!
  END
  SUBROUTINE DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
  use datatypes
  integer(i4)        INFO, K, LDA, LWORK, M, N
!     ..
!     .. Array Arguments ..
  real(dp)           A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DORGQR generates an M-by-N real matrix Q with orthonormal columns,
!  which is defined as the first N columns of a product of K elementary
!  reflectors of order M
!
!        Q  =  H(1) H(2) . . . H(k)
!
!  as returned by DGEQRF.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix Q. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix Q. M >= N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines the
!          matrix Q. N >= K >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the i-th column must contain the vector which
!          defines the elementary reflector H(i), for i = 1,2,...,k, as
!          returned by DGEQRF in the first k columns of its array
!          argument A.
!          On exit, the M-by-N matrix Q.
!
!  LDA     (input) INTEGER
!          The first dimension of the array A. LDA >= max(1,M).
!
!  TAU     (input) DOUBLE PRECISION array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGEQRF.
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK. LWORK >= max(1,N).
!          For optimum performance LWORK >= N*NB, where NB is the
!          optimal blocksize.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument has an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
  real(dp)           ZERO
  PARAMETER          ( ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
  LOGICAL            LQUERY
  integer(i4)        I, IB, IINFO, IWS, J, KI, KK, L, LDWORK, &
                     LWKOPT, NB, NBMIN, NX
!     ..
!     .. External Subroutines ..
  EXTERNAL           DLARFB, DLARFT, DORG2R, XERBLA
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          MAX, MIN
!     ..
!     .. External Functions ..
  integer(i4)        ILAENV
  EXTERNAL           ILAENV
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
  INFO = 0
  NB = ILAENV( 1_i4, 'DORGQR', ' ', M, N, K, -1_i4 )
  LWKOPT = MAX( 1, N )*NB
  WORK( 1 ) = LWKOPT
  LQUERY = ( LWORK.EQ.-1 )
  IF( M.LT.0 ) THEN
     INFO = -1
  ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
     INFO = -2
  ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
     INFO = -3
  ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
     INFO = -5
  ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
     INFO = -8
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DORGQR', -INFO )
     RETURN
  ELSE IF( LQUERY ) THEN
     RETURN
  END IF
!
!     Quick return if possible
!
  IF( N.LE.0 ) THEN
     WORK( 1 ) = 1
     RETURN
  END IF
!
  NBMIN = 2
  NX = 0
  IWS = N
  IF( NB.GT.1 .AND. NB.LT.K ) THEN
!
!        Determine when to cross over from blocked to unblocked code.
!
     NX = MAX( 0_i4, ILAENV( 3_i4, 'DORGQR', ' ', M, N, K, -1_i4 ) )
     IF( NX.LT.K ) THEN
!
!           Determine if workspace is large enough for blocked code.
!
        LDWORK = N
        IWS = LDWORK*NB
        IF( LWORK.LT.IWS ) THEN
!
!              Not enough workspace to use optimal NB:  reduce NB and
!              determine the minimum value of NB.
!
           NB = LWORK / LDWORK
           NBMIN = MAX( 2_i4, ILAENV( 2_i4, 'DORGQR', ' ', M, N, K, -1_i4 ) )
        END IF
     END IF
  END IF
!
  IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
!
!        Use blocked code after the last block.
!        The first kk columns are handled by the block method.
!
     KI = ( ( K-NX-1 ) / NB )*NB
     KK = MIN( K, KI+NB )
!
!        Set A(1:kk,kk+1:n) to zero.
!
     DO 20 J = KK + 1, N
        DO 10 I = 1, KK
           A( I, J ) = ZERO
10       CONTINUE
20    CONTINUE
  ELSE
     KK = 0
  END IF
!
!     Use unblocked code for the last or only block.
!
  IF( KK.LT.N ) &
     CALL DORG2R( M-KK, N-KK, K-KK, A( KK+1, KK+1 ), LDA, &
                  TAU( KK+1 ), WORK, IINFO )
!
  IF( KK.GT.0 ) THEN
!
!        Use blocked code
!
     DO 50 I = KI + 1, 1, -NB
        IB = MIN( NB, K-I+1 )
        IF( I+IB.LE.N ) THEN
!
!              Form the triangular factor of the block reflector
!              H = H(i) H(i+1) . . . H(i+ib-1)
!
           CALL DLARFT( 'Forward', 'Columnwise', M-I+1_i4, IB, &
                        A( I, I ), LDA, TAU( I ), WORK, LDWORK )
!
!              Apply H to A(i:m,i+ib:n) from the left
!
           CALL DLARFB( 'Left', 'No transpose', 'Forward', &
                        'Columnwise', M-I+1_i4, N-I-IB+1_i4, IB, &
                        A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ), &
                        LDA, WORK( IB+1 ), LDWORK )
        END IF
!
!           Apply H to rows i:m of current block
!
        CALL DORG2R( M-I+1_i4, IB, IB, A( I, I ), LDA, TAU( I ), WORK, &
                     IINFO )
!
!           Set rows 1:i-1 of current block to zero
!
        DO 40 J = I, I + IB - 1
           DO 30 L = 1, I - 1
              A( L, J ) = ZERO
30          CONTINUE
40       CONTINUE
50    CONTINUE
  END IF
!
  WORK( 1 ) = IWS
  RETURN
!
!     End of DORGQR
!
  END
  SUBROUTINE DSYTD2( UPLO, N, A, LDA, D, E, TAU, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
  use datatypes
  CHARACTER          UPLO
  integer(i4)        INFO, LDA, N
!     ..
!     .. Array Arguments ..
  real(dp)           A( LDA, * ), D( * ), E( * ), TAU( * )
!     ..
!
!  Purpose
!  =======
!
!  DSYTD2 reduces a real symmetric matrix A to symmetric tridiagonal
!  form T by an orthogonal similarity transformation: Q' * A * Q = T.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the upper or lower triangular part of the
!          symmetric matrix A is stored:
!          = 'U':  Upper triangular
!          = 'L':  Lower triangular
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!          n-by-n upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading n-by-n lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!          On exit, if UPLO = 'U', the diagonal and first superdiagonal
!          of A are overwritten by the corresponding elements of the
!          tridiagonal matrix T, and the elements above the first
!          superdiagonal, with the array TAU, represent the orthogonal
!          matrix Q as a product of elementary reflectors; if UPLO
!          = 'L', the diagonal and first subdiagonal of A are over-
!          written by the corresponding elements of the tridiagonal
!          matrix T, and the elements below the first subdiagonal, with
!          the array TAU, represent the orthogonal matrix Q as a product
!          of elementary reflectors. See Further Details.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  D       (output) DOUBLE PRECISION array, dimension (N)
!          The diagonal elements of the tridiagonal matrix T:
!          D(i) = A(i,i).
!
!  E       (output) DOUBLE PRECISION array, dimension (N-1)
!          The off-diagonal elements of the tridiagonal matrix T:
!          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
!
!  TAU     (output) DOUBLE PRECISION array, dimension (N-1)
!          The scalar factors of the elementary reflectors (see Further
!          Details).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!
!  Further Details
!  ===============
!
!  If UPLO = 'U', the matrix Q is represented as a product of elementary
!  reflectors
!
!     Q = H(n-1) . . . H(2) H(1).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
!  A(1:i-1,i+1), and tau in TAU(i).
!
!  If UPLO = 'L', the matrix Q is represented as a product of elementary
!  reflectors
!
!     Q = H(1) H(2) . . . H(n-1).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
!  and tau in TAU(i).
!
!  The contents of A on exit are illustrated by the following examples
!  with n = 5:
!
!  if UPLO = 'U':                       if UPLO = 'L':
!
!    (  d   e   v2  v3  v4 )              (  d                  )
!    (      d   e   v3  v4 )              (  e   d              )
!    (          d   e   v4 )              (  v1  e   d          )
!    (              d   e  )              (  v1  v2  e   d      )
!    (                  d  )              (  v1  v2  v3  e   d  )
!
!  where d and e denote diagonal and off-diagonal elements of T, and vi
!  denotes an element of the vector defining H(i).
!
!  =====================================================================
!
!     .. Parameters ..
  real(dp)           ONE, ZERO, HALF
  PARAMETER          ( ONE = 1.0D0, ZERO = 0.0D0, &
                     HALF = 1.0D0 / 2.0D0 )
!     ..
!     .. Local Scalars ..
  LOGICAL            UPPER
  integer(i4)        I
  real(dp)           ALPHA, TAUI
!     ..
!     .. External Subroutines ..
  EXTERNAL           DAXPY, DLARFG, DSYMV, DSYR2, XERBLA
!     ..
!     .. External Functions ..
  LOGICAL            LSAME
  real(dp)           DDOT
  EXTERNAL           LSAME, DDOT
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
  INFO = 0
  UPPER = LSAME( UPLO, 'U' )
  IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
     INFO = -1
  ELSE IF( N.LT.0 ) THEN
     INFO = -2
  ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
     INFO = -4
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DSYTD2', -INFO )
     RETURN
  END IF
!
!     Quick return if possible
!
  IF( N.LE.0 ) &
     RETURN
!
  IF( UPPER ) THEN
!
!        Reduce the upper triangle of A
!
     DO 10 I = N - 1, 1, -1
!
!           Generate elementary reflector H(i) = I - tau * v * v'
!           to annihilate A(1:i-1,i+1)
!
        CALL DLARFG( I, A( I, I+1 ), A( 1, I+1 ), 1_i4, TAUI )
        E( I ) = A( I, I+1 )
!
        IF( TAUI.NE.ZERO ) THEN
!
!              Apply H(i) from both sides to A(1:i,1:i)
!
           A( I, I+1 ) = ONE
!
!              Compute  x := tau * A * v  storing x in TAU(1:i)
!
           CALL DSYMV( UPLO, I, TAUI, A, LDA, A( 1, I+1 ), 1_i4, ZERO, &
                       TAU, 1_i4 )
!
!              Compute  w := x - 1/2 * tau * (x'*v) * v
!
           ALPHA = -HALF*TAUI*DDOT( I, TAU, 1_i4, A( 1, I+1 ), 1_i4 )
           CALL DAXPY( I, ALPHA, A( 1, I+1 ), 1_i4, TAU, 1_i4 )
!
!              Apply the transformation as a rank-2 update:
!                 A := A - v * w' - w * v'
!
           CALL DSYR2( UPLO, I, -ONE, A( 1, I+1 ), 1_i4, TAU, 1_i4, A, &
                       LDA )
!
           A( I, I+1 ) = E( I )
        END IF
        D( I+1 ) = A( I+1, I+1 )
        TAU( I ) = TAUI
10    CONTINUE
     D( 1 ) = A( 1, 1 )
  ELSE
!
!        Reduce the lower triangle of A
!
     DO 20 I = 1, N - 1
!
!           Generate elementary reflector H(i) = I - tau * v * v'
!           to annihilate A(i+2:n,i)
!
        CALL DLARFG( N-I, A( I+1, I ), A( MIN( I+2, N ), I ), 1_i4, &
                     TAUI )
        E( I ) = A( I+1, I )
!
        IF( TAUI.NE.ZERO ) THEN
!
!              Apply H(i) from both sides to A(i+1:n,i+1:n)
!
           A( I+1, I ) = ONE
!
!              Compute  x := tau * A * v  storing y in TAU(i:n-1)
!
           CALL DSYMV( UPLO, N-I, TAUI, A( I+1, I+1 ), LDA, &
                       A( I+1, I ), 1_i4, ZERO, TAU( I ), 1_i4 )
!
!              Compute  w := x - 1/2 * tau * (x'*v) * v
!
           ALPHA = -HALF*TAUI*DDOT( N-I, TAU( I ), 1_i4, A( I+1, I ), &
                   1_i4 )
           CALL DAXPY( N-I, ALPHA, A( I+1, I ), 1_i4, TAU( I ), 1_i4 )
!
!              Apply the transformation as a rank-2 update:
!                 A := A - v * w' - w * v'
!
           CALL DSYR2( UPLO, N-I, -ONE, A( I+1, I ), 1_i4, TAU( I ), 1_i4, &
                       A( I+1, I+1 ), LDA )
!
           A( I+1, I ) = E( I )
        END IF
        D( I ) = A( I, I )
        TAU( I ) = TAUI
20    CONTINUE
     D( N ) = A( N, N )
  END IF
!
  RETURN
!
!     End of DSYTD2
!
  END
  FUNCTION DLANSY( NORM, UPLO, N, A, LDA, WORK )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
  use datatypes
  CHARACTER          NORM, UPLO
  integer(i4)        LDA, N
!     ..
!     .. Array Arguments ..
  real(dp)           A( LDA, * ), WORK( * )
  real(dp)           DLANSY
!     ..
!
!  Purpose
!  =======
!
!  DLANSY  returns the value of the one norm,  or the Frobenius norm, or
!  the  infinity norm,  or the  element of  largest absolute value  of a
!  real symmetric matrix A.
!
!  Description
!  ===========
!
!  DLANSY returns the value
!
!     DLANSY = ( max(abs(A(i,j))), NORM = 'M' or 'm'
!              (
!              ( norm1(A),         NORM = '1', 'O' or 'o'
!              (
!              ( normI(A),         NORM = 'I' or 'i'
!              (
!              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
!
!  where  norm1  denotes the  one norm of a matrix (maximum column sum),
!  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
!  normF  denotes the  Frobenius norm of a matrix (square root of sum of
!  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.
!
!  Arguments
!  =========
!
!  NORM    (input) CHARACTER*1
!          Specifies the value to be returned in DLANSY as described
!          above.
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the upper or lower triangular part of the
!          symmetric matrix A is to be referenced.
!          = 'U':  Upper triangular part of A is referenced
!          = 'L':  Lower triangular part of A is referenced
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.  When N = 0, DLANSY is
!          set to zero.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
!          The symmetric matrix A.  If UPLO = 'U', the leading n by n
!          upper triangular part of A contains the upper triangular part
!          of the matrix A, and the strictly lower triangular part of A
!          is not referenced.  If UPLO = 'L', the leading n by n lower
!          triangular part of A contains the lower triangular part of
!          the matrix A, and the strictly upper triangular part of A is
!          not referenced.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(N,1).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK),
!          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,
!          WORK is not referenced.
!
! =====================================================================
!
!     .. Parameters ..
  real(dp)           ONE, ZERO
  PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
  integer(i4)        I, J
  real(dp)           ABSA, SCALE, SUM, VALUE
!     ..
!     .. External Subroutines ..
  EXTERNAL           DLASSQ
!     ..
!     .. External Functions ..
  LOGICAL            LSAME
  EXTERNAL           LSAME
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          ABS, MAX, SQRT
!     ..
!     .. Executable Statements ..
!
  IF( N.EQ.0 ) THEN
     VALUE = ZERO
  ELSE IF( LSAME( NORM, 'M' ) ) THEN
!
!        Find max(abs(A(i,j))).
!
     VALUE = ZERO
     IF( LSAME( UPLO, 'U' ) ) THEN
        DO 20 J = 1, N
           DO 10 I = 1, J
              VALUE = MAX( VALUE, ABS( A( I, J ) ) )
10          CONTINUE
20       CONTINUE
     ELSE
        DO 40 J = 1, N
           DO 30 I = J, N
              VALUE = MAX( VALUE, ABS( A( I, J ) ) )
30          CONTINUE
40       CONTINUE
     END IF
  ELSE IF( ( LSAME( NORM, 'I' ) ) .OR. ( LSAME( NORM, 'O' ) ) .OR. &
           ( NORM.EQ.'1' ) ) THEN
!
!        Find normI(A) ( = norm1(A), since A is symmetric).
!
     VALUE = ZERO
     IF( LSAME( UPLO, 'U' ) ) THEN
        DO 60 J = 1, N
           SUM = ZERO
           DO 50 I = 1, J - 1
              ABSA = ABS( A( I, J ) )
              SUM = SUM + ABSA
              WORK( I ) = WORK( I ) + ABSA
50          CONTINUE
           WORK( J ) = SUM + ABS( A( J, J ) )
60       CONTINUE
        DO 70 I = 1, N
           VALUE = MAX( VALUE, WORK( I ) )
70       CONTINUE
     ELSE
        DO 80 I = 1, N
           WORK( I ) = ZERO
80       CONTINUE
        DO 100 J = 1, N
           SUM = WORK( J ) + ABS( A( J, J ) )
           DO 90 I = J + 1, N
              ABSA = ABS( A( I, J ) )
              SUM = SUM + ABSA
              WORK( I ) = WORK( I ) + ABSA
90          CONTINUE
           VALUE = MAX( VALUE, SUM )
100       CONTINUE
     END IF
  ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
!
!        Find normF(A).
!
     SCALE = ZERO
     SUM = ONE
     IF( LSAME( UPLO, 'U' ) ) THEN
        DO 110 J = 2, N
           CALL DLASSQ( J-1_i4, A( 1, J ), 1_i4, SCALE, SUM )
110       CONTINUE
     ELSE
        DO 120 J = 1, N - 1
           CALL DLASSQ( N-J, A( J+1, J ), 1_i4, SCALE, SUM )
120       CONTINUE
     END IF
     SUM = 2*SUM
     CALL DLASSQ( N, A, LDA+1_i4, SCALE, SUM )
     VALUE = SCALE*SQRT( SUM )
  END IF
!
  DLANSY = VALUE
  RETURN
!
!     End of DLANSY
!
  END
  SUBROUTINE DORGTR( UPLO, N, A, LDA, TAU, WORK, LWORK, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
  use datatypes
  CHARACTER          UPLO
  integer(i4)        INFO, LDA, LWORK, N
!     ..
!     .. Array Arguments ..
  real(dp)           A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DORGTR generates a real orthogonal matrix Q which is defined as the
!  product of n-1 elementary reflectors of order N, as returned by
!  DSYTRD:
!
!  if UPLO = 'U', Q = H(n-1) . . . H(2) H(1),
!
!  if UPLO = 'L', Q = H(1) H(2) . . . H(n-1).
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U': Upper triangle of A contains elementary reflectors
!                 from DSYTRD;
!          = 'L': Lower triangle of A contains elementary reflectors
!                 from DSYTRD.
!
!  N       (input) INTEGER
!          The order of the matrix Q. N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the vectors which define the elementary reflectors,
!          as returned by DSYTRD.
!          On exit, the N-by-N orthogonal matrix Q.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A. LDA >= max(1,N).
!
!  TAU     (input) DOUBLE PRECISION array, dimension (N-1)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DSYTRD.
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK. LWORK >= max(1,N-1).
!          For optimum performance LWORK >= (N-1)*NB, where NB is
!          the optimal blocksize.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
  real(dp)           ZERO, ONE
  PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
  LOGICAL            LQUERY, UPPER
  integer(i4)        I, IINFO, J, LWKOPT, NB
!     ..
!     .. External Functions ..
  LOGICAL            LSAME
  integer(i4)        ILAENV
  EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
  EXTERNAL           DORGQL, DORGQR, XERBLA
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
  INFO = 0
  LQUERY = ( LWORK.EQ.-1 )
  UPPER = LSAME( UPLO, 'U' )
  IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
     INFO = -1
  ELSE IF( N.LT.0 ) THEN
     INFO = -2
  ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
     INFO = -4
  ELSE IF( LWORK.LT.MAX( 1, N-1 ) .AND. .NOT.LQUERY ) THEN
     INFO = -7
  END IF
!
  IF( INFO.EQ.0 ) THEN
     IF( UPPER ) THEN
        NB = ILAENV( 1_i4, 'DORGQL', ' ', N-1_i4, N-1_i4, N-1_i4, -1_i4 )
     ELSE
        NB = ILAENV( 1_i4, 'DORGQR', ' ', N-1_i4, N-1_i4, N-1_i4, -1_i4 )
     END IF
     LWKOPT = MAX( 1, N-1 )*NB
     WORK( 1 ) = LWKOPT
  END IF
!
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DORGTR', -INFO )
     RETURN
  ELSE IF( LQUERY ) THEN
     RETURN
  END IF
!
!     Quick return if possible
!
  IF( N.EQ.0 ) THEN
     WORK( 1 ) = 1
     RETURN
  END IF
!
  IF( UPPER ) THEN
!
!        Q was determined by a call to DSYTRD with UPLO = 'U'
!
!        Shift the vectors which define the elementary reflectors one
!        column to the left, and set the last row and column of Q to
!        those of the unit matrix
!
     DO 20 J = 1, N - 1
        DO 10 I = 1, J - 1
           A( I, J ) = A( I, J+1 )
10       CONTINUE
        A( N, J ) = ZERO
20    CONTINUE
     DO 30 I = 1, N - 1
        A( I, N ) = ZERO
30    CONTINUE
     A( N, N ) = ONE
!
!        Generate Q(1:n-1,1:n-1)
!
     CALL DORGQL( N-1_i4, N-1_i4, N-1_i4, A, LDA, TAU, WORK, LWORK, IINFO )
!
  ELSE
!
!        Q was determined by a call to DSYTRD with UPLO = 'L'.
!
!        Shift the vectors which define the elementary reflectors one
!        column to the right, and set the first row and column of Q to
!        those of the unit matrix
!
     DO 50 J = N, 2, -1
        A( 1, J ) = ZERO
        DO 40 I = J + 1, N
           A( I, J ) = A( I, J-1 )
40       CONTINUE
50    CONTINUE
     A( 1, 1 ) = ONE
     DO 60 I = 2, N
        A( I, 1 ) = ZERO
60    CONTINUE
     IF( N.GT.1 ) THEN
!
!           Generate Q(2:n,2:n)
!
        CALL DORGQR( N-1_i4, N-1_i4, N-1_i4, A( 2, 2 ), LDA, TAU, WORK, &
                     LWORK, IINFO )
     END IF
  END IF
  WORK( 1 ) = LWKOPT
  RETURN
!
!     End of DORGTR
!
  END
  SUBROUTINE DSYTRD( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
  use datatypes
  CHARACTER          UPLO
  integer(i4)        INFO, LDA, LWORK, N
!     ..
!     .. Array Arguments ..
  real(dp)           A( LDA, * ), D( * ), E( * ), TAU( * ), &
                     WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DSYTRD reduces a real symmetric matrix A to real symmetric
!  tridiagonal form T by an orthogonal similarity transformation:
!  Q**T * A * Q = T.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!          N-by-N upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading N-by-N lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!          On exit, if UPLO = 'U', the diagonal and first superdiagonal
!          of A are overwritten by the corresponding elements of the
!          tridiagonal matrix T, and the elements above the first
!          superdiagonal, with the array TAU, represent the orthogonal
!          matrix Q as a product of elementary reflectors; if UPLO
!          = 'L', the diagonal and first subdiagonal of A are over-
!          written by the corresponding elements of the tridiagonal
!          matrix T, and the elements below the first subdiagonal, with
!          the array TAU, represent the orthogonal matrix Q as a product
!          of elementary reflectors. See Further Details.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  D       (output) DOUBLE PRECISION array, dimension (N)
!          The diagonal elements of the tridiagonal matrix T:
!          D(i) = A(i,i).
!
!  E       (output) DOUBLE PRECISION array, dimension (N-1)
!          The off-diagonal elements of the tridiagonal matrix T:
!          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
!
!  TAU     (output) DOUBLE PRECISION array, dimension (N-1)
!          The scalar factors of the elementary reflectors (see Further
!          Details).
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.  LWORK >= 1.
!          For optimum performance LWORK >= N*NB, where NB is the
!          optimal blocksize.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  Further Details
!  ===============
!
!  If UPLO = 'U', the matrix Q is represented as a product of elementary
!  reflectors
!
!     Q = H(n-1) . . . H(2) H(1).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
!  A(1:i-1,i+1), and tau in TAU(i).
!
!  If UPLO = 'L', the matrix Q is represented as a product of elementary
!  reflectors
!
!     Q = H(1) H(2) . . . H(n-1).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
!  and tau in TAU(i).
!
!  The contents of A on exit are illustrated by the following examples
!  with n = 5:
!
!  if UPLO = 'U':                       if UPLO = 'L':
!
!    (  d   e   v2  v3  v4 )              (  d                  )
!    (      d   e   v3  v4 )              (  e   d              )
!    (          d   e   v4 )              (  v1  e   d          )
!    (              d   e  )              (  v1  v2  e   d      )
!    (                  d  )              (  v1  v2  v3  e   d  )
!
!  where d and e denote diagonal and off-diagonal elements of T, and vi
!  denotes an element of the vector defining H(i).
!
!  =====================================================================
!
!     .. Parameters ..
  real(dp)           ONE
  PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
  LOGICAL            LQUERY, UPPER
  integer(i4)        I, IINFO, IWS, J, KK, LDWORK, LWKOPT, NB, &
                     NBMIN, NX
!     ..
!     .. External Subroutines ..
  EXTERNAL           DLATRD, DSYR2K, DSYTD2, XERBLA
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          MAX
!     ..
!     .. External Functions ..
  LOGICAL            LSAME
  integer(i4)        ILAENV
  EXTERNAL           LSAME, ILAENV
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
  INFO = 0
  UPPER = LSAME( UPLO, 'U' )
  LQUERY = ( LWORK.EQ.-1 )
  IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
     INFO = -1
  ELSE IF( N.LT.0 ) THEN
     INFO = -2
  ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
     INFO = -4
  ELSE IF( LWORK.LT.1 .AND. .NOT.LQUERY ) THEN
     INFO = -9
  END IF
!
  IF( INFO.EQ.0 ) THEN
!
!        Determine the block size.
!
     NB = ILAENV( 1_i4, 'DSYTRD', UPLO, N, -1_i4, -1_i4, -1_i4 )
     LWKOPT = N*NB
     WORK( 1 ) = LWKOPT
  END IF
!
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DSYTRD', -INFO )
     RETURN
  ELSE IF( LQUERY ) THEN
     RETURN
  END IF
!
!     Quick return if possible
!
  IF( N.EQ.0 ) THEN
     WORK( 1 ) = 1
     RETURN
  END IF
!
  NX = N
  IWS = 1
  IF( NB.GT.1 .AND. NB.LT.N ) THEN
!
!        Determine when to cross over from blocked to unblocked code
!        (last block is always handled by unblocked code).
!
     NX = MAX( NB, ILAENV( 3_i4, 'DSYTRD', UPLO, N, -1_i4, -1_i4, -1_i4 ) )
     IF( NX.LT.N ) THEN
!
!           Determine if workspace is large enough for blocked code.
!
        LDWORK = N
        IWS = LDWORK*NB
        IF( LWORK.LT.IWS ) THEN
!
!              Not enough workspace to use optimal NB:  determine the
!              minimum value of NB, and reduce NB or force use of
!              unblocked code by setting NX = N.
!
           NB = MAX( LWORK / LDWORK, 1 )
           NBMIN = ILAENV( 2_i4, 'DSYTRD', UPLO, N, -1_i4, -1_i4, -1_i4 )
           IF( NB.LT.NBMIN ) &
              NX = N
        END IF
     ELSE
        NX = N
     END IF
  ELSE
     NB = 1
  END IF
!
  IF( UPPER ) THEN
!
!        Reduce the upper triangle of A.
!        Columns 1:kk are handled by the unblocked method.
!
     KK = N - ( ( N-NX+NB-1 ) / NB )*NB
     DO 20 I = N - NB + 1, KK + 1, -NB
!
!           Reduce columns i:i+nb-1 to tridiagonal form and form the
!           matrix W which is needed to update the unreduced part of
!           the matrix
!
        CALL DLATRD( UPLO, I+NB-1_i4, NB, A, LDA, E, TAU, WORK, &
                     LDWORK )
!
!           Update the unreduced submatrix A(1:i-1,1:i-1), using an
!           update of the form:  A := A - V*W' - W*V'
!
        CALL DSYR2K( UPLO, 'No transpose', I-1_i4, NB, -ONE, A( 1, I ), &
                     LDA, WORK, LDWORK, ONE, A, LDA )
!
!           Copy superdiagonal elements back into A, and diagonal
!           elements into D
!
        DO 10 J = I, I + NB - 1
           A( J-1, J ) = E( J-1 )
           D( J ) = A( J, J )
10       CONTINUE
20    CONTINUE
!
!        Use unblocked code to reduce the last or only block
!
     CALL DSYTD2( UPLO, KK, A, LDA, D, E, TAU, IINFO )
  ELSE
!
!        Reduce the lower triangle of A
!
     DO 40 I = 1, N - NX, NB
!
!           Reduce columns i:i+nb-1 to tridiagonal form and form the
!           matrix W which is needed to update the unreduced part of
!           the matrix
!
        CALL DLATRD( UPLO, N-I+1_i4, NB, A( I, I ), LDA, E( I ), &
                     TAU( I ), WORK, LDWORK )
!
!           Update the unreduced submatrix A(i+ib:n,i+ib:n), using
!           an update of the form:  A := A - V*W' - W*V'
!
        CALL DSYR2K( UPLO, 'No transpose', N-I-NB+1_i4, NB, -ONE, &
                     A( I+NB, I ), LDA, WORK( NB+1 ), LDWORK, ONE, &
                     A( I+NB, I+NB ), LDA )
!
!           Copy subdiagonal elements back into A, and diagonal
!           elements into D
!
        DO 30 J = I, I + NB - 1
           A( J+1, J ) = E( J )
           D( J ) = A( J, J )
30       CONTINUE
40    CONTINUE
!
!        Use unblocked code to reduce the last or only block
!
     CALL DSYTD2( UPLO, N-I+1_i4, A( I, I ), LDA, D( I ), E( I ), &
                  TAU( I ), IINFO )
  END IF
!
  WORK( 1 ) = LWKOPT
  RETURN
!
!     End of DSYTRD
!
  END
  SUBROUTINE DPOSV( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
!
!  -- LAPACK driver routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
  use datatypes
  CHARACTER          UPLO
  integer(i4)        INFO, LDA, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
  real(dp)           A( LDA, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  DPOSV computes the solution to a real system of linear equations
!     A * X = B,
!  where A is an N-by-N symmetric positive definite matrix and X and B
!  are N-by-NRHS matrices.
!
!  The Cholesky decomposition is used to factor A as
!     A = U**T* U,  if UPLO = 'U', or
!     A = L * L**T,  if UPLO = 'L',
!  where U is an upper triangular matrix and L is a lower triangular
!  matrix.  The factored form of A is then used to solve the system of
!  equations A * X = B.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The number of linear equations, i.e., the order of the
!          matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!          N-by-N upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading N-by-N lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!
!          On exit, if INFO = 0, the factor U or L from the Cholesky
!          factorization A = U**T*U or A = L*L**T.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!          On entry, the N-by-NRHS right hand side matrix B.
!          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, the leading minor of order i of A is not
!                positive definite, so the factorization could not be
!                completed, and the solution has not been computed.
!
!  =====================================================================
!
!     .. External Functions ..
  LOGICAL            LSAME
  EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
  EXTERNAL           DPOTRF, DPOTRS, XERBLA
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
  INFO = 0
  IF( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
     INFO = -1
  ELSE IF( N.LT.0 ) THEN
     INFO = -2
  ELSE IF( NRHS.LT.0 ) THEN
     INFO = -3
  ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
     INFO = -5
  ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
     INFO = -7
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DPOSV ', -INFO )
     RETURN
  END IF
!
!     Compute the Cholesky factorization A = U'*U or A = L*L'.
!
  CALL DPOTRF( UPLO, N, A, LDA, INFO )
  IF( INFO.EQ.0 ) THEN
!
!        Solve the system A*X = B, overwriting B with X.
!
     CALL DPOTRS( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
!
  END IF
  RETURN
!
!     End of DPOSV
!
  END
  SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
  use datatypes
  CHARACTER          UPLO
  integer(i4)        INFO, LDA, N
!     ..
!     .. Array Arguments ..
  real(dp)           A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DPOTRF computes the Cholesky factorization of a real symmetric
!  positive definite matrix A.
!
!  The factorization has the form
!     A = U**T * U,  if UPLO = 'U', or
!     A = L  * L**T,  if UPLO = 'L',
!  where U is an upper triangular matrix and L is lower triangular.
!
!  This is the block version of the algorithm, calling Level 3 BLAS.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!          N-by-N upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading N-by-N lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!
!          On exit, if INFO = 0, the factor U or L from the Cholesky
!          factorization A = U**T*U or A = L*L**T.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, the leading minor of order i is not
!                positive definite, and the factorization could not be
!                completed.
!
!  =====================================================================
!
!     .. Parameters ..
  real(dp)           ONE
  PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
  LOGICAL            UPPER
  integer(i4)        J, JB, NB
!     ..
!     .. External Functions ..
  LOGICAL            LSAME
  integer(i4)        ILAENV
  EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
  EXTERNAL           DGEMM, DPOTF2, DSYRK, DTRSM, XERBLA
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
  INFO = 0
  UPPER = LSAME( UPLO, 'U' )
  IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
     INFO = -1
  ELSE IF( N.LT.0 ) THEN
     INFO = -2
  ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
     INFO = -4
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DPOTRF', -INFO )
     RETURN
  END IF
!
!     Quick return if possible
!
  IF( N.EQ.0 ) &
     RETURN
!
!     Determine the block size for this environment.
!
  NB = ILAENV( 1_i4, 'DPOTRF', UPLO, N, -1_i4, -1_i4, -1_i4 )
  IF( NB.LE.1 .OR. NB.GE.N ) THEN
!
!        Use unblocked code.
!
     CALL DPOTF2( UPLO, N, A, LDA, INFO )
  ELSE
!
!        Use blocked code.
!
     IF( UPPER ) THEN
!
!           Compute the Cholesky factorization A = U'*U.
!
        DO 10 J = 1, N, NB
!
!              Update and factorize the current diagonal block and test
!              for non-positive-definiteness.
!
           JB = MIN( NB, N-J+1 )
           CALL DSYRK( 'Upper', 'Transpose', JB, J-1_i4, -ONE, &
                       A( 1, J ), LDA, ONE, A( J, J ), LDA )
           CALL DPOTF2( 'Upper', JB, A( J, J ), LDA, INFO )
           IF( INFO.NE.0 ) &
              GO TO 30
           IF( J+JB.LE.N ) THEN
!
!                 Compute the current block row.
!
              CALL DGEMM( 'Transpose', 'No transpose', JB, N-J-JB+1_i4, &
                          J-1_i4, -ONE, A( 1, J ), LDA, A( 1, J+JB ), &
                          LDA, ONE, A( J, J+JB ), LDA )
              CALL DTRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', &
                          JB, N-J-JB+1_i4, ONE, A( J, J ), LDA, &
                          A( J, J+JB ), LDA )
           END IF
10       CONTINUE
!
     ELSE
!
!           Compute the Cholesky factorization A = L*L'.
!
        DO 20 J = 1, N, NB
!
!              Update and factorize the current diagonal block and test
!              for non-positive-definiteness.
!
           JB = MIN( NB, N-J+1_i4 )
           CALL DSYRK( 'Lower', 'No transpose', JB, J-1_i4, -ONE, &
                       A( J, 1 ), LDA, ONE, A( J, J ), LDA )
           CALL DPOTF2( 'Lower', JB, A( J, J ), LDA, INFO )
           IF( INFO.NE.0 ) &
              GO TO 30
           IF( J+JB.LE.N ) THEN
!
!                 Compute the current block column.
!
              CALL DGEMM( 'No transpose', 'Transpose', N-J-JB+1_i4, JB, &
                          J-1_i4, -ONE, A( J+JB, 1 ), LDA, A( J, 1 ), &
                          LDA, ONE, A( J+JB, J ), LDA )
              CALL DTRSM( 'Right', 'Lower', 'Transpose', 'Non-unit', &
                          N-J-JB+1_i4, JB, ONE, A( J, J ), LDA, &
                          A( J+JB, J ), LDA )
           END IF
20       CONTINUE
     END IF
  END IF
  GO TO 40
!
30 CONTINUE
  INFO = INFO + J - 1
!
40 CONTINUE
  RETURN
!
!     End of DPOTRF
!
  END
  SUBROUTINE DPOTRS( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
  use datatypes
  CHARACTER          UPLO
  integer(i4)        INFO, LDA, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
  real(dp)           A( LDA, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  DPOTRS solves a system of linear equations A*X = B with a symmetric
!  positive definite matrix A using the Cholesky factorization
!  A = U**T*U or A = L*L**T computed by DPOTRF.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
!          The triangular factor U or L from the Cholesky factorization
!          A = U**T*U or A = L*L**T, as computed by DPOTRF.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!          On entry, the right hand side matrix B.
!          On exit, the solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
  real(dp)           ONE
  PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
  LOGICAL            UPPER
!     ..
!     .. External Functions ..
  LOGICAL            LSAME
  EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
  EXTERNAL           DTRSM, XERBLA
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
  INFO = 0
  UPPER = LSAME( UPLO, 'U' )
  IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
     INFO = -1
  ELSE IF( N.LT.0 ) THEN
     INFO = -2
  ELSE IF( NRHS.LT.0 ) THEN
     INFO = -3
  ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
     INFO = -5
  ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
     INFO = -7
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DPOTRS', -INFO )
     RETURN
  END IF
!
!     Quick return if possible
!
  IF( N.EQ.0 .OR. NRHS.EQ.0 ) &
     RETURN
!
  IF( UPPER ) THEN
!
!        Solve A*X = B where A = U'*U.
!
!        Solve U'*X = B, overwriting B with X.
!
     CALL DTRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', N, NRHS, &
                 ONE, A, LDA, B, LDB )
!
!        Solve U*X = B, overwriting B with X.
!
     CALL DTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N, &
                 NRHS, ONE, A, LDA, B, LDB )
  ELSE
!
!        Solve A*X = B where A = L*L'.
!
!        Solve L*X = B, overwriting B with X.
!
     CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Non-unit', N, &
                 NRHS, ONE, A, LDA, B, LDB )
!
!        Solve L'*X = B, overwriting B with X.
!
     CALL DTRSM( 'Left', 'Lower', 'Transpose', 'Non-unit', N, NRHS, &
                 ONE, A, LDA, B, LDB )
  END IF
!
  RETURN
!
!     End of DPOTRS
!
  END
  SUBROUTINE DPOTF2( UPLO, N, A, LDA, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
  use datatypes
  CHARACTER          UPLO
  integer(i4)        INFO, LDA, N
!     ..
!     .. Array Arguments ..
  real(dp)           A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DPOTF2 computes the Cholesky factorization of a real symmetric
!  positive definite matrix A.
!
!  The factorization has the form
!     A = U' * U ,  if UPLO = 'U', or
!     A = L  * L',  if UPLO = 'L',
!  where U is an upper triangular matrix and L is lower triangular.
!
!  This is the unblocked version of the algorithm, calling Level 2 BLAS.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the upper or lower triangular part of the
!          symmetric matrix A is stored.
!          = 'U':  Upper triangular
!          = 'L':  Lower triangular
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!          n by n upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading n by n lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!
!          On exit, if INFO = 0, the factor U or L from the Cholesky
!          factorization A = U'*U  or A = L*L'.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, the leading minor of order k is not
!               positive definite, and the factorization could not be
!               completed.
!
!  =====================================================================
!
!     .. Parameters ..
  real(dp)           ONE, ZERO
  PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
  LOGICAL            UPPER
  integer(i4)        J
  real(dp)           AJJ
!     ..
!     .. External Functions ..
  LOGICAL            LSAME
  real(dp)           DDOT
  EXTERNAL           LSAME, DDOT
!     ..
!     .. External Subroutines ..
  EXTERNAL           DGEMV, DSCAL, XERBLA
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          MAX, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
  INFO = 0
  UPPER = LSAME( UPLO, 'U' )
  IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
     INFO = -1
  ELSE IF( N.LT.0 ) THEN
     INFO = -2
  ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
     INFO = -4
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DPOTF2', -INFO )
     RETURN
  END IF
!
!     Quick return if possible
!
  IF( N.EQ.0 ) &
     RETURN
!
  IF( UPPER ) THEN
!
!        Compute the Cholesky factorization A = U'*U.
!
     DO 10 J = 1, N
!
!           Compute U(J,J) and test for non-positive-definiteness.
!
        AJJ = A( J, J ) - DDOT( J-1_i4, A( 1, J ), 1_i4, A( 1, J ), 1_i4 )
        IF( AJJ.LE.ZERO ) THEN
           A( J, J ) = AJJ
           GO TO 30
        END IF
        AJJ = SQRT( AJJ )
        A( J, J ) = AJJ
!
!           Compute elements J+1:N of row J.
!
        IF( J.LT.N ) THEN
           CALL DGEMV( 'Transpose', J-1_i4, N-J, -ONE, A( 1, J+1 ), &
                       LDA, A( 1, J ), 1_i4, ONE, A( J, J+1 ), LDA )
           CALL DSCAL( N-J, ONE / AJJ, A( J, J+1 ), LDA )
        END IF
10    CONTINUE
  ELSE
!
!        Compute the Cholesky factorization A = L*L'.
!
     DO 20 J = 1, N
!
!           Compute L(J,J) and test for non-positive-definiteness.
!
        AJJ = A( J, J ) - DDOT( J-1_i4, A( J, 1 ), LDA, A( J, 1 ), &
              LDA )
        IF( AJJ.LE.ZERO ) THEN
           A( J, J ) = AJJ
           GO TO 30
        END IF
        AJJ = SQRT( AJJ )
        A( J, J ) = AJJ
!
!           Compute elements J+1:N of column J.
!
        IF( J.LT.N ) THEN
           CALL DGEMV( 'No transpose', N-J, J-1_i4, -ONE, A( J+1, 1 ), &
                       LDA, A( J, 1 ), LDA, ONE, A( J+1, J ), 1_i4 )
           CALL DSCAL( N-J, ONE / AJJ, A( J+1, J ), 1_i4 )
        END IF
20    CONTINUE
  END IF
  GO TO 40
!
30 CONTINUE
  INFO = J
!
40 CONTINUE
  RETURN
!
!     End of DPOTF2
!
  END
  SUBROUTINE DSPTRS( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
  use datatypes
  CHARACTER          :: UPLO
  integer(i4)        :: INFO, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
  integer(i4)        :: IPIV( * )
  real(dp)           :: AP( * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  DSPTRS solves a system of linear equations A*X = B with a real
!  symmetric matrix A stored in packed format using the factorization
!  A = U*D*U**T or A = L*D*L**T computed by DSPTRF.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the details of the factorization are stored
!          as an upper or lower triangular matrix.
!          = 'U':  Upper triangular, form is A = U*D*U**T;
!          = 'L':  Lower triangular, form is A = L*D*L**T.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  AP      (input) DOUBLE PRECISION array, dimension (N*(N+1)/2)
!          The block diagonal matrix D and the multipliers used to
!          obtain the factor U or L as computed by DSPTRF, stored as a
!          packed triangular matrix.
!
!  IPIV    (input) INTEGER array, dimension (N)
!          Details of the interchanges and the block structure of D
!          as determined by DSPTRF.
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!          On entry, the right hand side matrix B.
!          On exit, the solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
  real(dp), parameter ::  ONE = 1.0_dp
!     ..
!     .. Local Scalars ..
  logical            :: UPPER
  integer(i4)        :: J, K, KC, KP
  real(dp)           :: AK, AKM1, AKM1K, BK, BKM1, DENOM
!     ..
!     .. External Functions ..
  LOGICAL            :: LSAME
  EXTERNAL           :: LSAME
!     ..
!     .. External Subroutines ..
  EXTERNAL           :: DGEMV, DGER, DSCAL, DSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          :: MAX
!     ..
!     .. Executable Statements ..
!
  INFO = 0
  UPPER = LSAME( UPLO, 'U' )
  IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
     INFO = -1
  ELSE IF( N.LT.0 ) THEN
     INFO = -2
  ELSE IF( NRHS.LT.0 ) THEN
     INFO = -3
  ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
     INFO = -7
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DSPTRS', -INFO )
     RETURN
  END IF
!
!     Quick return if possible
!
  IF( N.EQ.0 .OR. NRHS.EQ.0 ) RETURN
!
  IF( UPPER ) THEN
!
!        Solve A*X = B, where A = U*D*U'.
!
!        First solve U*D*X = B, overwriting B with X.
!
!        K is the main loop index, decreasing from N to 1 in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
    K = N
    KC = N*( N+1 ) / 2 + 1
   10    CONTINUE
!
!        If K < 1, exit from loop.
!
    IF( K.LT.1 ) GO TO 30
!
    KC = KC - K
    IF( IPIV( K ).GT.0 ) THEN
!
!           1 x 1 diagonal block
!
!           Interchange rows K and IPIV(K).
!
      KP = IPIV( K )
      IF ( KP.NE.K ) CALL DSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
!
!           Multiply by inv(U(K)), where U(K) is the transformation
!           stored in column K of A.
!
      CALL DGER( K-1_i4, NRHS, -ONE, AP( KC ), 1_i4, B( K, 1 ), LDB, B( 1, 1 ), LDB )
!
!           Multiply by the inverse of the diagonal block.
!
      CALL DSCAL( NRHS, ONE / AP( KC+K-1 ), B( K, 1 ), LDB )
      K = K - 1
    ELSE
!
!           2 x 2 diagonal block
!
!           Interchange rows K-1 and -IPIV(K).
!
      KP = -IPIV( K )
      IF( KP.NE.K-1 ) CALL DSWAP( NRHS, B( K-1, 1 ), LDB, B( KP, 1 ), LDB )
!
!           Multiply by inv(U(K)), where U(K) is the transformation
!           stored in columns K-1 and K of A.
!
      CALL DGER( K-2_i4, NRHS, -ONE, AP( KC ), 1_i4, B( K, 1 ), LDB, B( 1, 1 ), LDB )
      CALL DGER( K-2_i4, NRHS, -ONE, AP( KC-( K-1 ) ), 1_i4, B( K-1, 1 ), LDB, B( 1, 1 ), LDB )
!
!           Multiply by the inverse of the diagonal block.
!
      AKM1K = AP( KC+K-2 )
      AKM1 = AP( KC-1 ) / AKM1K
      AK = AP( KC+K-1 ) / AKM1K
      DENOM = AKM1*AK - ONE
      DO 20 J = 1, NRHS
         BKM1 = B( K-1, J ) / AKM1K
         BK = B( K, J ) / AKM1K
         B( K-1, J ) = ( AK*BKM1-BK ) / DENOM
         B( K, J ) = ( AKM1*BK-BKM1 ) / DENOM
   20 CONTINUE
      KC = KC - K + 1
      K = K - 2
    END IF
!
    GO TO 10
   30    CONTINUE
!
!        Next solve U'*X = B, overwriting B with X.
!
!        K is the main loop index, increasing from 1 to N in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
    K = 1
    KC = 1
   40    CONTINUE
!
!        If K > N, exit from loop.
!
    IF( K.GT.N ) GO TO 50
!
    IF( IPIV( K ).GT.0 ) THEN
!
!           1 x 1 diagonal block
!
!           Multiply by inv(U'(K)), where U(K) is the transformation
!           stored in column K of A.
!
    CALL DGEMV( 'Transpose', K-1_i4, NRHS, -ONE, B, LDB, AP( KC ), 1_i4, ONE, B( K, 1 ), LDB )
!
!           Interchange rows K and IPIV(K).
!
    KP = IPIV( K )
    IF( KP.NE.K ) CALL DSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
    KC = KC + K
    K = K + 1
  ELSE
!
!           2 x 2 diagonal block
!
!           Multiply by inv(U'(K+1)), where U(K+1) is the transformation
!           stored in columns K and K+1 of A.
!
    CALL DGEMV( 'Transpose', K-1_i4, NRHS, -ONE, B, LDB, AP( KC ), 1_i4, ONE, B( K, 1 ), LDB )
    CALL DGEMV( 'Transpose', K-1_i4, NRHS, -ONE, B, LDB, AP( KC+K ), 1_i4, ONE, B( K+1, 1 ), LDB )
!
!           Interchange rows K and -IPIV(K).
!
    KP = -IPIV( K )
    IF( KP.NE.K ) CALL DSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
    KC = KC + 2*K + 1
    K = K + 2
  END IF
!
  GO TO 40
   50    CONTINUE
!
  ELSE
!
!        Solve A*X = B, where A = L*D*L'.
!
!        First solve L*D*X = B, overwriting B with X.
!
!        K is the main loop index, increasing from 1 to N in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
    K = 1
    KC = 1
   60    CONTINUE
!
!        If K > N, exit from loop.
!
    IF( K.GT.N ) GO TO 80
!
    IF( IPIV( K ).GT.0 ) THEN
!
!           1 x 1 diagonal block
!
!           Interchange rows K and IPIV(K).
!
      KP = IPIV( K )
      IF( KP.NE.K ) CALL DSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
!
!           Multiply by inv(L(K)), where L(K) is the transformation
!           stored in column K of A.
!
      IF( K.LT.N ) CALL DGER( N-K, NRHS, -ONE, AP( KC+1 ), 1_i4, B( K, 1 ), LDB, B( K+1, 1 ), LDB )
!
!           Multiply by the inverse of the diagonal block.
!
      CALL DSCAL( NRHS, ONE / AP( KC ), B( K, 1 ), LDB )
      KC = KC + N - K + 1
      K = K + 1
    ELSE
!
!           2 x 2 diagonal block
!
!           Interchange rows K+1 and -IPIV(K).
!
      KP = -IPIV( K )
      IF( KP.NE.K+1 ) CALL DSWAP( NRHS, B( K+1, 1 ), LDB, B( KP, 1 ), LDB )
!
!           Multiply by inv(L(K)), where L(K) is the transformation
!           stored in columns K and K+1 of A.
!
      IF( K.LT.N-1 ) THEN
        CALL DGER( N-K-1_i4, NRHS, -ONE, AP( KC+2 ), 1_i4, B( K, 1 ), LDB, B( K+2, 1 ), LDB )
        CALL DGER( N-K-1_i4, NRHS, -ONE, AP( KC+N-K+2 ), 1_i4, B( K+1, 1 ), LDB, B( K+2, 1 ), LDB )
      END IF
!
!           Multiply by the inverse of the diagonal block.
!
      AKM1K = AP( KC+1 )
      AKM1 = AP( KC ) / AKM1K
      AK = AP( KC+N-K+1 ) / AKM1K
      DENOM = AKM1*AK - ONE
      DO 70 J = 1, NRHS
        BKM1 = B( K, J ) / AKM1K
        BK = B( K+1, J ) / AKM1K
        B( K, J ) = ( AK*BKM1-BK ) / DENOM
        B( K+1, J ) = ( AKM1*BK-BKM1 ) / DENOM
   70 CONTINUE
      KC = KC + 2*( N-K ) + 1
      K = K + 2
    END IF
!
    GO TO 60
   80    CONTINUE
!
!        Next solve L'*X = B, overwriting B with X.
!
!        K is the main loop index, decreasing from N to 1 in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
    K = N
    KC = N*( N+1 ) / 2 + 1
   90    CONTINUE
!
!        If K < 1, exit from loop.
!
    IF( K.LT.1 ) GO TO 100
!
    KC = KC - ( N-K+1 )
    IF( IPIV( K ).GT.0 ) THEN
!
!     1 x 1 diagonal block
!
!     Multiply by inv(L'(K)), where L(K) is the transformation
!     stored in column K of A.
!
      IF( K.LT.N ) CALL DGEMV( 'Transpose', N-K, NRHS, -ONE, B( K+1, 1 ), &
                           LDB, AP( KC+1 ), 1_i4, ONE, B( K, 1 ), LDB )
!
!     Interchange rows K and IPIV(K).
!
      KP = IPIV( K )
      IF( KP.NE.K ) CALL DSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
      K = K - 1
    ELSE
!
!     2 x 2 diagonal block
!
!     Multiply by inv(L'(K-1)), where L(K-1) is the transformation
!     stored in columns K-1 and K of A.
!
      IF( K.LT.N ) THEN
        CALL DGEMV( 'Transpose', N-K, NRHS, -ONE, B( K+1, 1 ), LDB, AP( KC+1 ), 1_i4, ONE, B( K, 1 ), LDB )
        CALL DGEMV( 'Transpose', N-K, NRHS, -ONE, B( K+1, 1 ), LDB, AP( KC-( N-K ) ), 1_i4, ONE, B( K-1, 1 ), LDB )
      END IF
!
!           Interchange rows K and -IPIV(K).
!
      KP = -IPIV( K )
      IF( KP.NE.K ) CALL DSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
      KC = KC - ( N-K+2 )
      K = K - 2
    ENDIF
!
    GO TO 90
  100    CONTINUE
  ENDIF
!
  RETURN
!
!     End of DSPTRS
!
  END
!
!> \brief <b> ZHEEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for HE matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download ZHEEV + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zheev.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zheev.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zheev.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK,
!                         INFO )
! 
!       .. Scalar Arguments ..
!       CHARACTER          JOBZ, UPLO
!       INTEGER            INFO, LDA, LWORK, N
!       ..
!       .. Array Arguments ..
!       real(dp) ::   RWORK( * ), W( * )
!       COMPLEX*16         A( LDA, * ), WORK( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZHEEV computes all eigenvalues and, optionally, eigenvectors of a
!> complex Hermitian matrix A.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOBZ
!> \verbatim
!>          JOBZ is CHARACTER*1
!>          = 'N':  Compute eigenvalues only;
!>          = 'V':  Compute eigenvalues and eigenvectors.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  Upper triangle of A is stored;
!>          = 'L':  Lower triangle of A is stored.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA, N)
!>          On entry, the Hermitian matrix A.  If UPLO = 'U', the
!>          leading N-by-N upper triangular part of A contains the
!>          upper triangular part of the matrix A.  If UPLO = 'L',
!>          the leading N-by-N lower triangular part of A contains
!>          the lower triangular part of the matrix A.
!>          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
!>          orthonormal eigenvectors of the matrix A.
!>          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
!>          or the upper triangle (if UPLO='U') of A, including the
!>          diagonal, is destroyed.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is real(dp) :: array, dimension (N)
!>          If INFO = 0, the eigenvalues in ascending order.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The length of the array WORK.  LWORK >= max(1,2*N-1).
!>          For optimal efficiency, LWORK >= (NB+1)*N,
!>          where NB is the blocksize for ZHETRD returned by ILAENV.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is real(dp) :: array, dimension (max(1, 3*N-2))
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  if INFO = i, the algorithm failed to converge; i
!>                off-diagonal elements of an intermediate tridiagonal
!>                form did not converge to zero.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup complex16HEeigen
!
!  =====================================================================
      SUBROUTINE ZHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, INFO )
!
!  -- LAPACK driver routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes
      CHARACTER    :: JOBZ, UPLO
      integer(i4)  :: INFO, LDA, LWORK, N
!     ..
!     .. Array Arguments ..
      real(dp)     :: RWORK( * ), W( * )
      complex(dpc) :: A( LDA, * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      real(dp)     :: ZERO, ONE
      PARAMETER ( ZERO = 0.0_dp, ONE = 1.0_dp )
      complex(dpc) :: CONE
      PARAMETER       ( CONE = ( 1.0D0, 0.0D0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL      :: LOWER, LQUERY, WANTZ
      INTEGER(i4)  :: IINFO, IMAX, INDE, INDTAU, INDWRK, ISCALE, LLWORK, LWKOPT, NB
      real(dp)     :: ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM
!     ..
!     .. External Functions ..
      LOGICAL      :: LSAME
      integer(i4)  :: ILAENV
      real(dp)     :: DLAMCH, ZLANHE
      EXTERNAL        LSAME, ILAENV, DLAMCH, ZLANHE
!     ..
!     .. External Subroutines ..
      EXTERNAL           DSCAL, DSTERF, XERBLA, ZHETRD, ZLASCL, ZSTEQR, ZUNGTR
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      WANTZ = LSAME( JOBZ, 'V' )
      LOWER = LSAME( UPLO, 'L' )
      LQUERY = ( LWORK.EQ.-1 )
!
      INFO = 0
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
!
      IF( INFO.EQ.0 ) THEN
         NB = ILAENV( 1, 'ZHETRD', UPLO, N, -1, -1, -1 )
         LWKOPT = MAX( 1, ( NB+1 )*N )
         WORK( 1 ) = LWKOPT
!
         IF( LWORK.LT.MAX( 1, 2*N-1 ) .AND. .NOT.LQUERY ) INFO = -8
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHEEV ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) THEN
         RETURN
      END IF
!
      IF( N.EQ.1 ) THEN
         W( 1 ) = A( 1, 1 )
         WORK( 1 ) = 1
         IF( WANTZ ) A( 1, 1 ) = CONE
         RETURN
      END IF
!
!     Get machine constants.
!
      SAFMIN = DLAMCH( 'Safe minimum' )
      EPS = DLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = SQRT( BIGNUM )
!
!     Scale matrix to allowable range, if necessary.
!
      ANRM = ZLANHE( 'M', UPLO, N, A, LDA, RWORK )
      ISCALE = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / ANRM
      ELSE IF( ANRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / ANRM
      END IF
      IF( ISCALE.EQ.1 ) CALL ZLASCL( UPLO, 0, 0, ONE, SIGMA, N, N, A, LDA, INFO )
!
!     Call ZHETRD to reduce Hermitian matrix to tridiagonal form.
!
      INDE = 1
      INDTAU = 1
      INDWRK = INDTAU + N
      LLWORK = LWORK - INDWRK + 1
      CALL ZHETRD( UPLO, N, A, LDA, W, RWORK( INDE ), WORK( INDTAU ), WORK( INDWRK ), LLWORK, IINFO )
!
!     For eigenvalues only, call DSTERF.  For eigenvectors, first call
!     ZUNGTR to generate the unitary matrix, then call ZSTEQR.
!
      IF( .NOT.WANTZ ) THEN
         CALL DSTERF( N, W, RWORK( INDE ), INFO )
      ELSE
         CALL ZUNGTR( UPLO, N, A, LDA, WORK( INDTAU ), WORK( INDWRK ), LLWORK, IINFO )
         INDWRK = INDE + N
         CALL ZSTEQR( JOBZ, N, W, RWORK( INDE ), A, LDA, RWORK( INDWRK ), INFO )
      END IF
!
!     If matrix was scaled, then rescale eigenvalues appropriately.
!
      IF( ISCALE.EQ.1 ) THEN
         IF( INFO.EQ.0 ) THEN
            IMAX = N
         ELSE
            IMAX = INFO - 1
         END IF
         CALL DSCAL( IMAX, ONE / SIGMA, W, 1 )
      END IF
!
!     Set WORK(1) to optimal complex workspace size.
!
      WORK( 1 ) = LWKOPT
!
      RETURN
!
!     End of ZHEEV
!
      END
!> \brief \b ZHETD2
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download ZHETD2 + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetd2.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetd2.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetd2.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZHETD2( UPLO, N, A, LDA, D, E, TAU, INFO )
! 
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       integer(i4) ::            INFO, LDA, N
!       ..
!       .. Array Arguments ..
!       real(dp) ::   D( * ), E( * )
!       COMPLEX*16         A( LDA, * ), TAU( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZHETD2 reduces a complex Hermitian matrix A to real symmetric
!> tridiagonal form T by a unitary similarity transformation:
!> Q**H * A * Q = T.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the upper or lower triangular part of the
!>          Hermitian matrix A is stored:
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is integer(i4) ::
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
!>          n-by-n upper triangular part of A contains the upper
!>          triangular part of the matrix A, and the strictly lower
!>          triangular part of A is not referenced.  If UPLO = 'L', the
!>          leading n-by-n lower triangular part of A contains the lower
!>          triangular part of the matrix A, and the strictly upper
!>          triangular part of A is not referenced.
!>          On exit, if UPLO = 'U', the diagonal and first superdiagonal
!>          of A are overwritten by the corresponding elements of the
!>          tridiagonal matrix T, and the elements above the first
!>          superdiagonal, with the array TAU, represent the unitary
!>          matrix Q as a product of elementary reflectors; if UPLO
!>          = 'L', the diagonal and first subdiagonal of A are over-
!>          written by the corresponding elements of the tridiagonal
!>          matrix T, and the elements below the first subdiagonal, with
!>          the array TAU, represent the unitary matrix Q as a product
!>          of elementary reflectors. See Further Details.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is integer(i4) ::
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] D
!> \verbatim
!>          D is real(dp) :: array, dimension (N)
!>          The diagonal elements of the tridiagonal matrix T:
!>          D(i) = A(i,i).
!> \endverbatim
!>
!> \param[out] E
!> \verbatim
!>          E is real(dp) :: array, dimension (N-1)
!>          The off-diagonal elements of the tridiagonal matrix T:
!>          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is COMPLEX*16 array, dimension (N-1)
!>          The scalar factors of the elementary reflectors (see Further
!>          Details).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is integer(i4) ::
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup complex16HEcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  If UPLO = 'U', the matrix Q is represented as a product of elementary
!>  reflectors
!>
!>     Q = H(n-1) . . . H(2) H(1).
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - tau * v * v**H
!>
!>  where tau is a complex scalar, and v is a complex vector with
!>  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
!>  A(1:i-1,i+1), and tau in TAU(i).
!>
!>  If UPLO = 'L', the matrix Q is represented as a product of elementary
!>  reflectors
!>
!>     Q = H(1) H(2) . . . H(n-1).
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - tau * v * v**H
!>
!>  where tau is a complex scalar, and v is a complex vector with
!>  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
!>  and tau in TAU(i).
!>
!>  The contents of A on exit are illustrated by the following examples
!>  with n = 5:
!>
!>  if UPLO = 'U':                       if UPLO = 'L':
!>
!>    (  d   e   v2  v3  v4 )              (  d                  )
!>    (      d   e   v3  v4 )              (  e   d              )
!>    (          d   e   v4 )              (  v1  e   d          )
!>    (              d   e  )              (  v1  v2  e   d      )
!>    (                  d  )              (  v1  v2  v3  e   d  )
!>
!>  where d and e denote diagonal and off-diagonal elements of T, and vi
!>  denotes an element of the vector defining H(i).
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZHETD2( UPLO, N, A, LDA, D, E, TAU, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes
      CHARACTER          UPLO
      integer(i4) ::            INFO, LDA, N
!     ..
!     .. Array Arguments ..
      real(dp) ::   D( * ), E( * )
      complex(dpc) ::     A( LDA, * ), TAU( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      complex(dpc) ::     ONE, ZERO, HALF
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ), &
                         ZERO = ( 0.0D+0, 0.0D+0 ), &
                         HALF = ( 0.5D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER
      integer(i4) ::            I
      complex(dpc) ::        ALPHA, TAUI
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZAXPY, ZHEMV, ZHER2, ZLARFG
!     ..
!     .. External Functions ..
      LOGICAL      ::    LSAME
      complex(dpc) ::    ZDOTC
      EXTERNAL           LSAME, ZDOTC
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U')
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHETD2', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.LE.0 ) RETURN
!
      IF( UPPER ) THEN
!
!        Reduce the upper triangle of A
!
         A( N, N ) = DBLE( A( N, N ) )
         DO 10 I = N - 1, 1, -1
!
!           Generate elementary reflector H(i) = I - tau * v * v**H
!           to annihilate A(1:i-1,i+1)
!
            ALPHA = A( I, I+1 )
            CALL ZLARFG( I, ALPHA, A( 1, I+1 ), 1, TAUI )
            E( I ) = ALPHA
!
            IF( TAUI.NE.ZERO ) THEN
!
!              Apply H(i) from both sides to A(1:i,1:i)
!
               A( I, I+1 ) = ONE
!
!              Compute  x := tau * A * v  storing x in TAU(1:i)
!
               CALL ZHEMV( UPLO, I, TAUI, A, LDA, A( 1, I+1 ), 1, ZERO, TAU, 1 )
!
!              Compute  w := x - 1/2 * tau * (x**H * v) * v
!
               ALPHA = -HALF*TAUI*ZDOTC( I, TAU, 1, A( 1, I+1 ), 1 )
               CALL ZAXPY( I, ALPHA, A( 1, I+1 ), 1, TAU, 1 )
!
!              Apply the transformation as a rank-2 update:
!                 A := A - v * w**H - w * v**H
!
               CALL ZHER2( UPLO, I, -ONE, A( 1, I+1 ), 1, TAU, 1, A, LDA )
!
            ELSE
               A( I, I ) = DBLE( A( I, I ) )
            END IF
            A( I, I+1 ) = E( I )
            D( I+1 ) = A( I+1, I+1 )
            TAU( I ) = TAUI
   10    CONTINUE
         D( 1 ) = A( 1, 1 )
      ELSE
!
!        Reduce the lower triangle of A
!
         A( 1, 1 ) = DBLE( A( 1, 1 ) )
         DO 20 I = 1, N - 1
!
!           Generate elementary reflector H(i) = I - tau * v * v**H
!           to annihilate A(i+2:n,i)
!
            ALPHA = A( I+1, I )
            CALL ZLARFG( N-I, ALPHA, A( MIN( I+2, N ), I ), 1, TAUI )
            E( I ) = ALPHA
!
            IF( TAUI.NE.ZERO ) THEN
!
!              Apply H(i) from both sides to A(i+1:n,i+1:n)
!
               A( I+1, I ) = ONE
!
!              Compute  x := tau * A * v  storing y in TAU(i:n-1)
!
               CALL ZHEMV( UPLO, N-I, TAUI, A( I+1, I+1 ), LDA, A( I+1, I ), 1, ZERO, TAU( I ), 1 )
!
!              Compute  w := x - 1/2 * tau * (x**H * v) * v
!
               ALPHA = -HALF*TAUI*ZDOTC( N-I, TAU( I ), 1, A( I+1, I ), 1 )
               CALL ZAXPY( N-I, ALPHA, A( I+1, I ), 1, TAU( I ), 1 )
!
!              Apply the transformation as a rank-2 update:
!                 A := A - v * w**H - w * v**H
!
               CALL ZHER2( UPLO, N-I, -ONE, A( I+1, I ), 1, TAU( I ), 1, A( I+1, I+1 ), LDA )
!
            ELSE
               A( I+1, I+1 ) = DBLE( A( I+1, I+1 ) )
            END IF
            A( I+1, I ) = E( I )
            D( I ) = A( I, I )
            TAU( I ) = TAUI
   20    CONTINUE
         D( N ) = A( N, N )
      END IF
!
      RETURN
!
!     End of ZHETD2
!
      END
!> \brief \b ZHETRD
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download ZHETRD + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhetrd.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhetrd.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhetrd.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZHETRD( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, INFO )
! 
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       integer(i4) ::            INFO, LDA, LWORK, N
!       ..
!       .. Array Arguments ..
!       real(dp) ::   D( * ), E( * )
!       COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZHETRD reduces a complex Hermitian matrix A to real symmetric
!> tridiagonal form T by a unitary similarity transformation:
!> Q**H * A * Q = T.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  Upper triangle of A is stored;
!>          = 'L':  Lower triangle of A is stored.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is integer(i4) ::
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
!>          N-by-N upper triangular part of A contains the upper
!>          triangular part of the matrix A, and the strictly lower
!>          triangular part of A is not referenced.  If UPLO = 'L', the
!>          leading N-by-N lower triangular part of A contains the lower
!>          triangular part of the matrix A, and the strictly upper
!>          triangular part of A is not referenced.
!>          On exit, if UPLO = 'U', the diagonal and first superdiagonal
!>          of A are overwritten by the corresponding elements of the
!>          tridiagonal matrix T, and the elements above the first
!>          superdiagonal, with the array TAU, represent the unitary
!>          matrix Q as a product of elementary reflectors; if UPLO
!>          = 'L', the diagonal and first subdiagonal of A are over-
!>          written by the corresponding elements of the tridiagonal
!>          matrix T, and the elements below the first subdiagonal, with
!>          the array TAU, represent the unitary matrix Q as a product
!>          of elementary reflectors. See Further Details.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is integer(i4) ::
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] D
!> \verbatim
!>          D is real(dp) :: array, dimension (N)
!>          The diagonal elements of the tridiagonal matrix T:
!>          D(i) = A(i,i).
!> \endverbatim
!>
!> \param[out] E
!> \verbatim
!>          E is real(dp) :: array, dimension (N-1)
!>          The off-diagonal elements of the tridiagonal matrix T:
!>          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is COMPLEX*16 array, dimension (N-1)
!>          The scalar factors of the elementary reflectors (see Further
!>          Details).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is integer(i4) ::
!>          The dimension of the array WORK.  LWORK >= 1.
!>          For optimum performance LWORK >= N*NB, where NB is the
!>          optimal blocksize.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is integer(i4) ::
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup complex16HEcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  If UPLO = 'U', the matrix Q is represented as a product of elementary
!>  reflectors
!>
!>     Q = H(n-1) . . . H(2) H(1).
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - tau * v * v**H
!>
!>  where tau is a complex scalar, and v is a complex vector with
!>  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
!>  A(1:i-1,i+1), and tau in TAU(i).
!>
!>  If UPLO = 'L', the matrix Q is represented as a product of elementary
!>  reflectors
!>
!>     Q = H(1) H(2) . . . H(n-1).
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - tau * v * v**H
!>
!>  where tau is a complex scalar, and v is a complex vector with
!>  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
!>  and tau in TAU(i).
!>
!>  The contents of A on exit are illustrated by the following examples
!>  with n = 5:
!>
!>  if UPLO = 'U':                       if UPLO = 'L':
!>
!>    (  d   e   v2  v3  v4 )              (  d                  )
!>    (      d   e   v3  v4 )              (  e   d              )
!>    (          d   e   v4 )              (  v1  e   d          )
!>    (              d   e  )              (  v1  v2  e   d      )
!>    (                  d  )              (  v1  v2  v3  e   d  )
!>
!>  where d and e denote diagonal and off-diagonal elements of T, and vi
!>  denotes an element of the vector defining H(i).
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZHETRD( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes
      CHARACTER          UPLO
      integer(i4) ::            INFO, LDA, LWORK, N
!     ..
!     .. Array Arguments ..
      real(dp) ::   D( * ), E( * )
      complex(dpc) ::        A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      real(dp) ::   ONE
      PARAMETER          ( ONE = 1.0D+0 )
      complex(dpc) ::        CONE
      PARAMETER          ( CONE = ( 1.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER
      integer(i4) ::            I, IINFO, IWS, J, KK, LDWORK, LWKOPT, NB, NBMIN, NX
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZHER2K, ZHETD2, ZLATRD
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      integer(i4) ::            ILAENV
      EXTERNAL           LSAME, ILAENV
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.1 .AND. .NOT.LQUERY ) THEN
         INFO = -9
      END IF
!
      IF( INFO.EQ.0 ) THEN
!
!        Determine the block size.
!
         NB = ILAENV( 1, 'ZHETRD', UPLO, N, -1, -1, -1 )
         LWKOPT = N*NB
         WORK( 1 ) = LWKOPT
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHETRD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!
      NX = N
      IWS = 1
      IF( NB.GT.1 .AND. NB.LT.N ) THEN
!
!        Determine when to cross over from blocked to unblocked code
!        (last block is always handled by unblocked code).
!
         NX = MAX( NB, ILAENV( 3, 'ZHETRD', UPLO, N, -1, -1, -1 ) )
         IF( NX.LT.N ) THEN
!
!           Determine if workspace is large enough for blocked code.
!
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
!
!              Not enough workspace to use optimal NB:  determine the
!              minimum value of NB, and reduce NB or force use of
!              unblocked code by setting NX = N.
!
               NB = MAX( LWORK / LDWORK, 1 )
               NBMIN = ILAENV( 2, 'ZHETRD', UPLO, N, -1, -1, -1 )
               IF( NB.LT.NBMIN ) NX = N
            END IF
         ELSE
            NX = N
         END IF
      ELSE
         NB = 1
      END IF
!
      IF( UPPER ) THEN
!
!        Reduce the upper triangle of A.
!        Columns 1:kk are handled by the unblocked method.
!
         KK = N - ( ( N-NX+NB-1 ) / NB )*NB
         DO 20 I = N - NB + 1, KK + 1, -NB
!
!           Reduce columns i:i+nb-1 to tridiagonal form and form the
!           matrix W which is needed to update the unreduced part of
!           the matrix
!
            CALL ZLATRD( UPLO, I+NB-1, NB, A, LDA, E, TAU, WORK, LDWORK )
!
!           Update the unreduced submatrix A(1:i-1,1:i-1), using an
!           update of the form:  A := A - V*W**H - W*V**H
!
            CALL ZHER2K( UPLO, 'No transpose', I-1, NB, -CONE, A( 1, I ), LDA, WORK, LDWORK, ONE, A, LDA )
!
!           Copy superdiagonal elements back into A, and diagonal
!           elements into D
!
            DO 10 J = I, I + NB - 1
               A( J-1, J ) = E( J-1 )
               D( J ) = A( J, J )
   10       CONTINUE
   20    CONTINUE
!
!        Use unblocked code to reduce the last or only block
!
         CALL ZHETD2( UPLO, KK, A, LDA, D, E, TAU, IINFO )
      ELSE
!
!        Reduce the lower triangle of A
!
         DO 40 I = 1, N - NX, NB
!
!           Reduce columns i:i+nb-1 to tridiagonal form and form the
!           matrix W which is needed to update the unreduced part of
!           the matrix
!
            CALL ZLATRD( UPLO, N-I+1, NB, A( I, I ), LDA, E( I ), TAU( I ), WORK, LDWORK )
!
!           Update the unreduced submatrix A(i+nb:n,i+nb:n), using
!           an update of the form:  A := A - V*W**H - W*V**H
!
            CALL ZHER2K( UPLO, 'No transpose', N-I-NB+1, NB, -CONE, &
                         A( I+NB, I ), LDA, WORK( NB+1 ), LDWORK, ONE, &
                         A( I+NB, I+NB ), LDA )
!
!           Copy subdiagonal elements back into A, and diagonal
!           elements into D
!
            DO 30 J = I, I + NB - 1
               A( J+1, J ) = E( J )
               D( J ) = A( J, J )
   30       CONTINUE
   40    CONTINUE
!
!        Use unblocked code to reduce the last or only block
!
         CALL ZHETD2( UPLO, N-I+1, A( I, I ), LDA, D( I ), E( I ), TAU( I ), IINFO )
      END IF
!
      WORK( 1 ) = LWKOPT
      RETURN
!
!     End of ZHETRD
!
      END
!> \brief \b ZLACGV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download ZLACGV + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlacgv.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlacgv.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlacgv.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLACGV( N, X, INCX )
! 
!       .. Scalar Arguments ..
!       integer(i4) ::            INCX, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         X( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLACGV conjugates a complex vector of length N.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is integer(i4) ::
!>          The length of the vector X.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension
!>                         (1+(N-1)*abs(INCX))
!>          On entry, the vector of length N to be conjugated.
!>          On exit, X is overwritten with conjg(X).
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is integer(i4) ::
!>          The spacing between successive elements of X.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup complex16OTHERauxiliary
!
!  =====================================================================
      SUBROUTINE ZLACGV( N, X, INCX )
!
!  -- LAPACK auxiliary routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes
      integer(i4) ::            INCX, N
!     ..
!     .. Array Arguments ..
      complex(dpc) :: X( * )
!     ..
!
! =====================================================================
!
!     .. Local Scalars ..
      integer(i4) ::            I, IOFF
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DCONJG
!     ..
!     .. Executable Statements ..
!
      IF( INCX.EQ.1 ) THEN
         DO 10 I = 1, N
            X( I ) = DCONJG( X( I ) )
   10    CONTINUE
      ELSE
         IOFF = 1
         IF( INCX.LT.0 ) IOFF = 1 - ( N-1 )*INCX
         DO 20 I = 1, N
            X( IOFF ) = DCONJG( X( IOFF ) )
            IOFF = IOFF + INCX
   20    CONTINUE
      END IF
      RETURN
!
!     End of ZLACGV
!
      END
!> \brief \b ZLADIV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download ZLADIV + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zladiv.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zladiv.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zladiv.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       COMPLEX*16     FUNCTION ZLADIV( X, Y )
! 
!       .. Scalar Arguments ..
!       COMPLEX*16         X, Y
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLADIV := X / Y, where X and Y are complex.  The computation of X / Y
!> will not overflow on an intermediary step unless the results
!> overflows.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] X
!> \verbatim
!>          X is COMPLEX*16
!> \endverbatim
!>
!> \param[in] Y
!> \verbatim
!>          Y is COMPLEX*16
!>          The complex scalars X and Y.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup complex16OTHERauxiliary
!
!  =====================================================================
      FUNCTION ZLADIV( X, Y )
!
!  -- LAPACK auxiliary routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes
      complex(dpc) :: X, Y
      complex(dpc) :: zladiv
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      real(dp) ::   ZI, ZR
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLADIV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, DIMAG
!     ..
!     .. Executable Statements ..
!
      CALL DLADIV( DBLE( X ), DIMAG( X ), DBLE( Y ), DIMAG( Y ), ZR, ZI )
      ZLADIV = DCMPLX( ZR, ZI )
!
      RETURN
!
!     End of ZLADIV
!
      END
!> \brief \b ZLANHE
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download ZLANHE + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlanhe.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlanhe.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlanhe.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       real(dp) :: FUNCTION ZLANHE( NORM, UPLO, N, A, LDA, WORK )
! 
!       .. Scalar Arguments ..
!       CHARACTER          NORM, UPLO
!       integer(i4) ::            LDA, N
!       ..
!       .. Array Arguments ..
!       real(dp) ::   WORK( * )
!       COMPLEX*16         A( LDA, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLANHE  returns the value of the one norm,  or the Frobenius norm, or
!> the  infinity norm,  or the  element of  largest absolute value  of a
!> complex hermitian matrix A.
!> \endverbatim
!>
!> \return ZLANHE
!> \verbatim
!>
!>    ZLANHE = ( max(abs(A(i,j))), NORM = 'M' or 'm'
!>             (
!>             ( norm1(A),         NORM = '1', 'O' or 'o'
!>             (
!>             ( normI(A),         NORM = 'I' or 'i'
!>             (
!>             ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
!>
!> where  norm1  denotes the  one norm of a matrix (maximum column sum),
!> normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
!> normF  denotes the  Frobenius norm of a matrix (square root of sum of
!> squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] NORM
!> \verbatim
!>          NORM is CHARACTER*1
!>          Specifies the value to be returned in ZLANHE as described
!>          above.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the upper or lower triangular part of the
!>          hermitian matrix A is to be referenced.
!>          = 'U':  Upper triangular part of A is referenced
!>          = 'L':  Lower triangular part of A is referenced
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is integer(i4) ::
!>          The order of the matrix A.  N >= 0.  When N = 0, ZLANHE is
!>          set to zero.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          The hermitian matrix A.  If UPLO = 'U', the leading n by n
!>          upper triangular part of A contains the upper triangular part
!>          of the matrix A, and the strictly lower triangular part of A
!>          is not referenced.  If UPLO = 'L', the leading n by n lower
!>          triangular part of A contains the lower triangular part of
!>          the matrix A, and the strictly upper triangular part of A is
!>          not referenced. Note that the imaginary parts of the diagonal
!>          elements need not be set and are assumed to be zero.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is integer(i4) ::
!>          The leading dimension of the array A.  LDA >= max(N,1).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is real(dp) :: array, dimension (MAX(1,LWORK)),
!>          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,
!>          WORK is not referenced.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup complex16HEauxiliary
!
!  =====================================================================
      FUNCTION ZLANHE( NORM, UPLO, N, A, LDA, WORK )
!
!  -- LAPACK auxiliary routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes
      CHARACTER    :: NORM, UPLO
      integer(i4)  :: LDA, N
      real(dp)     :: zlanhe
!     ..
!     .. Array Arguments ..
      real(dp)     :: WORK( * )
      complex(dpc) :: A( LDA, * )
!     ..
!
! =====================================================================
!
!     .. Parameters ..
      real(dp) ::   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      integer(i4) ::            I, J
      real(dp) ::   ABSA, SCALE, SUM, VALUE
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           ZLASSQ
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, SQRT
!     ..
!     .. Executable Statements ..
!
      IF( N.EQ.0 ) THEN
         VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
!
!        Find max(abs(A(i,j))).
!
         VALUE = ZERO
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 20 J = 1, N
               DO 10 I = 1, J - 1
                  VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   10          CONTINUE
               VALUE = MAX( VALUE, ABS( DBLE( A( J, J ) ) ) )
   20       CONTINUE
         ELSE
            DO 40 J = 1, N
               VALUE = MAX( VALUE, ABS( DBLE( A( J, J ) ) ) )
               DO 30 I = J + 1, N
                  VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   30          CONTINUE
   40       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'I' ) ) .OR. ( LSAME( NORM, 'O' ) ) .OR. ( NORM.EQ.'1' ) ) THEN
!
!        Find normI(A) ( = norm1(A), since A is hermitian).
!
         VALUE = ZERO
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 60 J = 1, N
               SUM = ZERO
               DO 50 I = 1, J - 1
                  ABSA = ABS( A( I, J ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
   50          CONTINUE
               WORK( J ) = SUM + ABS( DBLE( A( J, J ) ) )
   60       CONTINUE
            DO 70 I = 1, N
               VALUE = MAX( VALUE, WORK( I ) )
   70       CONTINUE
         ELSE
            DO 80 I = 1, N
               WORK( I ) = ZERO
   80       CONTINUE
            DO 100 J = 1, N
               SUM = WORK( J ) + ABS( DBLE( A( J, J ) ) )
               DO 90 I = J + 1, N
                  ABSA = ABS( A( I, J ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
   90          CONTINUE
               VALUE = MAX( VALUE, SUM )
  100       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
!
!        Find normF(A).
!
         SCALE = ZERO
         SUM = ONE
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 110 J = 2, N
               CALL ZLASSQ( J-1, A( 1, J ), 1, SCALE, SUM )
  110       CONTINUE
         ELSE
            DO 120 J = 1, N - 1
               CALL ZLASSQ( N-J, A( J+1, J ), 1, SCALE, SUM )
  120       CONTINUE
         END IF
         SUM = 2*SUM
         DO 130 I = 1, N
            IF( DBLE( A( I, I ) ).NE.ZERO ) THEN
               ABSA = ABS( DBLE( A( I, I ) ) )
               IF( SCALE.LT.ABSA ) THEN
                  SUM = ONE + SUM*( SCALE / ABSA )**2
                  SCALE = ABSA
               ELSE
                  SUM = SUM + ( ABSA / SCALE )**2
               END IF
            END IF
  130    CONTINUE
         VALUE = SCALE*SQRT( SUM )
      END IF
!
      ZLANHE = VALUE
      RETURN
!
!     End of ZLANHE
!
      END
!> \brief \b ZLARF
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download ZLARF + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlarf.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlarf.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlarf.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
! 
!       .. Scalar Arguments ..
!       CHARACTER          SIDE
!       integer(i4) ::            INCV, LDC, M, N
!       COMPLEX*16         TAU
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         C( LDC, * ), V( * ), WORK( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLARF applies a complex elementary reflector H to a complex M-by-N
!> matrix C, from either the left or the right. H is represented in the
!> form
!>
!>       H = I - tau * v * v**H
!>
!> where tau is a complex scalar and v is a complex vector.
!>
!> If tau = 0, then H is taken to be the unit matrix.
!>
!> To apply H**H, supply conjg(tau) instead
!> tau.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'L': form  H * C
!>          = 'R': form  C * H
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is integer(i4) ::
!>          The number of rows of the matrix C.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is integer(i4) ::
!>          The number of columns of the matrix C.
!> \endverbatim
!>
!> \param[in] V
!> \verbatim
!>          V is COMPLEX*16 array, dimension
!>                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
!>                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
!>          The vector v in the representation of H. V is not used if
!>          TAU = 0.
!> \endverbatim
!>
!> \param[in] INCV
!> \verbatim
!>          INCV is integer(i4) ::
!>          The increment between elements of v. INCV <> 0.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is COMPLEX*16
!>          The value tau in the representation of H.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is COMPLEX*16 array, dimension (LDC,N)
!>          On entry, the M-by-N matrix C.
!>          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
!>          or C * H if SIDE = 'R'.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is integer(i4) ::
!>          The leading dimension of the array C. LDC >= max(1,M).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension
!>                         (N) if SIDE = 'L'
!>                      or (M) if SIDE = 'R'
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup complex16OTHERauxiliary
!
!  =====================================================================
      SUBROUTINE ZLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
!
!  -- LAPACK auxiliary routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes
      CHARACTER    :: SIDE
      integer(i4)  :: INCV, LDC, M, N
      complex(dpc) :: TAU
!     ..
!     .. Array Arguments ..
      complex(dpc) :: C( LDC, * ), V( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      complex(dpc) :: ONE, ZERO
      PARAMETER ( ONE = ( 1.0D+0, 0.0D+0 ), ZERO = ( 0.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL     :: APPLYLEFT
      integer(i4) :: I, LASTV, LASTC
!     ..
!     .. External Subroutines ..
      EXTERNAL           ZGEMV, ZGERC
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      integer(i4) ::            ILAZLR, ILAZLC
      EXTERNAL           LSAME, ILAZLR, ILAZLC
!     ..
!     .. Executable Statements ..
!
      APPLYLEFT = LSAME( SIDE, 'L' )
      LASTV = 0
      LASTC = 0
      IF( TAU.NE.ZERO ) THEN
!     Set up variables for scanning V.  LASTV begins pointing to the end
!     of V.
         IF( APPLYLEFT ) THEN
            LASTV = M
         ELSE
            LASTV = N
         END IF
         IF( INCV.GT.0 ) THEN
            I = 1 + (LASTV-1) * INCV
         ELSE
            I = 1
         END IF
!     Look for the last non-zero row in V.
         DO WHILE( LASTV.GT.0 .AND. V( I ).EQ.ZERO )
            LASTV = LASTV - 1
            I = I - INCV
         END DO
         IF( APPLYLEFT ) THEN
!     Scan for the last non-zero column in C(1:lastv,:).
            LASTC = ILAZLC(LASTV, N, C, LDC)
         ELSE
!     Scan for the last non-zero row in C(:,1:lastv).
            LASTC = ILAZLR(M, LASTV, C, LDC)
         END IF
      END IF
!     Note that lastc.eq.0 renders the BLAS operations null; no special
!     case is needed at this level.
      IF( APPLYLEFT ) THEN
!
!        Form  H * C
!
         IF( LASTV.GT.0 ) THEN
!
!           w(1:lastc,1) := C(1:lastv,1:lastc)**H * v(1:lastv,1)
!
            CALL ZGEMV( 'Conjugate transpose', LASTV, LASTC, ONE, C, LDC, V, INCV, ZERO, WORK, 1 )
!
!           C(1:lastv,1:lastc) := C(...) - v(1:lastv,1) * w(1:lastc,1)**H
!
            CALL ZGERC( LASTV, LASTC, -TAU, V, INCV, WORK, 1, C, LDC )
         END IF
      ELSE
!
!        Form  C * H
!
         IF( LASTV.GT.0 ) THEN
!
!           w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1)
!
            CALL ZGEMV( 'No transpose', LASTC, LASTV, ONE, C, LDC, V, INCV, ZERO, WORK, 1 )
!
!           C(1:lastc,1:lastv) := C(...) - w(1:lastc,1) * v(1:lastv,1)**H
!
            CALL ZGERC( LASTC, LASTV, -TAU, WORK, 1, V, INCV, C, LDC )
         END IF
      END IF
      RETURN
!
!     End of ZLARF
!
      END
!> \brief \b ZLARFB
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download ZLARFB + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlarfb.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlarfb.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlarfb.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV,
!                          T, LDT, C, LDC, WORK, LDWORK )
! 
!       .. Scalar Arguments ..
!       CHARACTER          DIRECT, SIDE, STOREV, TRANS
!       integer(i4) ::            K, LDC, LDT, LDV, LDWORK, M, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         C( LDC, * ), T( LDT, * ), V( LDV, * ),
!      $                   WORK( LDWORK, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLARFB applies a complex block reflector H or its transpose H**H to a
!> complex M-by-N matrix C, from either the left or the right.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'L': apply H or H**H from the Left
!>          = 'R': apply H or H**H from the Right
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          = 'N': apply H (No transpose)
!>          = 'C': apply H**H (Conjugate transpose)
!> \endverbatim
!>
!> \param[in] DIRECT
!> \verbatim
!>          DIRECT is CHARACTER*1
!>          Indicates how H is formed from a product of elementary
!>          reflectors
!>          = 'F': H = H(1) H(2) . . . H(k) (Forward)
!>          = 'B': H = H(k) . . . H(2) H(1) (Backward)
!> \endverbatim
!>
!> \param[in] STOREV
!> \verbatim
!>          STOREV is CHARACTER*1
!>          Indicates how the vectors which define the elementary
!>          reflectors are stored:
!>          = 'C': Columnwise
!>          = 'R': Rowwise
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is integer(i4) ::
!>          The number of rows of the matrix C.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is integer(i4) ::
!>          The number of columns of the matrix C.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is integer(i4) ::
!>          The order of the matrix T (= the number of elementary
!>          reflectors whose product defines the block reflector).
!> \endverbatim
!>
!> \param[in] V
!> \verbatim
!>          V is COMPLEX*16 array, dimension
!>                                (LDV,K) if STOREV = 'C'
!>                                (LDV,M) if STOREV = 'R' and SIDE = 'L'
!>                                (LDV,N) if STOREV = 'R' and SIDE = 'R'
!>          See Further Details.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is integer(i4) ::
!>          The leading dimension of the array V.
!>          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);
!>          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);
!>          if STOREV = 'R', LDV >= K.
!> \endverbatim
!>
!> \param[in] T
!> \verbatim
!>          T is COMPLEX*16 array, dimension (LDT,K)
!>          The triangular K-by-K matrix T in the representation of the
!>          block reflector.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is integer(i4) ::
!>          The leading dimension of the array T. LDT >= K.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is COMPLEX*16 array, dimension (LDC,N)
!>          On entry, the M-by-N matrix C.
!>          On exit, C is overwritten by H*C or H**H*C or C*H or C*H**H.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is integer(i4) ::
!>          The leading dimension of the array C. LDC >= max(1,M).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (LDWORK,K)
!> \endverbatim
!>
!> \param[in] LDWORK
!> \verbatim
!>          LDWORK is integer(i4) ::
!>          The leading dimension of the array WORK.
!>          If SIDE = 'L', LDWORK >= max(1,N);
!>          if SIDE = 'R', LDWORK >= max(1,M).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup complex16OTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The shape of the matrix V and the storage of the vectors which define
!>  the H(i) is best illustrated by the following example with n = 5 and
!>  k = 3. The elements equal to 1 are not stored; the corresponding
!>  array elements are modified but restored on exit. The rest of the
!>  array is not used.
!>
!>  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
!>
!>               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
!>                   ( v1  1    )                     (     1 v2 v2 v2 )
!>                   ( v1 v2  1 )                     (        1 v3 v3 )
!>                   ( v1 v2 v3 )
!>                   ( v1 v2 v3 )
!>
!>  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
!>
!>               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
!>                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
!>                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
!>                   (     1 v3 )
!>                   (        1 )
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV, T, LDT, C, LDC, WORK, LDWORK )
!
!  -- LAPACK auxiliary routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes
      CHARACTER    :: DIRECT, SIDE, STOREV, TRANS
      integer(i4)  :: K, LDC, LDT, LDV, LDWORK, M, N
!     ..
!     .. Array Arguments ..
      complex(dpc) :: C( LDC, * ), T( LDT, * ), V( LDV, * ), WORK( LDWORK, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      complex(dpc) :: ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      CHARACTER   :: TRANST
      integer(i4) :: I, J, LASTV, LASTC
!     ..
!     .. External Functions ..
      LOGICAL     :: LSAME
      integer(i4) :: ILAZLR, ILAZLC
      EXTERNAL       LSAME, ILAZLR, ILAZLC
!     ..
!     .. External Subroutines ..
      EXTERNAL       ZCOPY, ZGEMM, ZLACGV, ZTRMM
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DCONJG
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF( M.LE.0 .OR. N.LE.0 ) RETURN
!
      IF( LSAME( TRANS, 'N' ) ) THEN
         TRANST = 'C'
      ELSE
         TRANST = 'N'
      END IF
!
      IF( LSAME( STOREV, 'C' ) ) THEN
!
         IF( LSAME( DIRECT, 'F' ) ) THEN
!
!           Let  V =  ( V1 )    (first K rows)
!                     ( V2 )
!           where  V1  is unit lower triangular.
!
            IF( LSAME( SIDE, 'L' ) ) THEN
!
!              Form  H * C  or  H**H * C  where  C = ( C1 )
!                                                    ( C2 )
!
               LASTV = MAX( K, ILAZLR( M, K, V, LDV ) )
               LASTC = ILAZLC( LASTV, N, C, LDC )
!
!              W := C**H * V  =  (C1**H * V1 + C2**H * V2)  (stored in WORK)
!
!              W := C1**H
!
               DO 10 J = 1, K
                  CALL ZCOPY( LASTC, C( J, 1 ), LDC, WORK( 1, J ), 1 )
                  CALL ZLACGV( LASTC, WORK( 1, J ), 1 )
   10          CONTINUE
!
!              W := W * V1
!
               CALL ZTRMM( 'Right', 'Lower', 'No transpose', 'Unit', LASTC, K, ONE, V, LDV, WORK, LDWORK )
               IF( LASTV.GT.K ) THEN
!
!                 W := W + C2**H *V2
!
                  CALL ZGEMM( 'Conjugate transpose', 'No transpose', &
                       LASTC, K, LASTV-K, ONE, C( K+1, 1 ), LDC, &
                       V( K+1, 1 ), LDV, ONE, WORK, LDWORK )
               END IF
!
!              W := W * T**H  or  W * T
!
               CALL ZTRMM( 'Right', 'Upper', TRANST, 'Non-unit', LASTC, K, ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V * W**H
!
               IF( M.GT.K ) THEN
!
!                 C2 := C2 - V2 * W**H
!
                  CALL ZGEMM( 'No transpose', 'Conjugate transpose', &
                       LASTV-K, LASTC, K, &
                       -ONE, V( K+1, 1 ), LDV, WORK, LDWORK, &
                       ONE, C( K+1, 1 ), LDC )
               END IF
!
!              W := W * V1**H
!
               CALL ZTRMM( 'Right', 'Lower', 'Conjugate transpose', 'Unit', LASTC, K, ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W**H
!
               DO 30 J = 1, K
                  DO 20 I = 1, LASTC
                     C( J, I ) = C( J, I ) - DCONJG( WORK( I, J ) )
   20             CONTINUE
   30          CONTINUE
!
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H**H  where  C = ( C1  C2 )
!
               LASTV = MAX( K, ILAZLR( N, K, V, LDV ) )
               LASTC = ILAZLR( M, LASTV, C, LDC )
!
!              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
!
!              W := C1
!
               DO 40 J = 1, K
                  CALL ZCOPY( LASTC, C( 1, J ), 1, WORK( 1, J ), 1 )
   40          CONTINUE
!
!              W := W * V1
!
               CALL ZTRMM( 'Right', 'Lower', 'No transpose', 'Unit', LASTC, K, ONE, V, LDV, WORK, LDWORK )
               IF( LASTV.GT.K ) THEN
!
!                 W := W + C2 * V2
!
                  CALL ZGEMM( 'No transpose', 'No transpose', LASTC, K, LASTV-K, &
                       ONE, C( 1, K+1 ), LDC, V( K+1, 1 ), LDV, ONE, WORK, LDWORK )
               END IF
!
!              W := W * T  or  W * T**H
!
               CALL ZTRMM( 'Right', 'Upper', TRANS, 'Non-unit', LASTC, K, ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V**H
!
               IF( LASTV.GT.K ) THEN
!
!                 C2 := C2 - W * V2**H
!
                  CALL ZGEMM( 'No transpose', 'Conjugate transpose', LASTC, LASTV-K, K, &
                       -ONE, WORK, LDWORK, V( K+1, 1 ), LDV, ONE, C( 1, K+1 ), LDC )
               END IF
!
!              W := W * V1**H
!
               CALL ZTRMM( 'Right', 'Lower', 'Conjugate transpose', 'Unit', LASTC, K, ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W
!
               DO 60 J = 1, K
                  DO 50 I = 1, LASTC
                     C( I, J ) = C( I, J ) - WORK( I, J )
   50             CONTINUE
   60          CONTINUE
            END IF
!
         ELSE
!
!           Let  V =  ( V1 )
!                     ( V2 )    (last K rows)
!           where  V2  is unit upper triangular.
!
            IF( LSAME( SIDE, 'L' ) ) THEN
!
!              Form  H * C  or  H**H * C  where  C = ( C1 )
!                                                    ( C2 )
!
               LASTV = MAX( K, ILAZLR( M, K, V, LDV ) )
               LASTC = ILAZLC( LASTV, N, C, LDC )
!
!              W := C**H * V  =  (C1**H * V1 + C2**H * V2)  (stored in WORK)
!
!              W := C2**H
!
               DO 70 J = 1, K
                  CALL ZCOPY( LASTC, C( LASTV-K+J, 1 ), LDC, WORK( 1, J ), 1 )
                  CALL ZLACGV( LASTC, WORK( 1, J ), 1 )
   70          CONTINUE
!
!              W := W * V2
!
               CALL ZTRMM( 'Right', 'Upper', 'No transpose', 'Unit', &
                    LASTC, K, ONE, V( LASTV-K+1, 1 ), LDV, WORK, LDWORK )
               IF( LASTV.GT.K ) THEN
!
!                 W := W + C1**H*V1
!
                  CALL ZGEMM( 'Conjugate transpose', 'No transpose', LASTC, K, LASTV-K, &
                       ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
!
!              W := W * T**H  or  W * T
!
               CALL ZTRMM( 'Right', 'Lower', TRANST, 'Non-unit', LASTC, K, ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V * W**H
!
               IF( LASTV.GT.K ) THEN
!
!                 C1 := C1 - V1 * W**H
!
                  CALL ZGEMM( 'No transpose', 'Conjugate transpose', LASTV-K, LASTC, K, &
                       -ONE, V, LDV, WORK, LDWORK, ONE, C, LDC )
               END IF
!
!              W := W * V2**H
!
               CALL ZTRMM( 'Right', 'Upper', 'Conjugate transpose', 'Unit', LASTC, K, ONE, V( LASTV-K+1, 1 ), LDV, &
                    WORK, LDWORK )
!
!              C2 := C2 - W**H
!
               DO 90 J = 1, K
                  DO 80 I = 1, LASTC
                     C( LASTV-K+J, I ) = C( LASTV-K+J, I ) - DCONJG( WORK( I, J ) )
   80             CONTINUE
   90          CONTINUE
!
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H**H  where  C = ( C1  C2 )
!
               LASTV = MAX( K, ILAZLR( N, K, V, LDV ) )
               LASTC = ILAZLR( M, LASTV, C, LDC )
!
!              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
!
!              W := C2
!
               DO 100 J = 1, K
                  CALL ZCOPY( LASTC, C( 1, LASTV-K+J ), 1, WORK( 1, J ), 1 )
  100          CONTINUE
!
!              W := W * V2
!
               CALL ZTRMM( 'Right', 'Upper', 'No transpose', 'Unit', LASTC, K, ONE, V( LASTV-K+1, 1 ), LDV, &
                    WORK, LDWORK )
               IF( LASTV.GT.K ) THEN
!
!                 W := W + C1 * V1
!
                  CALL ZGEMM( 'No transpose', 'No transpose', LASTC, K, LASTV-K, &
                       ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
!
!              W := W * T  or  W * T**H
!
               CALL ZTRMM( 'Right', 'Lower', TRANS, 'Non-unit', LASTC, K, ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V**H
!
               IF( LASTV.GT.K ) THEN
!
!                 C1 := C1 - W * V1**H
!
                  CALL ZGEMM( 'No transpose', 'Conjugate transpose', LASTC, LASTV-K, K, -ONE, WORK, LDWORK, V, LDV, &
                       ONE, C, LDC )
               END IF
!
!              W := W * V2**H
!
               CALL ZTRMM( 'Right', 'Upper', 'Conjugate transpose', 'Unit', LASTC, K, ONE, V( LASTV-K+1, 1 ), LDV, &
                    WORK, LDWORK )
!
!              C2 := C2 - W
!
               DO 120 J = 1, K
                  DO 110 I = 1, LASTC
                     C( I, LASTV-K+J ) = C( I, LASTV-K+J ) - WORK( I, J )
  110             CONTINUE
  120          CONTINUE
            END IF
         END IF
!
      ELSE IF( LSAME( STOREV, 'R' ) ) THEN
!
         IF( LSAME( DIRECT, 'F' ) ) THEN
!
!           Let  V =  ( V1  V2 )    (V1: first K columns)
!           where  V1  is unit upper triangular.
!
            IF( LSAME( SIDE, 'L' ) ) THEN
!
!              Form  H * C  or  H**H * C  where  C = ( C1 )
!                                                    ( C2 )
!
               LASTV = MAX( K, ILAZLC( K, M, V, LDV ) )
               LASTC = ILAZLC( LASTV, N, C, LDC )
!
!              W := C**H * V**H  =  (C1**H * V1**H + C2**H * V2**H) (stored in WORK)
!
!              W := C1**H
!
               DO 130 J = 1, K
                  CALL ZCOPY( LASTC, C( J, 1 ), LDC, WORK( 1, J ), 1 )
                  CALL ZLACGV( LASTC, WORK( 1, J ), 1 )
  130          CONTINUE
!
!              W := W * V1**H
!
               CALL ZTRMM( 'Right', 'Upper', 'Conjugate transpose', 'Unit', LASTC, K, ONE, V, LDV, WORK, LDWORK )
               IF( LASTV.GT.K ) THEN
!
!                 W := W + C2**H*V2**H
!
                  CALL ZGEMM( 'Conjugate transpose', 'Conjugate transpose', LASTC, K, LASTV-K, &
                       ONE, C( K+1, 1 ), LDC, V( 1, K+1 ), LDV, ONE, WORK, LDWORK )
               END IF
!
!              W := W * T**H  or  W * T
!
               CALL ZTRMM( 'Right', 'Upper', TRANST, 'Non-unit', LASTC, K, ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V**H * W**H
!
               IF( LASTV.GT.K ) THEN
!
!                 C2 := C2 - V2**H * W**H
!
                  CALL ZGEMM( 'Conjugate transpose', 'Conjugate transpose', LASTV-K, LASTC, K, &
                       -ONE, V( 1, K+1 ), LDV, WORK, LDWORK, ONE, C( K+1, 1 ), LDC )
               END IF
!
!              W := W * V1
!
               CALL ZTRMM( 'Right', 'Upper', 'No transpose', 'Unit', LASTC, K, ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W**H
!
               DO 150 J = 1, K
                  DO 140 I = 1, LASTC
                     C( J, I ) = C( J, I ) - DCONJG( WORK( I, J ) )
  140             CONTINUE
  150          CONTINUE
!
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H**H  where  C = ( C1  C2 )
!
               LASTV = MAX( K, ILAZLC( K, N, V, LDV ) )
               LASTC = ILAZLR( M, LASTV, C, LDC )
!
!              W := C * V**H  =  (C1*V1**H + C2*V2**H)  (stored in WORK)
!
!              W := C1
!
               DO 160 J = 1, K
                  CALL ZCOPY( LASTC, C( 1, J ), 1, WORK( 1, J ), 1 )
  160          CONTINUE
!
!              W := W * V1**H
!
               CALL ZTRMM( 'Right', 'Upper', 'Conjugate transpose', 'Unit', LASTC, K, ONE, V, LDV, WORK, LDWORK )
               IF( LASTV.GT.K ) THEN
!
!                 W := W + C2 * V2**H
!
                  CALL ZGEMM( 'No transpose', 'Conjugate transpose', LASTC, K, LASTV-K, ONE, C( 1, K+1 ), LDC, &
                       V( 1, K+1 ), LDV, ONE, WORK, LDWORK )
               END IF
!
!              W := W * T  or  W * T**H
!
               CALL ZTRMM( 'Right', 'Upper', TRANS, 'Non-unit', LASTC, K, ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V
!
               IF( LASTV.GT.K ) THEN
!
!                 C2 := C2 - W * V2
!
                  CALL ZGEMM( 'No transpose', 'No transpose', LASTC, LASTV-K, K, &
                       -ONE, WORK, LDWORK, V( 1, K+1 ), LDV, ONE, C( 1, K+1 ), LDC )
               END IF
!
!              W := W * V1
!
               CALL ZTRMM( 'Right', 'Upper', 'No transpose', 'Unit', LASTC, K, ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W
!
               DO 180 J = 1, K
                  DO 170 I = 1, LASTC
                     C( I, J ) = C( I, J ) - WORK( I, J )
  170             CONTINUE
  180          CONTINUE
!
            END IF
!
         ELSE
!
!           Let  V =  ( V1  V2 )    (V2: last K columns)
!           where  V2  is unit lower triangular.
!
            IF( LSAME( SIDE, 'L' ) ) THEN
!
!              Form  H * C  or  H**H * C  where  C = ( C1 )
!                                                    ( C2 )
!
               LASTV = MAX( K, ILAZLC( K, M, V, LDV ) )
               LASTC = ILAZLC( LASTV, N, C, LDC )
!
!              W := C**H * V**H  =  (C1**H * V1**H + C2**H * V2**H) (stored in WORK)
!
!              W := C2**H
!
               DO 190 J = 1, K
                  CALL ZCOPY( LASTC, C( LASTV-K+J, 1 ), LDC, WORK( 1, J ), 1 )
                  CALL ZLACGV( LASTC, WORK( 1, J ), 1 )
  190          CONTINUE
!
!              W := W * V2**H
!
               CALL ZTRMM( 'Right', 'Lower', 'Conjugate transpose', 'Unit', LASTC, K, ONE, V( 1, LASTV-K+1 ), LDV, &
                    WORK, LDWORK )
               IF( LASTV.GT.K ) THEN
!
!                 W := W + C1**H * V1**H
!
                  CALL ZGEMM( 'Conjugate transpose', 'Conjugate transpose', LASTC, K, LASTV-K, &
                       ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
!
!              W := W * T**H  or  W * T
!
               CALL ZTRMM( 'Right', 'Lower', TRANST, 'Non-unit', LASTC, K, ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V**H * W**H
!
               IF( LASTV.GT.K ) THEN
!
!                 C1 := C1 - V1**H * W**H
!
                  CALL ZGEMM( 'Conjugate transpose', 'Conjugate transpose', LASTV-K, LASTC, K, &
                       -ONE, V, LDV, WORK, LDWORK, ONE, C, LDC )
               END IF
!
!              W := W * V2
!
               CALL ZTRMM( 'Right', 'Lower', 'No transpose', 'Unit', LASTC, K, ONE, V( 1, LASTV-K+1 ), LDV, &
                    WORK, LDWORK )
!
!              C2 := C2 - W**H
!
               DO 210 J = 1, K
                  DO 200 I = 1, LASTC
                     C( LASTV-K+J, I ) = C( LASTV-K+J, I ) - DCONJG( WORK( I, J ) )
  200             CONTINUE
  210          CONTINUE
!
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H**H  where  C = ( C1  C2 )
!
               LASTV = MAX( K, ILAZLC( K, N, V, LDV ) )
               LASTC = ILAZLR( M, LASTV, C, LDC )
!
!              W := C * V**H  =  (C1*V1**H + C2*V2**H)  (stored in WORK)
!
!              W := C2
!
               DO 220 J = 1, K
                  CALL ZCOPY( LASTC, C( 1, LASTV-K+J ), 1, WORK( 1, J ), 1 )
  220          CONTINUE
!
!              W := W * V2**H
!
               CALL ZTRMM( 'Right', 'Lower', 'Conjugate transpose', 'Unit', LASTC, K, ONE, V( 1, LASTV-K+1 ), LDV, &
                    WORK, LDWORK )
               IF( LASTV.GT.K ) THEN
!
!                 W := W + C1 * V1**H
!
                  CALL ZGEMM( 'No transpose', 'Conjugate transpose', LASTC, K, LASTV-K, ONE, C, LDC, V, LDV, ONE, &
                       WORK, LDWORK )
               END IF
!
!              W := W * T  or  W * T**H
!
               CALL ZTRMM( 'Right', 'Lower', TRANS, 'Non-unit', LASTC, K, ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V
!
               IF( LASTV.GT.K ) THEN
!
!                 C1 := C1 - W * V1
!
                  CALL ZGEMM( 'No transpose', 'No transpose', LASTC, LASTV-K, K, -ONE, WORK, LDWORK, V, LDV, &
                       ONE, C, LDC )
               END IF
!
!              W := W * V2
!
               CALL ZTRMM( 'Right', 'Lower', 'No transpose', 'Unit', LASTC, K, ONE, V( 1, LASTV-K+1 ), LDV, &
                    WORK, LDWORK )
!
!              C1 := C1 - W
!
               DO 240 J = 1, K
                  DO 230 I = 1, LASTC
                     C( I, LASTV-K+J ) = C( I, LASTV-K+J ) - WORK( I, J )
  230             CONTINUE
  240          CONTINUE
!
            END IF
!
         END IF
      END IF
!
      RETURN
!
!     End of ZLARFB
!
      END
!> \brief \b ZLARFG
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download ZLARFG + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlarfg.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlarfg.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlarfg.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLARFG( N, ALPHA, X, INCX, TAU )
! 
!       .. Scalar Arguments ..
!       integer(i4) ::            INCX, N
!       COMPLEX*16         ALPHA, TAU
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         X( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLARFG generates a complex elementary reflector H of order n, such
!> that
!>
!>       H**H * ( alpha ) = ( beta ),   H**H * H = I.
!>              (   x   )   (   0  )
!>
!> where alpha and beta are scalars, with beta real, and x is an
!> (n-1)-element complex vector. H is represented in the form
!>
!>       H = I - tau * ( 1 ) * ( 1 v**H ) ,
!>                     ( v )
!>
!> where tau is a complex scalar and v is a complex (n-1)-element
!> vector. Note that H is not hermitian.
!>
!> If the elements of x are all zero and alpha is real, then tau = 0
!> and H is taken to be the unit matrix.
!>
!> Otherwise  1 <= real(tau) <= 2  and  abs(tau-1) <= 1 .
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is integer(i4) ::
!>          The order of the elementary reflector.
!> \endverbatim
!>
!> \param[in,out] ALPHA
!> \verbatim
!>          ALPHA is COMPLEX*16
!>          On entry, the value alpha.
!>          On exit, it is overwritten with the value beta.
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension
!>                         (1+(N-2)*abs(INCX))
!>          On entry, the vector x.
!>          On exit, it is overwritten with the vector v.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is integer(i4) ::
!>          The increment between elements of X. INCX > 0.
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is COMPLEX*16
!>          The value tau.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup complex16OTHERauxiliary
!
!  =====================================================================
      SUBROUTINE ZLARFG( N, ALPHA, X, INCX, TAU )
!
!  -- LAPACK auxiliary routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes
      integer(i4)  :: INCX, N
      complex(dpc) :: ALPHA, TAU
!     ..
!     .. Array Arguments ..
      complex(dpc) :: X( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      real(dp) ::   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      integer(i4)  :: J, KNT
      real(dp)     :: ALPHI, ALPHR, BETA, RSAFMN, SAFMIN, XNORM
!     ..
!     .. External Functions ..
      real(dp)     :: DLAMCH, DLAPY3, DZNRM2
      complex(dpc) :: ZLADIV
      EXTERNAL        DLAMCH, DLAPY3, DZNRM2, ZLADIV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC       ABS, DBLE, DCMPLX, DIMAG, SIGN
!     ..
!     .. External Subroutines ..
      EXTERNAL        ZDSCAL, ZSCAL
!     ..
!     .. Executable Statements ..
!
      IF( N.LE.0 ) THEN
         TAU = ZERO
         RETURN
      END IF
!
      XNORM = DZNRM2( N-1, X, INCX )
      ALPHR = DBLE( ALPHA )
      ALPHI = DIMAG( ALPHA )
!
      IF( XNORM.EQ.ZERO .AND. ALPHI.EQ.ZERO ) THEN
!
!        H  =  I
!
         TAU = ZERO
      ELSE
!
!        general case
!
         BETA = -SIGN( DLAPY3( ALPHR, ALPHI, XNORM ), ALPHR )
         SAFMIN = DLAMCH( 'S' ) / DLAMCH( 'E' )
         RSAFMN = ONE / SAFMIN
!
         KNT = 0
         IF( ABS( BETA ).LT.SAFMIN ) THEN
!
!           XNORM, BETA may be inaccurate; scale X and recompute them
!
   10       CONTINUE
            KNT = KNT + 1
            CALL ZDSCAL( N-1, RSAFMN, X, INCX )
            BETA = BETA*RSAFMN
            ALPHI = ALPHI*RSAFMN
            ALPHR = ALPHR*RSAFMN
            IF( ABS( BETA ).LT.SAFMIN ) GO TO 10
!
!           New BETA is at most 1, at least SAFMIN
!
            XNORM = DZNRM2( N-1, X, INCX )
            ALPHA = DCMPLX( ALPHR, ALPHI )
            BETA = -SIGN( DLAPY3( ALPHR, ALPHI, XNORM ), ALPHR )
         END IF
         TAU = DCMPLX( ( BETA-ALPHR ) / BETA, -ALPHI / BETA )
         ALPHA = ZLADIV( DCMPLX( ONE ), ALPHA-BETA )
         CALL ZSCAL( N-1, ALPHA, X, INCX )
!
!        If ALPHA is subnormal, it may lose relative accuracy
!
         DO 20 J = 1, KNT
            BETA = BETA*SAFMIN
 20      CONTINUE
         ALPHA = BETA
      END IF
!
      RETURN
!
!     End of ZLARFG
!
      END
!> \brief \b ZLARFT
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download ZLARFT + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlarft.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlarft.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlarft.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
! 
!       .. Scalar Arguments ..
!       CHARACTER          DIRECT, STOREV
!       integer(i4) ::            K, LDT, LDV, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         T( LDT, * ), TAU( * ), V( LDV, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLARFT forms the triangular factor T of a complex block reflector H
!> of order n, which is defined as a product of k elementary reflectors.
!>
!> If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
!>
!> If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
!>
!> If STOREV = 'C', the vector which defines the elementary reflector
!> H(i) is stored in the i-th column of the array V, and
!>
!>    H  =  I - V * T * V**H
!>
!> If STOREV = 'R', the vector which defines the elementary reflector
!> H(i) is stored in the i-th row of the array V, and
!>
!>    H  =  I - V**H * T * V
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] DIRECT
!> \verbatim
!>          DIRECT is CHARACTER*1
!>          Specifies the order in which the elementary reflectors are
!>          multiplied to form the block reflector:
!>          = 'F': H = H(1) H(2) . . . H(k) (Forward)
!>          = 'B': H = H(k) . . . H(2) H(1) (Backward)
!> \endverbatim
!>
!> \param[in] STOREV
!> \verbatim
!>          STOREV is CHARACTER*1
!>          Specifies how the vectors which define the elementary
!>          reflectors are stored (see also Further Details):
!>          = 'C': columnwise
!>          = 'R': rowwise
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is integer(i4) ::
!>          The order of the block reflector H. N >= 0.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is integer(i4) ::
!>          The order of the triangular factor T (= the number of
!>          elementary reflectors). K >= 1.
!> \endverbatim
!>
!> \param[in,out] V
!> \verbatim
!>          V is COMPLEX*16 array, dimension
!>                               (LDV,K) if STOREV = 'C'
!>                               (LDV,N) if STOREV = 'R'
!>          The matrix V. See further details.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is integer(i4) ::
!>          The leading dimension of the array V.
!>          If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is COMPLEX*16 array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i).
!> \endverbatim
!>
!> \param[out] T
!> \verbatim
!>          T is COMPLEX*16 array, dimension (LDT,K)
!>          The k by k triangular factor T of the block reflector.
!>          If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is
!>          lower triangular. The rest of the array is not used.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is integer(i4) ::
!>          The leading dimension of the array T. LDT >= K.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup complex16OTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The shape of the matrix V and the storage of the vectors which define
!>  the H(i) is best illustrated by the following example with n = 5 and
!>  k = 3. The elements equal to 1 are not stored; the corresponding
!>  array elements are modified but restored on exit. The rest of the
!>  array is not used.
!>
!>  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
!>
!>               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
!>                   ( v1  1    )                     (     1 v2 v2 v2 )
!>                   ( v1 v2  1 )                     (        1 v3 v3 )
!>                   ( v1 v2 v3 )
!>                   ( v1 v2 v3 )
!>
!>  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
!>
!>               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
!>                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
!>                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
!>                   (     1 v3 )
!>                   (        1 )
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
!
!  -- LAPACK auxiliary routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes
      CHARACTER    :: DIRECT, STOREV
      integer(i4)  :: K, LDT, LDV, N
!     ..
!     .. Array Arguments ..
      complex(dpc) :: T( LDT, * ), TAU( * ), V( LDV, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      complex(dpc) :: ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ), ZERO = ( 0.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      integer(i4)  :: I, J, PREVLASTV, LASTV
      complex(dpc) :: VII
!     ..
!     .. External Subroutines ..
      EXTERNAL           ZGEMV, ZLACGV, ZTRMV
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF( N.EQ.0 ) RETURN
!
      IF( LSAME( DIRECT, 'F' ) ) THEN
         PREVLASTV = N
         DO 20 I = 1, K
            PREVLASTV = MAX( PREVLASTV, I )
            IF( TAU( I ).EQ.ZERO ) THEN
!
!              H(i)  =  I
!
               DO 10 J = 1, I
                  T( J, I ) = ZERO
   10          CONTINUE
            ELSE
!
!              general case
!
               VII = V( I, I )
               V( I, I ) = ONE
               IF( LSAME( STOREV, 'C' ) ) THEN
!                 Skip any trailing zeros.
                  DO LASTV = N, I+1, -1
                     IF( V( LASTV, I ).NE.ZERO ) EXIT
                  END DO
                  J = MIN( LASTV, PREVLASTV )
!
!                 T(1:i-1,i) := - tau(i) * V(i:j,1:i-1)**H * V(i:j,i)
!
                  CALL ZGEMV( 'Conjugate transpose', J-I+1, I-1, -TAU( I ), V( I, 1 ), LDV, V( I, I ), 1, &
                              ZERO, T( 1, I ), 1 )
               ELSE
!                 Skip any trailing zeros.
                  DO LASTV = N, I+1, -1
                     IF( V( I, LASTV ).NE.ZERO ) EXIT
                  END DO
                  J = MIN( LASTV, PREVLASTV )
!
!                 T(1:i-1,i) := - tau(i) * V(1:i-1,i:j) * V(i,i:j)**H
!
                  IF( I.LT.J ) CALL ZLACGV( J-I, V( I, I+1 ), LDV )
                  CALL ZGEMV( 'No transpose', I-1, J-I+1, -TAU( I ), &
                              V( 1, I ), LDV, V( I, I ), LDV, ZERO, T( 1, I ), 1 )
                  IF( I.LT.J ) CALL ZLACGV( J-I, V( I, I+1 ), LDV )
               END IF
               V( I, I ) = VII
!
!              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)
!
               CALL ZTRMV( 'Upper', 'No transpose', 'Non-unit', I-1, T, LDT, T( 1, I ), 1 )
               T( I, I ) = TAU( I )
               IF( I.GT.1 ) THEN
                  PREVLASTV = MAX( PREVLASTV, LASTV )
               ELSE
                  PREVLASTV = LASTV
               END IF
             END IF
   20    CONTINUE
      ELSE
         PREVLASTV = 1
         DO 40 I = K, 1, -1
            IF( TAU( I ).EQ.ZERO ) THEN
!
!              H(i)  =  I
!
               DO 30 J = I, K
                  T( J, I ) = ZERO
   30          CONTINUE
            ELSE
!
!              general case
!
               IF( I.LT.K ) THEN
                  IF( LSAME( STOREV, 'C' ) ) THEN
                     VII = V( N-K+I, I )
                     V( N-K+I, I ) = ONE
!                    Skip any leading zeros.
                     DO LASTV = 1, I-1
                        IF( V( LASTV, I ).NE.ZERO ) EXIT
                     END DO
                     J = MAX( LASTV, PREVLASTV )
!
!                    T(i+1:k,i) :=
!                            - tau(i) * V(j:n-k+i,i+1:k)**H * V(j:n-k+i,i)
!
                     CALL ZGEMV( 'Conjugate transpose', N-K+I-J+1, K-I, -TAU( I ), V( J, I+1 ), LDV, V( J, I ), &
                                 1, ZERO, T( I+1, I ), 1 )
                     V( N-K+I, I ) = VII
                  ELSE
                     VII = V( I, N-K+I )
                     V( I, N-K+I ) = ONE
!                    Skip any leading zeros.
                     DO LASTV = 1, I-1
                        IF( V( I, LASTV ).NE.ZERO ) EXIT
                     END DO
                     J = MAX( LASTV, PREVLASTV )
!
!                    T(i+1:k,i) :=
!                            - tau(i) * V(i+1:k,j:n-k+i) * V(i,j:n-k+i)**H
!
                     CALL ZLACGV( N-K+I-1-J+1, V( I, J ), LDV )
                     CALL ZGEMV( 'No transpose', K-I, N-K+I-J+1, -TAU( I ), V( I+1, J ), LDV, V( I, J ), LDV, &
                          ZERO, T( I+1, I ), 1 )
                     CALL ZLACGV( N-K+I-1-J+1, V( I, J ), LDV )
                     V( I, N-K+I ) = VII
                  END IF
!
!                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)
!
                  CALL ZTRMV( 'Lower', 'No transpose', 'Non-unit', K-I, T( I+1, I+1 ), LDT, T( I+1, I ), 1 )
                  IF( I.GT.1 ) THEN
                     PREVLASTV = MIN( PREVLASTV, LASTV )
                  ELSE
                     PREVLASTV = LASTV
                  END IF
               END IF
               T( I, I ) = TAU( I )
            END IF
   40    CONTINUE
      END IF
      RETURN
!
!     End of ZLARFT
!
      END
!> \brief \b ZLASCL
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download ZLASCL + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlascl.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlascl.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlascl.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
! 
!       .. Scalar Arguments ..
!       CHARACTER          TYPE
!       integer(i4) ::            INFO, KL, KU, LDA, M, N
!       real(dp) ::   CFROM, CTO
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLASCL multiplies the M by N complex matrix A by the real scalar
!> CTO/CFROM.  This is done without over/underflow as long as the final
!> result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
!> A may be full, upper triangular, lower triangular, upper Hessenberg,
!> or banded.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TYPE
!> \verbatim
!>          TYPE is CHARACTER*1
!>          TYPE indices the storage type of the input matrix.
!>          = 'G':  A is a full matrix.
!>          = 'L':  A is a lower triangular matrix.
!>          = 'U':  A is an upper triangular matrix.
!>          = 'H':  A is an upper Hessenberg matrix.
!>          = 'B':  A is a symmetric band matrix with lower bandwidth KL
!>                  and upper bandwidth KU and with the only the lower
!>                  half stored.
!>          = 'Q':  A is a symmetric band matrix with lower bandwidth KL
!>                  and upper bandwidth KU and with the only the upper
!>                  half stored.
!>          = 'Z':  A is a band matrix with lower bandwidth KL and upper
!>                  bandwidth KU. See ZGBTRF for storage details.
!> \endverbatim
!>
!> \param[in] KL
!> \verbatim
!>          KL is integer(i4) ::
!>          The lower bandwidth of A.  Referenced only if TYPE = 'B',
!>          'Q' or 'Z'.
!> \endverbatim
!>
!> \param[in] KU
!> \verbatim
!>          KU is integer(i4) ::
!>          The upper bandwidth of A.  Referenced only if TYPE = 'B',
!>          'Q' or 'Z'.
!> \endverbatim
!>
!> \param[in] CFROM
!> \verbatim
!>          CFROM is real(dp) ::
!> \endverbatim
!>
!> \param[in] CTO
!> \verbatim
!>          CTO is real(dp) ::
!>
!>          The matrix A is multiplied by CTO/CFROM. A(I,J) is computed
!>          without over/underflow if the final result CTO*A(I,J)/CFROM
!>          can be represented without over/underflow.  CFROM must be
!>          nonzero.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is integer(i4) ::
!>          The number of rows of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is integer(i4) ::
!>          The number of columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          The matrix to be multiplied by CTO/CFROM.  See TYPE for the
!>          storage type.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is integer(i4) ::
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is integer(i4) ::
!>          0  - successful exit
!>          <0 - if INFO = -i, the i-th argument had an illegal value.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup complex16OTHERauxiliary
!
!  =====================================================================
      SUBROUTINE ZLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
!
!  -- LAPACK auxiliary routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes
      CHARACTER    :: TYPE
      integer(i4)  :: INFO, KL, KU, LDA, M, N
      real(dp)     :: CFROM, CTO
!     ..
!     .. Array Arguments ..
      complex(dpc) :: A( LDA, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      real(dp) ::   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
      LOGICAL      :: DONE
      integer(i4)  :: I, ITYPE, J, K1, K2, K3, K4
      real(dp)     :: BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM
!     ..
!     .. External Functions ..
      LOGICAL      :: LSAME, DISNAN
      real(dp)     :: DLAMCH
      EXTERNAL        LSAME, DLAMCH, DISNAN
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC       ABS, MAX, MIN
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
!
      IF( LSAME( TYPE, 'G' ) ) THEN
         ITYPE = 0
      ELSE IF( LSAME( TYPE, 'L' ) ) THEN
         ITYPE = 1
      ELSE IF( LSAME( TYPE, 'U' ) ) THEN
         ITYPE = 2
      ELSE IF( LSAME( TYPE, 'H' ) ) THEN
         ITYPE = 3
      ELSE IF( LSAME( TYPE, 'B' ) ) THEN
         ITYPE = 4
      ELSE IF( LSAME( TYPE, 'Q' ) ) THEN
         ITYPE = 5
      ELSE IF( LSAME( TYPE, 'Z' ) ) THEN
         ITYPE = 6
      ELSE
         ITYPE = -1
      END IF
!
      IF( ITYPE.EQ.-1 ) THEN
         INFO = -1
      ELSE IF( CFROM.EQ.ZERO .OR. DISNAN(CFROM) ) THEN
         INFO = -4
      ELSE IF( DISNAN(CTO) ) THEN
         INFO = -5
      ELSE IF( M.LT.0 ) THEN
         INFO = -6
      ELSE IF( N.LT.0 .OR. ( ITYPE.EQ.4 .AND. N.NE.M ) .OR. ( ITYPE.EQ.5 .AND. N.NE.M ) ) THEN
         INFO = -7
      ELSE IF( ITYPE.LE.3 .AND. LDA.LT.MAX( 1, M ) ) THEN
         INFO = -9
      ELSE IF( ITYPE.GE.4 ) THEN
         IF( KL.LT.0 .OR. KL.GT.MAX( M-1, 0 ) ) THEN
            INFO = -2
         ELSE IF( KU.LT.0 .OR. KU.GT.MAX( N-1, 0 ) .OR. &
                  ( ( ITYPE.EQ.4 .OR. ITYPE.EQ.5 ) .AND. KL.NE.KU ) ) THEN
            INFO = -3
         ELSE IF( ( ITYPE.EQ.4 .AND. LDA.LT.KL+1 ) .OR. &
                  ( ITYPE.EQ.5 .AND. LDA.LT.KU+1 ) .OR. &
                  ( ITYPE.EQ.6 .AND. LDA.LT.2*KL+KU+1 ) ) THEN
            INFO = -9
         END IF
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZLASCL', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 .OR. M.EQ.0 ) RETURN
!
!     Get machine parameters
!
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
!
      CFROMC = CFROM
      CTOC = CTO
!
   10 CONTINUE
      CFROM1 = CFROMC*SMLNUM
      IF( CFROM1.EQ.CFROMC ) THEN
!        CFROMC is an inf.  Multiply by a correctly signed zero for
!        finite CTOC, or a NaN if CTOC is infinite.
         MUL = CTOC / CFROMC
         DONE = .TRUE.
         CTO1 = CTOC
      ELSE
         CTO1 = CTOC / BIGNUM
         IF( CTO1.EQ.CTOC ) THEN
!           CTOC is either 0 or an inf.  In both cases, CTOC itself
!           serves as the correct multiplication factor.
            MUL = CTOC
            DONE = .TRUE.
            CFROMC = ONE
         ELSE IF( ABS( CFROM1 ).GT.ABS( CTOC ) .AND. CTOC.NE.ZERO ) THEN
            MUL = SMLNUM
            DONE = .FALSE.
            CFROMC = CFROM1
         ELSE IF( ABS( CTO1 ).GT.ABS( CFROMC ) ) THEN
            MUL = BIGNUM
            DONE = .FALSE.
            CTOC = CTO1
         ELSE
            MUL = CTOC / CFROMC
            DONE = .TRUE.
         END IF
      END IF
!
      IF( ITYPE.EQ.0 ) THEN
!
!        Full matrix
!
         DO 30 J = 1, N
            DO 20 I = 1, M
               A( I, J ) = A( I, J )*MUL
   20       CONTINUE
   30    CONTINUE
!
      ELSE IF( ITYPE.EQ.1 ) THEN
!
!        Lower triangular matrix
!
         DO 50 J = 1, N
            DO 40 I = J, M
               A( I, J ) = A( I, J )*MUL
   40       CONTINUE
   50    CONTINUE
!
      ELSE IF( ITYPE.EQ.2 ) THEN
!
!        Upper triangular matrix
!
         DO 70 J = 1, N
            DO 60 I = 1, MIN( J, M )
               A( I, J ) = A( I, J )*MUL
   60       CONTINUE
   70    CONTINUE
!
      ELSE IF( ITYPE.EQ.3 ) THEN
!
!        Upper Hessenberg matrix
!
         DO 90 J = 1, N
            DO 80 I = 1, MIN( J+1, M )
               A( I, J ) = A( I, J )*MUL
   80       CONTINUE
   90    CONTINUE
!
      ELSE IF( ITYPE.EQ.4 ) THEN
!
!        Lower half of a symmetric band matrix
!
         K3 = KL + 1
         K4 = N + 1
         DO 110 J = 1, N
            DO 100 I = 1, MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  100       CONTINUE
  110    CONTINUE
!
      ELSE IF( ITYPE.EQ.5 ) THEN
!
!        Upper half of a symmetric band matrix
!
         K1 = KU + 2
         K3 = KU + 1
         DO 130 J = 1, N
            DO 120 I = MAX( K1-J, 1 ), K3
               A( I, J ) = A( I, J )*MUL
  120       CONTINUE
  130    CONTINUE
!
      ELSE IF( ITYPE.EQ.6 ) THEN
!
!        Band matrix
!
         K1 = KL + KU + 2
         K2 = KL + 1
         K3 = 2*KL + KU + 1
         K4 = KL + KU + 1 + M
         DO 150 J = 1, N
            DO 140 I = MAX( K1-J, K2 ), MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  140       CONTINUE
  150    CONTINUE
!
      END IF
!
      IF( .NOT.DONE ) GO TO 10
!
      RETURN
!
!     End of ZLASCL
!
      END
!> \brief \b ZLASET
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download ZLASET + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaset.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaset.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaset.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
! 
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       integer(i4) ::            LDA, M, N
!       COMPLEX*16         ALPHA, BETA
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLASET initializes a 2-D array A to BETA on the diagonal and
!> ALPHA on the offdiagonals.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies the part of the matrix A to be set.
!>          = 'U':      Upper triangular part is set. The lower triangle
!>                      is unchanged.
!>          = 'L':      Lower triangular part is set. The upper triangle
!>                      is unchanged.
!>          Otherwise:  All of the matrix A is set.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is integer(i4) ::
!>          On entry, M specifies the number of rows of A.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is integer(i4) ::
!>          On entry, N specifies the number of columns of A.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is COMPLEX*16
!>          All the offdiagonal array elements are set to ALPHA.
!> \endverbatim
!>
!> \param[in] BETA
!> \verbatim
!>          BETA is COMPLEX*16
!>          All the diagonal array elements are set to BETA.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the m by n matrix A.
!>          On exit, A(i,j) = ALPHA, 1 <= i <= m, 1 <= j <= n, i.ne.j;
!>                   A(i,i) = BETA , 1 <= i <= min(m,n)
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is integer(i4) ::
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup complex16OTHERauxiliary
!
!  =====================================================================
      SUBROUTINE ZLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
!
!  -- LAPACK auxiliary routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes
      CHARACTER    :: UPLO
      integer(i4)  :: LDA, M, N
      complex(dpc) :: ALPHA, BETA
!     ..
!     .. Array Arguments ..
      complex(dpc) :: A( LDA, * )
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      integer(i4)  :: I, J
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MIN
!     ..
!     .. Executable Statements ..
!
      IF( LSAME( UPLO, 'U' ) ) THEN
!
!        Set the diagonal to BETA and the strictly upper triangular
!        part of the array to ALPHA.
!
         DO 20 J = 2, N
            DO 10 I = 1, MIN( J-1, M )
               A( I, J ) = ALPHA
   10       CONTINUE
   20    CONTINUE
         DO 30 I = 1, MIN( N, M )
            A( I, I ) = BETA
   30    CONTINUE
!
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
!
!        Set the diagonal to BETA and the strictly lower triangular
!        part of the array to ALPHA.
!
         DO 50 J = 1, MIN( M, N )
            DO 40 I = J + 1, M
               A( I, J ) = ALPHA
   40       CONTINUE
   50    CONTINUE
         DO 60 I = 1, MIN( N, M )
            A( I, I ) = BETA
   60    CONTINUE
!
      ELSE
!
!        Set the array to BETA on the diagonal and ALPHA on the
!        offdiagonal.
!
         DO 80 J = 1, N
            DO 70 I = 1, M
               A( I, J ) = ALPHA
   70       CONTINUE
   80    CONTINUE
         DO 90 I = 1, MIN( M, N )
            A( I, I ) = BETA
   90    CONTINUE
      END IF
!
      RETURN
!
!     End of ZLASET
!
      END
!> \brief \b ZLASR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download ZLASR + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlasr.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlasr.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlasr.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )
! 
!       .. Scalar Arguments ..
!       CHARACTER          DIRECT, PIVOT, SIDE
!       integer(i4) ::            LDA, M, N
!       ..
!       .. Array Arguments ..
!       real(dp) ::   C( * ), S( * )
!       COMPLEX*16         A( LDA, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLASR applies a sequence of real plane rotations to a complex matrix
!> A, from either the left or the right.
!>
!> When SIDE = 'L', the transformation takes the form
!>
!>    A := P*A
!>
!> and when SIDE = 'R', the transformation takes the form
!>
!>    A := A*P**T
!>
!> where P is an orthogonal matrix consisting of a sequence of z plane
!> rotations, with z = M when SIDE = 'L' and z = N when SIDE = 'R',
!> and P**T is the transpose of P.
!> 
!> When DIRECT = 'F' (Forward sequence), then
!> 
!>    P = P(z-1) * ... * P(2) * P(1)
!> 
!> and when DIRECT = 'B' (Backward sequence), then
!> 
!>    P = P(1) * P(2) * ... * P(z-1)
!> 
!> where P(k) is a plane rotation matrix defined by the 2-by-2 rotation
!> 
!>    R(k) = (  c(k)  s(k) )
!>         = ( -s(k)  c(k) ).
!> 
!> When PIVOT = 'V' (Variable pivot), the rotation is performed
!> for the plane (k,k+1), i.e., P(k) has the form
!> 
!>    P(k) = (  1                                            )
!>           (       ...                                     )
!>           (              1                                )
!>           (                   c(k)  s(k)                  )
!>           (                  -s(k)  c(k)                  )
!>           (                                1              )
!>           (                                     ...       )
!>           (                                            1  )
!> 
!> where R(k) appears as a rank-2 modification to the identity matrix in
!> rows and columns k and k+1.
!> 
!> When PIVOT = 'T' (Top pivot), the rotation is performed for the
!> plane (1,k+1), so P(k) has the form
!> 
!>    P(k) = (  c(k)                    s(k)                 )
!>           (         1                                     )
!>           (              ...                              )
!>           (                     1                         )
!>           ( -s(k)                    c(k)                 )
!>           (                                 1             )
!>           (                                      ...      )
!>           (                                             1 )
!> 
!> where R(k) appears in rows and columns 1 and k+1.
!> 
!> Similarly, when PIVOT = 'B' (Bottom pivot), the rotation is
!> performed for the plane (k,z), giving P(k) the form
!> 
!>    P(k) = ( 1                                             )
!>           (      ...                                      )
!>           (             1                                 )
!>           (                  c(k)                    s(k) )
!>           (                         1                     )
!>           (                              ...              )
!>           (                                     1         )
!>           (                 -s(k)                    c(k) )
!> 
!> where R(k) appears in rows and columns k and z.  The rotations are
!> performed without ever forming P(k) explicitly.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          Specifies whether the plane rotation matrix P is applied to
!>          A on the left or the right.
!>          = 'L':  Left, compute A := P*A
!>          = 'R':  Right, compute A:= A*P**T
!> \endverbatim
!>
!> \param[in] PIVOT
!> \verbatim
!>          PIVOT is CHARACTER*1
!>          Specifies the plane for which P(k) is a plane rotation
!>          matrix.
!>          = 'V':  Variable pivot, the plane (k,k+1)
!>          = 'T':  Top pivot, the plane (1,k+1)
!>          = 'B':  Bottom pivot, the plane (k,z)
!> \endverbatim
!>
!> \param[in] DIRECT
!> \verbatim
!>          DIRECT is CHARACTER*1
!>          Specifies whether P is a forward or backward sequence of
!>          plane rotations.
!>          = 'F':  Forward, P = P(z-1)*...*P(2)*P(1)
!>          = 'B':  Backward, P = P(1)*P(2)*...*P(z-1)
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is integer(i4) ::
!>          The number of rows of the matrix A.  If m <= 1, an immediate
!>          return is effected.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is integer(i4) ::
!>          The number of columns of the matrix A.  If n <= 1, an
!>          immediate return is effected.
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is real(dp) :: array, dimension
!>                  (M-1) if SIDE = 'L'
!>                  (N-1) if SIDE = 'R'
!>          The cosines c(k) of the plane rotations.
!> \endverbatim
!>
!> \param[in] S
!> \verbatim
!>          S is real(dp) :: array, dimension
!>                  (M-1) if SIDE = 'L'
!>                  (N-1) if SIDE = 'R'
!>          The sines s(k) of the plane rotations.  The 2-by-2 plane
!>          rotation part of the matrix P(k), R(k), has the form
!>          R(k) = (  c(k)  s(k) )
!>                 ( -s(k)  c(k) ).
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          The M-by-N matrix A.  On exit, A is overwritten by P*A if
!>          SIDE = 'R' or by A*P**T if SIDE = 'L'.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is integer(i4) ::
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup complex16OTHERauxiliary
!
!  =====================================================================
      SUBROUTINE ZLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )
!
!  -- LAPACK auxiliary routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes
      CHARACTER    :: DIRECT, PIVOT, SIDE
      integer(i4)  :: LDA, M, N
!     ..
!     .. Array Arguments ..
      real(dp)     :: C( * ), S( * )
      complex(dpc) :: A( LDA, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      real(dp)     :: ONE, ZERO
      PARAMETER       ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      integer(i4)  :: I, INFO, J
      real(dp)     :: CTEMP, STEMP
      complex(dpc) :: TEMP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      INFO = 0
      IF( .NOT.( LSAME( SIDE, 'L' ) .OR. LSAME( SIDE, 'R' ) ) ) THEN
         INFO = 1
      ELSE IF( .NOT.( LSAME( PIVOT, 'V' ) .OR. LSAME( PIVOT, 'T' ) .OR. LSAME( PIVOT, 'B' ) ) ) THEN
         INFO = 2
      ELSE IF( .NOT.( LSAME( DIRECT, 'F' ) .OR. LSAME( DIRECT, 'B' ) ) ) THEN
         INFO = 3
      ELSE IF( M.LT.0 ) THEN
         INFO = 4
      ELSE IF( N.LT.0 ) THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZLASR ', INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( ( M.EQ.0 ) .OR. ( N.EQ.0 ) ) RETURN
      IF( LSAME( SIDE, 'L' ) ) THEN
!
!        Form  P * A
!
         IF( LSAME( PIVOT, 'V' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 20 J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 10 I = 1, N
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
   10                CONTINUE
                  END IF
   20          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 40 J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 30 I = 1, N
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
   30                CONTINUE
                  END IF
   40          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 60 J = 2, M
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 50 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
   50                CONTINUE
                  END IF
   60          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 80 J = M, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 70 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
   70                CONTINUE
                  END IF
   80          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 100 J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 90 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
   90                CONTINUE
                  END IF
  100          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 120 J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 110 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
  110                CONTINUE
                  END IF
  120          CONTINUE
            END IF
         END IF
      ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!        Form A * P**T
!
         IF( LSAME( PIVOT, 'V' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 140 J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 130 I = 1, M
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  130                CONTINUE
                  END IF
  140          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 160 J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 150 I = 1, M
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  150                CONTINUE
                  END IF
  160          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 180 J = 2, N
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 170 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
  170                CONTINUE
                  END IF
  180          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 200 J = N, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 190 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
  190                CONTINUE
                  END IF
  200          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 220 J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 210 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
  210                CONTINUE
                  END IF
  220          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 240 J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 230 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
  230                CONTINUE
                  END IF
  240          CONTINUE
            END IF
         END IF
      END IF
!
      RETURN
!
!     End of ZLASR
!
      END
!> \brief \b ZLASSQ
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download ZLASSQ + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlassq.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlassq.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlassq.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLASSQ( N, X, INCX, SCALE, SUMSQ )
! 
!       .. Scalar Arguments ..
!       integer(i4) ::            INCX, N
!       real(dp) ::   SCALE, SUMSQ
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         X( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLASSQ returns the values scl and ssq such that
!>
!>    ( scl**2 )*ssq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
!>
!> where x( i ) = abs( X( 1 + ( i - 1 )*INCX ) ). The value of sumsq is
!> assumed to be at least unity and the value of ssq will then satisfy
!>
!>    1.0 .le. ssq .le. ( sumsq + 2*n ).
!>
!> scale is assumed to be non-negative and scl returns the value
!>
!>    scl = max( scale, abs( real( x( i ) ) ), abs( aimag( x( i ) ) ) ),
!>           i
!>
!> scale and sumsq must be supplied in SCALE and SUMSQ respectively.
!> SCALE and SUMSQ are overwritten by scl and ssq respectively.
!>
!> The routine makes only one pass through the vector X.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is integer(i4) ::
!>          The number of elements to be used from the vector X.
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension (N)
!>          The vector x as described above.
!>             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is integer(i4) ::
!>          The increment between successive values of the vector X.
!>          INCX > 0.
!> \endverbatim
!>
!> \param[in,out] SCALE
!> \verbatim
!>          SCALE is real(dp) ::
!>          On entry, the value  scale  in the equation above.
!>          On exit, SCALE is overwritten with the value  scl .
!> \endverbatim
!>
!> \param[in,out] SUMSQ
!> \verbatim
!>          SUMSQ is real(dp) ::
!>          On entry, the value  sumsq  in the equation above.
!>          On exit, SUMSQ is overwritten with the value  ssq .
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup complex16OTHERauxiliary
!
!  =====================================================================
      SUBROUTINE ZLASSQ( N, X, INCX, SCALE, SUMSQ )
!
!  -- LAPACK auxiliary routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes
      integer(i4)  :: INCX, N
      real(dp)     :: SCALE, SUMSQ
!     ..
!     .. Array Arguments ..
      complex(dpc) :: X( * )
!     ..
!
! =====================================================================
!
!     .. Parameters ..
      real(dp)     :: ZERO
      PARAMETER       ( ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      integer(i4)  :: IX
      real(dp)     :: TEMP1
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DIMAG
!     ..
!     .. Executable Statements ..
!
      IF( N.GT.0 ) THEN
         DO 10 IX = 1, 1 + ( N-1 )*INCX, INCX
            IF( DBLE( X( IX ) ).NE.ZERO ) THEN
               TEMP1 = ABS( DBLE( X( IX ) ) )
               IF( SCALE.LT.TEMP1 ) THEN
                  SUMSQ = 1 + SUMSQ*( SCALE / TEMP1 )**2
                  SCALE = TEMP1
               ELSE
                  SUMSQ = SUMSQ + ( TEMP1 / SCALE )**2
               END IF
            END IF
            IF( DIMAG( X( IX ) ).NE.ZERO ) THEN
               TEMP1 = ABS( DIMAG( X( IX ) ) )
               IF( SCALE.LT.TEMP1 ) THEN
                  SUMSQ = 1 + SUMSQ*( SCALE / TEMP1 )**2
                  SCALE = TEMP1
               ELSE
                  SUMSQ = SUMSQ + ( TEMP1 / SCALE )**2
               END IF
            END IF
   10    CONTINUE
      END IF
!
      RETURN
!
!     End of ZLASSQ
!
      END
!> \brief \b ZLATRD
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download ZLATRD + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlatrd.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlatrd.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlatrd.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLATRD( UPLO, N, NB, A, LDA, E, TAU, W, LDW )
! 
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       integer(i4) ::            LDA, LDW, N, NB
!       ..
!       .. Array Arguments ..
!       real(dp) ::   E( * )
!       COMPLEX*16         A( LDA, * ), TAU( * ), W( LDW, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLATRD reduces NB rows and columns of a complex Hermitian matrix A to
!> Hermitian tridiagonal form by a unitary similarity
!> transformation Q**H * A * Q, and returns the matrices V and W which are
!> needed to apply the transformation to the unreduced part of A.
!>
!> If UPLO = 'U', ZLATRD reduces the last NB rows and columns of a
!> matrix, of which the upper triangle is supplied;
!> if UPLO = 'L', ZLATRD reduces the first NB rows and columns of a
!> matrix, of which the lower triangle is supplied.
!>
!> This is an auxiliary routine called by ZHETRD.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the upper or lower triangular part of the
!>          Hermitian matrix A is stored:
!>          = 'U': Upper triangular
!>          = 'L': Lower triangular
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is integer(i4) ::
!>          The order of the matrix A.
!> \endverbatim
!>
!> \param[in] NB
!> \verbatim
!>          NB is integer(i4) ::
!>          The number of rows and columns to be reduced.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
!>          n-by-n upper triangular part of A contains the upper
!>          triangular part of the matrix A, and the strictly lower
!>          triangular part of A is not referenced.  If UPLO = 'L', the
!>          leading n-by-n lower triangular part of A contains the lower
!>          triangular part of the matrix A, and the strictly upper
!>          triangular part of A is not referenced.
!>          On exit:
!>          if UPLO = 'U', the last NB columns have been reduced to
!>            tridiagonal form, with the diagonal elements overwriting
!>            the diagonal elements of A; the elements above the diagonal
!>            with the array TAU, represent the unitary matrix Q as a
!>            product of elementary reflectors;
!>          if UPLO = 'L', the first NB columns have been reduced to
!>            tridiagonal form, with the diagonal elements overwriting
!>            the diagonal elements of A; the elements below the diagonal
!>            with the array TAU, represent the  unitary matrix Q as a
!>            product of elementary reflectors.
!>          See Further Details.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is integer(i4) ::
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] E
!> \verbatim
!>          E is real(dp) :: array, dimension (N-1)
!>          If UPLO = 'U', E(n-nb:n-1) contains the superdiagonal
!>          elements of the last NB columns of the reduced matrix;
!>          if UPLO = 'L', E(1:nb) contains the subdiagonal elements of
!>          the first NB columns of the reduced matrix.
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is COMPLEX*16 array, dimension (N-1)
!>          The scalar factors of the elementary reflectors, stored in
!>          TAU(n-nb:n-1) if UPLO = 'U', and in TAU(1:nb) if UPLO = 'L'.
!>          See Further Details.
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is COMPLEX*16 array, dimension (LDW,NB)
!>          The n-by-nb matrix W required to update the unreduced part
!>          of A.
!> \endverbatim
!>
!> \param[in] LDW
!> \verbatim
!>          LDW is integer(i4) ::
!>          The leading dimension of the array W. LDW >= max(1,N).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup complex16OTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  If UPLO = 'U', the matrix Q is represented as a product of elementary
!>  reflectors
!>
!>     Q = H(n) H(n-1) . . . H(n-nb+1).
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - tau * v * v**H
!>
!>  where tau is a complex scalar, and v is a complex vector with
!>  v(i:n) = 0 and v(i-1) = 1; v(1:i-1) is stored on exit in A(1:i-1,i),
!>  and tau in TAU(i-1).
!>
!>  If UPLO = 'L', the matrix Q is represented as a product of elementary
!>  reflectors
!>
!>     Q = H(1) H(2) . . . H(nb).
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - tau * v * v**H
!>
!>  where tau is a complex scalar, and v is a complex vector with
!>  v(1:i) = 0 and v(i+1) = 1; v(i+1:n) is stored on exit in A(i+1:n,i),
!>  and tau in TAU(i).
!>
!>  The elements of the vectors v together form the n-by-nb matrix V
!>  which is needed, with W, to apply the transformation to the unreduced
!>  part of the matrix, using a Hermitian rank-2k update of the form:
!>  A := A - V*W**H - W*V**H.
!>
!>  The contents of A on exit are illustrated by the following examples
!>  with n = 5 and nb = 2:
!>
!>  if UPLO = 'U':                       if UPLO = 'L':
!>
!>    (  a   a   a   v4  v5 )              (  d                  )
!>    (      a   a   v4  v5 )              (  1   d              )
!>    (          a   1   v5 )              (  v1  1   a          )
!>    (              d   1  )              (  v1  v2  a   a      )
!>    (                  d  )              (  v1  v2  a   a   a  )
!>
!>  where d denotes a diagonal element of the reduced matrix, a denotes
!>  an element of the original matrix that is unchanged, and vi denotes
!>  an element of the vector defining H(i).
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZLATRD( UPLO, N, NB, A, LDA, E, TAU, W, LDW )
!
!  -- LAPACK auxiliary routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes
      CHARACTER    :: UPLO
      integer(i4)  :: LDA, LDW, N, NB
!     ..
!     .. Array Arguments ..
      real(dp)     :: E( * )
      complex(dpc) :: A( LDA, * ), TAU( * ), W( LDW, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      complex(dpc) :: ZERO, ONE, HALF
      PARAMETER       ( ZERO = ( 0.0D+0, 0.0D+0 ), ONE = ( 1.0D+0, 0.0D+0 ), &
                         HALF = ( 0.5D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      integer(i4)  :: I, IW
      complex(dpc) :: ALPHA
!     ..
!     .. External Subroutines ..
      EXTERNAL           ZAXPY, ZGEMV, ZHEMV, ZLACGV, ZLARFG, ZSCAL
!     ..
!     .. External Functions ..
      LOGICAL      :: LSAME
      complex(dpc) :: ZDOTC
      EXTERNAL        LSAME, ZDOTC
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MIN
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF( N.LE.0 ) RETURN
!
      IF( LSAME( UPLO, 'U' ) ) THEN
!
!        Reduce last NB columns of upper triangle
!
         DO 10 I = N, N - NB + 1, -1
            IW = I - N + NB
            IF( I.LT.N ) THEN
!
!              Update A(1:i,i)
!
               A( I, I ) = DBLE( A( I, I ) )
               CALL ZLACGV( N-I, W( I, IW+1 ), LDW )
               CALL ZGEMV( 'No transpose', I, N-I, -ONE, A( 1, I+1 ), LDA, W( I, IW+1 ), LDW, ONE, A( 1, I ), 1 )
               CALL ZLACGV( N-I, W( I, IW+1 ), LDW )
               CALL ZLACGV( N-I, A( I, I+1 ), LDA )
               CALL ZGEMV( 'No transpose', I, N-I, -ONE, W( 1, IW+1 ), LDW, A( I, I+1 ), LDA, ONE, A( 1, I ), 1 )
               CALL ZLACGV( N-I, A( I, I+1 ), LDA )
               A( I, I ) = DBLE( A( I, I ) )
            END IF
            IF( I.GT.1 ) THEN
!
!              Generate elementary reflector H(i) to annihilate
!              A(1:i-2,i)
!
               ALPHA = A( I-1, I )
               CALL ZLARFG( I-1, ALPHA, A( 1, I ), 1, TAU( I-1 ) )
               E( I-1 ) = ALPHA
               A( I-1, I ) = ONE
!
!              Compute W(1:i-1,i)
!
               CALL ZHEMV( 'Upper', I-1, ONE, A, LDA, A( 1, I ), 1, ZERO, W( 1, IW ), 1 )
               IF( I.LT.N ) THEN
                  CALL ZGEMV( 'Conjugate transpose', I-1, N-I, ONE, W( 1, IW+1 ), LDW, A( 1, I ), 1, ZERO, &
                              W( I+1, IW ), 1 )
                  CALL ZGEMV( 'No transpose', I-1, N-I, -ONE, A( 1, I+1 ), LDA, W( I+1, IW ), 1, ONE, &
                              W( 1, IW ), 1 )
                  CALL ZGEMV( 'Conjugate transpose', I-1, N-I, ONE, A( 1, I+1 ), LDA, A( 1, I ), 1, ZERO, &
                              W( I+1, IW ), 1 )
                  CALL ZGEMV( 'No transpose', I-1, N-I, -ONE, W( 1, IW+1 ), LDW, W( I+1, IW ), 1, ONE, &
                              W( 1, IW ), 1 )
               END IF
               CALL ZSCAL( I-1, TAU( I-1 ), W( 1, IW ), 1 )
               ALPHA = -HALF*TAU( I-1 )*ZDOTC( I-1, W( 1, IW ), 1, A( 1, I ), 1 )
               CALL ZAXPY( I-1, ALPHA, A( 1, I ), 1, W( 1, IW ), 1 )
            END IF
!
   10    CONTINUE
      ELSE
!
!        Reduce first NB columns of lower triangle
!
         DO 20 I = 1, NB
!
!           Update A(i:n,i)
!
            A( I, I ) = DBLE( A( I, I ) )
            CALL ZLACGV( I-1, W( I, 1 ), LDW )
            CALL ZGEMV( 'No transpose', N-I+1, I-1, -ONE, A( I, 1 ), LDA, W( I, 1 ), LDW, ONE, A( I, I ), 1 )
            CALL ZLACGV( I-1, W( I, 1 ), LDW )
            CALL ZLACGV( I-1, A( I, 1 ), LDA )
            CALL ZGEMV( 'No transpose', N-I+1, I-1, -ONE, W( I, 1 ), LDW, A( I, 1 ), LDA, ONE, A( I, I ), 1 )
            CALL ZLACGV( I-1, A( I, 1 ), LDA )
            A( I, I ) = DBLE( A( I, I ) )
            IF( I.LT.N ) THEN
!
!              Generate elementary reflector H(i) to annihilate
!              A(i+2:n,i)
!
               ALPHA = A( I+1, I )
               CALL ZLARFG( N-I, ALPHA, A( MIN( I+2, N ), I ), 1, TAU( I ) )
               E( I ) = ALPHA
               A( I+1, I ) = ONE
!
!              Compute W(i+1:n,i)
!
               CALL ZHEMV( 'Lower', N-I, ONE, A( I+1, I+1 ), LDA, A( I+1, I ), 1, ZERO, W( I+1, I ), 1 )
               CALL ZGEMV( 'Conjugate transpose', N-I, I-1, ONE, W( I+1, 1 ), LDW, A( I+1, I ), 1, ZERO, &
                           W( 1, I ), 1 )
               CALL ZGEMV( 'No transpose', N-I, I-1, -ONE, A( I+1, 1 ), LDA, W( 1, I ), 1, ONE, W( I+1, I ), 1 )
               CALL ZGEMV( 'Conjugate transpose', N-I, I-1, ONE, A( I+1, 1 ), LDA, A( I+1, I ), 1, ZERO, &
                           W( 1, I ), 1 )
               CALL ZGEMV( 'No transpose', N-I, I-1, -ONE, W( I+1, 1 ), LDW, W( 1, I ), 1, ONE, W( I+1, I ), 1 )
               CALL ZSCAL( N-I, TAU( I ), W( I+1, I ), 1 )
               ALPHA = -HALF*TAU( I )*ZDOTC( N-I, W( I+1, I ), 1, A( I+1, I ), 1 )
               CALL ZAXPY( N-I, ALPHA, A( I+1, I ), 1, W( I+1, I ), 1 )
            END IF
!
   20    CONTINUE
      END IF
!
      RETURN
!
!     End of ZLATRD
!
      END
!> \brief \b ZSTEQR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download ZSTEQR + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsteqr.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsteqr.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsteqr.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
! 
!       .. Scalar Arguments ..
!       CHARACTER          COMPZ
!       integer(i4) ::            INFO, LDZ, N
!       ..
!       .. Array Arguments ..
!       real(dp) ::   D( * ), E( * ), WORK( * )
!       COMPLEX*16         Z( LDZ, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZSTEQR computes all eigenvalues and, optionally, eigenvectors of a
!> symmetric tridiagonal matrix using the implicit QL or QR method.
!> The eigenvectors of a full or band complex Hermitian matrix can also
!> be found if ZHETRD or ZHPTRD or ZHBTRD has been used to reduce this
!> matrix to tridiagonal form.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] COMPZ
!> \verbatim
!>          COMPZ is CHARACTER*1
!>          = 'N':  Compute eigenvalues only.
!>          = 'V':  Compute eigenvalues and eigenvectors of the original
!>                  Hermitian matrix.  On entry, Z must contain the
!>                  unitary matrix used to reduce the original matrix
!>                  to tridiagonal form.
!>          = 'I':  Compute eigenvalues and eigenvectors of the
!>                  tridiagonal matrix.  Z is initialized to the identity
!>                  matrix.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is integer(i4) ::
!>          The order of the matrix.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is real(dp) :: array, dimension (N)
!>          On entry, the diagonal elements of the tridiagonal matrix.
!>          On exit, if INFO = 0, the eigenvalues in ascending order.
!> \endverbatim
!>
!> \param[in,out] E
!> \verbatim
!>          E is real(dp) :: array, dimension (N-1)
!>          On entry, the (n-1) subdiagonal elements of the tridiagonal
!>          matrix.
!>          On exit, E has been destroyed.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is COMPLEX*16 array, dimension (LDZ, N)
!>          On entry, if  COMPZ = 'V', then Z contains the unitary
!>          matrix used in the reduction to tridiagonal form.
!>          On exit, if INFO = 0, then if COMPZ = 'V', Z contains the
!>          orthonormal eigenvectors of the original Hermitian matrix,
!>          and if COMPZ = 'I', Z contains the orthonormal eigenvectors
!>          of the symmetric tridiagonal matrix.
!>          If COMPZ = 'N', then Z is not referenced.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is integer(i4) ::
!>          The leading dimension of the array Z.  LDZ >= 1, and if
!>          eigenvectors are desired, then  LDZ >= max(1,N).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is real(dp) :: array, dimension (max(1,2*N-2))
!>          If COMPZ = 'N', then WORK is not referenced.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is integer(i4) ::
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  the algorithm has failed to find all the eigenvalues in
!>                a total of 30*N iterations; if INFO = i, then i
!>                elements of E have not converged to zero; on exit, D
!>                and E contain the elements of a symmetric tridiagonal
!>                matrix which is unitarily similar to the original
!>                matrix.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup complex16OTHERcomputational
!
!  =====================================================================
      SUBROUTINE ZSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes
      CHARACTER    :: COMPZ
      integer(i4)  :: INFO, LDZ, N
!     ..
!     .. Array Arguments ..
      real(dp)     :: D( * ), E( * ), WORK( * )
      complex(dpc) :: Z( LDZ, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      real(dp)     :: ZERO, ONE, TWO, THREE
      PARAMETER       ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, THREE = 3.0D0 )
      complex(dpc) :: CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D0, 0.0D0 ), CONE = ( 1.0D0, 0.0D0 ) )
      integer(i4)  :: MAXIT
      PARAMETER          ( MAXIT = 30 )
!     ..
!     .. Local Scalars ..
      integer(i4)  :: I, ICOMPZ, II, ISCALE, J, JTOT, K, L, L1, LEND, &
                      LENDM1, LENDP1, LENDSV, LM1, LSV, M, MM, MM1, NM1, NMAXIT
      real(dp)  :: ANORM, B, C, EPS, EPS2, F, G, P, R, RT1, RT2, S, SAFMAX, SAFMIN, SSFMAX, SSFMIN, TST
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      real(dp) ::   DLAMCH, DLANST, DLAPY2
      EXTERNAL           LSAME, DLAMCH, DLANST, DLAPY2
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLAE2, DLAEV2, DLARTG, DLASCL, DLASRT, XERBLA, ZLASET, ZLASR, ZSWAP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SIGN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
!
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ICOMPZ = 0
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ICOMPZ = 2
      ELSE
         ICOMPZ = -1
      END IF
      IF( ICOMPZ.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( ( LDZ.LT.1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1, N ) ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZSTEQR', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) RETURN
!
      IF( N.EQ.1 ) THEN
         IF( ICOMPZ.EQ.2 ) Z( 1, 1 ) = CONE
         RETURN
      END IF
!
!     Determine the unit roundoff and over/underflow thresholds.
!
      EPS = DLAMCH( 'E' )
      EPS2 = EPS**2
      SAFMIN = DLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      SSFMAX = SQRT( SAFMAX ) / THREE
      SSFMIN = SQRT( SAFMIN ) / EPS2
!
!     Compute the eigenvalues and eigenvectors of the tridiagonal
!     matrix.
!
      IF( ICOMPZ.EQ.2 ) CALL ZLASET( 'Full', N, N, CZERO, CONE, Z, LDZ )
!
      NMAXIT = N*MAXIT
      JTOT = 0
!
!     Determine where the matrix splits and choose QL or QR iteration
!     for each block, according to whether top or bottom diagonal
!     element is smaller.
!
      L1 = 1
      NM1 = N - 1
!
   10 CONTINUE
      IF( L1.GT.N ) GO TO 160
      IF( L1.GT.1 ) E( L1-1 ) = ZERO
      IF( L1.LE.NM1 ) THEN
         DO 20 M = L1, NM1
            TST = ABS( E( M ) )
            IF( TST.EQ.ZERO ) GO TO 30
            IF( TST.LE.( SQRT( ABS( D( M ) ) )*SQRT( ABS( D( M+1 ) ) ) )*EPS ) THEN
               E( M ) = ZERO
               GO TO 30
            END IF
   20    CONTINUE
      END IF
      M = N
!
   30 CONTINUE
      L = L1
      LSV = L
      LEND = M
      LENDSV = LEND
      L1 = M + 1
      IF( LEND.EQ.L ) GO TO 10
!
!     Scale submatrix in rows and columns L to LEND
!
      ANORM = DLANST( 'I', LEND-L+1, D( L ), E( L ) )
      ISCALE = 0
      IF( ANORM.EQ.ZERO ) GO TO 10
      IF( ANORM.GT.SSFMAX ) THEN
         ISCALE = 1
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L+1, 1, D( L ), N, INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L, 1, E( L ), N, INFO )
      ELSE IF( ANORM.LT.SSFMIN ) THEN
         ISCALE = 2
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L+1, 1, D( L ), N, INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L, 1, E( L ), N, INFO )
      END IF
!
!     Choose between QL and QR iteration
!
      IF( ABS( D( LEND ) ).LT.ABS( D( L ) ) ) THEN
         LEND = LSV
         L = LENDSV
      END IF
!
      IF( LEND.GT.L ) THEN
!
!        QL Iteration
!
!        Look for small subdiagonal element.
!
   40    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDM1 = LEND - 1
            DO 50 M = L, LENDM1
               TST = ABS( E( M ) )**2
               IF( TST.LE.( EPS2*ABS( D( M ) ) )*ABS( D( M+1 ) )+ SAFMIN )GO TO 60
   50       CONTINUE
         END IF
!
         M = LEND
!
   60    CONTINUE
         IF( M.LT.LEND ) E( M ) = ZERO
         P = D( L )
         IF( M.EQ.L ) GO TO 80
!
!        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
!        to compute its eigensystem.
!
         IF( M.EQ.L+1 ) THEN
            IF( ICOMPZ.GT.0 ) THEN
               CALL DLAEV2( D( L ), E( L ), D( L+1 ), RT1, RT2, C, S )
               WORK( L ) = C
               WORK( N-1+L ) = S
               CALL ZLASR( 'R', 'V', 'B', N, 2, WORK( L ), WORK( N-1+L ), Z( 1, L ), LDZ )
            ELSE
               CALL DLAE2( D( L ), E( L ), D( L+1 ), RT1, RT2 )
            END IF
            D( L ) = RT1
            D( L+1 ) = RT2
            E( L ) = ZERO
            L = L + 2
            IF( L.LE.LEND ) GO TO 40
            GO TO 140
         END IF
!
         IF( JTOT.EQ.NMAXIT ) GO TO 140
         JTOT = JTOT + 1
!
!        Form shift.
!
         G = ( D( L+1 )-P ) / ( TWO*E( L ) )
         R = DLAPY2( G, ONE )
         G = D( M ) - P + ( E( L ) / ( G+SIGN( R, G ) ) )
!
         S = ONE
         C = ONE
         P = ZERO
!
!        Inner loop
!
         MM1 = M - 1
         DO 70 I = MM1, L, -1
            F = S*E( I )
            B = C*E( I )
            CALL DLARTG( G, F, C, S, R )
            IF( I.NE.M-1 ) E( I+1 ) = R
            G = D( I+1 ) - P
            R = ( D( I )-G )*S + TWO*C*B
            P = S*R
            D( I+1 ) = G + P
            G = C*R - B
!
!           If eigenvectors are desired, then save rotations.
!
            IF( ICOMPZ.GT.0 ) THEN
               WORK( I ) = C
               WORK( N-1+I ) = -S
            END IF
!
   70    CONTINUE
!
!        If eigenvectors are desired, then apply saved rotations.
!
         IF( ICOMPZ.GT.0 ) THEN
            MM = M - L + 1
            CALL ZLASR( 'R', 'V', 'B', N, MM, WORK( L ), WORK( N-1+L ), Z( 1, L ), LDZ )
         END IF
!
         D( L ) = D( L ) - P
         E( L ) = G
         GO TO 40
!
!        Eigenvalue found.
!
   80    CONTINUE
         D( L ) = P
!
         L = L + 1
         IF( L.LE.LEND ) GO TO 40
         GO TO 140
!
      ELSE
!
!        QR Iteration
!
!        Look for small superdiagonal element.
!
   90    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDP1 = LEND + 1
            DO 100 M = L, LENDP1, -1
               TST = ABS( E( M-1 ) )**2
               IF( TST.LE.( EPS2*ABS( D( M ) ) )*ABS( D( M-1 ) )+ SAFMIN )GO TO 110
  100       CONTINUE
         END IF
!
         M = LEND
!
  110    CONTINUE
         IF( M.GT.LEND ) E( M-1 ) = ZERO
         P = D( L )
         IF( M.EQ.L ) GO TO 130
!
!        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
!        to compute its eigensystem.
!
         IF( M.EQ.L-1 ) THEN
            IF( ICOMPZ.GT.0 ) THEN
               CALL DLAEV2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2, C, S )
               WORK( M ) = C
               WORK( N-1+M ) = S
               CALL ZLASR( 'R', 'V', 'F', N, 2, WORK( M ), WORK( N-1+M ), Z( 1, L-1 ), LDZ )
            ELSE
               CALL DLAE2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2 )
            END IF
            D( L-1 ) = RT1
            D( L ) = RT2
            E( L-1 ) = ZERO
            L = L - 2
            IF( L.GE.LEND ) GO TO 90
            GO TO 140
         END IF
!
         IF( JTOT.EQ.NMAXIT ) GO TO 140
         JTOT = JTOT + 1
!
!        Form shift.
!
         G = ( D( L-1 )-P ) / ( TWO*E( L-1 ) )
         R = DLAPY2( G, ONE )
         G = D( M ) - P + ( E( L-1 ) / ( G+SIGN( R, G ) ) )
!
         S = ONE
         C = ONE
         P = ZERO
!
!        Inner loop
!
         LM1 = L - 1
         DO 120 I = M, LM1
            F = S*E( I )
            B = C*E( I )
            CALL DLARTG( G, F, C, S, R )
            IF( I.NE.M ) E( I-1 ) = R
            G = D( I ) - P
            R = ( D( I+1 )-G )*S + TWO*C*B
            P = S*R
            D( I ) = G + P
            G = C*R - B
!
!           If eigenvectors are desired, then save rotations.
!
            IF( ICOMPZ.GT.0 ) THEN
               WORK( I ) = C
               WORK( N-1+I ) = S
            END IF
!
  120    CONTINUE
!
!        If eigenvectors are desired, then apply saved rotations.
!
         IF( ICOMPZ.GT.0 ) THEN
            MM = L - M + 1
            CALL ZLASR( 'R', 'V', 'F', N, MM, WORK( M ), WORK( N-1+M ), Z( 1, M ), LDZ )
         END IF
!
         D( L ) = D( L ) - P
         E( LM1 ) = G
         GO TO 90
!
!        Eigenvalue found.
!
  130    CONTINUE
         D( L ) = P
!
         L = L - 1
         IF( L.GE.LEND ) GO TO 90
         GO TO 140
!
      END IF
!
!     Undo scaling if necessary
!
  140 CONTINUE
      IF( ISCALE.EQ.1 ) THEN
         CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV+1, 1, D( LSV ), N, INFO )
         CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV, 1, E( LSV ), N, INFO )
      ELSE IF( ISCALE.EQ.2 ) THEN
         CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV+1, 1, D( LSV ), N, INFO )
         CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV, 1, E( LSV ), N, INFO )
      END IF
!
!     Check for no convergence to an eigenvalue after a total
!     of N*MAXIT iterations.
!
      IF( JTOT.EQ.NMAXIT ) THEN
         DO 150 I = 1, N - 1
            IF( E( I ).NE.ZERO ) INFO = INFO + 1
  150    CONTINUE
         RETURN
      END IF
      GO TO 10
!
!     Order eigenvalues and eigenvectors.
!
  160 CONTINUE
      IF( ICOMPZ.EQ.0 ) THEN
!
!        Use Quick Sort
!
         CALL DLASRT( 'I', N, D, INFO )
!
      ELSE
!
!        Use Selection Sort to minimize swaps of eigenvectors
!
         DO 180 II = 2, N
            I = II - 1
            K = I
            P = D( I )
            DO 170 J = II, N
               IF( D( J ).LT.P ) THEN
                  K = J
                  P = D( J )
               END IF
  170       CONTINUE
            IF( K.NE.I ) THEN
               D( K ) = D( I )
               D( I ) = P
               CALL ZSWAP( N, Z( 1, I ), 1, Z( 1, K ), 1 )
            END IF
  180    CONTINUE
      END IF
      RETURN
!
!     End of ZSTEQR
!
      END
!> \brief \b ZUNG2L
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download ZUNG2L + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zung2l.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zung2l.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zung2l.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZUNG2L( M, N, K, A, LDA, TAU, WORK, INFO )
! 
!       .. Scalar Arguments ..
!       integer(i4) ::            INFO, K, LDA, M, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZUNG2L generates an m by n complex matrix Q with orthonormal columns,
!> which is defined as the last n columns of a product of k elementary
!> reflectors of order m
!>
!>       Q  =  H(k) . . . H(2) H(1)
!>
!> as returned by ZGEQLF.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is integer(i4) ::
!>          The number of rows of the matrix Q. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is integer(i4) ::
!>          The number of columns of the matrix Q. M >= N >= 0.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is integer(i4) ::
!>          The number of elementary reflectors whose product defines the
!>          matrix Q. N >= K >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the (n-k+i)-th column must contain the vector which
!>          defines the elementary reflector H(i), for i = 1,2,...,k, as
!>          returned by ZGEQLF in the last k columns of its array
!>          argument A.
!>          On exit, the m-by-n matrix Q.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is integer(i4) ::
!>          The first dimension of the array A. LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is COMPLEX*16 array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by ZGEQLF.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is integer(i4) ::
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument has an illegal value
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup complex16OTHERcomputational
!
!  =====================================================================
      SUBROUTINE ZUNG2L( M, N, K, A, LDA, TAU, WORK, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes
      integer(i4)  :: INFO, K, LDA, M, N
!     ..
!     .. Array Arguments ..
      complex(dpc) :: A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      complex(dpc) :: ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ), ZERO = ( 0.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      integer(i4)  :: I, II, J, L
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARF, ZSCAL
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNG2L', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.LE.0 ) RETURN
!
!     Initialise columns 1:n-k to columns of the unit matrix
!
      DO 20 J = 1, N - K
         DO 10 L = 1, M
            A( L, J ) = ZERO
   10    CONTINUE
         A( M-N+J, J ) = ONE
   20 CONTINUE
!
      DO 40 I = 1, K
         II = N - K + I
!
!        Apply H(i) to A(1:m-k+i,1:n-k+i) from the left
!
         A( M-N+II, II ) = ONE
         CALL ZLARF( 'Left', M-N+II, II-1, A( 1, II ), 1, TAU( I ), A, LDA, WORK )
         CALL ZSCAL( M-N+II-1, -TAU( I ), A( 1, II ), 1 )
         A( M-N+II, II ) = ONE - TAU( I )
!
!        Set A(m-k+i+1:m,n-k+i) to zero
!
         DO 30 L = M - N + II + 1, M
            A( L, II ) = ZERO
   30    CONTINUE
   40 CONTINUE
      RETURN
!
!     End of ZUNG2L
!
      END
!> \brief \b ZUNG2R
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download ZUNG2R + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zung2r.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zung2r.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zung2r.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZUNG2R( M, N, K, A, LDA, TAU, WORK, INFO )
! 
!       .. Scalar Arguments ..
!       integer(i4) ::            INFO, K, LDA, M, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZUNG2R generates an m by n complex matrix Q with orthonormal columns,
!> which is defined as the first n columns of a product of k elementary
!> reflectors of order m
!>
!>       Q  =  H(1) H(2) . . . H(k)
!>
!> as returned by ZGEQRF.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is integer(i4) ::
!>          The number of rows of the matrix Q. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is integer(i4) ::
!>          The number of columns of the matrix Q. M >= N >= 0.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is integer(i4) ::
!>          The number of elementary reflectors whose product defines the
!>          matrix Q. N >= K >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the i-th column must contain the vector which
!>          defines the elementary reflector H(i), for i = 1,2,...,k, as
!>          returned by ZGEQRF in the first k columns of its array
!>          argument A.
!>          On exit, the m by n matrix Q.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is integer(i4) ::
!>          The first dimension of the array A. LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is COMPLEX*16 array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by ZGEQRF.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is integer(i4) ::
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument has an illegal value
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup complex16OTHERcomputational
!
!  =====================================================================
      SUBROUTINE ZUNG2R( M, N, K, A, LDA, TAU, WORK, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes
      integer(i4)  :: INFO, K, LDA, M, N
!     ..
!     .. Array Arguments ..
      complex(dpc) :: A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      complex(dpc) :: ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ), ZERO = ( 0.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      integer(i4)  :: I, J, L
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARF, ZSCAL
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNG2R', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.LE.0 ) RETURN
!
!     Initialise columns k+1:n to columns of the unit matrix
!
      DO 20 J = K + 1, N
         DO 10 L = 1, M
            A( L, J ) = ZERO
   10    CONTINUE
         A( J, J ) = ONE
   20 CONTINUE
!
      DO 40 I = K, 1, -1
!
!        Apply H(i) to A(i:m,i:n) from the left
!
         IF( I.LT.N ) THEN
            A( I, I ) = ONE
            CALL ZLARF( 'Left', M-I+1, N-I, A( I, I ), 1, TAU( I ), A( I, I+1 ), LDA, WORK )
         END IF
         IF( I.LT.M ) CALL ZSCAL( M-I, -TAU( I ), A( I+1, I ), 1 )
         A( I, I ) = ONE - TAU( I )
!
!        Set A(1:i-1,i) to zero
!
         DO 30 L = 1, I - 1
            A( L, I ) = ZERO
   30    CONTINUE
   40 CONTINUE
      RETURN
!
!     End of ZUNG2R
!
      END
!> \brief \b ZUNGQL
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download ZUNGQL + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zungql.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zungql.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zungql.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZUNGQL( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
! 
!       .. Scalar Arguments ..
!       integer(i4) ::            INFO, K, LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZUNGQL generates an M-by-N complex matrix Q with orthonormal columns,
!> which is defined as the last N columns of a product of K elementary
!> reflectors of order M
!>
!>       Q  =  H(k) . . . H(2) H(1)
!>
!> as returned by ZGEQLF.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is integer(i4) ::
!>          The number of rows of the matrix Q. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is integer(i4) ::
!>          The number of columns of the matrix Q. M >= N >= 0.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is integer(i4) ::
!>          The number of elementary reflectors whose product defines the
!>          matrix Q. N >= K >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the (n-k+i)-th column must contain the vector which
!>          defines the elementary reflector H(i), for i = 1,2,...,k, as
!>          returned by ZGEQLF in the last k columns of its array
!>          argument A.
!>          On exit, the M-by-N matrix Q.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is integer(i4) ::
!>          The first dimension of the array A. LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is COMPLEX*16 array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by ZGEQLF.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is integer(i4) ::
!>          The dimension of the array WORK. LWORK >= max(1,N).
!>          For optimum performance LWORK >= N*NB, where NB is the
!>          optimal blocksize.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is integer(i4) ::
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument has an illegal value
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup complex16OTHERcomputational
!
!  =====================================================================
      SUBROUTINE ZUNGQL( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes
      integer(i4)  :: INFO, K, LDA, LWORK, M, N
!     ..
!     .. Array Arguments ..
      complex(dpc) :: A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      complex(dpc) :: ZERO
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY
      integer(i4) :: I, IB, IINFO, IWS, J, KK, L, LDWORK, LWKOPT, NB, NBMIN, NX
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARFB, ZLARFT, ZUNG2L
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. External Functions ..
      integer(i4) ::            ILAENV
      EXTERNAL           ILAENV
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
!
      IF( INFO.EQ.0 ) THEN
         IF( N.EQ.0 ) THEN
            LWKOPT = 1
         ELSE
            NB = ILAENV( 1, 'ZUNGQL', ' ', M, N, K, -1 )
            LWKOPT = N*NB
         END IF
         WORK( 1 ) = LWKOPT
!
         IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
            INFO = -8
         END IF
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNGQL', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.LE.0 ) THEN
         RETURN
      END IF
!
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
!
!        Determine when to cross over from blocked to unblocked code.
!
         NX = MAX( 0, ILAENV( 3, 'ZUNGQL', ' ', M, N, K, -1 ) )
         IF( NX.LT.K ) THEN
!
!           Determine if workspace is large enough for blocked code.
!
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
!
!              Not enough workspace to use optimal NB:  reduce NB and
!              determine the minimum value of NB.
!
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'ZUNGQL', ' ', M, N, K, -1 ) )
            END IF
         END IF
      END IF
!
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
!
!        Use blocked code after the first block.
!        The last kk columns are handled by the block method.
!
         KK = MIN( K, ( ( K-NX+NB-1 ) / NB )*NB )
!
!        Set A(m-kk+1:m,1:n-kk) to zero.
!
         DO 20 J = 1, N - KK
            DO 10 I = M - KK + 1, M
               A( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
      ELSE
         KK = 0
      END IF
!
!     Use unblocked code for the first or only block.
!
      CALL ZUNG2L( M-KK, N-KK, K-KK, A, LDA, TAU, WORK, IINFO )
!
      IF( KK.GT.0 ) THEN
!
!        Use blocked code
!
         DO 50 I = K - KK + 1, K, NB
            IB = MIN( NB, K-I+1 )
            IF( N-K+I.GT.1 ) THEN
!
!              Form the triangular factor of the block reflector
!              H = H(i+ib-1) . . . H(i+1) H(i)
!
               CALL ZLARFT( 'Backward', 'Columnwise', M-K+I+IB-1, IB, &
                            A( 1, N-K+I ), LDA, TAU( I ), WORK, LDWORK )
!
!              Apply H to A(1:m-k+i+ib-1,1:n-k+i-1) from the left
!
               CALL ZLARFB( 'Left', 'No transpose', 'Backward', &
                            'Columnwise', M-K+I+IB-1, N-K+I-1, IB, &
                            A( 1, N-K+I ), LDA, WORK, LDWORK, A, LDA, &
                            WORK( IB+1 ), LDWORK )
            END IF
!
!           Apply H to rows 1:m-k+i+ib-1 of current block
!
            CALL ZUNG2L( M-K+I+IB-1, IB, IB, A( 1, N-K+I ), LDA, TAU( I ), WORK, IINFO )
!
!           Set rows m-k+i+ib:m of current block to zero
!
            DO 40 J = N - K + I, N - K + I + IB - 1
               DO 30 L = M - K + I + IB, M
                  A( L, J ) = ZERO
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      END IF
!
      WORK( 1 ) = IWS
      RETURN
!
!     End of ZUNGQL
!
      END
!> \brief \b ZUNGQR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download ZUNGQR + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zungqr.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zungqr.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zungqr.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZUNGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
! 
!       .. Scalar Arguments ..
!       integer(i4) ::            INFO, K, LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZUNGQR generates an M-by-N complex matrix Q with orthonormal columns,
!> which is defined as the first N columns of a product of K elementary
!> reflectors of order M
!>
!>       Q  =  H(1) H(2) . . . H(k)
!>
!> as returned by ZGEQRF.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is integer(i4) ::
!>          The number of rows of the matrix Q. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is integer(i4) ::
!>          The number of columns of the matrix Q. M >= N >= 0.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is integer(i4) ::
!>          The number of elementary reflectors whose product defines the
!>          matrix Q. N >= K >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the i-th column must contain the vector which
!>          defines the elementary reflector H(i), for i = 1,2,...,k, as
!>          returned by ZGEQRF in the first k columns of its array
!>          argument A.
!>          On exit, the M-by-N matrix Q.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is integer(i4) ::
!>          The first dimension of the array A. LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is COMPLEX*16 array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by ZGEQRF.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is integer(i4) ::
!>          The dimension of the array WORK. LWORK >= max(1,N).
!>          For optimum performance LWORK >= N*NB, where NB is the
!>          optimal blocksize.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is integer(i4) ::
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument has an illegal value
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup complex16OTHERcomputational
!
!  =====================================================================
      SUBROUTINE ZUNGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes
      integer(i4)  :: INFO, K, LDA, LWORK, M, N
!     ..
!     .. Array Arguments ..
      complex(dpc) :: A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      complex(dpc) :: ZERO
      PARAMETER       ( ZERO = ( 0.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY
      integer(i4) :: I, IB, IINFO, IWS, J, KI, KK, L, LDWORK, LWKOPT, NB, NBMIN, NX
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARFB, ZLARFT, ZUNG2R
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. External Functions ..
      integer(i4)  :: ILAENV
      EXTERNAL        ILAENV
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      NB = ILAENV( 1, 'ZUNGQR', ' ', M, N, K, -1 )
      LWKOPT = MAX( 1, N )*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNGQR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.LE.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
!
!        Determine when to cross over from blocked to unblocked code.
!
         NX = MAX( 0, ILAENV( 3, 'ZUNGQR', ' ', M, N, K, -1 ) )
         IF( NX.LT.K ) THEN
!
!           Determine if workspace is large enough for blocked code.
!
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
!
!              Not enough workspace to use optimal NB:  reduce NB and
!              determine the minimum value of NB.
!
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'ZUNGQR', ' ', M, N, K, -1 ) )
            END IF
         END IF
      END IF
!
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
!
!        Use blocked code after the last block.
!        The first kk columns are handled by the block method.
!
         KI = ( ( K-NX-1 ) / NB )*NB
         KK = MIN( K, KI+NB )
!
!        Set A(1:kk,kk+1:n) to zero.
!
         DO 20 J = KK + 1, N
            DO 10 I = 1, KK
               A( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
      ELSE
         KK = 0
      END IF
!
!     Use unblocked code for the last or only block.
!
      IF( KK.LT.N ) CALL ZUNG2R( M-KK, N-KK, K-KK, A( KK+1, KK+1 ), LDA, TAU( KK+1 ), WORK, IINFO )
!
      IF( KK.GT.0 ) THEN
!
!        Use blocked code
!
         DO 50 I = KI + 1, 1, -NB
            IB = MIN( NB, K-I+1 )
            IF( I+IB.LE.N ) THEN
!
!              Form the triangular factor of the block reflector
!              H = H(i) H(i+1) . . . H(i+ib-1)
!
               CALL ZLARFT( 'Forward', 'Columnwise', M-I+1, IB, &
                            A( I, I ), LDA, TAU( I ), WORK, LDWORK )
!
!              Apply H to A(i:m,i+ib:n) from the left
!
               CALL ZLARFB( 'Left', 'No transpose', 'Forward', &
                            'Columnwise', M-I+1, N-I-IB+1, IB, &
                            A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ), &
                            LDA, WORK( IB+1 ), LDWORK )
            END IF
!
!           Apply H to rows i:m of current block
!
            CALL ZUNG2R( M-I+1, IB, IB, A( I, I ), LDA, TAU( I ), WORK, IINFO )
!
!           Set rows 1:i-1 of current block to zero
!
            DO 40 J = I, I + IB - 1
               DO 30 L = 1, I - 1
                  A( L, J ) = ZERO
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      END IF
!
      WORK( 1 ) = IWS
      RETURN
!
!     End of ZUNGQR
!
      END
!> \brief \b ZUNGTR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download ZUNGTR + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zungtr.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zungtr.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zungtr.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZUNGTR( UPLO, N, A, LDA, TAU, WORK, LWORK, INFO )
! 
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       integer(i4) ::            INFO, LDA, LWORK, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZUNGTR generates a complex unitary matrix Q which is defined as the
!> product of n-1 elementary reflectors of order N, as returned by
!> ZHETRD:
!>
!> if UPLO = 'U', Q = H(n-1) . . . H(2) H(1),
!>
!> if UPLO = 'L', Q = H(1) H(2) . . . H(n-1).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U': Upper triangle of A contains elementary reflectors
!>                 from ZHETRD;
!>          = 'L': Lower triangle of A contains elementary reflectors
!>                 from ZHETRD.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is integer(i4) ::
!>          The order of the matrix Q. N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the vectors which define the elementary reflectors,
!>          as returned by ZHETRD.
!>          On exit, the N-by-N unitary matrix Q.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is integer(i4) ::
!>          The leading dimension of the array A. LDA >= N.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is COMPLEX*16 array, dimension (N-1)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by ZHETRD.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is integer(i4) ::
!>          The dimension of the array WORK. LWORK >= N-1.
!>          For optimum performance LWORK >= (N-1)*NB, where NB is
!>          the optimal blocksize.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is integer(i4) ::
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup complex16OTHERcomputational
!
!  =====================================================================
      SUBROUTINE ZUNGTR( UPLO, N, A, LDA, TAU, WORK, LWORK, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes
      CHARACTER    :: UPLO
      integer(i4)  :: INFO, LDA, LWORK, N
!     ..
!     .. Array Arguments ..
      complex(dpc) :: A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      complex(dpc) :: ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ), ONE = ( 1.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL     ::     LQUERY, UPPER
      integer(i4) :: I, IINFO, J, LWKOPT, NB
!     ..
!     .. External Functions ..
      LOGICAL     :: LSAME
      integer(i4) :: ILAENV
      EXTERNAL       LSAME, ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZUNGQL, ZUNGQR
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.MAX( 1, N-1 ) .AND. .NOT.LQUERY ) THEN
         INFO = -7
      END IF
!
      IF( INFO.EQ.0 ) THEN
         IF( UPPER ) THEN
            NB = ILAENV( 1, 'ZUNGQL', ' ', N-1, N-1, N-1, -1 )
         ELSE
            NB = ILAENV( 1, 'ZUNGQR', ' ', N-1, N-1, N-1, -1 )
         END IF
         LWKOPT = MAX( 1, N-1 )*NB
         WORK( 1 ) = LWKOPT
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNGTR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!
      IF( UPPER ) THEN
!
!        Q was determined by a call to ZHETRD with UPLO = 'U'
!
!        Shift the vectors which define the elementary reflectors one
!        column to the left, and set the last row and column of Q to
!        those of the unit matrix
!
         DO 20 J = 1, N - 1
            DO 10 I = 1, J - 1
               A( I, J ) = A( I, J+1 )
   10       CONTINUE
            A( N, J ) = ZERO
   20    CONTINUE
         DO 30 I = 1, N - 1
            A( I, N ) = ZERO
   30    CONTINUE
         A( N, N ) = ONE
!
!        Generate Q(1:n-1,1:n-1)
!
         CALL ZUNGQL( N-1, N-1, N-1, A, LDA, TAU, WORK, LWORK, IINFO )
!
      ELSE
!
!        Q was determined by a call to ZHETRD with UPLO = 'L'.
!
!        Shift the vectors which define the elementary reflectors one
!        column to the right, and set the first row and column of Q to
!        those of the unit matrix
!
         DO 50 J = N, 2, -1
            A( 1, J ) = ZERO
            DO 40 I = J + 1, N
               A( I, J ) = A( I, J-1 )
   40       CONTINUE
   50    CONTINUE
         A( 1, 1 ) = ONE
         DO 60 I = 2, N
            A( I, 1 ) = ZERO
   60    CONTINUE
         IF( N.GT.1 ) THEN
!
!           Generate Q(2:n,2:n)
!
            CALL ZUNGQR( N-1, N-1, N-1, A( 2, 2 ), LDA, TAU, WORK, LWORK, IINFO )
         END IF
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
!
!     End of ZUNGTR
!
      END
!
!> \brief \b DLADIV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DLADIV + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dladiv.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dladiv.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dladiv.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLADIV( A, B, C, D, P, Q )
! 
!       .. Scalar Arguments ..
!       real(dp)     ::    A, B, C, D, P, Q
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLADIV performs complex division in  real arithmetic
!>
!>                       a + i*b
!>            p + i*q = ---------
!>                       c + i*d
!>
!> The algorithm is due to Robert L. Smith and can be found
!> in D. Knuth, The art of Computer Programming, Vol.2, p.195
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] A
!> \verbatim
!>          A is real(dp)     :: 
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is real(dp)     :: 
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is real(dp)     :: 
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is real(dp)     :: 
!>          The scalars a, b, c, and d in the above expression.
!> \endverbatim
!>
!> \param[out] P
!> \verbatim
!>          P is real(dp)     :: 
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is real(dp)     :: 
!>          The scalars p and q in the above expression.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup auxOTHERauxiliary
!
!  =====================================================================
      SUBROUTINE DLADIV( A, B, C, D, P, Q )
!
!  -- LAPACK auxiliary routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes
      real(dp)     ::    A, B, C, D, P, Q
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      real(dp)     ::    E, F
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS
!     ..
!     .. Executable Statements ..
!
      IF( ABS( D ).LT.ABS( C ) ) THEN
         E = D / C
         F = C + D*E
         P = ( A+B*E ) / F
         Q = ( B-A*E ) / F
      ELSE
         E = C / D
         F = D + C*E
         P = ( B+A*E ) / F
         Q = ( -A+B*E ) / F
      END IF
!
      RETURN
!
!     End of DLADIV
!
      END
!> \brief \b DLAPY3
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DLAPY3 + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlapy3.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlapy3.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlapy3.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       real(dp)     ::  FUNCTION DLAPY3( X, Y, Z )
! 
!       .. Scalar Arguments ..
!       real(dp)     ::    X, Y, Z
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAPY3 returns sqrt(x**2+y**2+z**2), taking care not to cause
!> unnecessary overflow.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] X
!> \verbatim
!>          X is real(dp)     :: 
!> \endverbatim
!>
!> \param[in] Y
!> \verbatim
!>          Y is real(dp)     :: 
!> \endverbatim
!>
!> \param[in] Z
!> \verbatim
!>          Z is real(dp)     :: 
!>          X, Y and Z specify the values x, y and z.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup auxOTHERauxiliary
!
!  =====================================================================
      FUNCTION DLAPY3( X, Y, Z )
!
!  -- LAPACK auxiliary routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes
      real(dp)     :: X, Y, Z
      real(dp)     :: DLAPY3
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      real(dp)     ::    ZERO
      PARAMETER          ( ZERO = 0.0D0 )
!     ..
!     .. Local Scalars ..
      real(dp)     ::    W, XABS, YABS, ZABS
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
!     ..
!     .. Executable Statements ..
!
      XABS = ABS( X )
      YABS = ABS( Y )
      ZABS = ABS( Z )
      W = MAX( XABS, YABS, ZABS )
      IF( W.EQ.ZERO ) THEN
!     W can be zero for max(0,nan,0)
!     adding all three entries together will make sure
!     NaN will not disappear.
         DLAPY3 =  XABS + YABS + ZABS
      ELSE
         DLAPY3 = W*SQRT( ( XABS / W )**2+( YABS / W )**2+ ( ZABS / W )**2 )
      END IF
      RETURN
!
!     End of DLAPY3
!
      END
!> \brief \b DISNAN
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DISNAN + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/disnan.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/disnan.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/disnan.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       LOGICAL FUNCTION DISNAN( DIN )
! 
!       .. Scalar Arguments ..
!       real(dp)     ::    DIN
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DISNAN returns .TRUE. if its argument is NaN, and .FALSE.
!> otherwise.  To be replaced by the Fortran 2003 intrinsic in the
!> future.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] DIN
!> \verbatim
!>          DIN is real(dp)     :: 
!>          Input to test for NaN.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup auxOTHERauxiliary
!
!  =====================================================================
      LOGICAL FUNCTION DISNAN( DIN )
!
!  -- LAPACK auxiliary routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes
      real(dp)     ::    DIN
!     ..
!
!  =====================================================================
!
!  .. External Functions ..
      LOGICAL DLAISNAN
      EXTERNAL DLAISNAN
!  ..
!  .. Executable Statements ..
      DISNAN = DLAISNAN(DIN,DIN)
      RETURN
      END
!> \brief \b ILAZLC
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download ILAZLC + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilazlc.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilazlc.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilazlc.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       integer(i4)  ::  FUNCTION ILAZLC( M, N, A, LDA )
! 
!       .. Scalar Arguments ..
!       integer(i4)  ::             M, N, LDA
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ILAZLC scans A for its last non-zero column.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is integer(i4)  :: 
!>          The number of rows of the matrix A.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is integer(i4)  :: 
!>          The number of columns of the matrix A.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          The m by n matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is integer(i4)  :: 
!>          The leading dimension of the array A. LDA >= max(1,M).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup complex16OTHERauxiliary
!
!  =====================================================================
      FUNCTION ILAZLC( M, N, A, LDA )
!
!  -- LAPACK auxiliary routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes
      integer(i4)  :: M, N, LDA
      integer(i4)  :: ILAZLC
!     ..
!     .. Array Arguments ..
      complex(dpc) :: A( LDA, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      complex(dpc) :: ZERO
      PARAMETER ( ZERO = (0.0D+0, 0.0D+0) )
!     ..
!     .. Local Scalars ..
      integer(i4)  ::  I
!     ..
!     .. Executable Statements ..
!
!     Quick test for the common case where one corner is non-zero.
      IF( N.EQ.0 ) THEN
         ILAZLC = N
      ELSE IF( A(1, N).NE.ZERO .OR. A(M, N).NE.ZERO ) THEN
         ILAZLC = N
      ELSE
!     Now scan each column from the end, returning with the first non-zero.
         DO ILAZLC = N, 1, -1
            DO I = 1, M
               IF( A(I, ILAZLC).NE.ZERO ) RETURN
            END DO
         END DO
      END IF
      RETURN
      END
!> \brief \b ILAZLR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download ILAZLR + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilazlr.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilazlr.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilazlr.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       integer(i4)  ::  FUNCTION ILAZLR( M, N, A, LDA )
! 
!       .. Scalar Arguments ..
!       integer(i4)  ::             M, N, LDA
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ILAZLR scans A for its last non-zero row.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is integer(i4)  :: 
!>          The number of rows of the matrix A.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is integer(i4)  :: 
!>          The number of columns of the matrix A.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          The m by n matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is integer(i4)  :: 
!>          The leading dimension of the array A. LDA >= max(1,M).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup complex16OTHERauxiliary
!
!  =====================================================================
      FUNCTION ILAZLR( M, N, A, LDA )
!
!  -- LAPACK auxiliary routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes
      integer(i4)  :: M, N, LDA
      integer(i4)  :: ILAZLR
!     ..
!     .. Array Arguments ..
      complex(dpc) :: A( LDA, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      complex(dpc) :: ZERO
      PARAMETER ( ZERO = (0.0D+0, 0.0D+0) )
!     ..
!     .. Local Scalars ..
      integer(i4)  ::  I, J
!     ..
!     .. Executable Statements ..
!
!     Quick test for the common case where one corner is non-zero.
      IF( M.EQ.0 ) THEN
         ILAZLR = M
      ELSE IF( A(M, 1).NE.ZERO .OR. A(M, N).NE.ZERO ) THEN
         ILAZLR = M
      ELSE
!     Scan up each column tracking the last zero row seen.
         ILAZLR = 0
         DO J = 1, N
            I=M
            DO WHILE ((A(I,J).NE.ZERO).AND.(I.GE.1))
              I=I-1
              IF (I.EQ.0) THEN
                 EXIT
              END IF
            ENDDO         
            ILAZLR = MAX( ILAZLR, I )
         END DO
      END IF
      RETURN
      END
!> \brief \b DLAISNAN
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DLAISNAN + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaisnan.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaisnan.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaisnan.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       LOGICAL FUNCTION DLAISNAN( DIN1, DIN2 )
! 
!       .. Scalar Arguments ..
!       DOUBLE PRECISION   DIN1, DIN2
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> This routine is not for general use.  It exists solely to avoid
!> over-optimization in DISNAN.
!>
!> DLAISNAN checks for NaNs by comparing its two arguments for
!> inequality.  NaN is the only floating-point value where NaN != NaN
!> returns .TRUE.  To check for NaNs, pass the same variable as both
!> arguments.
!>
!> A compiler must assume that the two arguments are
!> not the same variable, and the test will not be optimized away.
!> Interprocedural or whole-program optimization may delete this
!> test.  The ISNAN functions will be replaced by the correct
!> Fortran 03 intrinsic once the intrinsic is widely available.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] DIN1
!> \verbatim
!>          DIN1 is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] DIN2
!> \verbatim
!>          DIN2 is DOUBLE PRECISION
!>          Two numbers to compare for inequality.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup auxOTHERauxiliary
!
!  =====================================================================
      LOGICAL FUNCTION DLAISNAN( DIN1, DIN2 )
!
!  -- LAPACK auxiliary routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes
      real(dp) :: DIN1, DIN2
!     ..
!
!  =====================================================================
!
!  .. Executable Statements ..
      DLAISNAN = (DIN1.NE.DIN2)
      RETURN
      END
!
!  DIVIDE AND CONQUER RELATED ROUTINES AFTER HERE
!
!> \brief <b> ZHEEVD computes the eigenvalues and, optionally, the left and/or right eigenvectors for HE matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download ZHEEVD + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zheevd.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zheevd.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zheevd.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZHEEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK,
!                          LRWORK, IWORK, LIWORK, INFO )
! 
!       .. Scalar Arguments ..
!       CHARACTER          JOBZ, UPLO
!       INTEGER            INFO, LDA, LIWORK, LRWORK, LWORK, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
!       DOUBLE PRECISION   RWORK( * ), W( * )
!       COMPLEX*16         A( LDA, * ), WORK( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZHEEVD computes all eigenvalues and, optionally, eigenvectors of a
!> complex Hermitian matrix A.  If eigenvectors are desired, it uses a
!> divide and conquer algorithm.
!>
!> The divide and conquer algorithm makes very mild assumptions about
!> floating point arithmetic. It will work on machines with a guard
!> digit in add/subtract, or on those binary machines without guard
!> digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
!> Cray-2. It could conceivably fail on hexadecimal or decimal machines
!> without guard digits, but we know of none.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOBZ
!> \verbatim
!>          JOBZ is CHARACTER*1
!>          = 'N':  Compute eigenvalues only;
!>          = 'V':  Compute eigenvalues and eigenvectors.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  Upper triangle of A is stored;
!>          = 'L':  Lower triangle of A is stored.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA, N)
!>          On entry, the Hermitian matrix A.  If UPLO = 'U', the
!>          leading N-by-N upper triangular part of A contains the
!>          upper triangular part of the matrix A.  If UPLO = 'L',
!>          the leading N-by-N lower triangular part of A contains
!>          the lower triangular part of the matrix A.
!>          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
!>          orthonormal eigenvectors of the matrix A.
!>          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
!>          or the upper triangle (if UPLO='U') of A, including the
!>          diagonal, is destroyed.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is DOUBLE PRECISION array, dimension (N)
!>          If INFO = 0, the eigenvalues in ascending order.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The length of the array WORK.
!>          If N <= 1,                LWORK must be at least 1.
!>          If JOBZ  = 'N' and N > 1, LWORK must be at least N + 1.
!>          If JOBZ  = 'V' and N > 1, LWORK must be at least 2*N + N**2.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal sizes of the WORK, RWORK and
!>          IWORK arrays, returns these values as the first entries of
!>          the WORK, RWORK and IWORK arrays, and no error message
!>          related to LWORK or LRWORK or LIWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array,
!>                                         dimension (LRWORK)
!>          On exit, if INFO = 0, RWORK(1) returns the optimal LRWORK.
!> \endverbatim
!>
!> \param[in] LRWORK
!> \verbatim
!>          LRWORK is INTEGER
!>          The dimension of the array RWORK.
!>          If N <= 1,                LRWORK must be at least 1.
!>          If JOBZ  = 'N' and N > 1, LRWORK must be at least N.
!>          If JOBZ  = 'V' and N > 1, LRWORK must be at least
!>                         1 + 5*N + 2*N**2.
!>
!>          If LRWORK = -1, then a workspace query is assumed; the
!>          routine only calculates the optimal sizes of the WORK, RWORK
!>          and IWORK arrays, returns these values as the first entries
!>          of the WORK, RWORK and IWORK arrays, and no error message
!>          related to LWORK or LRWORK or LIWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (MAX(1,LIWORK))
!>          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
!> \endverbatim
!>
!> \param[in] LIWORK
!> \verbatim
!>          LIWORK is INTEGER
!>          The dimension of the array IWORK.
!>          If N <= 1,                LIWORK must be at least 1.
!>          If JOBZ  = 'N' and N > 1, LIWORK must be at least 1.
!>          If JOBZ  = 'V' and N > 1, LIWORK must be at least 3 + 5*N.
!>
!>          If LIWORK = -1, then a workspace query is assumed; the
!>          routine only calculates the optimal sizes of the WORK, RWORK
!>          and IWORK arrays, returns these values as the first entries
!>          of the WORK, RWORK and IWORK arrays, and no error message
!>          related to LWORK or LRWORK or LIWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  if INFO = i and JOBZ = 'N', then the algorithm failed
!>                to converge; i off-diagonal elements of an intermediate
!>                tridiagonal form did not converge to zero;
!>                if INFO = i and JOBZ = 'V', then the algorithm failed
!>                to compute an eigenvalue while working on the submatrix
!>                lying in rows and columns INFO/(N+1) through
!>                mod(INFO,N+1).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup complex16HEeigen
!
!> \par Further Details:
!  =====================
!>
!>  Modified description of INFO. Sven, 16 Feb 05.
!
!> \par Contributors:
!  ==================
!>
!> Jeff Rutter, Computer Science Division, University of California
!> at Berkeley, USA
!>
!  =====================================================================
      SUBROUTINE ZHEEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )
!
!  -- LAPACK driver routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes
      CHARACTER          JOBZ, UPLO
      INTEGER(i4)        INFO, LDA, LIWORK, LRWORK, LWORK, N
!     ..
!     .. Array Arguments ..
      INTEGER(i4)        IWORK( * )
      real(dp)           RWORK( * ), W( * )
      COMPLEX(dpc)       A( LDA, * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      real(dp)           ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
      complex(dpc)       CONE
      PARAMETER          ( CONE = ( 1.0D0, 0.0D0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            LOWER, LQUERY, WANTZ
      INTEGER(i4)        IINFO, IMAX, INDE, INDRWK, INDTAU, INDWK2, &
                         INDWRK, ISCALE, LIOPT, LIWMIN, LLRWK, LLWORK, &
                         LLWRK2, LOPT, LROPT, LRWMIN, LWMIN
      real(dp)           ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER(i4)        ILAENV
      real(dp)           DLAMCH, ZLANHE
      EXTERNAL           LSAME, ILAENV, DLAMCH, ZLANHE
!     ..
!     .. External Subroutines ..
      EXTERNAL           DSCAL, DSTERF, XERBLA, ZHETRD, ZLACPY, ZLASCL, ZSTEDC, ZUNMTR
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      WANTZ = LSAME( JOBZ, 'V' )
      LOWER = LSAME( UPLO, 'L' )
      LQUERY = ( LWORK.EQ.-1 .OR. LRWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
!
      INFO = 0
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
!
      IF( INFO.EQ.0 ) THEN
         IF( N.LE.1 ) THEN
            LWMIN = 1
            LRWMIN = 1
            LIWMIN = 1
            LOPT = LWMIN
            LROPT = LRWMIN
            LIOPT = LIWMIN
         ELSE
            IF( WANTZ ) THEN
               LWMIN = 2*N + N*N
               LRWMIN = 1 + 5*N + 2*N**2
               LIWMIN = 3 + 5*N
            ELSE
               LWMIN = N + 1
               LRWMIN = N
               LIWMIN = 1
            END IF
            LOPT = MAX( LWMIN, N + ILAENV( 1, 'ZHETRD', UPLO, N, -1, -1, -1 ) )
            LROPT = LRWMIN
            LIOPT = LIWMIN
         END IF
         WORK( 1 ) = LOPT
         RWORK( 1 ) = LROPT
         IWORK( 1 ) = LIOPT
!
         IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -8
         ELSE IF( LRWORK.LT.LRWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -10
         ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -12
         END IF
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHEEVD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) RETURN
!
      IF( N.EQ.1 ) THEN
         W( 1 ) = A( 1, 1 )
         IF( WANTZ ) A( 1, 1 ) = CONE
         RETURN
      END IF
!
!     Get machine constants.
!
      SAFMIN = DLAMCH( 'Safe minimum' )
      EPS = DLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = SQRT( BIGNUM )
!
!     Scale matrix to allowable range, if necessary.
!
      ANRM = ZLANHE( 'M', UPLO, N, A, LDA, RWORK )
      ISCALE = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / ANRM
      ELSE IF( ANRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / ANRM
      END IF
      IF( ISCALE.EQ.1 ) CALL ZLASCL( UPLO, 0, 0, ONE, SIGMA, N, N, A, LDA, INFO )
!
!     Call ZHETRD to reduce Hermitian matrix to tridiagonal form.
!
      INDE = 1
      INDTAU = 1
      INDWRK = INDTAU + N
      INDRWK = INDE + N
      INDWK2 = INDWRK + N*N
      LLWORK = LWORK - INDWRK + 1
      LLWRK2 = LWORK - INDWK2 + 1
      LLRWK = LRWORK - INDRWK + 1
      CALL ZHETRD( UPLO, N, A, LDA, W, RWORK( INDE ), WORK( INDTAU ), WORK( INDWRK ), LLWORK, IINFO )
!
!     For eigenvalues only, call DSTERF.  For eigenvectors, first call
!     ZSTEDC to generate the eigenvector matrix, WORK(INDWRK), of the
!     tridiagonal matrix, then call ZUNMTR to multiply it to the
!     Householder transformations represented as Householder vectors in
!     A.
!
      IF( .NOT.WANTZ ) THEN
         CALL DSTERF( N, W, RWORK( INDE ), INFO )
      ELSE
         CALL ZSTEDC( 'I', N, W, RWORK( INDE ), WORK( INDWRK ), N, &
                      WORK( INDWK2 ), LLWRK2, RWORK( INDRWK ), LLRWK, &
                      IWORK, LIWORK, INFO )
         CALL ZUNMTR( 'L', UPLO, 'N', N, N, A, LDA, WORK( INDTAU ), &
                      WORK( INDWRK ), N, WORK( INDWK2 ), LLWRK2, IINFO )
         CALL ZLACPY( 'A', N, N, WORK( INDWRK ), N, A, LDA )
      END IF
!
!     If matrix was scaled, then rescale eigenvalues appropriately.
!
      IF( ISCALE.EQ.1 ) THEN
         IF( INFO.EQ.0 ) THEN
            IMAX = N
         ELSE
            IMAX = INFO - 1
         END IF
         CALL DSCAL( IMAX, ONE / SIGMA, W, 1 )
      END IF
!
      WORK( 1 ) = LOPT
      RWORK( 1 ) = LROPT
      IWORK( 1 ) = LIOPT
!
      RETURN
!
!     End of ZHEEVD
!
      END
!> \brief \b ZSTEDC
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download ZSTEDC + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zstedc.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zstedc.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zstedc.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZSTEDC( COMPZ, N, D, E, Z, LDZ, WORK, LWORK, RWORK,
!                          LRWORK, IWORK, LIWORK, INFO )
! 
!       .. Scalar Arguments ..
!       CHARACTER          COMPZ
!       INTEGER(i4)        INFO, LDZ, LIWORK, LRWORK, LWORK, N
!       ..
!       .. Array Arguments ..
!       INTEGER(i4)        IWORK( * )
!       DOUBLE PRECISION   D( * ), E( * ), RWORK( * )
!       COMPLEX*16         WORK( * ), Z( LDZ, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZSTEDC computes all eigenvalues and, optionally, eigenvectors of a
!> symmetric tridiagonal matrix using the divide and conquer method.
!> The eigenvectors of a full or band complex Hermitian matrix can also
!> be found if ZHETRD or ZHPTRD or ZHBTRD has been used to reduce this
!> matrix to tridiagonal form.
!>
!> This code makes very mild assumptions about floating point
!> arithmetic. It will work on machines with a guard digit in
!> add/subtract, or on those binary machines without guard digits
!> which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
!> It could conceivably fail on hexadecimal or decimal machines
!> without guard digits, but we know of none.  See DLAED3 for details.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] COMPZ
!> \verbatim
!>          COMPZ is CHARACTER*1
!>          = 'N':  Compute eigenvalues only.
!>          = 'I':  Compute eigenvectors of tridiagonal matrix also.
!>          = 'V':  Compute eigenvectors of original Hermitian matrix
!>                  also.  On entry, Z contains the unitary matrix used
!>                  to reduce the original matrix to tridiagonal form.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The dimension of the symmetric tridiagonal matrix.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>          On entry, the diagonal elements of the tridiagonal matrix.
!>          On exit, if INFO = 0, the eigenvalues in ascending order.
!> \endverbatim
!>
!> \param[in,out] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension (N-1)
!>          On entry, the subdiagonal elements of the tridiagonal matrix.
!>          On exit, E has been destroyed.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is COMPLEX*16 array, dimension (LDZ,N)
!>          On entry, if COMPZ = 'V', then Z contains the unitary
!>          matrix used in the reduction to tridiagonal form.
!>          On exit, if INFO = 0, then if COMPZ = 'V', Z contains the
!>          orthonormal eigenvectors of the original Hermitian matrix,
!>          and if COMPZ = 'I', Z contains the orthonormal eigenvectors
!>          of the symmetric tridiagonal matrix.
!>          If  COMPZ = 'N', then Z is not referenced.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of the array Z.  LDZ >= 1.
!>          If eigenvectors are desired, then LDZ >= max(1,N).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.
!>          If COMPZ = 'N' or 'I', or N <= 1, LWORK must be at least 1.
!>          If COMPZ = 'V' and N > 1, LWORK must be at least N*N.
!>          Note that for COMPZ = 'V', then if N is less than or
!>          equal to the minimum divide size, usually 25, then LWORK need
!>          only be 1.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal sizes of the WORK, RWORK and
!>          IWORK arrays, returns these values as the first entries of
!>          the WORK, RWORK and IWORK arrays, and no error message
!>          related to LWORK or LRWORK or LIWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array,
!>                                         dimension (LRWORK)
!>          On exit, if INFO = 0, RWORK(1) returns the optimal LRWORK.
!> \endverbatim
!>
!> \param[in] LRWORK
!> \verbatim
!>          LRWORK is INTEGER
!>          The dimension of the array RWORK.
!>          If COMPZ = 'N' or N <= 1, LRWORK must be at least 1.
!>          If COMPZ = 'V' and N > 1, LRWORK must be at least
!>                         1 + 3*N + 2*N*lg N + 4*N**2 ,
!>                         where lg( N ) = smallest integer k such
!>                         that 2**k >= N.
!>          If COMPZ = 'I' and N > 1, LRWORK must be at least
!>                         1 + 4*N + 2*N**2 .
!>          Note that for COMPZ = 'I' or 'V', then if N is less than or
!>          equal to the minimum divide size, usually 25, then LRWORK
!>          need only be max(1,2*(N-1)).
!>
!>          If LRWORK = -1, then a workspace query is assumed; the
!>          routine only calculates the optimal sizes of the WORK, RWORK
!>          and IWORK arrays, returns these values as the first entries
!>          of the WORK, RWORK and IWORK arrays, and no error message
!>          related to LWORK or LRWORK or LIWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (MAX(1,LIWORK))
!>          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
!> \endverbatim
!>
!> \param[in] LIWORK
!> \verbatim
!>          LIWORK is INTEGER
!>          The dimension of the array IWORK.
!>          If COMPZ = 'N' or N <= 1, LIWORK must be at least 1.
!>          If COMPZ = 'V' or N > 1,  LIWORK must be at least
!>                                    6 + 6*N + 5*N*lg N.
!>          If COMPZ = 'I' or N > 1,  LIWORK must be at least
!>                                    3 + 5*N .
!>          Note that for COMPZ = 'I' or 'V', then if N is less than or
!>          equal to the minimum divide size, usually 25, then LIWORK
!>          need only be 1.
!>
!>          If LIWORK = -1, then a workspace query is assumed; the
!>          routine only calculates the optimal sizes of the WORK, RWORK
!>          and IWORK arrays, returns these values as the first entries
!>          of the WORK, RWORK and IWORK arrays, and no error message
!>          related to LWORK or LRWORK or LIWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  The algorithm failed to compute an eigenvalue while
!>                working on the submatrix lying in rows and columns
!>                INFO/(N+1) through mod(INFO,N+1).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup complex16OTHERcomputational
!
!> \par Contributors:
!  ==================
!>
!> Jeff Rutter, Computer Science Division, University of California
!> at Berkeley, USA
!
!  =====================================================================
      SUBROUTINE ZSTEDC( COMPZ, N, D, E, Z, LDZ, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes
      CHARACTER          COMPZ
      INTEGER(i4)        INFO, LDZ, LIWORK, LRWORK, LWORK, N
!     ..
!     .. Array Arguments ..
      INTEGER(i4)        IWORK( * )
      real(dp)           D( * ), E( * ), RWORK( * )
      complex(dpc)       WORK( * ), Z( LDZ, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      real(dp)    ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER(i4)        FINISH, I, ICOMPZ, II, J, K, LGN, LIWMIN, LL, LRWMIN, LWMIN, M, SMLSIZ, START
      real(dp)    EPS, ORGNRM, P, TINY
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER(i4)        ILAENV
      real(dp)           DLAMCH, DLANST
      EXTERNAL           LSAME, ILAENV, DLAMCH, DLANST
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLASCL, DLASET, DSTEDC, DSTEQR, DSTERF, XERBLA, &
                         ZLACPY, ZLACRM, ZLAED0, ZSTEQR, ZSWAP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, INT, LOG, MAX, MOD, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 .OR. LRWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
!
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ICOMPZ = 0
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ICOMPZ = 2
      ELSE
         ICOMPZ = -1
      END IF
      IF( ICOMPZ.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( ( LDZ.LT.1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1, N ) ) ) THEN
         INFO = -6
      END IF
!
      IF( INFO.EQ.0 ) THEN
!
!        Compute the workspace requirements
!
         SMLSIZ = ILAENV( 9, 'ZSTEDC', ' ', 0, 0, 0, 0 )
         IF( N.LE.1 .OR. ICOMPZ.EQ.0 ) THEN
            LWMIN = 1
            LIWMIN = 1
            LRWMIN = 1
         ELSE IF( N.LE.SMLSIZ ) THEN
            LWMIN = 1
            LIWMIN = 1
            LRWMIN = 2*( N - 1 )
         ELSE IF( ICOMPZ.EQ.1 ) THEN
            LGN = INT( LOG( DBLE( N ) ) / LOG( TWO ) )
            IF( 2**LGN.LT.N ) LGN = LGN + 1
            IF( 2**LGN.LT.N ) LGN = LGN + 1
            LWMIN = N*N
            LRWMIN = 1 + 3*N + 2*N*LGN + 4*N**2
            LIWMIN = 6 + 6*N + 5*N*LGN
         ELSE IF( ICOMPZ.EQ.2 ) THEN
            LWMIN = 1
            LRWMIN = 1 + 4*N + 2*N**2
            LIWMIN = 3 + 5*N
         END IF
         WORK( 1 ) = LWMIN
         RWORK( 1 ) = LRWMIN
         IWORK( 1 ) = LIWMIN
!
         IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -8
         ELSE IF( LRWORK.LT.LRWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -10
         ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -12
         END IF
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZSTEDC', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) RETURN
      IF( N.EQ.1 ) THEN
         IF( ICOMPZ.NE.0 ) Z( 1, 1 ) = ONE
         RETURN
      END IF
!
!     If the following conditional clause is removed, then the routine
!     will use the Divide and Conquer routine to compute only the
!     eigenvalues, which requires (3N + 3N**2) real workspace and
!     (2 + 5N + 2N lg(N)) integer workspace.
!     Since on many architectures DSTERF is much faster than any other
!     algorithm for finding eigenvalues only, it is used here
!     as the default. If the conditional clause is removed, then
!     information on the size of workspace needs to be changed.
!
!     If COMPZ = 'N', use DSTERF to compute the eigenvalues.
!
      IF( ICOMPZ.EQ.0 ) THEN
         CALL DSTERF( N, D, E, INFO )
         GO TO 70
      END IF
!
!     If N is smaller than the minimum divide size (SMLSIZ+1), then
!     solve the problem with another solver.
!
      IF( N.LE.SMLSIZ ) THEN
!
         CALL ZSTEQR( COMPZ, N, D, E, Z, LDZ, RWORK, INFO )
!
      ELSE
!
!        If COMPZ = 'I', we simply call DSTEDC instead.
!
         IF( ICOMPZ.EQ.2 ) THEN
            CALL DLASET( 'Full', N, N, ZERO, ONE, RWORK, N )
            LL = N*N + 1
            CALL DSTEDC( 'I', N, D, E, RWORK, N, RWORK( LL ), LRWORK-LL+1, IWORK, LIWORK, INFO )
            DO 20 J = 1, N
               DO 10 I = 1, N
                  Z( I, J ) = RWORK( ( J-1 )*N+I )
   10          CONTINUE
   20       CONTINUE
            GO TO 70
         END IF
!
!        From now on, only option left to be handled is COMPZ = 'V',
!        i.e. ICOMPZ = 1.
!
!        Scale.
!
         ORGNRM = DLANST( 'M', N, D, E )
         IF( ORGNRM.EQ.ZERO ) GO TO 70
!
         EPS = DLAMCH( 'Epsilon' )
!
         START = 1
!
!        while ( START <= N )
!
   30    CONTINUE
         IF( START.LE.N ) THEN
!
!           Let FINISH be the position of the next subdiagonal entry
!           such that E( FINISH ) <= TINY or FINISH = N if no such
!           subdiagonal exists.  The matrix identified by the elements
!           between START and FINISH constitutes an independent
!           sub-problem.
!
            FINISH = START
   40       CONTINUE
            IF( FINISH.LT.N ) THEN
               TINY = EPS*SQRT( ABS( D( FINISH ) ) )* SQRT( ABS( D( FINISH+1 ) ) )
               IF( ABS( E( FINISH ) ).GT.TINY ) THEN
                  FINISH = FINISH + 1
                  GO TO 40
               END IF
            END IF
!
!           (Sub) Problem determined.  Compute its size and solve it.
!
            M = FINISH - START + 1
            IF( M.GT.SMLSIZ ) THEN
!
!              Scale.
!
               ORGNRM = DLANST( 'M', M, D( START ), E( START ) )
               CALL DLASCL( 'G', 0, 0, ORGNRM, ONE, M, 1, D( START ), M, INFO )
               CALL DLASCL( 'G', 0, 0, ORGNRM, ONE, M-1, 1, E( START ), M-1, INFO )
!
               CALL ZLAED0( N, M, D( START ), E( START ), Z( 1, START ), LDZ, WORK, N, RWORK, IWORK, INFO )
               IF( INFO.GT.0 ) THEN
                  INFO = ( INFO / ( M+1 )+START-1 )*( N+1 ) + MOD( INFO, ( M+1 ) ) + START - 1
                  GO TO 70
               END IF
!
!              Scale back.
!
               CALL DLASCL( 'G', 0, 0, ONE, ORGNRM, M, 1, D( START ), M, INFO )
!
            ELSE
               CALL DSTEQR( 'I', M, D( START ), E( START ), RWORK, M, RWORK( M*M+1 ), INFO )
               CALL ZLACRM( N, M, Z( 1, START ), LDZ, RWORK, M, WORK, N, RWORK( M*M+1 ) )
               CALL ZLACPY( 'A', N, M, WORK, N, Z( 1, START ), LDZ )
               IF( INFO.GT.0 ) THEN
                  INFO = START*( N+1 ) + FINISH
                  GO TO 70
               END IF
            END IF
!
            START = FINISH + 1
            GO TO 30
         END IF
!
!        endwhile
!
!        If the problem split any number of times, then the eigenvalues
!        will not be properly ordered.  Here we permute the eigenvalues
!        (and the associated eigenvectors) into ascending order.
!
         IF( M.NE.N ) THEN
!
!           Use Selection Sort to minimize swaps of eigenvectors
!
            DO 60 II = 2, N
               I = II - 1
               K = I
               P = D( I )
               DO 50 J = II, N
                  IF( D( J ).LT.P ) THEN
                     K = J
                     P = D( J )
                  END IF
   50          CONTINUE
               IF( K.NE.I ) THEN
                  D( K ) = D( I )
                  D( I ) = P
                  CALL ZSWAP( N, Z( 1, I ), 1, Z( 1, K ), 1 )
               END IF
   60       CONTINUE
         END IF
      END IF
!
   70 CONTINUE
      WORK( 1 ) = LWMIN
      RWORK( 1 ) = LRWMIN
      IWORK( 1 ) = LIWMIN
!
      RETURN
!
!     End of ZSTEDC
!
      END
!> \brief \b ZUNMTR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download ZUNMTR + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zunmtr.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zunmtr.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zunmtr.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZUNMTR( SIDE, UPLO, TRANS, M, N, A, LDA, TAU, C, LDC,
!                          WORK, LWORK, INFO )
! 
!       .. Scalar Arguments ..
!       CHARACTER          SIDE, TRANS, UPLO
!       INTEGER            INFO, LDA, LDC, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZUNMTR overwrites the general complex M-by-N matrix C with
!>
!>                 SIDE = 'L'     SIDE = 'R'
!> TRANS = 'N':      Q * C          C * Q
!> TRANS = 'C':      Q**H * C       C * Q**H
!>
!> where Q is a complex unitary matrix of order nq, with nq = m if
!> SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of
!> nq-1 elementary reflectors, as returned by ZHETRD:
!>
!> if UPLO = 'U', Q = H(nq-1) . . . H(2) H(1);
!>
!> if UPLO = 'L', Q = H(1) H(2) . . . H(nq-1).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'L': apply Q or Q**H from the Left;
!>          = 'R': apply Q or Q**H from the Right.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U': Upper triangle of A contains elementary reflectors
!>                 from ZHETRD;
!>          = 'L': Lower triangle of A contains elementary reflectors
!>                 from ZHETRD.
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          = 'N':  No transpose, apply Q;
!>          = 'C':  Conjugate transpose, apply Q**H.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix C. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix C. N >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension
!>                               (LDA,M) if SIDE = 'L'
!>                               (LDA,N) if SIDE = 'R'
!>          The vectors which define the elementary reflectors, as
!>          returned by ZHETRD.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.
!>          LDA >= max(1,M) if SIDE = 'L'; LDA >= max(1,N) if SIDE = 'R'.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is COMPLEX*16 array, dimension
!>                               (M-1) if SIDE = 'L'
!>                               (N-1) if SIDE = 'R'
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by ZHETRD.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is COMPLEX*16 array, dimension (LDC,N)
!>          On entry, the M-by-N matrix C.
!>          On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C. LDC >= max(1,M).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.
!>          If SIDE = 'L', LWORK >= max(1,N);
!>          if SIDE = 'R', LWORK >= max(1,M).
!>          For optimum performance LWORK >= N*NB if SIDE = 'L', and
!>          LWORK >=M*NB if SIDE = 'R', where NB is the optimal
!>          blocksize.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup complex16OTHERcomputational
!
!  =====================================================================
      SUBROUTINE ZUNMTR( SIDE, UPLO, TRANS, M, N, A, LDA, TAU, C, LDC, WORK, LWORK, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes, only : i4, dpc
      CHARACTER     :: SIDE, TRANS, UPLO
      INTEGER(i4)   :: INFO, LDA, LDC, LWORK, M, N
!     ..
!     .. Array Arguments ..
      complex(dpc)  :: A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL      ::    LEFT, LQUERY, UPPER
      INTEGER(i4)  ::    I1, I2, IINFO, LWKOPT, MI, NB, NI, NQ, NW
!     ..
!     .. External Functions ..
      LOGICAL      ::    LSAME
      INTEGER(i4)  ::    ILAENV
      EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZUNMQL, ZUNMQR
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 )
!
!     NQ is the order of Q and NW is the minimum dimension of WORK
!
      IF( LEFT ) THEN
         NQ = M
         NW = N
      ELSE
         NQ = N
         NW = M
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.LSAME( TRANS, 'N' ) .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF( LWORK.LT.MAX( 1, NW ) .AND. .NOT.LQUERY ) THEN
         INFO = -12
      END IF
!
      IF( INFO.EQ.0 ) THEN
         IF( UPPER ) THEN
            IF( LEFT ) THEN
               NB = ILAENV( 1, 'ZUNMQL', SIDE // TRANS, M-1, N, M-1, -1 )
            ELSE
               NB = ILAENV( 1, 'ZUNMQL', SIDE // TRANS, M, N-1, N-1, -1 )
            END IF
         ELSE
            IF( LEFT ) THEN
               NB = ILAENV( 1, 'ZUNMQR', SIDE // TRANS, M-1, N, M-1, -1 )
            ELSE
               NB = ILAENV( 1, 'ZUNMQR', SIDE // TRANS, M, N-1, N-1, -1 )
            END IF
         END IF
         LWKOPT = MAX( 1, NW )*NB
         WORK( 1 ) = LWKOPT
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNMTR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 .OR. NQ.EQ.1 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!
      IF( LEFT ) THEN
         MI = M - 1
         NI = N
      ELSE
         MI = M
         NI = N - 1
      END IF
!
      IF( UPPER ) THEN
!
!        Q was determined by a call to ZHETRD with UPLO = 'U'
!
         CALL ZUNMQL( SIDE, TRANS, MI, NI, NQ-1, A( 1, 2 ), LDA, TAU, C, LDC, WORK, LWORK, IINFO )
      ELSE
!
!        Q was determined by a call to ZHETRD with UPLO = 'L'
!
         IF( LEFT ) THEN
            I1 = 2
            I2 = 1
         ELSE
            I1 = 1
            I2 = 2
         END IF
         CALL ZUNMQR( SIDE, TRANS, MI, NI, NQ-1, A( 2, 1 ), LDA, TAU, C( I1, I2 ), LDC, WORK, LWORK, IINFO )
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
!
!     End of ZUNMTR
!
      END
!> \brief \b ZLAED0
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download ZLAED0 + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaed0.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaed0.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaed0.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLAED0( QSIZ, N, D, E, Q, LDQ, QSTORE, LDQS, RWORK,
!                          IWORK, INFO )
! 
!       .. Scalar Arguments ..
!       integer(i4)  ::            INFO, LDQ, LDQS, N, QSIZ
!       ..
!       .. Array Arguments ..
!       integer(i4)  ::            IWORK( * )
!       real(dp) ::   D( * ), E( * ), RWORK( * )
!       COMPLEX*16         Q( LDQ, * ), QSTORE( LDQS, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> Using the divide and conquer method, ZLAED0 computes all eigenvalues
!> of a symmetric tridiagonal matrix which is one diagonal block of
!> those from reducing a dense or band Hermitian matrix and
!> corresponding eigenvectors of the dense or band matrix.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] QSIZ
!> \verbatim
!>          QSIZ is integer(i4)  ::
!>         The dimension of the unitary matrix used to reduce
!>         the full matrix to tridiagonal form.  QSIZ >= N if ICOMPQ = 1.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is integer(i4)  ::
!>         The dimension of the symmetric tridiagonal matrix.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is real(dp) :: array, dimension (N)
!>         On entry, the diagonal elements of the tridiagonal matrix.
!>         On exit, the eigenvalues in ascending order.
!> \endverbatim
!>
!> \param[in,out] E
!> \verbatim
!>          E is real(dp) :: array, dimension (N-1)
!>         On entry, the off-diagonal elements of the tridiagonal matrix.
!>         On exit, E has been destroyed.
!> \endverbatim
!>
!> \param[in,out] Q
!> \verbatim
!>          Q is COMPLEX*16 array, dimension (LDQ,N)
!>         On entry, Q must contain an QSIZ x N matrix whose columns
!>         unitarily orthonormal. It is a part of the unitary matrix
!>         that reduces the full dense Hermitian matrix to a
!>         (reducible) symmetric tridiagonal matrix.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is integer(i4)  ::
!>         The leading dimension of the array Q.  LDQ >= max(1,N).
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is integer(i4)  :: array,
!>         the dimension of IWORK must be at least
!>                      6 + 6*N + 5*N*lg N
!>                      ( lg( N ) = smallest integer k
!>                                  such that 2^k >= N )
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is real(dp) :: array,
!>                               dimension (1 + 3*N + 2*N*lg N + 3*N**2)
!>                        ( lg( N ) = smallest integer k
!>                                    such that 2^k >= N )
!> \endverbatim
!>
!> \param[out] QSTORE
!> \verbatim
!>          QSTORE is COMPLEX*16 array, dimension (LDQS, N)
!>         Used to store parts of
!>         the eigenvector matrix when the updating matrix multiplies
!>         take place.
!> \endverbatim
!>
!> \param[in] LDQS
!> \verbatim
!>          LDQS is integer(i4)  ::
!>         The leading dimension of the array QSTORE.
!>         LDQS >= max(1,N).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is integer(i4)  ::
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  The algorithm failed to compute an eigenvalue while
!>                working on the submatrix lying in rows and columns
!>                INFO/(N+1) through mod(INFO,N+1).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup complex16OTHERcomputational
!
!  =====================================================================
      SUBROUTINE ZLAED0( QSIZ, N, D, E, Q, LDQ, QSTORE, LDQS, RWORK, IWORK, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes
      integer(i4)  ::            INFO, LDQ, LDQS, N, QSIZ
!     ..
!     .. Array Arguments ..
      integer(i4)  ::            IWORK( * )
      real(dp) ::   D( * ), E( * ), RWORK( * )
      COMPLEX(dpc)       Q( LDQ, * ), QSTORE( LDQS, * )
!     ..
!
!  =====================================================================
!
!  Warning:      N could be as big as QSIZ!
!
!     .. Parameters ..
      real(dp) ::   TWO
      PARAMETER          ( TWO = 2.D+0 )
!     ..
!     .. Local Scalars ..
      integer(i4)  ::            CURLVL, CURPRB, CURR, I, IGIVCL, IGIVNM, &
                         IGIVPT, INDXQ, IPERM, IPRMPT, IQ, IQPTR, IWREM, &
                         J, K, LGN, LL, MATSIZ, MSD2, SMLSIZ, SMM1, &
                         SPM1, SPM2, SUBMAT, SUBPBS, TLVLS
      real(dp) ::   TEMP
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPY, DSTEQR, XERBLA, ZCOPY, ZLACRM, ZLAED7
!     ..
!     .. External Functions ..
      integer(i4)  ::            ILAENV
      EXTERNAL           ILAENV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, INT, LOG, MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
!
!     IF( ICOMPQ .LT. 0 .OR. ICOMPQ .GT. 2 ) THEN
!        INFO = -1
!     ELSE IF( ( ICOMPQ .EQ. 1 ) .AND. ( QSIZ .LT. MAX( 0, N ) ) )
!    $        THEN
      IF( QSIZ.LT.MAX( 0, N ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDQS.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZLAED0', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) RETURN
!
      SMLSIZ = ILAENV( 9, 'ZLAED0', ' ', 0, 0, 0, 0 )
!
!     Determine the size and placement of the submatrices, and save in
!     the leading elements of IWORK.
!
      IWORK( 1 ) = N
      SUBPBS = 1
      TLVLS = 0
   10 CONTINUE
      IF( IWORK( SUBPBS ).GT.SMLSIZ ) THEN
         DO 20 J = SUBPBS, 1, -1
            IWORK( 2*J ) = ( IWORK( J )+1 ) / 2
            IWORK( 2*J-1 ) = IWORK( J ) / 2
   20    CONTINUE
         TLVLS = TLVLS + 1
         SUBPBS = 2*SUBPBS
         GO TO 10
      END IF
      DO 30 J = 2, SUBPBS
         IWORK( J ) = IWORK( J ) + IWORK( J-1 )
   30 CONTINUE
!
!     Divide the matrix into SUBPBS submatrices of size at most SMLSIZ+1
!     using rank-1 modifications (cuts).
!
      SPM1 = SUBPBS - 1
      DO 40 I = 1, SPM1
         SUBMAT = IWORK( I ) + 1
         SMM1 = SUBMAT - 1
         D( SMM1 ) = D( SMM1 ) - ABS( E( SMM1 ) )
         D( SUBMAT ) = D( SUBMAT ) - ABS( E( SMM1 ) )
   40 CONTINUE
!
      INDXQ = 4*N + 3
!
!     Set up workspaces for eigenvalues only/accumulate new vectors
!     routine
!
      TEMP = LOG( DBLE( N ) ) / LOG( TWO )
      LGN = INT( TEMP )
      IF( 2**LGN.LT.N ) LGN = LGN + 1
      IF( 2**LGN.LT.N ) LGN = LGN + 1
      IPRMPT = INDXQ + N + 1
      IPERM = IPRMPT + N*LGN
      IQPTR = IPERM + N*LGN
      IGIVPT = IQPTR + N + 2
      IGIVCL = IGIVPT + N*LGN
!
      IGIVNM = 1
      IQ = IGIVNM + 2*N*LGN
      IWREM = IQ + N**2 + 1
!     Initialize pointers
      DO 50 I = 0, SUBPBS
         IWORK( IPRMPT+I ) = 1
         IWORK( IGIVPT+I ) = 1
   50 CONTINUE
      IWORK( IQPTR ) = 1
!
!     Solve each submatrix eigenproblem at the bottom of the divide and
!     conquer tree.
!
      CURR = 0
      DO 70 I = 0, SPM1
         IF( I.EQ.0 ) THEN
            SUBMAT = 1
            MATSIZ = IWORK( 1 )
         ELSE
            SUBMAT = IWORK( I ) + 1
            MATSIZ = IWORK( I+1 ) - IWORK( I )
         END IF
         LL = IQ - 1 + IWORK( IQPTR+CURR )
         CALL DSTEQR( 'I', MATSIZ, D( SUBMAT ), E( SUBMAT ), RWORK( LL ), MATSIZ, RWORK, INFO )
         CALL ZLACRM( QSIZ, MATSIZ, Q( 1, SUBMAT ), LDQ, RWORK( LL ), MATSIZ, QSTORE( 1, SUBMAT ), LDQS, RWORK( IWREM ) )
         IWORK( IQPTR+CURR+1 ) = IWORK( IQPTR+CURR ) + MATSIZ**2
         CURR = CURR + 1
         IF( INFO.GT.0 ) THEN
            INFO = SUBMAT*( N+1 ) + SUBMAT + MATSIZ - 1
            RETURN
         END IF
         K = 1
         DO 60 J = SUBMAT, IWORK( I+1 )
            IWORK( INDXQ+J ) = K
            K = K + 1
   60    CONTINUE
   70 CONTINUE
!
!     Successively merge eigensystems of adjacent submatrices
!     into eigensystem for the corresponding larger matrix.
!
!     while ( SUBPBS > 1 )
!
      CURLVL = 1
   80 CONTINUE
      IF( SUBPBS.GT.1 ) THEN
         SPM2 = SUBPBS - 2
         DO 90 I = 0, SPM2, 2
            IF( I.EQ.0 ) THEN
               SUBMAT = 1
               MATSIZ = IWORK( 2 )
               MSD2 = IWORK( 1 )
               CURPRB = 0
            ELSE
               SUBMAT = IWORK( I ) + 1
               MATSIZ = IWORK( I+2 ) - IWORK( I )
               MSD2 = MATSIZ / 2
               CURPRB = CURPRB + 1
            END IF
!
!     Merge lower order eigensystems (of size MSD2 and MATSIZ - MSD2)
!     into an eigensystem of size MATSIZ.  ZLAED7 handles the case
!     when the eigenvectors of a full or band Hermitian matrix (which
!     was reduced to tridiagonal form) are desired.
!
!     I am free to use Q as a valuable working space until Loop 150.
!
            CALL ZLAED7( MATSIZ, MSD2, QSIZ, TLVLS, CURLVL, CURPRB, &
                         D( SUBMAT ), QSTORE( 1, SUBMAT ), LDQS, &
                         E( SUBMAT+MSD2-1 ), IWORK( INDXQ+SUBMAT ), &
                         RWORK( IQ ), IWORK( IQPTR ), IWORK( IPRMPT ), &
                         IWORK( IPERM ), IWORK( IGIVPT ), &
                         IWORK( IGIVCL ), RWORK( IGIVNM ), &
                         Q( 1, SUBMAT ), RWORK( IWREM ), IWORK( SUBPBS+1 ), INFO )
            IF( INFO.GT.0 ) THEN
               INFO = SUBMAT*( N+1 ) + SUBMAT + MATSIZ - 1
               RETURN
            END IF
            IWORK( I / 2+1 ) = IWORK( I+2 )
   90    CONTINUE
         SUBPBS = SUBPBS / 2
         CURLVL = CURLVL + 1
         GO TO 80
      END IF
!
!     end while
!
!     Re-merge the eigenvalues/vectors which were deflated at the final
!     merge step.
!
      DO 100 I = 1, N
         J = IWORK( INDXQ+I )
         RWORK( I ) = D( J )
         CALL ZCOPY( QSIZ, QSTORE( 1, J ), 1, Q( 1, I ), 1 )
  100 CONTINUE
      CALL DCOPY( N, RWORK, 1, D, 1 )
!
      RETURN
!
!     End of ZLAED0
!
      END
!> \brief \b ZLACRM
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download ZLACRM + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlacrm.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlacrm.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlacrm.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLACRM( M, N, A, LDA, B, LDB, C, LDC, RWORK )
! 
!       .. Scalar Arguments ..
!       integer(i4)  ::            LDA, LDB, LDC, M, N
!       ..
!       .. Array Arguments ..
!       real(dp) ::   B( LDB, * ), RWORK( * )
!       COMPLEX*16         A( LDA, * ), C( LDC, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLACRM performs a very simple matrix-matrix multiplication:
!>          C := A * B,
!> where A is M by N and complex; B is N by N and real;
!> C is M by N and complex.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is integer(i4)  ::
!>          The number of rows of the matrix A and of the matrix C.
!>          M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is integer(i4)  ::
!>          The number of columns and rows of the matrix B and
!>          the number of columns of the matrix C.
!>          N >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA, N)
!>          A contains the M by N matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is integer(i4)  ::
!>          The leading dimension of the array A. LDA >=max(1,M).
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is real(dp) :: array, dimension (LDB, N)
!>          B contains the N by N matrix B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is integer(i4)  ::
!>          The leading dimension of the array B. LDB >=max(1,N).
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is COMPLEX*16 array, dimension (LDC, N)
!>          C contains the M by N matrix C.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is integer(i4)  ::
!>          The leading dimension of the array C. LDC >=max(1,N).
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is real(dp) :: array, dimension (2*M*N)
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup complex16OTHERauxiliary
!
!  =====================================================================
      SUBROUTINE ZLACRM( M, N, A, LDA, B, LDB, C, LDC, RWORK )
!
!  -- LAPACK auxiliary routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes
      integer(i4)  ::            LDA, LDB, LDC, M, N
!     ..
!     .. Array Arguments ..
      real(dp) ::   B( LDB, * ), RWORK( * )
      COMPLEX(dpc)       A( LDA, * ), C( LDC, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      real(dp) ::   ONE, ZERO
      PARAMETER          ( ONE = 1.0D0, ZERO = 0.0D0 )
!     ..
!     .. Local Scalars ..
      integer(i4)  ::            I, J, L
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, DIMAG
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGEMM
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible.
!
      IF( ( M.EQ.0 ) .OR. ( N.EQ.0 ) ) RETURN
!
      DO 20 J = 1, N
         DO 10 I = 1, M
            RWORK( ( J-1 )*M+I ) = DBLE( A( I, J ) )
   10    CONTINUE
   20 CONTINUE
!
      L = M*N + 1
      CALL DGEMM( 'N', 'N', M, N, N, ONE, RWORK, M, B, LDB, ZERO, RWORK( L ), M )
      DO 40 J = 1, N
         DO 30 I = 1, M
            C( I, J ) = RWORK( L+( J-1 )*M+I-1 )
   30    CONTINUE
   40 CONTINUE
!
      DO 60 J = 1, N
         DO 50 I = 1, M
            RWORK( ( J-1 )*M+I ) = DIMAG( A( I, J ) )
   50    CONTINUE
   60 CONTINUE
      CALL DGEMM( 'N', 'N', M, N, N, ONE, RWORK, M, B, LDB, ZERO, RWORK( L ), M )
      DO 80 J = 1, N
         DO 70 I = 1, M
            C( I, J ) = DCMPLX( DBLE( C( I, J ) ), RWORK( L+( J-1 )*M+I-1 ) )
   70    CONTINUE
   80 CONTINUE
!
      RETURN
!
!     End of ZLACRM
!
      END
!> \brief \b ZLACPY
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download ZLACPY + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlacpy.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlacpy.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlacpy.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLACPY( UPLO, M, N, A, LDA, B, LDB )
! 
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       integer(i4)  ::            LDA, LDB, M, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), B( LDB, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLACPY copies all or part of a two-dimensional matrix A to another
!> matrix B.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies the part of the matrix A to be copied to B.
!>          = 'U':      Upper triangular part
!>          = 'L':      Lower triangular part
!>          Otherwise:  All of the matrix A
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is integer(i4)  ::
!>          The number of rows of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is integer(i4)  ::
!>          The number of columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          The m by n matrix A.  If UPLO = 'U', only the upper trapezium
!>          is accessed; if UPLO = 'L', only the lower trapezium is
!>          accessed.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is integer(i4)  ::
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension (LDB,N)
!>          On exit, B = A in the locations specified by UPLO.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is integer(i4)  ::
!>          The leading dimension of the array B.  LDB >= max(1,M).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup complex16OTHERauxiliary
!
!  =====================================================================
      SUBROUTINE ZLACPY( UPLO, M, N, A, LDA, B, LDB )
!
!  -- LAPACK auxiliary routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes
      CHARACTER          UPLO
      integer(i4)  ::            LDA, LDB, M, N
!     ..
!     .. Array Arguments ..
      COMPLEX(dpc)       A( LDA, * ), B( LDB, * )
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      integer(i4)  ::            I, J
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MIN
!     ..
!     .. Executable Statements ..
!
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 20 J = 1, N
            DO 10 I = 1, MIN( J, M )
               B( I, J ) = A( I, J )
   10       CONTINUE
   20    CONTINUE
!
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
         DO 40 J = 1, N
            DO 30 I = J, M
               B( I, J ) = A( I, J )
   30       CONTINUE
   40    CONTINUE
!
      ELSE
         DO 60 J = 1, N
            DO 50 I = 1, M
               B( I, J ) = A( I, J )
   50       CONTINUE
   60    CONTINUE
      END IF
!
      RETURN
!
!     End of ZLACPY
!
      END
!> \brief \b DSTEBZ
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DSTEDC + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dstedc.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dstedc.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dstedc.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSTEDC( COMPZ, N, D, E, Z, LDZ, WORK, LWORK, IWORK,
!                          LIWORK, INFO )
! 
!       .. Scalar Arguments ..
!       CHARACTER          COMPZ
!       integer(i4)  ::            INFO, LDZ, LIWORK, LWORK, N
!       ..
!       .. Array Arguments ..
!       integer(i4)  ::            IWORK( * )
!       real(dp) ::   D( * ), E( * ), WORK( * ), Z( LDZ, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DSTEDC computes all eigenvalues and, optionally, eigenvectors of a
!> symmetric tridiagonal matrix using the divide and conquer method.
!> The eigenvectors of a full or band real symmetric matrix can also be
!> found if DSYTRD or DSPTRD or DSBTRD has been used to reduce this
!> matrix to tridiagonal form.
!>
!> This code makes very mild assumptions about floating point
!> arithmetic. It will work on machines with a guard digit in
!> add/subtract, or on those binary machines without guard digits
!> which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
!> It could conceivably fail on hexadecimal or decimal machines
!> without guard digits, but we know of none.  See DLAED3 for details.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] COMPZ
!> \verbatim
!>          COMPZ is CHARACTER*1
!>          = 'N':  Compute eigenvalues only.
!>          = 'I':  Compute eigenvectors of tridiagonal matrix also.
!>          = 'V':  Compute eigenvectors of original dense symmetric
!>                  matrix also.  On entry, Z contains the orthogonal
!>                  matrix used to reduce the original matrix to
!>                  tridiagonal form.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is integer(i4)  ::
!>          The dimension of the symmetric tridiagonal matrix.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is real(dp) :: array, dimension (N)
!>          On entry, the diagonal elements of the tridiagonal matrix.
!>          On exit, if INFO = 0, the eigenvalues in ascending order.
!> \endverbatim
!>
!> \param[in,out] E
!> \verbatim
!>          E is real(dp) :: array, dimension (N-1)
!>          On entry, the subdiagonal elements of the tridiagonal matrix.
!>          On exit, E has been destroyed.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is real(dp) :: array, dimension (LDZ,N)
!>          On entry, if COMPZ = 'V', then Z contains the orthogonal
!>          matrix used in the reduction to tridiagonal form.
!>          On exit, if INFO = 0, then if COMPZ = 'V', Z contains the
!>          orthonormal eigenvectors of the original symmetric matrix,
!>          and if COMPZ = 'I', Z contains the orthonormal eigenvectors
!>          of the symmetric tridiagonal matrix.
!>          If  COMPZ = 'N', then Z is not referenced.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is integer(i4)  ::
!>          The leading dimension of the array Z.  LDZ >= 1.
!>          If eigenvectors are desired, then LDZ >= max(1,N).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is real(dp) :: array,
!>                                         dimension (LWORK)
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is integer(i4)  ::
!>          The dimension of the array WORK.
!>          If COMPZ = 'N' or N <= 1 then LWORK must be at least 1.
!>          If COMPZ = 'V' and N > 1 then LWORK must be at least
!>                         ( 1 + 3*N + 2*N*lg N + 4*N**2 ),
!>                         where lg( N ) = smallest integer k such
!>                         that 2**k >= N.
!>          If COMPZ = 'I' and N > 1 then LWORK must be at least
!>                         ( 1 + 4*N + N**2 ).
!>          Note that for COMPZ = 'I' or 'V', then if N is less than or
!>          equal to the minimum divide size, usually 25, then LWORK need
!>          only be max(1,2*(N-1)).
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is integer(i4)  :: array, dimension (MAX(1,LIWORK))
!>          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
!> \endverbatim
!>
!> \param[in] LIWORK
!> \verbatim
!>          LIWORK is integer(i4)  ::
!>          The dimension of the array IWORK.
!>          If COMPZ = 'N' or N <= 1 then LIWORK must be at least 1.
!>          If COMPZ = 'V' and N > 1 then LIWORK must be at least
!>                         ( 6 + 6*N + 5*N*lg N ).
!>          If COMPZ = 'I' and N > 1 then LIWORK must be at least
!>                         ( 3 + 5*N ).
!>          Note that for COMPZ = 'I' or 'V', then if N is less than or
!>          equal to the minimum divide size, usually 25, then LIWORK
!>          need only be 1.
!>
!>          If LIWORK = -1, then a workspace query is assumed; the
!>          routine only calculates the optimal size of the IWORK array,
!>          returns this value as the first entry of the IWORK array, and
!>          no error message related to LIWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is integer(i4)  ::
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  The algorithm failed to compute an eigenvalue while
!>                working on the submatrix lying in rows and columns
!>                INFO/(N+1) through mod(INFO,N+1).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup auxOTHERcomputational
!
!> \par Contributors:
!  ==================
!>
!> Jeff Rutter, Computer Science Division, University of California
!> at Berkeley, USA \n
!>  Modified by Francoise Tisseur, University of Tennessee
!>
!  =====================================================================
      SUBROUTINE DSTEDC( COMPZ, N, D, E, Z, LDZ, WORK, LWORK, IWORK, LIWORK, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes
      CHARACTER          COMPZ
      integer(i4)  ::            INFO, LDZ, LIWORK, LWORK, N
!     ..
!     .. Array Arguments ..
      integer(i4)  ::            IWORK( * )
      real(dp) ::   D( * ), E( * ), WORK( * ), Z( LDZ, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      real(dp) ::   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY
      integer(i4)  ::            FINISH, I, ICOMPZ, II, J, K, LGN, LIWMIN, &
                         LWMIN, M, SMLSIZ, START, STOREZ, STRTRW
      real(dp) ::   EPS, ORGNRM, P, TINY
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      integer(i4)  ::            ILAENV
      real(dp) ::   DLAMCH, DLANST
      EXTERNAL           LSAME, ILAENV, DLAMCH, DLANST
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGEMM, DLACPY, DLAED0, DLASCL, DLASET, DLASRT, &
                         DSTEQR, DSTERF, DSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, INT, LOG, MAX, MOD, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
!
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ICOMPZ = 0
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ICOMPZ = 2
      ELSE
         ICOMPZ = -1
      END IF
      IF( ICOMPZ.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( ( LDZ.LT.1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1, N ) ) ) THEN
         INFO = -6
      END IF
!
      IF( INFO.EQ.0 ) THEN
!
!        Compute the workspace requirements
!
         SMLSIZ = ILAENV( 9, 'DSTEDC', ' ', 0, 0, 0, 0 )
         IF( N.LE.1 .OR. ICOMPZ.EQ.0 ) THEN
            LIWMIN = 1
            LWMIN = 1
         ELSE IF( N.LE.SMLSIZ ) THEN
            LIWMIN = 1
            LWMIN = 2*( N - 1 )
         ELSE
            LGN = INT( LOG( DBLE( N ) )/LOG( TWO ) )
            IF( 2**LGN.LT.N ) LGN = LGN + 1
            IF( 2**LGN.LT.N ) LGN = LGN + 1
            IF( ICOMPZ.EQ.1 ) THEN
               LWMIN = 1 + 3*N + 2*N*LGN + 4*N**2
               LIWMIN = 6 + 6*N + 5*N*LGN
            ELSE IF( ICOMPZ.EQ.2 ) THEN
               LWMIN = 1 + 4*N + N**2
               LIWMIN = 3 + 5*N
            END IF
         END IF
         WORK( 1 ) = LWMIN
         IWORK( 1 ) = LIWMIN
!
         IF( LWORK.LT.LWMIN .AND. .NOT. LQUERY ) THEN
            INFO = -8
         ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT. LQUERY ) THEN
            INFO = -10
         END IF
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSTEDC', -INFO )
         RETURN
      ELSE IF (LQUERY) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) RETURN
      IF( N.EQ.1 ) THEN
         IF( ICOMPZ.NE.0 ) Z( 1, 1 ) = ONE
         RETURN
      END IF
!
!     If the following conditional clause is removed, then the routine
!     will use the Divide and Conquer routine to compute only the
!     eigenvalues, which requires (3N + 3N**2) real workspace and
!     (2 + 5N + 2N lg(N)) integer workspace.
!     Since on many architectures DSTERF is much faster than any other
!     algorithm for finding eigenvalues only, it is used here
!     as the default. If the conditional clause is removed, then
!     information on the size of workspace needs to be changed.
!
!     If COMPZ = 'N', use DSTERF to compute the eigenvalues.
!
      IF( ICOMPZ.EQ.0 ) THEN
         CALL DSTERF( N, D, E, INFO )
         GO TO 50
      END IF
!
!     If N is smaller than the minimum divide size (SMLSIZ+1), then
!     solve the problem with another solver.
!
      IF( N.LE.SMLSIZ ) THEN
!
         CALL DSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
!
      ELSE
!
!        If COMPZ = 'V', the Z matrix must be stored elsewhere for later
!        use.
!
         IF( ICOMPZ.EQ.1 ) THEN
            STOREZ = 1 + N*N
         ELSE
            STOREZ = 1
         END IF
!
         IF( ICOMPZ.EQ.2 ) THEN
            CALL DLASET( 'Full', N, N, ZERO, ONE, Z, LDZ )
         END IF
!
!        Scale.
!
         ORGNRM = DLANST( 'M', N, D, E )
         IF( ORGNRM.EQ.ZERO ) GO TO 50
!
         EPS = DLAMCH( 'Epsilon' )
!
         START = 1
!
!        while ( START <= N )
!
   10    CONTINUE
         IF( START.LE.N ) THEN
!
!           Let FINISH be the position of the next subdiagonal entry
!           such that E( FINISH ) <= TINY or FINISH = N if no such
!           subdiagonal exists.  The matrix identified by the elements
!           between START and FINISH constitutes an independent
!           sub-problem.
!
            FINISH = START
   20       CONTINUE
            IF( FINISH.LT.N ) THEN
               TINY = EPS*SQRT( ABS( D( FINISH ) ) )* SQRT( ABS( D( FINISH+1 ) ) )
               IF( ABS( E( FINISH ) ).GT.TINY ) THEN
                  FINISH = FINISH + 1
                  GO TO 20
               END IF
            END IF
!
!           (Sub) Problem determined.  Compute its size and solve it.
!
            M = FINISH - START + 1
            IF( M.EQ.1 ) THEN
               START = FINISH + 1
               GO TO 10
            END IF
            IF( M.GT.SMLSIZ ) THEN
!
!              Scale.
!
               ORGNRM = DLANST( 'M', M, D( START ), E( START ) )
               CALL DLASCL( 'G', 0, 0, ORGNRM, ONE, M, 1, D( START ), M, INFO )
               CALL DLASCL( 'G', 0, 0, ORGNRM, ONE, M-1, 1, E( START ), M-1, INFO )
!
               IF( ICOMPZ.EQ.1 ) THEN
                  STRTRW = 1
               ELSE
                  STRTRW = START
               END IF
               CALL DLAED0( ICOMPZ, N, M, D( START ), E( START ), &
                            Z( STRTRW, START ), LDZ, WORK( 1 ), N, &
                            WORK( STOREZ ), IWORK, INFO )
               IF( INFO.NE.0 ) THEN
                  INFO = ( INFO / ( M+1 )+START-1 )*( N+1 ) + &
                         MOD( INFO, ( M+1 ) ) + START - 1
                  GO TO 50
               END IF
!
!              Scale back.
!
               CALL DLASCL( 'G', 0, 0, ONE, ORGNRM, M, 1, D( START ), M, INFO )
!
            ELSE
               IF( ICOMPZ.EQ.1 ) THEN
!
!                 Since QR won't update a Z matrix which is larger than
!                 the length of D, we must solve the sub-problem in a
!                 workspace and then multiply back into Z.
!
                  CALL DSTEQR( 'I', M, D( START ), E( START ), WORK, M, WORK( M*M+1 ), INFO )
                  CALL DLACPY( 'A', N, M, Z( 1, START ), LDZ, WORK( STOREZ ), N )
                  CALL DGEMM( 'N', 'N', N, M, M, ONE, WORK( STOREZ ), N, WORK, M, ZERO, Z( 1, START ), LDZ )
               ELSE IF( ICOMPZ.EQ.2 ) THEN
                  CALL DSTEQR( 'I', M, D( START ), E( START ), Z( START, START ), LDZ, WORK, INFO )
               ELSE
                  CALL DSTERF( M, D( START ), E( START ), INFO )
               END IF
               IF( INFO.NE.0 ) THEN
                  INFO = START*( N+1 ) + FINISH
                  GO TO 50
               END IF
            END IF
!
            START = FINISH + 1
            GO TO 10
         END IF
!
!        endwhile
!
!        If the problem split any number of times, then the eigenvalues
!        will not be properly ordered.  Here we permute the eigenvalues
!        (and the associated eigenvectors) into ascending order.
!
         IF( M.NE.N ) THEN
            IF( ICOMPZ.EQ.0 ) THEN
!
!              Use Quick Sort
!
               CALL DLASRT( 'I', N, D, INFO )
!
            ELSE
!
!              Use Selection Sort to minimize swaps of eigenvectors
!
               DO 40 II = 2, N
                  I = II - 1
                  K = I
                  P = D( I )
                  DO 30 J = II, N
                     IF( D( J ).LT.P ) THEN
                        K = J
                        P = D( J )
                     END IF
   30             CONTINUE
                  IF( K.NE.I ) THEN
                     D( K ) = D( I )
                     D( I ) = P
                     CALL DSWAP( N, Z( 1, I ), 1, Z( 1, K ), 1 )
                  END IF
   40          CONTINUE
            END IF
         END IF
      END IF
!
   50 CONTINUE
      WORK( 1 ) = LWMIN
      IWORK( 1 ) = LIWMIN
!
      RETURN
!
!     End of DSTEDC
!
      END
!> \brief \b ZUNMQL
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download ZUNMQL + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zunmql.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zunmql.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zunmql.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZUNMQL( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
!                          WORK, LWORK, INFO )
! 
!       .. Scalar Arguments ..
!       CHARACTER          SIDE, TRANS
!       integer(i4)  ::            INFO, K, LDA, LDC, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZUNMQL overwrites the general complex M-by-N matrix C with
!>
!>                 SIDE = 'L'     SIDE = 'R'
!> TRANS = 'N':      Q * C          C * Q
!> TRANS = 'C':      Q**H * C       C * Q**H
!>
!> where Q is a complex unitary matrix defined as the product of k
!> elementary reflectors
!>
!>       Q = H(k) . . . H(2) H(1)
!>
!> as returned by ZGEQLF. Q is of order M if SIDE = 'L' and of order N
!> if SIDE = 'R'.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'L': apply Q or Q**H from the Left;
!>          = 'R': apply Q or Q**H from the Right.
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          = 'N':  No transpose, apply Q;
!>          = 'C':  Transpose, apply Q**H.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is integer(i4)  ::
!>          The number of rows of the matrix C. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is integer(i4)  ::
!>          The number of columns of the matrix C. N >= 0.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is integer(i4)  ::
!>          The number of elementary reflectors whose product defines
!>          the matrix Q.
!>          If SIDE = 'L', M >= K >= 0;
!>          if SIDE = 'R', N >= K >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,K)
!>          The i-th column must contain the vector which defines the
!>          elementary reflector H(i), for i = 1,2,...,k, as returned by
!>          ZGEQLF in the last k columns of its array argument A.
!>          A is modified by the routine but restored on exit.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is integer(i4)  ::
!>          The leading dimension of the array A.
!>          If SIDE = 'L', LDA >= max(1,M);
!>          if SIDE = 'R', LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is COMPLEX*16 array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by ZGEQLF.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is COMPLEX*16 array, dimension (LDC,N)
!>          On entry, the M-by-N matrix C.
!>          On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is integer(i4)  ::
!>          The leading dimension of the array C. LDC >= max(1,M).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is integer(i4)  ::
!>          The dimension of the array WORK.
!>          If SIDE = 'L', LWORK >= max(1,N);
!>          if SIDE = 'R', LWORK >= max(1,M).
!>          For optimum performance LWORK >= N*NB if SIDE = 'L', and
!>          LWORK >= M*NB if SIDE = 'R', where NB is the optimal
!>          blocksize.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is integer(i4)  ::
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup complex16OTHERcomputational
!
!  =====================================================================
      SUBROUTINE ZUNMQL( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes, only : i4, dpc
      CHARACTER          SIDE, TRANS
      integer(i4)  ::            INFO, K, LDA, LDC, LWORK, M, N
!     ..
!     .. Array Arguments ..
      COMPLEX(dpc)       A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      integer(i4)  ::            NBMAX, LDT
      PARAMETER          ( NBMAX = 64, LDT = NBMAX+1 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LEFT, LQUERY, NOTRAN
      integer(i4)  ::            I, I1, I2, I3, IB, IINFO, IWS, LDWORK, LWKOPT, MI, NB, NBMIN, NI, NQ, NW
!     ..
!     .. Local Arrays ..
      COMPLEX(dpc)       T( LDT, NBMAX )
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      integer(i4)  ::            ILAENV
      EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARFB, ZLARFT, ZUNM2L
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
      LQUERY = ( LWORK.EQ.-1 )
!
!     NQ is the order of Q and NW is the minimum dimension of WORK
!
      IF( LEFT ) THEN
         NQ = M
         NW = MAX( 1, N )
      ELSE
         NQ = N
         NW = MAX( 1, M )
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      END IF
!
      IF( INFO.EQ.0 ) THEN
         IF( M.EQ.0 .OR. N.EQ.0 ) THEN
            LWKOPT = 1
         ELSE
!
!           Determine the block size.  NB may be at most NBMAX, where
!           NBMAX is used to define the local array T.
!
            NB = MIN( NBMAX, ILAENV( 1, 'ZUNMQL', SIDE // TRANS, M, N, K, -1 ) )
            LWKOPT = NW*NB
         END IF
         WORK( 1 ) = LWKOPT
!
         IF( LWORK.LT.NW .AND. .NOT.LQUERY ) THEN
            INFO = -12
         END IF
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNMQL', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 ) THEN
         RETURN
      END IF
!
      NBMIN = 2
      LDWORK = NW
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
         IWS = NW*NB
         IF( LWORK.LT.IWS ) THEN
            NB = LWORK / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'ZUNMQL', SIDE // TRANS, M, N, K, -1 ) )
         END IF
      ELSE
         IWS = NW
      END IF
!
      IF( NB.LT.NBMIN .OR. NB.GE.K ) THEN
!
!        Use unblocked code
!
         CALL ZUNM2L( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, IINFO )
      ELSE
!
!        Use blocked code
!
         IF( ( LEFT .AND. NOTRAN ) .OR. ( .NOT.LEFT .AND. .NOT.NOTRAN ) ) THEN
            I1 = 1
            I2 = K
            I3 = NB
         ELSE
            I1 = ( ( K-1 ) / NB )*NB + 1
            I2 = 1
            I3 = -NB
         END IF
!
         IF( LEFT ) THEN
            NI = N
         ELSE
            MI = M
         END IF
!
         DO 10 I = I1, I2, I3
            IB = MIN( NB, K-I+1 )
!
!           Form the triangular factor of the block reflector
!           H = H(i+ib-1) . . . H(i+1) H(i)
!
            CALL ZLARFT( 'Backward', 'Columnwise', NQ-K+I+IB-1, IB, A( 1, I ), LDA, TAU( I ), T, LDT )
            IF( LEFT ) THEN
!
!              H or H**H is applied to C(1:m-k+i+ib-1,1:n)
!
               MI = M - K + I + IB - 1
            ELSE
!
!              H or H**H is applied to C(1:m,1:n-k+i+ib-1)
!
               NI = N - K + I + IB - 1
            END IF
!
!           Apply H or H**H
!
            CALL ZLARFB( SIDE, TRANS, 'Backward', 'Columnwise', MI, NI, &
                         IB, A( 1, I ), LDA, T, LDT, C, LDC, WORK, LDWORK )
   10    CONTINUE
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
!
!     End of ZUNMQL
!
      END
!> \brief \b ZUNMQR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download ZUNMQR + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zunmqr.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zunmqr.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zunmqr.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZUNMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
!                          WORK, LWORK, INFO )
! 
!       .. Scalar Arguments ..
!       CHARACTER          SIDE, TRANS
!       integer(i4)  ::            INFO, K, LDA, LDC, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZUNMQR overwrites the general complex M-by-N matrix C with
!>
!>                 SIDE = 'L'     SIDE = 'R'
!> TRANS = 'N':      Q * C          C * Q
!> TRANS = 'C':      Q**H * C       C * Q**H
!>
!> where Q is a complex unitary matrix defined as the product of k
!> elementary reflectors
!>
!>       Q = H(1) H(2) . . . H(k)
!>
!> as returned by ZGEQRF. Q is of order M if SIDE = 'L' and of order N
!> if SIDE = 'R'.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'L': apply Q or Q**H from the Left;
!>          = 'R': apply Q or Q**H from the Right.
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          = 'N':  No transpose, apply Q;
!>          = 'C':  Conjugate transpose, apply Q**H.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is integer(i4)  ::
!>          The number of rows of the matrix C. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is integer(i4)  ::
!>          The number of columns of the matrix C. N >= 0.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is integer(i4)  ::
!>          The number of elementary reflectors whose product defines
!>          the matrix Q.
!>          If SIDE = 'L', M >= K >= 0;
!>          if SIDE = 'R', N >= K >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,K)
!>          The i-th column must contain the vector which defines the
!>          elementary reflector H(i), for i = 1,2,...,k, as returned by
!>          ZGEQRF in the first k columns of its array argument A.
!>          A is modified by the routine but restored on exit.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is integer(i4)  ::
!>          The leading dimension of the array A.
!>          If SIDE = 'L', LDA >= max(1,M);
!>          if SIDE = 'R', LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is COMPLEX*16 array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by ZGEQRF.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is COMPLEX*16 array, dimension (LDC,N)
!>          On entry, the M-by-N matrix C.
!>          On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is integer(i4)  ::
!>          The leading dimension of the array C. LDC >= max(1,M).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is integer(i4)  ::
!>          The dimension of the array WORK.
!>          If SIDE = 'L', LWORK >= max(1,N);
!>          if SIDE = 'R', LWORK >= max(1,M).
!>          For optimum performance LWORK >= N*NB if SIDE = 'L', and
!>          LWORK >= M*NB if SIDE = 'R', where NB is the optimal
!>          blocksize.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is integer(i4)  ::
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup complex16OTHERcomputational
!
!  =====================================================================
      SUBROUTINE ZUNMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes, only : i4, dpc
      CHARACTER          SIDE, TRANS
      integer(i4)  ::            INFO, K, LDA, LDC, LWORK, M, N
!     ..
!     .. Array Arguments ..
      COMPLEX(dpc)       A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      integer(i4)  ::            NBMAX, LDT
      PARAMETER          ( NBMAX = 64, LDT = NBMAX+1 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LEFT, LQUERY, NOTRAN
      integer(i4)  ::            I, I1, I2, I3, IB, IC, IINFO, IWS, JC, LDWORK, LWKOPT, MI, NB, NBMIN, NI, NQ, NW
!     ..
!     .. Local Arrays ..
      COMPLEX(dpc)       T( LDT, NBMAX )
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      integer(i4)  ::            ILAENV
      EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARFB, ZLARFT, ZUNM2R
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
      LQUERY = ( LWORK.EQ.-1 )
!
!     NQ is the order of Q and NW is the minimum dimension of WORK
!
      IF( LEFT ) THEN
         NQ = M
         NW = N
      ELSE
         NQ = N
         NW = M
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF( LWORK.LT.MAX( 1, NW ) .AND. .NOT.LQUERY ) THEN
         INFO = -12
      END IF
!
      IF( INFO.EQ.0 ) THEN
!
!        Determine the block size.  NB may be at most NBMAX, where NBMAX
!        is used to define the local array T.
!
         NB = MIN( NBMAX, ILAENV( 1, 'ZUNMQR', SIDE // TRANS, M, N, K, -1 ) )
         LWKOPT = MAX( 1, NW )*NB
         WORK( 1 ) = LWKOPT
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNMQR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!
      NBMIN = 2
      LDWORK = NW
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
         IWS = NW*NB
         IF( LWORK.LT.IWS ) THEN
            NB = LWORK / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'ZUNMQR', SIDE // TRANS, M, N, K, -1 ) )
         END IF
      ELSE
         IWS = NW
      END IF
!
      IF( NB.LT.NBMIN .OR. NB.GE.K ) THEN
!
!        Use unblocked code
!
         CALL ZUNM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, IINFO )
      ELSE
!
!        Use blocked code
!
         IF( ( LEFT .AND. .NOT.NOTRAN ) .OR. ( .NOT.LEFT .AND. NOTRAN ) ) THEN
            I1 = 1
            I2 = K
            I3 = NB
         ELSE
            I1 = ( ( K-1 ) / NB )*NB + 1
            I2 = 1
            I3 = -NB
         END IF
!
         IF( LEFT ) THEN
            NI = N
            JC = 1
         ELSE
            MI = M
            IC = 1
         END IF
!
         DO 10 I = I1, I2, I3
            IB = MIN( NB, K-I+1 )
!
!           Form the triangular factor of the block reflector
!           H = H(i) H(i+1) . . . H(i+ib-1)
!
            CALL ZLARFT( 'Forward', 'Columnwise', NQ-I+1, IB, A( I, I ), LDA, TAU( I ), T, LDT )
            IF( LEFT ) THEN
!
!              H or H**H is applied to C(i:m,1:n)
!
               MI = M - I + 1
               IC = I
            ELSE
!
!              H or H**H is applied to C(1:m,i:n)
!
               NI = N - I + 1
               JC = I
            END IF
!
!           Apply H or H**H
!
            CALL ZLARFB( SIDE, TRANS, 'Forward', 'Columnwise', MI, NI, &
                         IB, A( I, I ), LDA, T, LDT, C( IC, JC ), LDC, WORK, LDWORK )
   10    CONTINUE
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
!
!     End of ZUNMQR
!
      END
!> \brief \b ZLAED7
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download ZLAED7 + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaed7.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaed7.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaed7.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLAED7( N, CUTPNT, QSIZ, TLVLS, CURLVL, CURPBM, D, Q,
!                          LDQ, RHO, INDXQ, QSTORE, QPTR, PRMPTR, PERM,
!                          GIVPTR, GIVCOL, GIVNUM, WORK, RWORK, IWORK,
!                          INFO )
! 
!       .. Scalar Arguments ..
!       integer(i4)  ::            CURLVL, CURPBM, CUTPNT, INFO, LDQ, N, QSIZ,
!      $                   TLVLS
!       real(dp)    ::   RHO
!       ..
!       .. Array Arguments ..
!       integer(i4)  ::            GIVCOL( 2, * ), GIVPTR( * ), INDXQ( * ),
!      $                   IWORK( * ), PERM( * ), PRMPTR( * ), QPTR( * )
!       real(dp)    ::   D( * ), GIVNUM( 2, * ), QSTORE( * ), RWORK( * )
!       COMPLEX*16         Q( LDQ, * ), WORK( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLAED7 computes the updated eigensystem of a diagonal
!> matrix after modification by a rank-one symmetric matrix. This
!> routine is used only for the eigenproblem which requires all
!> eigenvalues and optionally eigenvectors of a dense or banded
!> Hermitian matrix that has been reduced to tridiagonal form.
!>
!>   T = Q(in) ( D(in) + RHO * Z*Z**H ) Q**H(in) = Q(out) * D(out) * Q**H(out)
!>
!>   where Z = Q**Hu, u is a vector of length N with ones in the
!>   CUTPNT and CUTPNT + 1 th elements and zeros elsewhere.
!>
!>    The eigenvectors of the original matrix are stored in Q, and the
!>    eigenvalues are in D.  The algorithm consists of three stages:
!>
!>       The first stage consists of deflating the size of the problem
!>       when there are multiple eigenvalues or if there is a zero in
!>       the Z vector.  For each such occurence the dimension of the
!>       secular equation problem is reduced by one.  This stage is
!>       performed by the routine DLAED2.
!>
!>       The second stage consists of calculating the updated
!>       eigenvalues. This is done by finding the roots of the secular
!>       equation via the routine DLAED4 (as called by SLAED3).
!>       This routine also calculates the eigenvectors of the current
!>       problem.
!>
!>       The final stage consists of computing the updated eigenvectors
!>       directly using the updated eigenvalues.  The eigenvectors for
!>       the current problem are multiplied with the eigenvectors from
!>       the overall problem.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is integer(i4)  ::
!>         The dimension of the symmetric tridiagonal matrix.  N >= 0.
!> \endverbatim
!>
!> \param[in] CUTPNT
!> \verbatim
!>          CUTPNT is integer(i4)  ::
!>         Contains the location of the last eigenvalue in the leading
!>         sub-matrix.  min(1,N) <= CUTPNT <= N.
!> \endverbatim
!>
!> \param[in] QSIZ
!> \verbatim
!>          QSIZ is integer(i4)  ::
!>         The dimension of the unitary matrix used to reduce
!>         the full matrix to tridiagonal form.  QSIZ >= N.
!> \endverbatim
!>
!> \param[in] TLVLS
!> \verbatim
!>          TLVLS is integer(i4)  ::
!>         The total number of merging levels in the overall divide and
!>         conquer tree.
!> \endverbatim
!>
!> \param[in] CURLVL
!> \verbatim
!>          CURLVL is integer(i4)  ::
!>         The current level in the overall merge routine,
!>         0 <= curlvl <= tlvls.
!> \endverbatim
!>
!> \param[in] CURPBM
!> \verbatim
!>          CURPBM is integer(i4)  ::
!>         The current problem in the current level in the overall
!>         merge routine (counting from upper left to lower right).
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is real(dp)    :: array, dimension (N)
!>         On entry, the eigenvalues of the rank-1-perturbed matrix.
!>         On exit, the eigenvalues of the repaired matrix.
!> \endverbatim
!>
!> \param[in,out] Q
!> \verbatim
!>          Q is COMPLEX*16 array, dimension (LDQ,N)
!>         On entry, the eigenvectors of the rank-1-perturbed matrix.
!>         On exit, the eigenvectors of the repaired tridiagonal matrix.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is integer(i4)  ::
!>         The leading dimension of the array Q.  LDQ >= max(1,N).
!> \endverbatim
!>
!> \param[in] RHO
!> \verbatim
!>          RHO is real(dp)    ::
!>         Contains the subdiagonal element used to create the rank-1
!>         modification.
!> \endverbatim
!>
!> \param[out] INDXQ
!> \verbatim
!>          INDXQ is integer(i4)  :: array, dimension (N)
!>         This contains the permutation which will reintegrate the
!>         subproblem just solved back into sorted order,
!>         ie. D( INDXQ( I = 1, N ) ) will be in ascending order.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is integer(i4)  :: array, dimension (4*N)
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is real(dp)    :: array,
!>                                 dimension (3*N+2*QSIZ*N)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (QSIZ*N)
!> \endverbatim
!>
!> \param[in,out] QSTORE
!> \verbatim
!>          QSTORE is real(dp)    :: array, dimension (N**2+1)
!>         Stores eigenvectors of submatrices encountered during
!>         divide and conquer, packed together. QPTR points to
!>         beginning of the submatrices.
!> \endverbatim
!>
!> \param[in,out] QPTR
!> \verbatim
!>          QPTR is integer(i4)  :: array, dimension (N+2)
!>         List of indices pointing to beginning of submatrices stored
!>         in QSTORE. The submatrices are numbered starting at the
!>         bottom left of the divide and conquer tree, from left to
!>         right and bottom to top.
!> \endverbatim
!>
!> \param[in] PRMPTR
!> \verbatim
!>          PRMPTR is integer(i4)  :: array, dimension (N lg N)
!>         Contains a list of pointers which indicate where in PERM a
!>         level's permutation is stored.  PRMPTR(i+1) - PRMPTR(i)
!>         indicates the size of the permutation and also the size of
!>         the full, non-deflated problem.
!> \endverbatim
!>
!> \param[in] PERM
!> \verbatim
!>          PERM is integer(i4)  :: array, dimension (N lg N)
!>         Contains the permutations (from deflation and sorting) to be
!>         applied to each eigenblock.
!> \endverbatim
!>
!> \param[in] GIVPTR
!> \verbatim
!>          GIVPTR is integer(i4)  :: array, dimension (N lg N)
!>         Contains a list of pointers which indicate where in GIVCOL a
!>         level's Givens rotations are stored.  GIVPTR(i+1) - GIVPTR(i)
!>         indicates the number of Givens rotations.
!> \endverbatim
!>
!> \param[in] GIVCOL
!> \verbatim
!>          GIVCOL is integer(i4)  :: array, dimension (2, N lg N)
!>         Each pair of numbers indicates a pair of columns to take place
!>         in a Givens rotation.
!> \endverbatim
!>
!> \param[in] GIVNUM
!> \verbatim
!>          GIVNUM is real(dp)    :: array, dimension (2, N lg N)
!>         Each number indicates the S value to be used in the
!>         corresponding Givens rotation.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is integer(i4)  ::
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  if INFO = 1, an eigenvalue did not converge
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup complex16OTHERcomputational
!
!  =====================================================================
      SUBROUTINE ZLAED7( N, CUTPNT, QSIZ, TLVLS, CURLVL, CURPBM, D, Q, &
                         LDQ, RHO, INDXQ, QSTORE, QPTR, PRMPTR, PERM, &
                         GIVPTR, GIVCOL, GIVNUM, WORK, RWORK, IWORK, &
                         INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes, only : i4, dp, dpc
      integer(i4)  ::            CURLVL, CURPBM, CUTPNT, INFO, LDQ, N, QSIZ, TLVLS
      real(dp)    ::   RHO
!     ..
!     .. Array Arguments ..
      integer(i4)  ::            GIVCOL( 2, * ), GIVPTR( * ), INDXQ( * ), IWORK( * ), PERM( * ), PRMPTR( * ), QPTR( * )
      real(dp)    ::   D( * ), GIVNUM( 2, * ), QSTORE( * ), RWORK( * )
      COMPLEX(dpc)       Q( LDQ, * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      integer(i4)  ::            COLTYP, CURR, I, IDLMDA, INDX, INDXC, INDXP, IQ, IW, IZ, K, N1, N2, PTR
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLAED9, DLAEDA, DLAMRG, XERBLA, ZLACRM, ZLAED8
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
!
!     IF( ICOMPQ.LT.0 .OR. ICOMPQ.GT.1 ) THEN
!        INFO = -1
!     ELSE IF( N.LT.0 ) THEN
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( MIN( 1, N ).GT.CUTPNT .OR. N.LT.CUTPNT ) THEN
         INFO = -2
      ELSE IF( QSIZ.LT.N ) THEN
         INFO = -3
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZLAED7', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) RETURN
!
!     The following values are for bookkeeping purposes only.  They are
!     integer pointers which indicate the portion of the workspace
!     used by a particular array in DLAED2 and SLAED3.
!
      IZ = 1
      IDLMDA = IZ + N
      IW = IDLMDA + N
      IQ = IW + N
!
      INDX = 1
      INDXC = INDX + N
      COLTYP = INDXC + N
      INDXP = COLTYP + N
!
!     Form the z-vector which consists of the last row of Q_1 and the
!     first row of Q_2.
!
      PTR = 1 + 2**TLVLS
      DO 10 I = 1, CURLVL - 1
         PTR = PTR + 2**( TLVLS-I )
   10 CONTINUE
      CURR = PTR + CURPBM
      CALL DLAEDA( N, TLVLS, CURLVL, CURPBM, PRMPTR, PERM, GIVPTR, GIVCOL, GIVNUM, QSTORE, QPTR, RWORK( IZ ), &
                   RWORK( IZ+N ), INFO )
!
!     When solving the final problem, we no longer need the stored data,
!     so we will overwrite the data from this level onto the previously
!     used storage space.
!
      IF( CURLVL.EQ.TLVLS ) THEN
         QPTR( CURR ) = 1
         PRMPTR( CURR ) = 1
         GIVPTR( CURR ) = 1
      END IF
!
!     Sort and Deflate eigenvalues.
!
      CALL ZLAED8( K, N, QSIZ, Q, LDQ, D, RHO, CUTPNT, RWORK( IZ ), RWORK( IDLMDA ), WORK, QSIZ, RWORK( IW ), &
                   IWORK( INDXP ), IWORK( INDX ), INDXQ, PERM( PRMPTR( CURR ) ), GIVPTR( CURR+1 ), &
                   GIVCOL( 1, GIVPTR( CURR ) ), GIVNUM( 1, GIVPTR( CURR ) ), INFO )
      PRMPTR( CURR+1 ) = PRMPTR( CURR ) + N
      GIVPTR( CURR+1 ) = GIVPTR( CURR+1 ) + GIVPTR( CURR )
!
!     Solve Secular Equation.
!
      IF( K.NE.0 ) THEN
         CALL DLAED9( K, 1, K, N, D, RWORK( IQ ), K, RHO, RWORK( IDLMDA ), RWORK( IW ), QSTORE( QPTR( CURR ) ), K, INFO )
         CALL ZLACRM( QSIZ, K, WORK, QSIZ, QSTORE( QPTR( CURR ) ), K, Q, LDQ, RWORK( IQ ) )
         QPTR( CURR+1 ) = QPTR( CURR ) + K**2
         IF( INFO.NE.0 ) THEN
            RETURN
         END IF
!
!     Prepare the INDXQ sorting premutation.
!
         N1 = K
         N2 = N - K
         CALL DLAMRG( N1, N2, D, 1, -1, INDXQ )
      ELSE
         QPTR( CURR+1 ) = QPTR( CURR )
         DO 20 I = 1, N
            INDXQ( I ) = I
   20    CONTINUE
      END IF
!
      RETURN
!
!     End of ZLAED7
!
      END
!> \brief \b DLAED0
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DLAED0 + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaed0.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaed0.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaed0.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAED0( ICOMPQ, QSIZ, N, D, E, Q, LDQ, QSTORE, LDQS,
!                          WORK, IWORK, INFO )
! 
!       .. Scalar Arguments ..
!       integer(i4)  ::            ICOMPQ, INFO, LDQ, LDQS, N, QSIZ
!       ..
!       .. Array Arguments ..
!       integer(i4)  ::            IWORK( * )
!       real(dp)    ::   D( * ), E( * ), Q( LDQ, * ), QSTORE( LDQS, * ),
!      $                   WORK( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAED0 computes all eigenvalues and corresponding eigenvectors of a
!> symmetric tridiagonal matrix using the divide and conquer method.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ICOMPQ
!> \verbatim
!>          ICOMPQ is integer(i4)  ::
!>          = 0:  Compute eigenvalues only.
!>          = 1:  Compute eigenvectors of original dense symmetric matrix
!>                also.  On entry, Q contains the orthogonal matrix used
!>                to reduce the original matrix to tridiagonal form.
!>          = 2:  Compute eigenvalues and eigenvectors of tridiagonal
!>                matrix.
!> \endverbatim
!>
!> \param[in] QSIZ
!> \verbatim
!>          QSIZ is integer(i4)  ::
!>         The dimension of the orthogonal matrix used to reduce
!>         the full matrix to tridiagonal form.  QSIZ >= N if ICOMPQ = 1.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is integer(i4)  ::
!>         The dimension of the symmetric tridiagonal matrix.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is real(dp)    :: array, dimension (N)
!>         On entry, the main diagonal of the tridiagonal matrix.
!>         On exit, its eigenvalues.
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is real(dp)    :: array, dimension (N-1)
!>         The off-diagonal elements of the tridiagonal matrix.
!>         On exit, E has been destroyed.
!> \endverbatim
!>
!> \param[in,out] Q
!> \verbatim
!>          Q is real(dp)    :: array, dimension (LDQ, N)
!>         On entry, Q must contain an N-by-N orthogonal matrix.
!>         If ICOMPQ = 0    Q is not referenced.
!>         If ICOMPQ = 1    On entry, Q is a subset of the columns of the
!>                          orthogonal matrix used to reduce the full
!>                          matrix to tridiagonal form corresponding to
!>                          the subset of the full matrix which is being
!>                          decomposed at this time.
!>         If ICOMPQ = 2    On entry, Q will be the identity matrix.
!>                          On exit, Q contains the eigenvectors of the
!>                          tridiagonal matrix.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is integer(i4)  ::
!>         The leading dimension of the array Q.  If eigenvectors are
!>         desired, then  LDQ >= max(1,N).  In any case,  LDQ >= 1.
!> \endverbatim
!>
!> \param[out] QSTORE
!> \verbatim
!>          QSTORE is real(dp)    :: array, dimension (LDQS, N)
!>         Referenced only when ICOMPQ = 1.  Used to store parts of
!>         the eigenvector matrix when the updating matrix multiplies
!>         take place.
!> \endverbatim
!>
!> \param[in] LDQS
!> \verbatim
!>          LDQS is integer(i4)  ::
!>         The leading dimension of the array QSTORE.  If ICOMPQ = 1,
!>         then  LDQS >= max(1,N).  In any case,  LDQS >= 1.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is real(dp)    :: array,
!>         If ICOMPQ = 0 or 1, the dimension of WORK must be at least
!>                     1 + 3*N + 2*N*lg N + 3*N**2
!>                     ( lg( N ) = smallest integer k
!>                                 such that 2^k >= N )
!>         If ICOMPQ = 2, the dimension of WORK must be at least
!>                     4*N + N**2.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is integer(i4)  :: array,
!>         If ICOMPQ = 0 or 1, the dimension of IWORK must be at least
!>                        6 + 6*N + 5*N*lg N.
!>                        ( lg( N ) = smallest integer k
!>                                    such that 2^k >= N )
!>         If ICOMPQ = 2, the dimension of IWORK must be at least
!>                        3 + 5*N.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is integer(i4)  ::
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  The algorithm failed to compute an eigenvalue while
!>                working on the submatrix lying in rows and columns
!>                INFO/(N+1) through mod(INFO,N+1).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup auxOTHERcomputational
!
!> \par Contributors:
!  ==================
!>
!> Jeff Rutter, Computer Science Division, University of California
!> at Berkeley, USA
!
!  =====================================================================
      SUBROUTINE DLAED0( ICOMPQ, QSIZ, N, D, E, Q, LDQ, QSTORE, LDQS, WORK, IWORK, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes, only : i4, dp, dpc
      integer(i4)  ::            ICOMPQ, INFO, LDQ, LDQS, N, QSIZ
!     ..
!     .. Array Arguments ..
      integer(i4)  ::            IWORK( * )
      real(dp)    ::   D( * ), E( * ), Q( LDQ, * ), QSTORE( LDQS, * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      real(dp)    ::   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.D0, ONE = 1.D0, TWO = 2.D0 )
!     ..
!     .. Local Scalars ..
      integer(i4)  ::            CURLVL, CURPRB, CURR, I, IGIVCL, IGIVNM, IGIVPT, INDXQ, IPERM, IPRMPT, IQ, IQPTR, IWREM, &
                         J, K, LGN, MATSIZ, MSD2, SMLSIZ, SMM1, SPM1, SPM2, SUBMAT, SUBPBS, TLVLS
      real(dp)    ::   TEMP
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMM, DLACPY, DLAED1, DLAED7, DSTEQR, XERBLA
!     ..
!     .. External Functions ..
      integer(i4)  ::            ILAENV
      EXTERNAL           ILAENV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, INT, LOG, MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
!
      IF( ICOMPQ.LT.0 .OR. ICOMPQ.GT.2 ) THEN
         INFO = -1
      ELSE IF( ( ICOMPQ.EQ.1 ) .AND. ( QSIZ.LT.MAX( 0, N ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDQS.LT.MAX( 1, N ) ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLAED0', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) RETURN
!
      SMLSIZ = ILAENV( 9, 'DLAED0', ' ', 0, 0, 0, 0 )
!
!     Determine the size and placement of the submatrices, and save in
!     the leading elements of IWORK.
!
      IWORK( 1 ) = N
      SUBPBS = 1
      TLVLS = 0
   10 CONTINUE
      IF( IWORK( SUBPBS ).GT.SMLSIZ ) THEN
         DO 20 J = SUBPBS, 1, -1
            IWORK( 2*J ) = ( IWORK( J )+1 ) / 2
            IWORK( 2*J-1 ) = IWORK( J ) / 2
   20    CONTINUE
         TLVLS = TLVLS + 1
         SUBPBS = 2*SUBPBS
         GO TO 10
      END IF
      DO 30 J = 2, SUBPBS
         IWORK( J ) = IWORK( J ) + IWORK( J-1 )
   30 CONTINUE
!
!     Divide the matrix into SUBPBS submatrices of size at most SMLSIZ+1
!     using rank-1 modifications (cuts).
!
      SPM1 = SUBPBS - 1
      DO 40 I = 1, SPM1
         SUBMAT = IWORK( I ) + 1
         SMM1 = SUBMAT - 1
         D( SMM1 ) = D( SMM1 ) - ABS( E( SMM1 ) )
         D( SUBMAT ) = D( SUBMAT ) - ABS( E( SMM1 ) )
   40 CONTINUE
!
      INDXQ = 4*N + 3
      IF( ICOMPQ.NE.2 ) THEN
!
!        Set up workspaces for eigenvalues only/accumulate new vectors
!        routine
!
         TEMP = LOG( DBLE( N ) ) / LOG( TWO )
         LGN = INT( TEMP )
         IF( 2**LGN.LT.N ) LGN = LGN + 1
         IF( 2**LGN.LT.N ) LGN = LGN + 1
         IPRMPT = INDXQ + N + 1
         IPERM = IPRMPT + N*LGN
         IQPTR = IPERM + N*LGN
         IGIVPT = IQPTR + N + 2
         IGIVCL = IGIVPT + N*LGN
!
         IGIVNM = 1
         IQ = IGIVNM + 2*N*LGN
         IWREM = IQ + N**2 + 1
!
!        Initialize pointers
!
         DO 50 I = 0, SUBPBS
            IWORK( IPRMPT+I ) = 1
            IWORK( IGIVPT+I ) = 1
   50    CONTINUE
         IWORK( IQPTR ) = 1
      END IF
!
!     Solve each submatrix eigenproblem at the bottom of the divide and
!     conquer tree.
!
      CURR = 0
      DO 70 I = 0, SPM1
         IF( I.EQ.0 ) THEN
            SUBMAT = 1
            MATSIZ = IWORK( 1 )
         ELSE
            SUBMAT = IWORK( I ) + 1
            MATSIZ = IWORK( I+1 ) - IWORK( I )
         END IF
         IF( ICOMPQ.EQ.2 ) THEN
            CALL DSTEQR( 'I', MATSIZ, D( SUBMAT ), E( SUBMAT ), Q( SUBMAT, SUBMAT ), LDQ, WORK, INFO )
            IF( INFO.NE.0 ) GO TO 130
         ELSE
            CALL DSTEQR( 'I', MATSIZ, D( SUBMAT ), E( SUBMAT ), &
                         WORK( IQ-1+IWORK( IQPTR+CURR ) ), MATSIZ, WORK, INFO )
            IF( INFO.NE.0 ) GO TO 130
            IF( ICOMPQ.EQ.1 ) THEN
               CALL DGEMM( 'N', 'N', QSIZ, MATSIZ, MATSIZ, ONE, Q( 1, SUBMAT ), LDQ, WORK( IQ-1+IWORK( IQPTR+ &
                           CURR ) ), MATSIZ, ZERO, QSTORE( 1, SUBMAT ), LDQS )
            END IF
            IWORK( IQPTR+CURR+1 ) = IWORK( IQPTR+CURR ) + MATSIZ**2
            CURR = CURR + 1
         END IF
         K = 1
         DO 60 J = SUBMAT, IWORK( I+1 )
            IWORK( INDXQ+J ) = K
            K = K + 1
   60    CONTINUE
   70 CONTINUE
!
!     Successively merge eigensystems of adjacent submatrices
!     into eigensystem for the corresponding larger matrix.
!
!     while ( SUBPBS > 1 )
!
      CURLVL = 1
   80 CONTINUE
      IF( SUBPBS.GT.1 ) THEN
         SPM2 = SUBPBS - 2
         DO 90 I = 0, SPM2, 2
            IF( I.EQ.0 ) THEN
               SUBMAT = 1
               MATSIZ = IWORK( 2 )
               MSD2 = IWORK( 1 )
               CURPRB = 0
            ELSE
               SUBMAT = IWORK( I ) + 1
               MATSIZ = IWORK( I+2 ) - IWORK( I )
               MSD2 = MATSIZ / 2
               CURPRB = CURPRB + 1
            END IF
!
!     Merge lower order eigensystems (of size MSD2 and MATSIZ - MSD2)
!     into an eigensystem of size MATSIZ.
!     DLAED1 is used only for the full eigensystem of a tridiagonal
!     matrix.
!     DLAED7 handles the cases in which eigenvalues only or eigenvalues
!     and eigenvectors of a full symmetric matrix (which was reduced to
!     tridiagonal form) are desired.
!
            IF( ICOMPQ.EQ.2 ) THEN
               CALL DLAED1( MATSIZ, D( SUBMAT ), Q( SUBMAT, SUBMAT ), LDQ, IWORK( INDXQ+SUBMAT ), &
                            E( SUBMAT+MSD2-1 ), MSD2, WORK, IWORK( SUBPBS+1 ), INFO )
            ELSE
               CALL DLAED7( ICOMPQ, MATSIZ, QSIZ, TLVLS, CURLVL, CURPRB, D( SUBMAT ), QSTORE( 1, SUBMAT ), LDQS, &
                            IWORK( INDXQ+SUBMAT ), E( SUBMAT+MSD2-1 ), MSD2, WORK( IQ ), IWORK( IQPTR ), &
                            IWORK( IPRMPT ), IWORK( IPERM ), IWORK( IGIVPT ), IWORK( IGIVCL ), &
                            WORK( IGIVNM ), WORK( IWREM ), IWORK( SUBPBS+1 ), INFO )
            END IF
            IF( INFO.NE.0 ) GO TO 130
            IWORK( I / 2+1 ) = IWORK( I+2 )
   90    CONTINUE
         SUBPBS = SUBPBS / 2
         CURLVL = CURLVL + 1
         GO TO 80
      END IF
!
!     end while
!
!     Re-merge the eigenvalues/vectors which were deflated at the final
!     merge step.
!
      IF( ICOMPQ.EQ.1 ) THEN
         DO 100 I = 1, N
            J = IWORK( INDXQ+I )
            WORK( I ) = D( J )
            CALL DCOPY( QSIZ, QSTORE( 1, J ), 1, Q( 1, I ), 1 )
  100    CONTINUE
         CALL DCOPY( N, WORK, 1, D, 1 )
      ELSE IF( ICOMPQ.EQ.2 ) THEN
         DO 110 I = 1, N
            J = IWORK( INDXQ+I )
            WORK( I ) = D( J )
            CALL DCOPY( N, Q( 1, J ), 1, WORK( N*I+1 ), 1 )
  110    CONTINUE
         CALL DCOPY( N, WORK, 1, D, 1 )
         CALL DLACPY( 'A', N, N, WORK( N+1 ), N, Q, LDQ )
      ELSE
         DO 120 I = 1, N
            J = IWORK( INDXQ+I )
            WORK( I ) = D( J )
  120    CONTINUE
         CALL DCOPY( N, WORK, 1, D, 1 )
      END IF
      GO TO 140
!
  130 CONTINUE
      INFO = SUBMAT*( N+1 ) + SUBMAT + MATSIZ - 1
!
  140 CONTINUE
      RETURN
!
!     End of DLAED0
!
      END
!> \brief \b DLACPY
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DLACPY + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlacpy.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlacpy.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlacpy.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLACPY( UPLO, M, N, A, LDA, B, LDB )
! 
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       integer(i4)  ::            LDA, LDB, M, N
!       ..
!       .. Array Arguments ..
!       real(dp)    ::   A( LDA, * ), B( LDB, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLACPY copies all or part of a two-dimensional matrix A to another
!> matrix B.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies the part of the matrix A to be copied to B.
!>          = 'U':      Upper triangular part
!>          = 'L':      Lower triangular part
!>          Otherwise:  All of the matrix A
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is integer(i4)  ::
!>          The number of rows of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is integer(i4)  ::
!>          The number of columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is real(dp)    :: array, dimension (LDA,N)
!>          The m by n matrix A.  If UPLO = 'U', only the upper triangle
!>          or trapezoid is accessed; if UPLO = 'L', only the lower
!>          triangle or trapezoid is accessed.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is integer(i4)  ::
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is real(dp)    :: array, dimension (LDB,N)
!>          On exit, B = A in the locations specified by UPLO.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is integer(i4)  ::
!>          The leading dimension of the array B.  LDB >= max(1,M).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup auxOTHERauxiliary
!
!  =====================================================================
      SUBROUTINE DLACPY( UPLO, M, N, A, LDA, B, LDB )
!
!  -- LAPACK auxiliary routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes, only : i4, dp, dpc
      CHARACTER          UPLO
      integer(i4)  ::            LDA, LDB, M, N
!     ..
!     .. Array Arguments ..
      real(dp)    ::   A( LDA, * ), B( LDB, * )
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      integer(i4)  ::            I, J
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MIN
!     ..
!     .. Executable Statements ..
!
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 20 J = 1, N
            DO 10 I = 1, MIN( J, M )
               B( I, J ) = A( I, J )
   10       CONTINUE
   20    CONTINUE
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
         DO 40 J = 1, N
            DO 30 I = J, M
               B( I, J ) = A( I, J )
   30       CONTINUE
   40    CONTINUE
      ELSE
         DO 60 J = 1, N
            DO 50 I = 1, M
               B( I, J ) = A( I, J )
   50       CONTINUE
   60    CONTINUE
      END IF
      RETURN
!
!     End of DLACPY
!
      END
!> \brief \b ZUNM2L
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download ZUNM2L + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zunm2l.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zunm2l.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zunm2l.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZUNM2L( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
!                          WORK, INFO )
! 
!       .. Scalar Arguments ..
!       CHARACTER          SIDE, TRANS
!       integer(i4)  ::            INFO, K, LDA, LDC, M, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZUNM2L overwrites the general complex m-by-n matrix C with
!>
!>       Q * C  if SIDE = 'L' and TRANS = 'N', or
!>
!>       Q**H* C  if SIDE = 'L' and TRANS = 'C', or
!>
!>       C * Q  if SIDE = 'R' and TRANS = 'N', or
!>
!>       C * Q**H if SIDE = 'R' and TRANS = 'C',
!>
!> where Q is a complex unitary matrix defined as the product of k
!> elementary reflectors
!>
!>       Q = H(k) . . . H(2) H(1)
!>
!> as returned by ZGEQLF. Q is of order m if SIDE = 'L' and of order n
!> if SIDE = 'R'.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'L': apply Q or Q**H from the Left
!>          = 'R': apply Q or Q**H from the Right
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          = 'N': apply Q  (No transpose)
!>          = 'C': apply Q**H (Conjugate transpose)
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is integer(i4)  ::
!>          The number of rows of the matrix C. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is integer(i4)  ::
!>          The number of columns of the matrix C. N >= 0.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is integer(i4)  ::
!>          The number of elementary reflectors whose product defines
!>          the matrix Q.
!>          If SIDE = 'L', M >= K >= 0;
!>          if SIDE = 'R', N >= K >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,K)
!>          The i-th column must contain the vector which defines the
!>          elementary reflector H(i), for i = 1,2,...,k, as returned by
!>          ZGEQLF in the last k columns of its array argument A.
!>          A is modified by the routine but restored on exit.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is integer(i4)  ::
!>          The leading dimension of the array A.
!>          If SIDE = 'L', LDA >= max(1,M);
!>          if SIDE = 'R', LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is COMPLEX*16 array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by ZGEQLF.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is COMPLEX*16 array, dimension (LDC,N)
!>          On entry, the m-by-n matrix C.
!>          On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is integer(i4)  ::
!>          The leading dimension of the array C. LDC >= max(1,M).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension
!>                                   (N) if SIDE = 'L',
!>                                   (M) if SIDE = 'R'
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is integer(i4)  ::
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup complex16OTHERcomputational
!
!  =====================================================================
      SUBROUTINE ZUNM2L( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes, only : i4, dp, dpc
      CHARACTER          SIDE, TRANS
      integer(i4)  ::            INFO, K, LDA, LDC, M, N
!     ..
!     .. Array Arguments ..
      COMPLEX(dpc)       A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX(dpc)       ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            LEFT, NOTRAN
      integer(i4)  ::            I, I1, I2, I3, MI, NI, NQ
      COMPLEX(dpc)       AII, TAUI
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARF
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
!
!     NQ is the order of Q
!
      IF( LEFT ) THEN
         NQ = M
      ELSE
         NQ = N
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNM2L', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) RETURN
!
      IF( ( LEFT .AND. NOTRAN .OR. .NOT.LEFT .AND. .NOT.NOTRAN ) ) THEN
         I1 = 1
         I2 = K
         I3 = 1
      ELSE
         I1 = K
         I2 = 1
         I3 = -1
      END IF
!
      IF( LEFT ) THEN
         NI = N
      ELSE
         MI = M
      END IF
!
      DO 10 I = I1, I2, I3
         IF( LEFT ) THEN
!
!           H(i) or H(i)**H is applied to C(1:m-k+i,1:n)
!
            MI = M - K + I
         ELSE
!
!           H(i) or H(i)**H is applied to C(1:m,1:n-k+i)
!
            NI = N - K + I
         END IF
!
!        Apply H(i) or H(i)**H
!
         IF( NOTRAN ) THEN
            TAUI = TAU( I )
         ELSE
            TAUI = DCONJG( TAU( I ) )
         END IF
         AII = A( NQ-K+I, I )
         A( NQ-K+I, I ) = ONE
         CALL ZLARF( SIDE, MI, NI, A( 1, I ), 1, TAUI, C, LDC, WORK )
         A( NQ-K+I, I ) = AII
   10 CONTINUE
      RETURN
!
!     End of ZUNM2L
!
      END
!> \brief \b ZUNM2R
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download ZUNM2R + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zunm2r.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zunm2r.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zunm2r.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZUNM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
!                          WORK, INFO )
! 
!       .. Scalar Arguments ..
!       CHARACTER          SIDE, TRANS
!       integer(i4)  ::            INFO, K, LDA, LDC, M, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZUNM2R overwrites the general complex m-by-n matrix C with
!>
!>       Q * C  if SIDE = 'L' and TRANS = 'N', or
!>
!>       Q**H* C  if SIDE = 'L' and TRANS = 'C', or
!>
!>       C * Q  if SIDE = 'R' and TRANS = 'N', or
!>
!>       C * Q**H if SIDE = 'R' and TRANS = 'C',
!>
!> where Q is a complex unitary matrix defined as the product of k
!> elementary reflectors
!>
!>       Q = H(1) H(2) . . . H(k)
!>
!> as returned by ZGEQRF. Q is of order m if SIDE = 'L' and of order n
!> if SIDE = 'R'.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'L': apply Q or Q**H from the Left
!>          = 'R': apply Q or Q**H from the Right
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          = 'N': apply Q  (No transpose)
!>          = 'C': apply Q**H (Conjugate transpose)
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is integer(i4)  ::
!>          The number of rows of the matrix C. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is integer(i4)  ::
!>          The number of columns of the matrix C. N >= 0.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is integer(i4)  ::
!>          The number of elementary reflectors whose product defines
!>          the matrix Q.
!>          If SIDE = 'L', M >= K >= 0;
!>          if SIDE = 'R', N >= K >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,K)
!>          The i-th column must contain the vector which defines the
!>          elementary reflector H(i), for i = 1,2,...,k, as returned by
!>          ZGEQRF in the first k columns of its array argument A.
!>          A is modified by the routine but restored on exit.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is integer(i4)  ::
!>          The leading dimension of the array A.
!>          If SIDE = 'L', LDA >= max(1,M);
!>          if SIDE = 'R', LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is COMPLEX*16 array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by ZGEQRF.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is COMPLEX*16 array, dimension (LDC,N)
!>          On entry, the m-by-n matrix C.
!>          On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is integer(i4)  ::
!>          The leading dimension of the array C. LDC >= max(1,M).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension
!>                                   (N) if SIDE = 'L',
!>                                   (M) if SIDE = 'R'
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is integer(i4)  ::
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup complex16OTHERcomputational
!
!  =====================================================================
      SUBROUTINE ZUNM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes, only : i4, dp, dpc
      CHARACTER          SIDE, TRANS
      integer(i4)  ::            INFO, K, LDA, LDC, M, N
!     ..
!     .. Array Arguments ..
      COMPLEX(dpc)       A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX(dpc)       ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            LEFT, NOTRAN
      integer(i4)  ::            I, I1, I2, I3, IC, JC, MI, NI, NQ
      COMPLEX(dpc)       AII, TAUI
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARF
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
!
!     NQ is the order of Q
!
      IF( LEFT ) THEN
         NQ = M
      ELSE
         NQ = N
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNM2R', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) RETURN
!
      IF( ( LEFT .AND. .NOT.NOTRAN .OR. .NOT.LEFT .AND. NOTRAN ) ) THEN
         I1 = 1
         I2 = K
         I3 = 1
      ELSE
         I1 = K
         I2 = 1
         I3 = -1
      END IF
!
      IF( LEFT ) THEN
         NI = N
         JC = 1
      ELSE
         MI = M
         IC = 1
      END IF
!
      DO 10 I = I1, I2, I3
         IF( LEFT ) THEN
!
!           H(i) or H(i)**H is applied to C(i:m,1:n)
!
            MI = M - I + 1
            IC = I
         ELSE
!
!           H(i) or H(i)**H is applied to C(1:m,i:n)
!
            NI = N - I + 1
            JC = I
         END IF
!
!        Apply H(i) or H(i)**H
!
         IF( NOTRAN ) THEN
            TAUI = TAU( I )
         ELSE
            TAUI = DCONJG( TAU( I ) )
         END IF
         AII = A( I, I )
         A( I, I ) = ONE
         CALL ZLARF( SIDE, MI, NI, A( I, I ), 1, TAUI, C( IC, JC ), LDC, WORK )
         A( I, I ) = AII
   10 CONTINUE
      RETURN
!
!     End of ZUNM2R
!
      END
!> \brief \b ZLAED8
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download ZLAED8 + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaed8.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaed8.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaed8.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLAED8( K, N, QSIZ, Q, LDQ, D, RHO, CUTPNT, Z, DLAMDA,
!                          Q2, LDQ2, W, INDXP, INDX, INDXQ, PERM, GIVPTR,
!                          GIVCOL, GIVNUM, INFO )
! 
!       .. Scalar Arguments ..
!       integer(i4)  ::            CUTPNT, GIVPTR, INFO, K, LDQ, LDQ2, N, QSIZ
!       real(dp)    ::   RHO
!       ..
!       .. Array Arguments ..
!       integer(i4)  ::            GIVCOL( 2, * ), INDX( * ), INDXP( * ),
!      $                   INDXQ( * ), PERM( * )
!       real(dp)    ::   D( * ), DLAMDA( * ), GIVNUM( 2, * ), W( * ),
!      $                   Z( * )
!       COMPLEX*16         Q( LDQ, * ), Q2( LDQ2, * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLAED8 merges the two sets of eigenvalues together into a single
!> sorted set.  Then it tries to deflate the size of the problem.
!> There are two ways in which deflation can occur:  when two or more
!> eigenvalues are close together or if there is a tiny element in the
!> Z vector.  For each such occurrence the order of the related secular
!> equation problem is reduced by one.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[out] K
!> \verbatim
!>          K is integer(i4)  ::
!>         Contains the number of non-deflated eigenvalues.
!>         This is the order of the related secular equation.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is integer(i4)  ::
!>         The dimension of the symmetric tridiagonal matrix.  N >= 0.
!> \endverbatim
!>
!> \param[in] QSIZ
!> \verbatim
!>          QSIZ is integer(i4)  ::
!>         The dimension of the unitary matrix used to reduce
!>         the dense or band matrix to tridiagonal form.
!>         QSIZ >= N if ICOMPQ = 1.
!> \endverbatim
!>
!> \param[in,out] Q
!> \verbatim
!>          Q is COMPLEX*16 array, dimension (LDQ,N)
!>         On entry, Q contains the eigenvectors of the partially solved
!>         system which has been previously updated in matrix
!>         multiplies with other partially solved eigensystems.
!>         On exit, Q contains the trailing (N-K) updated eigenvectors
!>         (those which were deflated) in its last N-K columns.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is integer(i4)  ::
!>         The leading dimension of the array Q.  LDQ >= max( 1, N ).
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is real(dp)    :: array, dimension (N)
!>         On entry, D contains the eigenvalues of the two submatrices to
!>         be combined.  On exit, D contains the trailing (N-K) updated
!>         eigenvalues (those which were deflated) sorted into increasing
!>         order.
!> \endverbatim
!>
!> \param[in,out] RHO
!> \verbatim
!>          RHO is real(dp)    ::
!>         Contains the off diagonal element associated with the rank-1
!>         cut which originally split the two submatrices which are now
!>         being recombined. RHO is modified during the computation to
!>         the value required by DLAED3.
!> \endverbatim
!>
!> \param[in] CUTPNT
!> \verbatim
!>          CUTPNT is integer(i4)  ::
!>         Contains the location of the last eigenvalue in the leading
!>         sub-matrix.  MIN(1,N) <= CUTPNT <= N.
!> \endverbatim
!>
!> \param[in] Z
!> \verbatim
!>          Z is real(dp)    :: array, dimension (N)
!>         On input this vector contains the updating vector (the last
!>         row of the first sub-eigenvector matrix and the first row of
!>         the second sub-eigenvector matrix).  The contents of Z are
!>         destroyed during the updating process.
!> \endverbatim
!>
!> \param[out] DLAMDA
!> \verbatim
!>          DLAMDA is real(dp)    :: array, dimension (N)
!>         Contains a copy of the first K eigenvalues which will be used
!>         by DLAED3 to form the secular equation.
!> \endverbatim
!>
!> \param[out] Q2
!> \verbatim
!>          Q2 is COMPLEX*16 array, dimension (LDQ2,N)
!>         If ICOMPQ = 0, Q2 is not referenced.  Otherwise,
!>         Contains a copy of the first K eigenvectors which will be used
!>         by DLAED7 in a matrix multiply (DGEMM) to update the new
!>         eigenvectors.
!> \endverbatim
!>
!> \param[in] LDQ2
!> \verbatim
!>          LDQ2 is integer(i4)  ::
!>         The leading dimension of the array Q2.  LDQ2 >= max( 1, N ).
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is real(dp)    :: array, dimension (N)
!>         This will hold the first k values of the final
!>         deflation-altered z-vector and will be passed to DLAED3.
!> \endverbatim
!>
!> \param[out] INDXP
!> \verbatim
!>          INDXP is integer(i4)  :: array, dimension (N)
!>         This will contain the permutation used to place deflated
!>         values of D at the end of the array. On output INDXP(1:K)
!>         points to the nondeflated D-values and INDXP(K+1:N)
!>         points to the deflated eigenvalues.
!> \endverbatim
!>
!> \param[out] INDX
!> \verbatim
!>          INDX is integer(i4)  :: array, dimension (N)
!>         This will contain the permutation used to sort the contents of
!>         D into ascending order.
!> \endverbatim
!>
!> \param[in] INDXQ
!> \verbatim
!>          INDXQ is integer(i4)  :: array, dimension (N)
!>         This contains the permutation which separately sorts the two
!>         sub-problems in D into ascending order.  Note that elements in
!>         the second half of this permutation must first have CUTPNT
!>         added to their values in order to be accurate.
!> \endverbatim
!>
!> \param[out] PERM
!> \verbatim
!>          PERM is integer(i4)  :: array, dimension (N)
!>         Contains the permutations (from deflation and sorting) to be
!>         applied to each eigenblock.
!> \endverbatim
!>
!> \param[out] GIVPTR
!> \verbatim
!>          GIVPTR is integer(i4)  ::
!>         Contains the number of Givens rotations which took place in
!>         this subproblem.
!> \endverbatim
!>
!> \param[out] GIVCOL
!> \verbatim
!>          GIVCOL is integer(i4)  :: array, dimension (2, N)
!>         Each pair of numbers indicates a pair of columns to take place
!>         in a Givens rotation.
!> \endverbatim
!>
!> \param[out] GIVNUM
!> \verbatim
!>          GIVNUM is real(dp)    :: array, dimension (2, N)
!>         Each number indicates the S value to be used in the
!>         corresponding Givens rotation.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is integer(i4)  ::
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup complex16OTHERcomputational
!
!  =====================================================================
      SUBROUTINE ZLAED8( K, N, QSIZ, Q, LDQ, D, RHO, CUTPNT, Z, DLAMDA, &
                         Q2, LDQ2, W, INDXP, INDX, INDXQ, PERM, GIVPTR, &
                         GIVCOL, GIVNUM, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes, only : i4, dp, dpc
      integer(i4) :: CUTPNT, GIVPTR, INFO, K, LDQ, LDQ2, N, QSIZ
      real(dp)    ::   RHO
!     ..
!     .. Array Arguments ..
      integer(i4)  ::            GIVCOL( 2, * ), INDX( * ), INDXP( * ), INDXQ( * ), PERM( * )
      real(dp)    ::   D( * ), DLAMDA( * ), GIVNUM( 2, * ), W( * ), Z( * )
      COMPLEX(dpc)       Q( LDQ, * ), Q2( LDQ2, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      real(dp)    ::   MONE, ZERO, ONE, TWO, EIGHT
      PARAMETER          ( MONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, EIGHT = 8.0D0 )
!     ..
!     .. Local Scalars ..
      integer(i4)  ::            I, IMAX, J, JLAM, JMAX, JP, K2, N1, N1P1, N2
      real(dp)    ::   C, EPS, S, T, TAU, TOL
!     ..
!     .. External Functions ..
      integer(i4)  ::            IDAMAX
      real(dp)    ::   DLAMCH, DLAPY2
      EXTERNAL           IDAMAX, DLAMCH, DLAPY2
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPY, DLAMRG, DSCAL, XERBLA, ZCOPY, ZDROT, ZLACPY
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
!
      IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( QSIZ.LT.N ) THEN
         INFO = -3
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( CUTPNT.LT.MIN( 1, N ) .OR. CUTPNT.GT.N ) THEN
         INFO = -8
      ELSE IF( LDQ2.LT.MAX( 1, N ) ) THEN
         INFO = -12
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZLAED8', -INFO )
         RETURN
      END IF
!
!     Need to initialize GIVPTR to O here in case of quick exit
!     to prevent an unspecified code behavior (usually sigfault) 
!     when IWORK array on entry to *stedc is not zeroed 
!     (or at least some IWORK entries which used in *laed7 for GIVPTR).
!
      GIVPTR = 0
!
!     Quick return if possible
!
      IF( N.EQ.0 ) RETURN
!
      N1 = CUTPNT
      N2 = N - N1
      N1P1 = N1 + 1
!
      IF( RHO.LT.ZERO ) THEN
         CALL DSCAL( N2, MONE, Z( N1P1 ), 1 )
      END IF
!
!     Normalize z so that norm(z) = 1
!
      T = ONE / SQRT( TWO )
      DO 10 J = 1, N
         INDX( J ) = J
   10 CONTINUE
      CALL DSCAL( N, T, Z, 1 )
      RHO = ABS( TWO*RHO )
!
!     Sort the eigenvalues into increasing order
!
      DO 20 I = CUTPNT + 1, N
         INDXQ( I ) = INDXQ( I ) + CUTPNT
   20 CONTINUE
      DO 30 I = 1, N
         DLAMDA( I ) = D( INDXQ( I ) )
         W( I ) = Z( INDXQ( I ) )
   30 CONTINUE
      I = 1
      J = CUTPNT + 1
      CALL DLAMRG( N1, N2, DLAMDA, 1, 1, INDX )
      DO 40 I = 1, N
         D( I ) = DLAMDA( INDX( I ) )
         Z( I ) = W( INDX( I ) )
   40 CONTINUE
!
!     Calculate the allowable deflation tolerance
!
      IMAX = IDAMAX( N, Z, 1 )
      JMAX = IDAMAX( N, D, 1 )
      EPS = DLAMCH( 'Epsilon' )
      TOL = EIGHT*EPS*ABS( D( JMAX ) )
!
!     If the rank-1 modifier is small enough, no more needs to be done
!     -- except to reorganize Q so that its columns correspond with the
!     elements in D.
!
      IF( RHO*ABS( Z( IMAX ) ).LE.TOL ) THEN
         K = 0
         DO 50 J = 1, N
            PERM( J ) = INDXQ( INDX( J ) )
            CALL ZCOPY( QSIZ, Q( 1, PERM( J ) ), 1, Q2( 1, J ), 1 )
   50    CONTINUE
         CALL ZLACPY( 'A', QSIZ, N, Q2( 1, 1 ), LDQ2, Q( 1, 1 ), LDQ )
         RETURN
      END IF
!
!     If there are multiple eigenvalues then the problem deflates.  Here
!     the number of equal eigenvalues are found.  As each equal
!     eigenvalue is found, an elementary reflector is computed to rotate
!     the corresponding eigensubspace so that the corresponding
!     components of Z are zero in this new basis.
!
      K = 0
      K2 = N + 1
      DO 60 J = 1, N
         IF( RHO*ABS( Z( J ) ).LE.TOL ) THEN
!
!           Deflate due to small z component.
!
            K2 = K2 - 1
            INDXP( K2 ) = J
            IF( J.EQ.N ) GO TO 100
         ELSE
            JLAM = J
            GO TO 70
         END IF
   60 CONTINUE
   70 CONTINUE
      J = J + 1
      IF( J.GT.N ) GO TO 90
      IF( RHO*ABS( Z( J ) ).LE.TOL ) THEN
!
!        Deflate due to small z component.
!
         K2 = K2 - 1
         INDXP( K2 ) = J
      ELSE
!
!        Check if eigenvalues are close enough to allow deflation.
!
         S = Z( JLAM )
         C = Z( J )
!
!        Find sqrt(a**2+b**2) without overflow or
!        destructive underflow.
!
         TAU = DLAPY2( C, S )
         T = D( J ) - D( JLAM )
         C = C / TAU
         S = -S / TAU
         IF( ABS( T*C*S ).LE.TOL ) THEN
!
!           Deflation is possible.
!
            Z( J ) = TAU
            Z( JLAM ) = ZERO
!
!           Record the appropriate Givens rotation
!
            GIVPTR = GIVPTR + 1
            GIVCOL( 1, GIVPTR ) = INDXQ( INDX( JLAM ) )
            GIVCOL( 2, GIVPTR ) = INDXQ( INDX( J ) )
            GIVNUM( 1, GIVPTR ) = C
            GIVNUM( 2, GIVPTR ) = S
            CALL ZDROT( QSIZ, Q( 1, INDXQ( INDX( JLAM ) ) ), 1, Q( 1, INDXQ( INDX( J ) ) ), 1, C, S )
            T = D( JLAM )*C*C + D( J )*S*S
            D( J ) = D( JLAM )*S*S + D( J )*C*C
            D( JLAM ) = T
            K2 = K2 - 1
            I = 1
   80       CONTINUE
            IF( K2+I.LE.N ) THEN
               IF( D( JLAM ).LT.D( INDXP( K2+I ) ) ) THEN
                  INDXP( K2+I-1 ) = INDXP( K2+I )
                  INDXP( K2+I ) = JLAM
                  I = I + 1
                  GO TO 80
               ELSE
                  INDXP( K2+I-1 ) = JLAM
               END IF
            ELSE
               INDXP( K2+I-1 ) = JLAM
            END IF
            JLAM = J
         ELSE
            K = K + 1
            W( K ) = Z( JLAM )
            DLAMDA( K ) = D( JLAM )
            INDXP( K ) = JLAM
            JLAM = J
         END IF
      END IF
      GO TO 70
   90 CONTINUE
!
!     Record the last eigenvalue.
!
      K = K + 1
      W( K ) = Z( JLAM )
      DLAMDA( K ) = D( JLAM )
      INDXP( K ) = JLAM
!
  100 CONTINUE
!
!     Sort the eigenvalues and corresponding eigenvectors into DLAMDA
!     and Q2 respectively.  The eigenvalues/vectors which were not
!     deflated go into the first K slots of DLAMDA and Q2 respectively,
!     while those which were deflated go into the last N - K slots.
!
      DO 110 J = 1, N
         JP = INDXP( J )
         DLAMDA( J ) = D( JP )
         PERM( J ) = INDXQ( INDX( JP ) )
         CALL ZCOPY( QSIZ, Q( 1, PERM( J ) ), 1, Q2( 1, J ), 1 )
  110 CONTINUE
!
!     The deflated eigenvalues and their corresponding vectors go back
!     into the last N - K slots of D and Q respectively.
!
      IF( K.LT.N ) THEN
         CALL DCOPY( N-K, DLAMDA( K+1 ), 1, D( K+1 ), 1 )
         CALL ZLACPY( 'A', QSIZ, N-K, Q2( 1, K+1 ), LDQ2, Q( 1, K+1 ), LDQ )
      END IF
!
      RETURN
!
!     End of ZLAED8
!
      END
!> \brief \b DLAED1
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DLAED1 + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaed1.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaed1.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaed1.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAED1( N, D, Q, LDQ, INDXQ, RHO, CUTPNT, WORK, IWORK,
!                          INFO )
! 
!       .. Scalar Arguments ..
!       integer(i4)  ::            CUTPNT, INFO, LDQ, N
!       real(dp)    ::   RHO
!       ..
!       .. Array Arguments ..
!       integer(i4)  ::            INDXQ( * ), IWORK( * )
!       real(dp)    ::   D( * ), Q( LDQ, * ), WORK( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAED1 computes the updated eigensystem of a diagonal
!> matrix after modification by a rank-one symmetric matrix.  This
!> routine is used only for the eigenproblem which requires all
!> eigenvalues and eigenvectors of a tridiagonal matrix.  DLAED7 handles
!> the case in which eigenvalues only or eigenvalues and eigenvectors
!> of a full symmetric matrix (which was reduced to tridiagonal form)
!> are desired.
!>
!>   T = Q(in) ( D(in) + RHO * Z*Z**T ) Q**T(in) = Q(out) * D(out) * Q**T(out)
!>
!>    where Z = Q**T*u, u is a vector of length N with ones in the
!>    CUTPNT and CUTPNT + 1 th elements and zeros elsewhere.
!>
!>    The eigenvectors of the original matrix are stored in Q, and the
!>    eigenvalues are in D.  The algorithm consists of three stages:
!>
!>       The first stage consists of deflating the size of the problem
!>       when there are multiple eigenvalues or if there is a zero in
!>       the Z vector.  For each such occurence the dimension of the
!>       secular equation problem is reduced by one.  This stage is
!>       performed by the routine DLAED2.
!>
!>       The second stage consists of calculating the updated
!>       eigenvalues. This is done by finding the roots of the secular
!>       equation via the routine DLAED4 (as called by DLAED3).
!>       This routine also calculates the eigenvectors of the current
!>       problem.
!>
!>       The final stage consists of computing the updated eigenvectors
!>       directly using the updated eigenvalues.  The eigenvectors for
!>       the current problem are multiplied with the eigenvectors from
!>       the overall problem.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is integer(i4)  ::
!>         The dimension of the symmetric tridiagonal matrix.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is real(dp)    :: array, dimension (N)
!>         On entry, the eigenvalues of the rank-1-perturbed matrix.
!>         On exit, the eigenvalues of the repaired matrix.
!> \endverbatim
!>
!> \param[in,out] Q
!> \verbatim
!>          Q is real(dp)    :: array, dimension (LDQ,N)
!>         On entry, the eigenvectors of the rank-1-perturbed matrix.
!>         On exit, the eigenvectors of the repaired tridiagonal matrix.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is integer(i4)  ::
!>         The leading dimension of the array Q.  LDQ >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] INDXQ
!> \verbatim
!>          INDXQ is integer(i4)  :: array, dimension (N)
!>         On entry, the permutation which separately sorts the two
!>         subproblems in D into ascending order.
!>         On exit, the permutation which will reintegrate the
!>         subproblems back into sorted order,
!>         i.e. D( INDXQ( I = 1, N ) ) will be in ascending order.
!> \endverbatim
!>
!> \param[in] RHO
!> \verbatim
!>          RHO is real(dp)    ::
!>         The subdiagonal entry used to create the rank-1 modification.
!> \endverbatim
!>
!> \param[in] CUTPNT
!> \verbatim
!>          CUTPNT is integer(i4)  ::
!>         The location of the last eigenvalue in the leading sub-matrix.
!>         min(1,N) <= CUTPNT <= N/2.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is real(dp)    :: array, dimension (4*N + N**2)
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is integer(i4)  :: array, dimension (4*N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is integer(i4)  ::
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  if INFO = 1, an eigenvalue did not converge
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup auxOTHERcomputational
!
!> \par Contributors:
!  ==================
!>
!> Jeff Rutter, Computer Science Division, University of California
!> at Berkeley, USA \n
!>  Modified by Francoise Tisseur, University of Tennessee
!>
!  =====================================================================
      SUBROUTINE DLAED1( N, D, Q, LDQ, INDXQ, RHO, CUTPNT, WORK, IWORK, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes, only : i4, dp, dpc
      integer(i4)  ::            CUTPNT, INFO, LDQ, N
      real(dp)    ::   RHO
!     ..
!     .. Array Arguments ..
      integer(i4)  ::            INDXQ( * ), IWORK( * )
      real(dp)    ::   D( * ), Q( LDQ, * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      integer(i4)  ::            COLTYP, I, IDLMDA, INDX, INDXC, INDXP, IQ2, IS, IW, IZ, K, N1, N2, ZPP1
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPY, DLAED2, DLAED3, DLAMRG, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
!
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( MIN( 1, N / 2 ).GT.CUTPNT .OR. ( N / 2 ).LT.CUTPNT ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLAED1', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) RETURN
!
!     The following values are integer pointers which indicate
!     the portion of the workspace
!     used by a particular array in DLAED2 and DLAED3.
!
      IZ = 1
      IDLMDA = IZ + N
      IW = IDLMDA + N
      IQ2 = IW + N
!
      INDX = 1
      INDXC = INDX + N
      COLTYP = INDXC + N
      INDXP = COLTYP + N
!
!
!     Form the z-vector which consists of the last row of Q_1 and the
!     first row of Q_2.
!
      CALL DCOPY( CUTPNT, Q( CUTPNT, 1 ), LDQ, WORK( IZ ), 1 )
      ZPP1 = CUTPNT + 1
      CALL DCOPY( N-CUTPNT, Q( ZPP1, ZPP1 ), LDQ, WORK( IZ+CUTPNT ), 1 )
!
!     Deflate eigenvalues.
!
      CALL DLAED2( K, N, CUTPNT, D, Q, LDQ, INDXQ, RHO, WORK( IZ ), &
                   WORK( IDLMDA ), WORK( IW ), WORK( IQ2 ), &
                   IWORK( INDX ), IWORK( INDXC ), IWORK( INDXP ), IWORK( COLTYP ), INFO )
!
      IF( INFO.NE.0 ) GO TO 20
!
!     Solve Secular Equation.
!
      IF( K.NE.0 ) THEN
         IS = ( IWORK( COLTYP )+IWORK( COLTYP+1 ) )*CUTPNT + ( IWORK( COLTYP+1 )+IWORK( COLTYP+2 ) )*( N-CUTPNT ) + IQ2
         CALL DLAED3( K, N, CUTPNT, D, Q, LDQ, RHO, WORK( IDLMDA ), WORK( IQ2 ), IWORK( INDXC ), IWORK( COLTYP ), &
                      WORK( IW ), WORK( IS ), INFO )
         IF( INFO.NE.0 ) GO TO 20
!
!     Prepare the INDXQ sorting permutation.
!
         N1 = K
         N2 = N - K
         CALL DLAMRG( N1, N2, D, 1, -1, INDXQ )
      ELSE
         DO 10 I = 1, N
            INDXQ( I ) = I
   10    CONTINUE
      END IF
!
   20 CONTINUE
      RETURN
!
!     End of DLAED1
!
      END
!> \brief \b DLAED7
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DLAED7 + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaed7.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaed7.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaed7.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAED7( ICOMPQ, N, QSIZ, TLVLS, CURLVL, CURPBM, D, Q,
!                          LDQ, INDXQ, RHO, CUTPNT, QSTORE, QPTR, PRMPTR,
!                          PERM, GIVPTR, GIVCOL, GIVNUM, WORK, IWORK,
!                          INFO )
! 
!       .. Scalar Arguments ..
!       integer(i4)  ::            CURLVL, CURPBM, CUTPNT, ICOMPQ, INFO, LDQ, N,
!      $                   QSIZ, TLVLS
!       real(dp)    ::   RHO
!       ..
!       .. Array Arguments ..
!       integer(i4)  ::            GIVCOL( 2, * ), GIVPTR( * ), INDXQ( * ),
!      $                   IWORK( * ), PERM( * ), PRMPTR( * ), QPTR( * )
!       real(dp)    ::   D( * ), GIVNUM( 2, * ), Q( LDQ, * ),
!      $                   QSTORE( * ), WORK( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAED7 computes the updated eigensystem of a diagonal
!> matrix after modification by a rank-one symmetric matrix. This
!> routine is used only for the eigenproblem which requires all
!> eigenvalues and optionally eigenvectors of a dense symmetric matrix
!> that has been reduced to tridiagonal form.  DLAED1 handles
!> the case in which all eigenvalues and eigenvectors of a symmetric
!> tridiagonal matrix are desired.
!>
!>   T = Q(in) ( D(in) + RHO * Z*Z**T ) Q**T(in) = Q(out) * D(out) * Q**T(out)
!>
!>    where Z = Q**Tu, u is a vector of length N with ones in the
!>    CUTPNT and CUTPNT + 1 th elements and zeros elsewhere.
!>
!>    The eigenvectors of the original matrix are stored in Q, and the
!>    eigenvalues are in D.  The algorithm consists of three stages:
!>
!>       The first stage consists of deflating the size of the problem
!>       when there are multiple eigenvalues or if there is a zero in
!>       the Z vector.  For each such occurence the dimension of the
!>       secular equation problem is reduced by one.  This stage is
!>       performed by the routine DLAED8.
!>
!>       The second stage consists of calculating the updated
!>       eigenvalues. This is done by finding the roots of the secular
!>       equation via the routine DLAED4 (as called by DLAED9).
!>       This routine also calculates the eigenvectors of the current
!>       problem.
!>
!>       The final stage consists of computing the updated eigenvectors
!>       directly using the updated eigenvalues.  The eigenvectors for
!>       the current problem are multiplied with the eigenvectors from
!>       the overall problem.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ICOMPQ
!> \verbatim
!>          ICOMPQ is integer(i4)  ::
!>          = 0:  Compute eigenvalues only.
!>          = 1:  Compute eigenvectors of original dense symmetric matrix
!>                also.  On entry, Q contains the orthogonal matrix used
!>                to reduce the original matrix to tridiagonal form.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is integer(i4)  ::
!>         The dimension of the symmetric tridiagonal matrix.  N >= 0.
!> \endverbatim
!>
!> \param[in] QSIZ
!> \verbatim
!>          QSIZ is integer(i4)  ::
!>         The dimension of the orthogonal matrix used to reduce
!>         the full matrix to tridiagonal form.  QSIZ >= N if ICOMPQ = 1.
!> \endverbatim
!>
!> \param[in] TLVLS
!> \verbatim
!>          TLVLS is integer(i4)  ::
!>         The total number of merging levels in the overall divide and
!>         conquer tree.
!> \endverbatim
!>
!> \param[in] CURLVL
!> \verbatim
!>          CURLVL is integer(i4)  ::
!>         The current level in the overall merge routine,
!>         0 <= CURLVL <= TLVLS.
!> \endverbatim
!>
!> \param[in] CURPBM
!> \verbatim
!>          CURPBM is integer(i4)  ::
!>         The current problem in the current level in the overall
!>         merge routine (counting from upper left to lower right).
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is real(dp)    :: array, dimension (N)
!>         On entry, the eigenvalues of the rank-1-perturbed matrix.
!>         On exit, the eigenvalues of the repaired matrix.
!> \endverbatim
!>
!> \param[in,out] Q
!> \verbatim
!>          Q is real(dp)    :: array, dimension (LDQ, N)
!>         On entry, the eigenvectors of the rank-1-perturbed matrix.
!>         On exit, the eigenvectors of the repaired tridiagonal matrix.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is integer(i4)  ::
!>         The leading dimension of the array Q.  LDQ >= max(1,N).
!> \endverbatim
!>
!> \param[out] INDXQ
!> \verbatim
!>          INDXQ is integer(i4)  :: array, dimension (N)
!>         The permutation which will reintegrate the subproblem just
!>         solved back into sorted order, i.e., D( INDXQ( I = 1, N ) )
!>         will be in ascending order.
!> \endverbatim
!>
!> \param[in] RHO
!> \verbatim
!>          RHO is real(dp)    ::
!>         The subdiagonal element used to create the rank-1
!>         modification.
!> \endverbatim
!>
!> \param[in] CUTPNT
!> \verbatim
!>          CUTPNT is integer(i4)  ::
!>         Contains the location of the last eigenvalue in the leading
!>         sub-matrix.  min(1,N) <= CUTPNT <= N.
!> \endverbatim
!>
!> \param[in,out] QSTORE
!> \verbatim
!>          QSTORE is real(dp)    :: array, dimension (N**2+1)
!>         Stores eigenvectors of submatrices encountered during
!>         divide and conquer, packed together. QPTR points to
!>         beginning of the submatrices.
!> \endverbatim
!>
!> \param[in,out] QPTR
!> \verbatim
!>          QPTR is integer(i4)  :: array, dimension (N+2)
!>         List of indices pointing to beginning of submatrices stored
!>         in QSTORE. The submatrices are numbered starting at the
!>         bottom left of the divide and conquer tree, from left to
!>         right and bottom to top.
!> \endverbatim
!>
!> \param[in] PRMPTR
!> \verbatim
!>          PRMPTR is integer(i4)  :: array, dimension (N lg N)
!>         Contains a list of pointers which indicate where in PERM a
!>         level's permutation is stored.  PRMPTR(i+1) - PRMPTR(i)
!>         indicates the size of the permutation and also the size of
!>         the full, non-deflated problem.
!> \endverbatim
!>
!> \param[in] PERM
!> \verbatim
!>          PERM is integer(i4)  :: array, dimension (N lg N)
!>         Contains the permutations (from deflation and sorting) to be
!>         applied to each eigenblock.
!> \endverbatim
!>
!> \param[in] GIVPTR
!> \verbatim
!>          GIVPTR is integer(i4)  :: array, dimension (N lg N)
!>         Contains a list of pointers which indicate where in GIVCOL a
!>         level's Givens rotations are stored.  GIVPTR(i+1) - GIVPTR(i)
!>         indicates the number of Givens rotations.
!> \endverbatim
!>
!> \param[in] GIVCOL
!> \verbatim
!>          GIVCOL is integer(i4)  :: array, dimension (2, N lg N)
!>         Each pair of numbers indicates a pair of columns to take place
!>         in a Givens rotation.
!> \endverbatim
!>
!> \param[in] GIVNUM
!> \verbatim
!>          GIVNUM is real(dp)    :: array, dimension (2, N lg N)
!>         Each number indicates the S value to be used in the
!>         corresponding Givens rotation.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is real(dp)    :: array, dimension (3*N+2*QSIZ*N)
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is integer(i4)  :: array, dimension (4*N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is integer(i4)  ::
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  if INFO = 1, an eigenvalue did not converge
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup auxOTHERcomputational
!
!> \par Contributors:
!  ==================
!>
!> Jeff Rutter, Computer Science Division, University of California
!> at Berkeley, USA
!
!  =====================================================================
      SUBROUTINE DLAED7( ICOMPQ, N, QSIZ, TLVLS, CURLVL, CURPBM, D, Q, &
                         LDQ, INDXQ, RHO, CUTPNT, QSTORE, QPTR, PRMPTR, &
                         PERM, GIVPTR, GIVCOL, GIVNUM, WORK, IWORK, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes, only : i4, dp, dpc
      integer(i4)  ::            CURLVL, CURPBM, CUTPNT, ICOMPQ, INFO, LDQ, N, QSIZ, TLVLS
      real(dp)    ::   RHO
!     ..
!     .. Array Arguments ..
      integer(i4)  ::            GIVCOL( 2, * ), GIVPTR( * ), INDXQ( * ), IWORK( * ), PERM( * ), PRMPTR( * ), QPTR( * )
      real(dp)    ::   D( * ), GIVNUM( 2, * ), Q( LDQ, * ), QSTORE( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      real(dp)    ::   ONE, ZERO
      PARAMETER          ( ONE = 1.0D0, ZERO = 0.0D0 )
!     ..
!     .. Local Scalars ..
      integer(i4)  ::            COLTYP, CURR, I, IDLMDA, INDX, INDXC, INDXP, IQ2, IS, IW, IZ, K, LDQ2, N1, N2, PTR
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGEMM, DLAED8, DLAED9, DLAEDA, DLAMRG, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
!
      IF( ICOMPQ.LT.0 .OR. ICOMPQ.GT.1 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( ICOMPQ.EQ.1 .AND. QSIZ.LT.N ) THEN
         INFO = -4
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( MIN( 1, N ).GT.CUTPNT .OR. N.LT.CUTPNT ) THEN
         INFO = -12
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLAED7', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) RETURN
!
!     The following values are for bookkeeping purposes only.  They are
!     integer pointers which indicate the portion of the workspace
!     used by a particular array in DLAED8 and DLAED9.
!
      IF( ICOMPQ.EQ.1 ) THEN
         LDQ2 = QSIZ
      ELSE
         LDQ2 = N
      END IF
!
      IZ = 1
      IDLMDA = IZ + N
      IW = IDLMDA + N
      IQ2 = IW + N
      IS = IQ2 + N*LDQ2
!
      INDX = 1
      INDXC = INDX + N
      COLTYP = INDXC + N
      INDXP = COLTYP + N
!
!     Form the z-vector which consists of the last row of Q_1 and the
!     first row of Q_2.
!
      PTR = 1 + 2**TLVLS
      DO 10 I = 1, CURLVL - 1
         PTR = PTR + 2**( TLVLS-I )
   10 CONTINUE
      CURR = PTR + CURPBM
      CALL DLAEDA( N, TLVLS, CURLVL, CURPBM, PRMPTR, PERM, GIVPTR, &
                   GIVCOL, GIVNUM, QSTORE, QPTR, WORK( IZ ), WORK( IZ+N ), INFO )
!
!     When solving the final problem, we no longer need the stored data,
!     so we will overwrite the data from this level onto the previously
!     used storage space.
!
      IF( CURLVL.EQ.TLVLS ) THEN
         QPTR( CURR ) = 1
         PRMPTR( CURR ) = 1
         GIVPTR( CURR ) = 1
      END IF
!
!     Sort and Deflate eigenvalues.
!
      CALL DLAED8( ICOMPQ, K, N, QSIZ, D, Q, LDQ, INDXQ, RHO, CUTPNT, &
                   WORK( IZ ), WORK( IDLMDA ), WORK( IQ2 ), LDQ2, &
                   WORK( IW ), PERM( PRMPTR( CURR ) ), GIVPTR( CURR+1 ), &
                   GIVCOL( 1, GIVPTR( CURR ) ), GIVNUM( 1, GIVPTR( CURR ) ), IWORK( INDXP ), IWORK( INDX ), INFO )
      PRMPTR( CURR+1 ) = PRMPTR( CURR ) + N
      GIVPTR( CURR+1 ) = GIVPTR( CURR+1 ) + GIVPTR( CURR )
!
!     Solve Secular Equation.
!
      IF( K.NE.0 ) THEN
         CALL DLAED9( K, 1, K, N, D, WORK( IS ), K, RHO, WORK( IDLMDA ), WORK( IW ), QSTORE( QPTR( CURR ) ), K, INFO )
         IF( INFO.NE.0 ) GO TO 30
         IF( ICOMPQ.EQ.1 ) THEN
            CALL DGEMM( 'N', 'N', QSIZ, K, K, ONE, WORK( IQ2 ), LDQ2, QSTORE( QPTR( CURR ) ), K, ZERO, Q, LDQ )
         END IF
         QPTR( CURR+1 ) = QPTR( CURR ) + K**2
!
!     Prepare the INDXQ sorting permutation.
!
         N1 = K
         N2 = N - K
         CALL DLAMRG( N1, N2, D, 1, -1, INDXQ )
      ELSE
         QPTR( CURR+1 ) = QPTR( CURR )
         DO 20 I = 1, N
            INDXQ( I ) = I
   20    CONTINUE
      END IF
!
   30 CONTINUE
      RETURN
!
!     End of DLAED7
!
      END
!> \brief \b DLAED9
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DLAED9 + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaed9.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaed9.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaed9.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAED9( K, KSTART, KSTOP, N, D, Q, LDQ, RHO, DLAMDA, W,
!                          S, LDS, INFO )
! 
!       .. Scalar Arguments ..
!       integer(i4)  ::            INFO, K, KSTART, KSTOP, LDQ, LDS, N
!       real(dp)    ::   RHO
!       ..
!       .. Array Arguments ..
!       real(dp)    ::   D( * ), DLAMDA( * ), Q( LDQ, * ), S( LDS, * ),
!      $                   W( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAED9 finds the roots of the secular equation, as defined by the
!> values in D, Z, and RHO, between KSTART and KSTOP.  It makes the
!> appropriate calls to DLAED4 and then stores the new matrix of
!> eigenvectors for use in calculating the next level of Z vectors.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] K
!> \verbatim
!>          K is integer(i4)  ::
!>          The number of terms in the rational function to be solved by
!>          DLAED4.  K >= 0.
!> \endverbatim
!>
!> \param[in] KSTART
!> \verbatim
!>          KSTART is integer(i4)  ::
!> \endverbatim
!>
!> \param[in] KSTOP
!> \verbatim
!>          KSTOP is integer(i4)  ::
!>          The updated eigenvalues Lambda(I), KSTART <= I <= KSTOP
!>          are to be computed.  1 <= KSTART <= KSTOP <= K.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is integer(i4)  ::
!>          The number of rows and columns in the Q matrix.
!>          N >= K (delation may result in N > K).
!> \endverbatim
!>
!> \param[out] D
!> \verbatim
!>          D is real(dp)    :: array, dimension (N)
!>          D(I) contains the updated eigenvalues
!>          for KSTART <= I <= KSTOP.
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is real(dp)    :: array, dimension (LDQ,N)
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is integer(i4)  ::
!>          The leading dimension of the array Q.  LDQ >= max( 1, N ).
!> \endverbatim
!>
!> \param[in] RHO
!> \verbatim
!>          RHO is real(dp)    ::
!>          The value of the parameter in the rank one update equation.
!>          RHO >= 0 required.
!> \endverbatim
!>
!> \param[in] DLAMDA
!> \verbatim
!>          DLAMDA is real(dp)    :: array, dimension (K)
!>          The first K elements of this array contain the old roots
!>          of the deflated updating problem.  These are the poles
!>          of the secular equation.
!> \endverbatim
!>
!> \param[in] W
!> \verbatim
!>          W is real(dp)    :: array, dimension (K)
!>          The first K elements of this array contain the components
!>          of the deflation-adjusted updating vector.
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is real(dp)    :: array, dimension (LDS, K)
!>          Will contain the eigenvectors of the repaired matrix which
!>          will be stored for subsequent Z vector calculation and
!>          multiplied by the previously accumulated eigenvectors
!>          to update the system.
!> \endverbatim
!>
!> \param[in] LDS
!> \verbatim
!>          LDS is integer(i4)  ::
!>          The leading dimension of S.  LDS >= max( 1, K ).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is integer(i4)  ::
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  if INFO = 1, an eigenvalue did not converge
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup auxOTHERcomputational
!
!> \par Contributors:
!  ==================
!>
!> Jeff Rutter, Computer Science Division, University of California
!> at Berkeley, USA
!
!  =====================================================================
      SUBROUTINE DLAED9( K, KSTART, KSTOP, N, D, Q, LDQ, RHO, DLAMDA, W, S, LDS, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes, only : i4, dp, dpc
      integer(i4)  ::            INFO, K, KSTART, KSTOP, LDQ, LDS, N
      real(dp)    ::   RHO
!     ..
!     .. Array Arguments ..
      real(dp)    ::   D( * ), DLAMDA( * ), Q( LDQ, * ), S( LDS, * ), W( * )
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      integer(i4)  ::            I, J
      real(dp)    ::   TEMP
!     ..
!     .. External Functions ..
      real(dp)    ::   DLAMC3, DNRM2
      EXTERNAL           DLAMC3, DNRM2
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPY, DLAED4, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, SIGN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
!
      IF( K.LT.0 ) THEN
         INFO = -1
      ELSE IF( KSTART.LT.1 .OR. KSTART.GT.MAX( 1, K ) ) THEN
         INFO = -2
      ELSE IF( MAX( 1, KSTOP ).LT.KSTART .OR. KSTOP.GT.MAX( 1, K ) ) THEN
         INFO = -3
      ELSE IF( N.LT.K ) THEN
         INFO = -4
      ELSE IF( LDQ.LT.MAX( 1, K ) ) THEN
         INFO = -7
      ELSE IF( LDS.LT.MAX( 1, K ) ) THEN
         INFO = -12
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLAED9', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( K.EQ.0 ) RETURN
!
!     Modify values DLAMDA(i) to make sure all DLAMDA(i)-DLAMDA(j) can
!     be computed with high relative accuracy (barring over/underflow).
!     This is a problem on machines without a guard digit in
!     add/subtract (Cray XMP, Cray YMP, Cray C 90 and Cray 2).
!     The following code replaces DLAMDA(I) by 2*DLAMDA(I)-DLAMDA(I),
!     which on any of these machines zeros out the bottommost
!     bit of DLAMDA(I) if it is 1; this makes the subsequent
!     subtractions DLAMDA(I)-DLAMDA(J) unproblematic when cancellation
!     occurs. On binary machines with a guard digit (almost all
!     machines) it does not change DLAMDA(I) at all. On hexadecimal
!     and decimal machines with a guard digit, it slightly
!     changes the bottommost bits of DLAMDA(I). It does not account
!     for hexadecimal or decimal machines without guard digits
!     (we know of none). We use a subroutine call to compute
!     2*DLAMBDA(I) to prevent optimizing compilers from eliminating
!     this code.
!
      DO 10 I = 1, N
         DLAMDA( I ) = DLAMC3( DLAMDA( I ), DLAMDA( I ) ) - DLAMDA( I )
   10 CONTINUE
!
      DO 20 J = KSTART, KSTOP
         CALL DLAED4( K, J, DLAMDA, W, Q( 1, J ), RHO, D( J ), INFO )
!
!        If the zero finder fails, the computation is terminated.
!
         IF( INFO.NE.0 ) GO TO 120
   20 CONTINUE
!
      IF( K.EQ.1 .OR. K.EQ.2 ) THEN
         DO 40 I = 1, K
            DO 30 J = 1, K
               S( J, I ) = Q( J, I )
   30       CONTINUE
   40    CONTINUE
         GO TO 120
      END IF
!
!     Compute updated W.
!
      CALL DCOPY( K, W, 1, S, 1 )
!
!     Initialize W(I) = Q(I,I)
!
      CALL DCOPY( K, Q, LDQ+1, W, 1 )
      DO 70 J = 1, K
         DO 50 I = 1, J - 1
            W( I ) = W( I )*( Q( I, J ) / ( DLAMDA( I )-DLAMDA( J ) ) )
   50    CONTINUE
         DO 60 I = J + 1, K
            W( I ) = W( I )*( Q( I, J ) / ( DLAMDA( I )-DLAMDA( J ) ) )
   60    CONTINUE
   70 CONTINUE
      DO 80 I = 1, K
         W( I ) = SIGN( SQRT( -W( I ) ), S( I, 1 ) )
   80 CONTINUE
!
!     Compute eigenvectors of the modified rank-1 modification.
!
      DO 110 J = 1, K
         DO 90 I = 1, K
            Q( I, J ) = W( I ) / Q( I, J )
   90    CONTINUE
         TEMP = DNRM2( K, Q( 1, J ), 1 )
         DO 100 I = 1, K
            S( I, J ) = Q( I, J ) / TEMP
  100    CONTINUE
  110 CONTINUE
!
  120 CONTINUE
      RETURN
!
!     End of DLAED9
!
      END
!> \brief \b DLAEDA
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DLAEDA + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaeda.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaeda.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaeda.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAEDA( N, TLVLS, CURLVL, CURPBM, PRMPTR, PERM, GIVPTR,
!                          GIVCOL, GIVNUM, Q, QPTR, Z, ZTEMP, INFO )
! 
!       .. Scalar Arguments ..
!       integer(i4)  ::            CURLVL, CURPBM, INFO, N, TLVLS
!       ..
!       .. Array Arguments ..
!       integer(i4)  ::            GIVCOL( 2, * ), GIVPTR( * ), PERM( * ),
!      $                   PRMPTR( * ), QPTR( * )
!       real(dp)    ::   GIVNUM( 2, * ), Q( * ), Z( * ), ZTEMP( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAEDA computes the Z vector corresponding to the merge step in the
!> CURLVLth step of the merge process with TLVLS steps for the CURPBMth
!> problem.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is integer(i4)  ::
!>         The dimension of the symmetric tridiagonal matrix.  N >= 0.
!> \endverbatim
!>
!> \param[in] TLVLS
!> \verbatim
!>          TLVLS is integer(i4)  ::
!>         The total number of merging levels in the overall divide and
!>         conquer tree.
!> \endverbatim
!>
!> \param[in] CURLVL
!> \verbatim
!>          CURLVL is integer(i4)  ::
!>         The current level in the overall merge routine,
!>         0 <= curlvl <= tlvls.
!> \endverbatim
!>
!> \param[in] CURPBM
!> \verbatim
!>          CURPBM is integer(i4)  ::
!>         The current problem in the current level in the overall
!>         merge routine (counting from upper left to lower right).
!> \endverbatim
!>
!> \param[in] PRMPTR
!> \verbatim
!>          PRMPTR is integer(i4)  :: array, dimension (N lg N)
!>         Contains a list of pointers which indicate where in PERM a
!>         level's permutation is stored.  PRMPTR(i+1) - PRMPTR(i)
!>         indicates the size of the permutation and incidentally the
!>         size of the full, non-deflated problem.
!> \endverbatim
!>
!> \param[in] PERM
!> \verbatim
!>          PERM is integer(i4)  :: array, dimension (N lg N)
!>         Contains the permutations (from deflation and sorting) to be
!>         applied to each eigenblock.
!> \endverbatim
!>
!> \param[in] GIVPTR
!> \verbatim
!>          GIVPTR is integer(i4)  :: array, dimension (N lg N)
!>         Contains a list of pointers which indicate where in GIVCOL a
!>         level's Givens rotations are stored.  GIVPTR(i+1) - GIVPTR(i)
!>         indicates the number of Givens rotations.
!> \endverbatim
!>
!> \param[in] GIVCOL
!> \verbatim
!>          GIVCOL is integer(i4)  :: array, dimension (2, N lg N)
!>         Each pair of numbers indicates a pair of columns to take place
!>         in a Givens rotation.
!> \endverbatim
!>
!> \param[in] GIVNUM
!> \verbatim
!>          GIVNUM is real(dp)    :: array, dimension (2, N lg N)
!>         Each number indicates the S value to be used in the
!>         corresponding Givens rotation.
!> \endverbatim
!>
!> \param[in] Q
!> \verbatim
!>          Q is real(dp)    :: array, dimension (N**2)
!>         Contains the square eigenblocks from previous levels, the
!>         starting positions for blocks are given by QPTR.
!> \endverbatim
!>
!> \param[in] QPTR
!> \verbatim
!>          QPTR is integer(i4)  :: array, dimension (N+2)
!>         Contains a list of pointers which indicate where in Q an
!>         eigenblock is stored.  SQRT( QPTR(i+1) - QPTR(i) ) indicates
!>         the size of the block.
!> \endverbatim
!>
!> \param[out] Z
!> \verbatim
!>          Z is real(dp)    :: array, dimension (N)
!>         On output this vector contains the updating vector (the last
!>         row of the first sub-eigenvector matrix and the first row of
!>         the second sub-eigenvector matrix).
!> \endverbatim
!>
!> \param[out] ZTEMP
!> \verbatim
!>          ZTEMP is real(dp)    :: array, dimension (N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is integer(i4)  ::
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup auxOTHERcomputational
!
!> \par Contributors:
!  ==================
!>
!> Jeff Rutter, Computer Science Division, University of California
!> at Berkeley, USA
!
!  =====================================================================
      SUBROUTINE DLAEDA( N, TLVLS, CURLVL, CURPBM, PRMPTR, PERM, GIVPTR, GIVCOL, GIVNUM, Q, QPTR, Z, ZTEMP, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes, only : i4, dp, dpc
      integer(i4)  ::            CURLVL, CURPBM, INFO, N, TLVLS
!     ..
!     .. Array Arguments ..
      integer(i4)  ::            GIVCOL( 2, * ), GIVPTR( * ), PERM( * ), PRMPTR( * ), QPTR( * )
      real(dp)    ::   GIVNUM( 2, * ), Q( * ), Z( * ), ZTEMP( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      real(dp)    ::   ZERO, HALF, ONE
      PARAMETER          ( ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
      integer(i4)  ::            BSIZ1, BSIZ2, CURR, I, K, MID, PSIZ1, PSIZ2, PTR, ZPTR1
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMV, DROT, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DBLE, INT, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
!
      IF( N.LT.0 ) THEN
         INFO = -1
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLAEDA', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) RETURN
!
!     Determine location of first number in second half.
!
      MID = N / 2 + 1
!
!     Gather last/first rows of appropriate eigenblocks into center of Z
!
      PTR = 1
!
!     Determine location of lowest level subproblem in the full storage
!     scheme
!
      CURR = PTR + CURPBM*2**CURLVL + 2**( CURLVL-1 ) - 1
!
!     Determine size of these matrices.  We add HALF to the value of
!     the SQRT in case the machine underestimates one of these square
!     roots.
!
      BSIZ1 = INT( HALF+SQRT( DBLE( QPTR( CURR+1 )-QPTR( CURR ) ) ) )
      BSIZ2 = INT( HALF+SQRT( DBLE( QPTR( CURR+2 )-QPTR( CURR+1 ) ) ) )
      DO 10 K = 1, MID - BSIZ1 - 1
         Z( K ) = ZERO
   10 CONTINUE
      CALL DCOPY( BSIZ1, Q( QPTR( CURR )+BSIZ1-1 ), BSIZ1, Z( MID-BSIZ1 ), 1 )
      CALL DCOPY( BSIZ2, Q( QPTR( CURR+1 ) ), BSIZ2, Z( MID ), 1 )
      DO 20 K = MID + BSIZ2, N
         Z( K ) = ZERO
   20 CONTINUE
!
!     Loop through remaining levels 1 -> CURLVL applying the Givens
!     rotations and permutation and then multiplying the center matrices
!     against the current Z.
!
      PTR = 2**TLVLS + 1
      DO 70 K = 1, CURLVL - 1
         CURR = PTR + CURPBM*2**( CURLVL-K ) + 2**( CURLVL-K-1 ) - 1
         PSIZ1 = PRMPTR( CURR+1 ) - PRMPTR( CURR )
         PSIZ2 = PRMPTR( CURR+2 ) - PRMPTR( CURR+1 )
         ZPTR1 = MID - PSIZ1
!
!       Apply Givens at CURR and CURR+1
!
         DO 30 I = GIVPTR( CURR ), GIVPTR( CURR+1 ) - 1
            CALL DROT( 1, Z( ZPTR1+GIVCOL( 1, I )-1 ), 1, Z( ZPTR1+GIVCOL( 2, I )-1 ), 1, GIVNUM( 1, I ), GIVNUM( 2, I ) )
   30    CONTINUE
         DO 40 I = GIVPTR( CURR+1 ), GIVPTR( CURR+2 ) - 1
            CALL DROT( 1, Z( MID-1+GIVCOL( 1, I ) ), 1, Z( MID-1+GIVCOL( 2, I ) ), 1, GIVNUM( 1, I ), GIVNUM( 2, I ) )
   40    CONTINUE
         PSIZ1 = PRMPTR( CURR+1 ) - PRMPTR( CURR )
         PSIZ2 = PRMPTR( CURR+2 ) - PRMPTR( CURR+1 )
         DO 50 I = 0, PSIZ1 - 1
            ZTEMP( I+1 ) = Z( ZPTR1+PERM( PRMPTR( CURR )+I )-1 )
   50    CONTINUE
         DO 60 I = 0, PSIZ2 - 1
            ZTEMP( PSIZ1+I+1 ) = Z( MID+PERM( PRMPTR( CURR+1 )+I )-1 )
   60    CONTINUE
!
!        Multiply Blocks at CURR and CURR+1
!
!        Determine size of these matrices.  We add HALF to the value of
!        the SQRT in case the machine underestimates one of these
!        square roots.
!
         BSIZ1 = INT( HALF+SQRT( DBLE( QPTR( CURR+1 )-QPTR( CURR ) ) ) )
         BSIZ2 = INT( HALF+SQRT( DBLE( QPTR( CURR+2 )-QPTR( CURR+ 1 ) ) ) )
         IF( BSIZ1.GT.0 ) THEN
            CALL DGEMV( 'T', BSIZ1, BSIZ1, ONE, Q( QPTR( CURR ) ), BSIZ1, ZTEMP( 1 ), 1, ZERO, Z( ZPTR1 ), 1 )
         END IF
         CALL DCOPY( PSIZ1-BSIZ1, ZTEMP( BSIZ1+1 ), 1, Z( ZPTR1+BSIZ1 ), 1 )
         IF( BSIZ2.GT.0 ) THEN
            CALL DGEMV( 'T', BSIZ2, BSIZ2, ONE, Q( QPTR( CURR+1 ) ), BSIZ2, ZTEMP( PSIZ1+1 ), 1, ZERO, Z( MID ), 1 )
         END IF
         CALL DCOPY( PSIZ2-BSIZ2, ZTEMP( PSIZ1+BSIZ2+1 ), 1, Z( MID+BSIZ2 ), 1 )
!
         PTR = PTR + 2**( TLVLS-K )
   70 CONTINUE
!
      RETURN
!
!     End of DLAEDA
!
      END
!> \brief \b DLAMRG
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DLAMRG + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlamrg.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlamrg.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlamrg.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAMRG( N1, N2, A, DTRD1, DTRD2, INDEX )
! 
!       .. Scalar Arguments ..
!       integer(i4)  ::            DTRD1, DTRD2, N1, N2
!       ..
!       .. Array Arguments ..
!       integer(i4)  ::            INDEX( * )
!       real(dp)    ::   A( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAMRG will create a permutation list which will merge the elements
!> of A (which is composed of two independently sorted sets) into a
!> single set which is sorted in ascending order.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N1
!> \verbatim
!>          N1 is integer(i4)  ::
!> \endverbatim
!>
!> \param[in] N2
!> \verbatim
!>          N2 is integer(i4)  ::
!>         These arguements contain the respective lengths of the two
!>         sorted lists to be merged.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is real(dp)    :: array, dimension (N1+N2)
!>         The first N1 elements of A contain a list of numbers which
!>         are sorted in either ascending or descending order.  Likewise
!>         for the final N2 elements.
!> \endverbatim
!>
!> \param[in] DTRD1
!> \verbatim
!>          DTRD1 is integer(i4)  ::
!> \endverbatim
!>
!> \param[in] DTRD2
!> \verbatim
!>          DTRD2 is integer(i4)  ::
!>         These are the strides to be taken through the array A.
!>         Allowable strides are 1 and -1.  They indicate whether a
!>         subset of A is sorted in ascending (DTRDx = 1) or descending
!>         (DTRDx = -1) order.
!> \endverbatim
!>
!> \param[out] INDEX
!> \verbatim
!>          INDEX is integer(i4)  :: array, dimension (N1+N2)
!>         On exit this array will contain a permutation such that
!>         if B( I ) = A( INDEX( I ) ) for I=1,N1+N2, then B will be
!>         sorted in ascending order.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup auxOTHERcomputational
!
!  =====================================================================
      SUBROUTINE DLAMRG( N1, N2, A, DTRD1, DTRD2, INDEX )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes, only : i4, dp, dpc
      integer(i4)  ::            DTRD1, DTRD2, N1, N2
!     ..
!     .. Array Arguments ..
      integer(i4)  ::            INDEX( * )
      real(dp)    ::   A( * )
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      integer(i4)  ::            I, IND1, IND2, N1SV, N2SV
!     ..
!     .. Executable Statements ..
!
      N1SV = N1
      N2SV = N2
      IF( DTRD1.GT.0 ) THEN
         IND1 = 1
      ELSE
         IND1 = N1
      END IF
      IF( DTRD2.GT.0 ) THEN
         IND2 = 1 + N1
      ELSE
         IND2 = N1 + N2
      END IF
      I = 1
!     while ( (N1SV > 0) & (N2SV > 0) )
   10 CONTINUE
      IF( N1SV.GT.0 .AND. N2SV.GT.0 ) THEN
         IF( A( IND1 ).LE.A( IND2 ) ) THEN
            INDEX( I ) = IND1
            I = I + 1
            IND1 = IND1 + DTRD1
            N1SV = N1SV - 1
         ELSE
            INDEX( I ) = IND2
            I = I + 1
            IND2 = IND2 + DTRD2
            N2SV = N2SV - 1
         END IF
         GO TO 10
      END IF
!     end while
      IF( N1SV.EQ.0 ) THEN
         DO 20 N1SV = 1, N2SV
            INDEX( I ) = IND2
            I = I + 1
            IND2 = IND2 + DTRD2
   20    CONTINUE
      ELSE
!     N2SV .EQ. 0
         DO 30 N2SV = 1, N1SV
            INDEX( I ) = IND1
            I = I + 1
            IND1 = IND1 + DTRD1
   30    CONTINUE
      END IF
!
      RETURN
!
!     End of DLAMRG
!
      END

!> \brief \b DLAED2
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DLAED2 + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaed2.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaed2.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaed2.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAED2( K, N, N1, D, Q, LDQ, INDXQ, RHO, Z, DLAMDA, W,
!                          Q2, INDX, INDXC, INDXP, COLTYP, INFO )
! 
!       .. Scalar Arguments ..
!       integer(i4)  ::          INFO, K, LDQ, N, N1
!       real(dp)    ::    RHO
!       ..
!       .. Array Arguments ..
!       integer(i4)  ::          COLTYP( * ), INDX( * ), INDXC( * ), INDXP( * ),
!      $                   INDXQ( * )
!       real(dp)    ::    D( * ), DLAMDA( * ), Q( LDQ, * ), Q2( * ),
!      $                   W( * ), Z( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAED2 merges the two sets of eigenvalues together into a single
!> sorted set.  Then it tries to deflate the size of the problem.
!> There are two ways in which deflation can occur:  when two or more
!> eigenvalues are close together or if there is a tiny entry in the
!> Z vector.  For each such occurrence the order of the related secular
!> equation problem is reduced by one.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[out] K
!> \verbatim
!>          K is INTEGER
!>         The number of non-deflated eigenvalues, and the order of the
!>         related secular equation. 0 <= K <=N.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         The dimension of the symmetric tridiagonal matrix.  N >= 0.
!> \endverbatim
!>
!> \param[in] N1
!> \verbatim
!>          N1 is INTEGER
!>         The location of the last eigenvalue in the leading sub-matrix.
!>         min(1,N) <= N1 <= N/2.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is real(dp)    ::  array, dimension (N)
!>         On entry, D contains the eigenvalues of the two submatrices to
!>         be combined.
!>         On exit, D contains the trailing (N-K) updated eigenvalues
!>         (those which were deflated) sorted into increasing order.
!> \endverbatim
!>
!> \param[in,out] Q
!> \verbatim
!>          Q is real(dp)    ::  array, dimension (LDQ, N)
!>         On entry, Q contains the eigenvectors of two submatrices in
!>         the two square blocks with corners at (1,1), (N1,N1)
!>         and (N1+1, N1+1), (N,N).
!>         On exit, Q contains the trailing (N-K) updated eigenvectors
!>         (those which were deflated) in its last N-K columns.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>         The leading dimension of the array Q.  LDQ >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] INDXQ
!> \verbatim
!>          INDXQ is INTEGER array, dimension (N)
!>         The permutation which separately sorts the two sub-problems
!>         in D into ascending order.  Note that elements in the second
!>         half of this permutation must first have N1 added to their
!>         values. Destroyed on exit.
!> \endverbatim
!>
!> \param[in,out] RHO
!> \verbatim
!>          RHO is real(dp)    :: 
!>         On entry, the off-diagonal element associated with the rank-1
!>         cut which originally split the two submatrices which are now
!>         being recombined.
!>         On exit, RHO has been modified to the value required by
!>         DLAED3.
!> \endverbatim
!>
!> \param[in] Z
!> \verbatim
!>          Z is real(dp)    ::  array, dimension (N)
!>         On entry, Z contains the updating vector (the last
!>         row of the first sub-eigenvector matrix and the first row of
!>         the second sub-eigenvector matrix).
!>         On exit, the contents of Z have been destroyed by the updating
!>         process.
!> \endverbatim
!>
!> \param[out] DLAMDA
!> \verbatim
!>          DLAMDA is real(dp)    ::  array, dimension (N)
!>         A copy of the first K eigenvalues which will be used by
!>         DLAED3 to form the secular equation.
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is real(dp)    ::  array, dimension (N)
!>         The first k values of the final deflation-altered z-vector
!>         which will be passed to DLAED3.
!> \endverbatim
!>
!> \param[out] Q2
!> \verbatim
!>          Q2 is real(dp)    ::  array, dimension (N1**2+(N-N1)**2)
!>         A copy of the first K eigenvectors which will be used by
!>         DLAED3 in a matrix multiply (DGEMM) to solve for the new
!>         eigenvectors.
!> \endverbatim
!>
!> \param[out] INDX
!> \verbatim
!>          INDX is INTEGER array, dimension (N)
!>         The permutation used to sort the contents of DLAMDA into
!>         ascending order.
!> \endverbatim
!>
!> \param[out] INDXC
!> \verbatim
!>          INDXC is INTEGER array, dimension (N)
!>         The permutation used to arrange the columns of the deflated
!>         Q matrix into three groups:  the first group contains non-zero
!>         elements only at and above N1, the second contains
!>         non-zero elements only below N1, and the third is dense.
!> \endverbatim
!>
!> \param[out] INDXP
!> \verbatim
!>          INDXP is INTEGER array, dimension (N)
!>         The permutation used to place deflated values of D at the end
!>         of the array.  INDXP(1:K) points to the nondeflated D-values
!>         and INDXP(K+1:N) points to the deflated eigenvalues.
!> \endverbatim
!>
!> \param[out] COLTYP
!> \verbatim
!>          COLTYP is INTEGER array, dimension (N)
!>         During execution, a label which will indicate which of the
!>         following types a column in the Q2 matrix is:
!>         1 : non-zero in the upper half only;
!>         2 : dense;
!>         3 : non-zero in the lower half only;
!>         4 : deflated.
!>         On exit, COLTYP(i) is the number of columns of type i,
!>         for i=1 to 4 only.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup auxOTHERcomputational
!
!> \par Contributors:
!  ==================
!>
!> Jeff Rutter, Computer Science Division, University of California
!> at Berkeley, USA \n
!>  Modified by Francoise Tisseur, University of Tennessee
!>
!  =====================================================================
      SUBROUTINE DLAED2( K, N, N1, D, Q, LDQ, INDXQ, RHO, Z, DLAMDA, W, Q2, INDX, INDXC, INDXP, COLTYP, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes
      integer(i4)  ::          INFO, K, LDQ, N, N1
      real(dp)    ::    RHO
!     ..
!     .. Array Arguments ..
      integer(i4)  ::          COLTYP( * ), INDX( * ), INDXC( * ), INDXP( * ), INDXQ( * )
      real(dp)    ::    D( * ), DLAMDA( * ), Q( LDQ, * ), Q2( * ), W( * ), Z( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      real(dp)    ::    MONE, ZERO, ONE, TWO, EIGHT
      PARAMETER          ( MONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, EIGHT = 8.0D0 )
!     ..
!     .. Local Arrays ..
      integer(i4)  ::          CTOT( 4 ), PSM( 4 )
!     ..
!     .. Local Scalars ..
      integer(i4)  ::          CT, I, IMAX, IQ1, IQ2, J, JMAX, JS, K2, N1P1, N2, NJ, PJ
      real(dp)    ::    C, EPS, S, T, TAU, TOL
!     ..
!     .. External Functions ..
      integer(i4)  ::          IDAMAX
      real(dp)    ::    DLAMCH, DLAPY2
      EXTERNAL           IDAMAX, DLAMCH, DLAPY2
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPY, DLACPY, DLAMRG, DROT, DSCAL, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
!
      IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( MIN( 1, ( N / 2 ) ).GT.N1 .OR. ( N / 2 ).LT.N1 ) THEN
         INFO = -3
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLAED2', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) RETURN
!
      N2 = N - N1
      N1P1 = N1 + 1
!
      IF( RHO.LT.ZERO ) THEN
         CALL DSCAL( N2, MONE, Z( N1P1 ), 1 )
      END IF
!
!     Normalize z so that norm(z) = 1.  Since z is the concatenation of
!     two normalized vectors, norm2(z) = sqrt(2).
!
      T = ONE / SQRT( TWO )
      CALL DSCAL( N, T, Z, 1 )
!
!     RHO = ABS( norm(z)**2 * RHO )
!
      RHO = ABS( TWO*RHO )
!
!     Sort the eigenvalues into increasing order
!
      DO 10 I = N1P1, N
         INDXQ( I ) = INDXQ( I ) + N1
   10 CONTINUE
!
!     re-integrate the deflated parts from the last pass
!
      DO 20 I = 1, N
         DLAMDA( I ) = D( INDXQ( I ) )
   20 CONTINUE
      CALL DLAMRG( N1, N2, DLAMDA, 1, 1, INDXC )
      DO 30 I = 1, N
         INDX( I ) = INDXQ( INDXC( I ) )
   30 CONTINUE
!
!     Calculate the allowable deflation tolerance
!
      IMAX = IDAMAX( N, Z, 1 )
      JMAX = IDAMAX( N, D, 1 )
      EPS = DLAMCH( 'Epsilon' )
      TOL = EIGHT*EPS*MAX( ABS( D( JMAX ) ), ABS( Z( IMAX ) ) )
!
!     If the rank-1 modifier is small enough, no more needs to be done
!     except to reorganize Q so that its columns correspond with the
!     elements in D.
!
      IF( RHO*ABS( Z( IMAX ) ).LE.TOL ) THEN
         K = 0
         IQ2 = 1
         DO 40 J = 1, N
            I = INDX( J )
            CALL DCOPY( N, Q( 1, I ), 1, Q2( IQ2 ), 1 )
            DLAMDA( J ) = D( I )
            IQ2 = IQ2 + N
   40    CONTINUE
         CALL DLACPY( 'A', N, N, Q2, N, Q, LDQ )
         CALL DCOPY( N, DLAMDA, 1, D, 1 )
         GO TO 190
      END IF
!
!     If there are multiple eigenvalues then the problem deflates.  Here
!     the number of equal eigenvalues are found.  As each equal
!     eigenvalue is found, an elementary reflector is computed to rotate
!     the corresponding eigensubspace so that the corresponding
!     components of Z are zero in this new basis.
!
      DO 50 I = 1, N1
         COLTYP( I ) = 1
   50 CONTINUE
      DO 60 I = N1P1, N
         COLTYP( I ) = 3
   60 CONTINUE
!
!
      K = 0
      K2 = N + 1
      DO 70 J = 1, N
         NJ = INDX( J )
         IF( RHO*ABS( Z( NJ ) ).LE.TOL ) THEN
!
!           Deflate due to small z component.
!
            K2 = K2 - 1
            COLTYP( NJ ) = 4
            INDXP( K2 ) = NJ
            IF( J.EQ.N ) GO TO 100
         ELSE
            PJ = NJ
            GO TO 80
         END IF
   70 CONTINUE
   80 CONTINUE
      J = J + 1
      NJ = INDX( J )
      IF( J.GT.N ) GO TO 100
      IF( RHO*ABS( Z( NJ ) ).LE.TOL ) THEN
!
!        Deflate due to small z component.
!
         K2 = K2 - 1
         COLTYP( NJ ) = 4
         INDXP( K2 ) = NJ
      ELSE
!
!        Check if eigenvalues are close enough to allow deflation.
!
         S = Z( PJ )
         C = Z( NJ )
!
!        Find sqrt(a**2+b**2) without overflow or
!        destructive underflow.
!
         TAU = DLAPY2( C, S )
         T = D( NJ ) - D( PJ )
         C = C / TAU
         S = -S / TAU
         IF( ABS( T*C*S ).LE.TOL ) THEN
!
!           Deflation is possible.
!
            Z( NJ ) = TAU
            Z( PJ ) = ZERO
            IF( COLTYP( NJ ).NE.COLTYP( PJ ) ) COLTYP( NJ ) = 2
            COLTYP( PJ ) = 4
            CALL DROT( N, Q( 1, PJ ), 1, Q( 1, NJ ), 1, C, S )
            T = D( PJ )*C**2 + D( NJ )*S**2
            D( NJ ) = D( PJ )*S**2 + D( NJ )*C**2
            D( PJ ) = T
            K2 = K2 - 1
            I = 1
   90       CONTINUE
            IF( K2+I.LE.N ) THEN
               IF( D( PJ ).LT.D( INDXP( K2+I ) ) ) THEN
                  INDXP( K2+I-1 ) = INDXP( K2+I )
                  INDXP( K2+I ) = PJ
                  I = I + 1
                  GO TO 90
               ELSE
                  INDXP( K2+I-1 ) = PJ
               END IF
            ELSE
               INDXP( K2+I-1 ) = PJ
            END IF
            PJ = NJ
         ELSE
            K = K + 1
            DLAMDA( K ) = D( PJ )
            W( K ) = Z( PJ )
            INDXP( K ) = PJ
            PJ = NJ
         END IF
      END IF
      GO TO 80
  100 CONTINUE
!
!     Record the last eigenvalue.
!
      K = K + 1
      DLAMDA( K ) = D( PJ )
      W( K ) = Z( PJ )
      INDXP( K ) = PJ
!
!     Count up the total number of the various types of columns, then
!     form a permutation which positions the four column types into
!     four uniform groups (although one or more of these groups may be
!     empty).
!
      DO 110 J = 1, 4
         CTOT( J ) = 0
  110 CONTINUE
      DO 120 J = 1, N
         CT = COLTYP( J )
         CTOT( CT ) = CTOT( CT ) + 1
  120 CONTINUE
!
!     PSM(*) = Position in SubMatrix (of types 1 through 4)
!
      PSM( 1 ) = 1
      PSM( 2 ) = 1 + CTOT( 1 )
      PSM( 3 ) = PSM( 2 ) + CTOT( 2 )
      PSM( 4 ) = PSM( 3 ) + CTOT( 3 )
      K = N - CTOT( 4 )
!
!     Fill out the INDXC array so that the permutation which it induces
!     will place all type-1 columns first, all type-2 columns next,
!     then all type-3's, and finally all type-4's.
!
      DO 130 J = 1, N
         JS = INDXP( J )
         CT = COLTYP( JS )
         INDX( PSM( CT ) ) = JS
         INDXC( PSM( CT ) ) = J
         PSM( CT ) = PSM( CT ) + 1
  130 CONTINUE
!
!     Sort the eigenvalues and corresponding eigenvectors into DLAMDA
!     and Q2 respectively.  The eigenvalues/vectors which were not
!     deflated go into the first K slots of DLAMDA and Q2 respectively,
!     while those which were deflated go into the last N - K slots.
!
      I = 1
      IQ1 = 1
      IQ2 = 1 + ( CTOT( 1 )+CTOT( 2 ) )*N1
      DO 140 J = 1, CTOT( 1 )
         JS = INDX( I )
         CALL DCOPY( N1, Q( 1, JS ), 1, Q2( IQ1 ), 1 )
         Z( I ) = D( JS )
         I = I + 1
         IQ1 = IQ1 + N1
  140 CONTINUE
!
      DO 150 J = 1, CTOT( 2 )
         JS = INDX( I )
         CALL DCOPY( N1, Q( 1, JS ), 1, Q2( IQ1 ), 1 )
         CALL DCOPY( N2, Q( N1+1, JS ), 1, Q2( IQ2 ), 1 )
         Z( I ) = D( JS )
         I = I + 1
         IQ1 = IQ1 + N1
         IQ2 = IQ2 + N2
  150 CONTINUE
!
      DO 160 J = 1, CTOT( 3 )
         JS = INDX( I )
         CALL DCOPY( N2, Q( N1+1, JS ), 1, Q2( IQ2 ), 1 )
         Z( I ) = D( JS )
         I = I + 1
         IQ2 = IQ2 + N2
  160 CONTINUE
!
      IQ1 = IQ2
      DO 170 J = 1, CTOT( 4 )
         JS = INDX( I )
         CALL DCOPY( N, Q( 1, JS ), 1, Q2( IQ2 ), 1 )
         IQ2 = IQ2 + N
         Z( I ) = D( JS )
         I = I + 1
  170 CONTINUE
!
!     The deflated eigenvalues and their corresponding vectors go back
!     into the last N - K slots of D and Q respectively.
!
      IF( K.LT.N ) THEN
         CALL DLACPY( 'A', N, CTOT( 4 ), Q2( IQ1 ), N, Q( 1, K+1 ), LDQ )
         CALL DCOPY( N-K, Z( K+1 ), 1, D( K+1 ), 1 )
      END IF         
!
!     Copy CTOT into COLTYP for referencing in DLAED3.
!
      DO 180 J = 1, 4
         COLTYP( J ) = CTOT( J )
  180 CONTINUE
!
  190 CONTINUE
      RETURN
!
!     End of DLAED2
!
      END
!> \brief \b DLAED3
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DLAED3 + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaed3.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaed3.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaed3.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAED3( K, N, N1, D, Q, LDQ, RHO, DLAMDA, Q2, INDX,
!                          CTOT, W, S, INFO )
! 
!       .. Scalar Arguments ..
!       integer(i4)  ::          INFO, K, LDQ, N, N1
!       real(dp)    ::    RHO
!       ..
!       .. Array Arguments ..
!       integer(i4)  ::          CTOT( * ), INDX( * )
!       real(dp)    ::    D( * ), DLAMDA( * ), Q( LDQ, * ), Q2( * ),
!      $                   S( * ), W( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAED3 finds the roots of the secular equation, as defined by the
!> values in D, W, and RHO, between 1 and K.  It makes the
!> appropriate calls to DLAED4 and then updates the eigenvectors by
!> multiplying the matrix of eigenvectors of the pair of eigensystems
!> being combined by the matrix of eigenvectors of the K-by-K system
!> which is solved here.
!>
!> This code makes very mild assumptions about floating point
!> arithmetic. It will work on machines with a guard digit in
!> add/subtract, or on those binary machines without guard digits
!> which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
!> It could conceivably fail on hexadecimal or decimal machines
!> without guard digits, but we know of none.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of terms in the rational function to be solved by
!>          DLAED4.  K >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of rows and columns in the Q matrix.
!>          N >= K (deflation may result in N>K).
!> \endverbatim
!>
!> \param[in] N1
!> \verbatim
!>          N1 is INTEGER
!>          The location of the last eigenvalue in the leading submatrix.
!>          min(1,N) <= N1 <= N/2.
!> \endverbatim
!>
!> \param[out] D
!> \verbatim
!>          D is real(dp)    ::  array, dimension (N)
!>          D(I) contains the updated eigenvalues for
!>          1 <= I <= K.
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is real(dp)    ::  array, dimension (LDQ,N)
!>          Initially the first K columns are used as workspace.
!>          On output the columns 1 to K contain
!>          the updated eigenvectors.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of the array Q.  LDQ >= max(1,N).
!> \endverbatim
!>
!> \param[in] RHO
!> \verbatim
!>          RHO is real(dp)    :: 
!>          The value of the parameter in the rank one update equation.
!>          RHO >= 0 required.
!> \endverbatim
!>
!> \param[in,out] DLAMDA
!> \verbatim
!>          DLAMDA is real(dp)    ::  array, dimension (K)
!>          The first K elements of this array contain the old roots
!>          of the deflated updating problem.  These are the poles
!>          of the secular equation. May be changed on output by
!>          having lowest order bit set to zero on Cray X-MP, Cray Y-MP,
!>          Cray-2, or Cray C-90, as described above.
!> \endverbatim
!>
!> \param[in] Q2
!> \verbatim
!>          Q2 is real(dp)    ::  array, dimension (LDQ2, N)
!>          The first K columns of this matrix contain the non-deflated
!>          eigenvectors for the split problem.
!> \endverbatim
!>
!> \param[in] INDX
!> \verbatim
!>          INDX is INTEGER array, dimension (N)
!>          The permutation used to arrange the columns of the deflated
!>          Q matrix into three groups (see DLAED2).
!>          The rows of the eigenvectors found by DLAED4 must be likewise
!>          permuted before the matrix multiply can take place.
!> \endverbatim
!>
!> \param[in] CTOT
!> \verbatim
!>          CTOT is INTEGER array, dimension (4)
!>          A count of the total number of the various types of columns
!>          in Q, as described in INDX.  The fourth column type is any
!>          column which has been deflated.
!> \endverbatim
!>
!> \param[in,out] W
!> \verbatim
!>          W is real(dp)    ::  array, dimension (K)
!>          The first K elements of this array contain the components
!>          of the deflation-adjusted updating vector. Destroyed on
!>          output.
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is real(dp)    ::  array, dimension (N1 + 1)*K
!>          Will contain the eigenvectors of the repaired matrix which
!>          will be multiplied by the previously accumulated eigenvectors
!>          to update the system.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  if INFO = 1, an eigenvalue did not converge
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup auxOTHERcomputational
!
!> \par Contributors:
!  ==================
!>
!> Jeff Rutter, Computer Science Division, University of California
!> at Berkeley, USA \n
!>  Modified by Francoise Tisseur, University of Tennessee
!>
!  =====================================================================
      SUBROUTINE DLAED3( K, N, N1, D, Q, LDQ, RHO, DLAMDA, Q2, INDX, CTOT, W, S, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes
      integer(i4)  ::          INFO, K, LDQ, N, N1
      real(dp)    ::    RHO
!     ..
!     .. Array Arguments ..
      integer(i4)  ::          CTOT( * ), INDX( * )
      real(dp)    ::    D( * ), DLAMDA( * ), Q( LDQ, * ), Q2( * ), S( * ), W( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      real(dp)    ::    ONE, ZERO
      PARAMETER          ( ONE = 1.0D0, ZERO = 0.0D0 )
!     ..
!     .. Local Scalars ..
      integer(i4)  ::          I, II, IQ2, J, N12, N2, N23
      real(dp)    ::    TEMP
!     ..
!     .. External Functions ..
      real(dp)    ::    DLAMC3, DNRM2
      EXTERNAL           DLAMC3, DNRM2
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMM, DLACPY, DLAED4, DLASET, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, SIGN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
!
      IF( K.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.K ) THEN
         INFO = -2
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLAED3', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( K.EQ.0 ) RETURN
!
!     Modify values DLAMDA(i) to make sure all DLAMDA(i)-DLAMDA(j) can
!     be computed with high relative accuracy (barring over/underflow).
!     This is a problem on machines without a guard digit in
!     add/subtract (Cray XMP, Cray YMP, Cray C 90 and Cray 2).
!     The following code replaces DLAMDA(I) by 2*DLAMDA(I)-DLAMDA(I),
!     which on any of these machines zeros out the bottommost
!     bit of DLAMDA(I) if it is 1; this makes the subsequent
!     subtractions DLAMDA(I)-DLAMDA(J) unproblematic when cancellation
!     occurs. On binary machines with a guard digit (almost all
!     machines) it does not change DLAMDA(I) at all. On hexadecimal
!     and decimal machines with a guard digit, it slightly
!     changes the bottommost bits of DLAMDA(I). It does not account
!     for hexadecimal or decimal machines without guard digits
!     (we know of none). We use a subroutine call to compute
!     2*DLAMBDA(I) to prevent optimizing compilers from eliminating
!     this code.
!
      DO 10 I = 1, K
         DLAMDA( I ) = DLAMC3( DLAMDA( I ), DLAMDA( I ) ) - DLAMDA( I )
   10 CONTINUE
!
      DO 20 J = 1, K
         CALL DLAED4( K, J, DLAMDA, W, Q( 1, J ), RHO, D( J ), INFO )
!
!        If the zero finder fails, the computation is terminated.
!
         IF( INFO.NE.0 ) GO TO 120
   20 CONTINUE
!
      IF( K.EQ.1 ) GO TO 110
      IF( K.EQ.2 ) THEN
         DO 30 J = 1, K
            W( 1 ) = Q( 1, J )
            W( 2 ) = Q( 2, J )
            II = INDX( 1 )
            Q( 1, J ) = W( II )
            II = INDX( 2 )
            Q( 2, J ) = W( II )
   30    CONTINUE
         GO TO 110
      END IF
!
!     Compute updated W.
!
      CALL DCOPY( K, W, 1, S, 1 )
!
!     Initialize W(I) = Q(I,I)
!
      CALL DCOPY( K, Q, LDQ+1, W, 1 )
      DO 60 J = 1, K
         DO 40 I = 1, J - 1
            W( I ) = W( I )*( Q( I, J ) / ( DLAMDA( I )-DLAMDA( J ) ) )
   40    CONTINUE
         DO 50 I = J + 1, K
            W( I ) = W( I )*( Q( I, J ) / ( DLAMDA( I )-DLAMDA( J ) ) )
   50    CONTINUE
   60 CONTINUE
      DO 70 I = 1, K
         W( I ) = SIGN( SQRT( -W( I ) ), S( I ) )
   70 CONTINUE
!
!     Compute eigenvectors of the modified rank-1 modification.
!
      DO 100 J = 1, K
         DO 80 I = 1, K
            S( I ) = W( I ) / Q( I, J )
   80    CONTINUE
         TEMP = DNRM2( K, S, 1 )
         DO 90 I = 1, K
            II = INDX( I )
            Q( I, J ) = S( II ) / TEMP
   90    CONTINUE
  100 CONTINUE
!
!     Compute the updated eigenvectors.
!
  110 CONTINUE
!
      N2 = N - N1
      N12 = CTOT( 1 ) + CTOT( 2 )
      N23 = CTOT( 2 ) + CTOT( 3 )
!
      CALL DLACPY( 'A', N23, K, Q( CTOT( 1 )+1, 1 ), LDQ, S, N23 )
      IQ2 = N1*N12 + 1
      IF( N23.NE.0 ) THEN
         CALL DGEMM( 'N', 'N', N2, K, N23, ONE, Q2( IQ2 ), N2, S, N23, ZERO, Q( N1+1, 1 ), LDQ )
      ELSE
         CALL DLASET( 'A', N2, K, ZERO, ZERO, Q( N1+1, 1 ), LDQ )
      END IF
!
      CALL DLACPY( 'A', N12, K, Q, LDQ, S, N12 )
      IF( N12.NE.0 ) THEN
         CALL DGEMM( 'N', 'N', N1, K, N12, ONE, Q2, N1, S, N12, ZERO, Q, LDQ )
      ELSE
         CALL DLASET( 'A', N1, K, ZERO, ZERO, Q( 1, 1 ), LDQ )
      END IF
!
!
  120 CONTINUE
      RETURN
!
!     End of DLAED3
!
      END
!> \brief \b DLAED4
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DLAED4 + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaed4.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaed4.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaed4.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAED4( N, I, D, Z, DELTA, RHO, DLAM, INFO )
! 
!       .. Scalar Arguments ..
!       integer(i4)  ::          I, INFO, N
!       real(dp)    ::    DLAM, RHO
!       ..
!       .. Array Arguments ..
!       real(dp)    ::    D( * ), DELTA( * ), Z( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> This subroutine computes the I-th updated eigenvalue of a symmetric
!> rank-one modification to a diagonal matrix whose elements are
!> given in the array d, and that
!>
!>            D(i) < D(j)  for  i < j
!>
!> and that RHO > 0.  This is arranged by the calling routine, and is
!> no loss in generality.  The rank-one modified system is thus
!>
!>            diag( D )  +  RHO * Z * Z_transpose.
!>
!> where we assume the Euclidean norm of Z is 1.
!>
!> The method consists of approximating the rational functions in the
!> secular equation by simpler interpolating rational functions.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         The length of all arrays.
!> \endverbatim
!>
!> \param[in] I
!> \verbatim
!>          I is INTEGER
!>         The index of the eigenvalue to be computed.  1 <= I <= N.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is real(dp)    ::  array, dimension (N)
!>         The original eigenvalues.  It is assumed that they are in
!>         order, D(I) < D(J)  for I < J.
!> \endverbatim
!>
!> \param[in] Z
!> \verbatim
!>          Z is real(dp)    ::  array, dimension (N)
!>         The components of the updating vector.
!> \endverbatim
!>
!> \param[out] DELTA
!> \verbatim
!>          DELTA is real(dp)    ::  array, dimension (N)
!>         If N .GT. 2, DELTA contains (D(j) - lambda_I) in its  j-th
!>         component.  If N = 1, then DELTA(1) = 1. If N = 2, see DLAED5
!>         for detail. The vector DELTA contains the information necessary
!>         to construct the eigenvectors by DLAED3 and DLAED9.
!> \endverbatim
!>
!> \param[in] RHO
!> \verbatim
!>          RHO is real(dp)    :: 
!>         The scalar in the symmetric updating formula.
!> \endverbatim
!>
!> \param[out] DLAM
!> \verbatim
!>          DLAM is real(dp)    :: 
!>         The computed lambda_I, the I-th updated eigenvalue.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>         = 0:  successful exit
!>         > 0:  if INFO = 1, the updating process failed.
!> \endverbatim
!
!> \par Internal Parameters:
!  =========================
!>
!> \verbatim
!>  Logical variable ORGATI (origin-at-i?) is used for distinguishing
!>  whether D(i) or D(i+1) is treated as the origin.
!>
!>            ORGATI = .true.    origin at i
!>            ORGATI = .false.   origin at i+1
!>
!>   Logical variable SWTCH3 (switch-for-3-poles?) is for noting
!>   if we are working with THREE poles!
!>
!>   MAXIT is the maximum number of iterations allowed for each
!>   eigenvalue.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup auxOTHERcomputational
!
!> \par Contributors:
!  ==================
!>
!>     Ren-Cang Li, Computer Science Division, University of California
!>     at Berkeley, USA
!>
!  =====================================================================
      SUBROUTINE DLAED4( N, I, D, Z, DELTA, RHO, DLAM, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes
      integer(i4)  ::          I, INFO, N
      real(dp)    ::    DLAM, RHO
!     ..
!     .. Array Arguments ..
      real(dp)    ::    D( * ), DELTA( * ), Z( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      integer(i4)  ::          MAXIT
      PARAMETER          ( MAXIT = 30 )
      real(dp)    ::    ZERO, ONE, TWO, THREE, FOUR, EIGHT, TEN
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, THREE = 3.0D0, FOUR = 4.0D0, EIGHT = 8.0D0, TEN = 10.0D0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            ORGATI, SWTCH, SWTCH3
      integer(i4)  ::          II, IIM1, IIP1, IP1, ITER, J, NITER
      real(dp)    ::    A, B, C, DEL, DLTLB, DLTUB, DPHI, DPSI, DW, EPS, ERRETM, ETA, MIDPT, PHI, PREW, PSI, &
                         RHOINV, TAU, TEMP, TEMP1, W
!     ..
!     .. Local Arrays ..
      real(dp)    ::    ZZ( 3 )
!     ..
!     .. External Functions ..
      real(dp)    ::    DLAMCH
      EXTERNAL           DLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLAED5, DLAED6
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Since this routine is called in an inner loop, we do no argument
!     checking.
!
!     Quick return for N=1 and 2.
!
      INFO = 0
      IF( N.EQ.1 ) THEN
!
!         Presumably, I=1 upon entry
!
         DLAM = D( 1 ) + RHO*Z( 1 )*Z( 1 )
         DELTA( 1 ) = ONE
         RETURN
      END IF
      IF( N.EQ.2 ) THEN
         CALL DLAED5( I, D, Z, DELTA, RHO, DLAM )
         RETURN
      END IF
!
!     Compute machine epsilon
!
      EPS = DLAMCH( 'Epsilon' )
      RHOINV = ONE / RHO
!
!     The case I = N
!
      IF( I.EQ.N ) THEN
!
!        Initialize some basic variables
!
         II = N - 1
         NITER = 1
!
!        Calculate initial guess
!
         MIDPT = RHO / TWO
!
!        If ||Z||_2 is not one, then TEMP should be set to
!        RHO * ||Z||_2^2 / TWO
!
         DO 10 J = 1, N
            DELTA( J ) = ( D( J )-D( I ) ) - MIDPT
   10    CONTINUE
!
         PSI = ZERO
         DO 20 J = 1, N - 2
            PSI = PSI + Z( J )*Z( J ) / DELTA( J )
   20    CONTINUE
!
         C = RHOINV + PSI
         W = C + Z( II )*Z( II ) / DELTA( II ) + Z( N )*Z( N ) / DELTA( N )
!
         IF( W.LE.ZERO ) THEN
            TEMP = Z( N-1 )*Z( N-1 ) / ( D( N )-D( N-1 )+RHO ) + Z( N )*Z( N ) / RHO
            IF( C.LE.TEMP ) THEN
               TAU = RHO
            ELSE
               DEL = D( N ) - D( N-1 )
               A = -C*DEL + Z( N-1 )*Z( N-1 ) + Z( N )*Z( N )
               B = Z( N )*Z( N )*DEL
               IF( A.LT.ZERO ) THEN
                  TAU = TWO*B / ( SQRT( A*A+FOUR*B*C )-A )
               ELSE
                  TAU = ( A+SQRT( A*A+FOUR*B*C ) ) / ( TWO*C )
               END IF
            END IF
!
!           It can be proved that
!               D(N)+RHO/2 <= LAMBDA(N) < D(N)+TAU <= D(N)+RHO
!
            DLTLB = MIDPT
            DLTUB = RHO
         ELSE
            DEL = D( N ) - D( N-1 )
            A = -C*DEL + Z( N-1 )*Z( N-1 ) + Z( N )*Z( N )
            B = Z( N )*Z( N )*DEL
            IF( A.LT.ZERO ) THEN
               TAU = TWO*B / ( SQRT( A*A+FOUR*B*C )-A )
            ELSE
               TAU = ( A+SQRT( A*A+FOUR*B*C ) ) / ( TWO*C )
            END IF
!
!           It can be proved that
!               D(N) < D(N)+TAU < LAMBDA(N) < D(N)+RHO/2
!
            DLTLB = ZERO
            DLTUB = MIDPT
         END IF
!
         DO 30 J = 1, N
            DELTA( J ) = ( D( J )-D( I ) ) - TAU
   30    CONTINUE
!
!        Evaluate PSI and the derivative DPSI
!
         DPSI = ZERO
         PSI = ZERO
         ERRETM = ZERO
         DO 40 J = 1, II
            TEMP = Z( J ) / DELTA( J )
            PSI = PSI + Z( J )*TEMP
            DPSI = DPSI + TEMP*TEMP
            ERRETM = ERRETM + PSI
   40    CONTINUE
         ERRETM = ABS( ERRETM )
!
!        Evaluate PHI and the derivative DPHI
!
         TEMP = Z( N ) / DELTA( N )
         PHI = Z( N )*TEMP
         DPHI = TEMP*TEMP
         ERRETM = EIGHT*( -PHI-PSI ) + ERRETM - PHI + RHOINV + ABS( TAU )*( DPSI+DPHI )
!
         W = RHOINV + PHI + PSI
!
!        Test for convergence
!
         IF( ABS( W ).LE.EPS*ERRETM ) THEN
            DLAM = D( I ) + TAU
            GO TO 250
         END IF
!
         IF( W.LE.ZERO ) THEN
            DLTLB = MAX( DLTLB, TAU )
         ELSE
            DLTUB = MIN( DLTUB, TAU )
         END IF
!
!        Calculate the new step
!
         NITER = NITER + 1
         C = W - DELTA( N-1 )*DPSI - DELTA( N )*DPHI
         A = ( DELTA( N-1 )+DELTA( N ) )*W - DELTA( N-1 )*DELTA( N )*( DPSI+DPHI )
         B = DELTA( N-1 )*DELTA( N )*W
         IF( C.LT.ZERO ) C = ABS( C )
         IF( C.EQ.ZERO ) THEN
!      h   ETA = B/A
!           ETA = RHO - TAU
            ETA = DLTUB - TAU
         ELSE IF( A.GE.ZERO ) THEN
            ETA = ( A+SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
         ELSE
            ETA = TWO*B / ( A-SQRT( ABS( A*A-FOUR*B*C ) ) )
         END IF
!
!        Note, eta should be positive if w is negative, and
!        eta should be negative otherwise. However,
!        if for some reason caused by roundoff, eta*w > 0,
!        we simply use one Newton step instead. This way
!        will guarantee eta*w < 0.
!
         IF( W*ETA.GT.ZERO ) ETA = -W / ( DPSI+DPHI )
         TEMP = TAU + ETA
         IF( TEMP.GT.DLTUB .OR. TEMP.LT.DLTLB ) THEN
            IF( W.LT.ZERO ) THEN
               ETA = ( DLTUB-TAU ) / TWO
            ELSE
               ETA = ( DLTLB-TAU ) / TWO
            END IF
         END IF
         DO 50 J = 1, N
            DELTA( J ) = DELTA( J ) - ETA
   50    CONTINUE
!
         TAU = TAU + ETA
!
!        Evaluate PSI and the derivative DPSI
!
         DPSI = ZERO
         PSI = ZERO
         ERRETM = ZERO
         DO 60 J = 1, II
            TEMP = Z( J ) / DELTA( J )
            PSI = PSI + Z( J )*TEMP
            DPSI = DPSI + TEMP*TEMP
            ERRETM = ERRETM + PSI
   60    CONTINUE
         ERRETM = ABS( ERRETM )
!
!        Evaluate PHI and the derivative DPHI
!
         TEMP = Z( N ) / DELTA( N )
         PHI = Z( N )*TEMP
         DPHI = TEMP*TEMP
         ERRETM = EIGHT*( -PHI-PSI ) + ERRETM - PHI + RHOINV + ABS( TAU )*( DPSI+DPHI )
!
         W = RHOINV + PHI + PSI
!
!        Main loop to update the values of the array   DELTA
!
         ITER = NITER + 1
!
         DO 90 NITER = ITER, MAXIT
!
!           Test for convergence
!
            IF( ABS( W ).LE.EPS*ERRETM ) THEN
               DLAM = D( I ) + TAU
               GO TO 250
            END IF
!
            IF( W.LE.ZERO ) THEN
               DLTLB = MAX( DLTLB, TAU )
            ELSE
               DLTUB = MIN( DLTUB, TAU )
            END IF
!
!           Calculate the new step
!
            C = W - DELTA( N-1 )*DPSI - DELTA( N )*DPHI
            A = ( DELTA( N-1 )+DELTA( N ) )*W - DELTA( N-1 )*DELTA( N )*( DPSI+DPHI )
            B = DELTA( N-1 )*DELTA( N )*W
            IF( A.GE.ZERO ) THEN
               ETA = ( A+SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
            ELSE
               ETA = TWO*B / ( A-SQRT( ABS( A*A-FOUR*B*C ) ) )
            END IF
!
!           Note, eta should be positive if w is negative, and
!           eta should be negative otherwise. However,
!           if for some reason caused by roundoff, eta*w > 0,
!           we simply use one Newton step instead. This way
!           will guarantee eta*w < 0.
!
            IF( W*ETA.GT.ZERO ) ETA = -W / ( DPSI+DPHI )
            TEMP = TAU + ETA
            IF( TEMP.GT.DLTUB .OR. TEMP.LT.DLTLB ) THEN
               IF( W.LT.ZERO ) THEN
                  ETA = ( DLTUB-TAU ) / TWO
               ELSE
                  ETA = ( DLTLB-TAU ) / TWO
               END IF
            END IF
            DO 70 J = 1, N
               DELTA( J ) = DELTA( J ) - ETA
   70       CONTINUE
!
            TAU = TAU + ETA
!
!           Evaluate PSI and the derivative DPSI
!
            DPSI = ZERO
            PSI = ZERO
            ERRETM = ZERO
            DO 80 J = 1, II
               TEMP = Z( J ) / DELTA( J )
               PSI = PSI + Z( J )*TEMP
               DPSI = DPSI + TEMP*TEMP
               ERRETM = ERRETM + PSI
   80       CONTINUE
            ERRETM = ABS( ERRETM )
!
!           Evaluate PHI and the derivative DPHI
!
            TEMP = Z( N ) / DELTA( N )
            PHI = Z( N )*TEMP
            DPHI = TEMP*TEMP
            ERRETM = EIGHT*( -PHI-PSI ) + ERRETM - PHI + RHOINV + ABS( TAU )*( DPSI+DPHI )
!
            W = RHOINV + PHI + PSI
   90    CONTINUE
!
!        Return with INFO = 1, NITER = MAXIT and not converged
!
         INFO = 1
         DLAM = D( I ) + TAU
         GO TO 250
!
!        End for the case I = N
!
      ELSE
!
!        The case for I < N
!
         NITER = 1
         IP1 = I + 1
!
!        Calculate initial guess
!
         DEL = D( IP1 ) - D( I )
         MIDPT = DEL / TWO
         DO 100 J = 1, N
            DELTA( J ) = ( D( J )-D( I ) ) - MIDPT
  100    CONTINUE
!
         PSI = ZERO
         DO 110 J = 1, I - 1
            PSI = PSI + Z( J )*Z( J ) / DELTA( J )
  110    CONTINUE
!
         PHI = ZERO
         DO 120 J = N, I + 2, -1
            PHI = PHI + Z( J )*Z( J ) / DELTA( J )
  120    CONTINUE
         C = RHOINV + PSI + PHI
         W = C + Z( I )*Z( I ) / DELTA( I ) + Z( IP1 )*Z( IP1 ) / DELTA( IP1 )
!
         IF( W.GT.ZERO ) THEN
!
!           d(i)< the ith eigenvalue < (d(i)+d(i+1))/2
!
!           We choose d(i) as origin.
!
            ORGATI = .TRUE.
            A = C*DEL + Z( I )*Z( I ) + Z( IP1 )*Z( IP1 )
            B = Z( I )*Z( I )*DEL
            IF( A.GT.ZERO ) THEN
               TAU = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
            ELSE
               TAU = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
            END IF
            DLTLB = ZERO
            DLTUB = MIDPT
         ELSE
!
!           (d(i)+d(i+1))/2 <= the ith eigenvalue < d(i+1)
!
!           We choose d(i+1) as origin.
!
            ORGATI = .FALSE.
            A = C*DEL - Z( I )*Z( I ) - Z( IP1 )*Z( IP1 )
            B = Z( IP1 )*Z( IP1 )*DEL
            IF( A.LT.ZERO ) THEN
               TAU = TWO*B / ( A-SQRT( ABS( A*A+FOUR*B*C ) ) )
            ELSE
               TAU = -( A+SQRT( ABS( A*A+FOUR*B*C ) ) ) / ( TWO*C )
            END IF
            DLTLB = -MIDPT
            DLTUB = ZERO
         END IF
!
         IF( ORGATI ) THEN
            DO 130 J = 1, N
               DELTA( J ) = ( D( J )-D( I ) ) - TAU
  130       CONTINUE
         ELSE
            DO 140 J = 1, N
               DELTA( J ) = ( D( J )-D( IP1 ) ) - TAU
  140       CONTINUE
         END IF
         IF( ORGATI ) THEN
            II = I
         ELSE
            II = I + 1
         END IF
         IIM1 = II - 1
         IIP1 = II + 1
!
!        Evaluate PSI and the derivative DPSI
!
         DPSI = ZERO
         PSI = ZERO
         ERRETM = ZERO
         DO 150 J = 1, IIM1
            TEMP = Z( J ) / DELTA( J )
            PSI = PSI + Z( J )*TEMP
            DPSI = DPSI + TEMP*TEMP
            ERRETM = ERRETM + PSI
  150    CONTINUE
         ERRETM = ABS( ERRETM )
!
!        Evaluate PHI and the derivative DPHI
!
         DPHI = ZERO
         PHI = ZERO
         DO 160 J = N, IIP1, -1
            TEMP = Z( J ) / DELTA( J )
            PHI = PHI + Z( J )*TEMP
            DPHI = DPHI + TEMP*TEMP
            ERRETM = ERRETM + PHI
  160    CONTINUE
!
         W = RHOINV + PHI + PSI
!
!        W is the value of the secular function with
!        its ii-th element removed.
!
         SWTCH3 = .FALSE.
         IF( ORGATI ) THEN
            IF( W.LT.ZERO ) SWTCH3 = .TRUE.
         ELSE
            IF( W.GT.ZERO ) SWTCH3 = .TRUE.
         END IF
         IF( II.EQ.1 .OR. II.EQ.N ) SWTCH3 = .FALSE.
!
         TEMP = Z( II ) / DELTA( II )
         DW = DPSI + DPHI + TEMP*TEMP
         TEMP = Z( II )*TEMP
         W = W + TEMP
         ERRETM = EIGHT*( PHI-PSI ) + ERRETM + TWO*RHOINV + THREE*ABS( TEMP ) + ABS( TAU )*DW
!
!        Test for convergence
!
         IF( ABS( W ).LE.EPS*ERRETM ) THEN
            IF( ORGATI ) THEN
               DLAM = D( I ) + TAU
            ELSE
               DLAM = D( IP1 ) + TAU
            END IF
            GO TO 250
         END IF
!
         IF( W.LE.ZERO ) THEN
            DLTLB = MAX( DLTLB, TAU )
         ELSE
            DLTUB = MIN( DLTUB, TAU )
         END IF
!
!        Calculate the new step
!
         NITER = NITER + 1
         IF( .NOT.SWTCH3 ) THEN
            IF( ORGATI ) THEN
               C = W - DELTA( IP1 )*DW - ( D( I )-D( IP1 ) )* ( Z( I ) / DELTA( I ) )**2
            ELSE
               C = W - DELTA( I )*DW - ( D( IP1 )-D( I ) )* ( Z( IP1 ) / DELTA( IP1 ) )**2
            END IF
            A = ( DELTA( I )+DELTA( IP1 ) )*W - DELTA( I )*DELTA( IP1 )*DW
            B = DELTA( I )*DELTA( IP1 )*W
            IF( C.EQ.ZERO ) THEN
               IF( A.EQ.ZERO ) THEN
                  IF( ORGATI ) THEN
                     A = Z( I )*Z( I ) + DELTA( IP1 )*DELTA( IP1 )* ( DPSI+DPHI )
                  ELSE
                     A = Z( IP1 )*Z( IP1 ) + DELTA( I )*DELTA( I )* ( DPSI+DPHI )
                  END IF
               END IF
               ETA = B / A
            ELSE IF( A.LE.ZERO ) THEN
               ETA = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
            ELSE
               ETA = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
            END IF
         ELSE
!
!           Interpolation using THREE most relevant poles
!
            TEMP = RHOINV + PSI + PHI
            IF( ORGATI ) THEN
               TEMP1 = Z( IIM1 ) / DELTA( IIM1 )
               TEMP1 = TEMP1*TEMP1
               C = TEMP - DELTA( IIP1 )*( DPSI+DPHI ) - ( D( IIM1 )-D( IIP1 ) )*TEMP1
               ZZ( 1 ) = Z( IIM1 )*Z( IIM1 )
               ZZ( 3 ) = DELTA( IIP1 )*DELTA( IIP1 )* ( ( DPSI-TEMP1 )+DPHI )
            ELSE
               TEMP1 = Z( IIP1 ) / DELTA( IIP1 )
               TEMP1 = TEMP1*TEMP1
               C = TEMP - DELTA( IIM1 )*( DPSI+DPHI ) - ( D( IIP1 )-D( IIM1 ) )*TEMP1
               ZZ( 1 ) = DELTA( IIM1 )*DELTA( IIM1 )* ( DPSI+( DPHI-TEMP1 ) )
               ZZ( 3 ) = Z( IIP1 )*Z( IIP1 )
            END IF
            ZZ( 2 ) = Z( II )*Z( II )
            CALL DLAED6( NITER, ORGATI, C, DELTA( IIM1 ), ZZ, W, ETA, INFO )
            IF( INFO.NE.0 ) GO TO 250
         END IF
!
!        Note, eta should be positive if w is negative, and
!        eta should be negative otherwise. However,
!        if for some reason caused by roundoff, eta*w > 0,
!        we simply use one Newton step instead. This way
!        will guarantee eta*w < 0.
!
         IF( W*ETA.GE.ZERO ) ETA = -W / DW
         TEMP = TAU + ETA
         IF( TEMP.GT.DLTUB .OR. TEMP.LT.DLTLB ) THEN
            IF( W.LT.ZERO ) THEN
               ETA = ( DLTUB-TAU ) / TWO
            ELSE
               ETA = ( DLTLB-TAU ) / TWO
            END IF
         END IF
!
         PREW = W
!
         DO 180 J = 1, N
            DELTA( J ) = DELTA( J ) - ETA
  180    CONTINUE
!
!        Evaluate PSI and the derivative DPSI
!
         DPSI = ZERO
         PSI = ZERO
         ERRETM = ZERO
         DO 190 J = 1, IIM1
            TEMP = Z( J ) / DELTA( J )
            PSI = PSI + Z( J )*TEMP
            DPSI = DPSI + TEMP*TEMP
            ERRETM = ERRETM + PSI
  190    CONTINUE
         ERRETM = ABS( ERRETM )
!
!        Evaluate PHI and the derivative DPHI
!
         DPHI = ZERO
         PHI = ZERO
         DO 200 J = N, IIP1, -1
            TEMP = Z( J ) / DELTA( J )
            PHI = PHI + Z( J )*TEMP
            DPHI = DPHI + TEMP*TEMP
            ERRETM = ERRETM + PHI
  200    CONTINUE
!
         TEMP = Z( II ) / DELTA( II )
         DW = DPSI + DPHI + TEMP*TEMP
         TEMP = Z( II )*TEMP
         W = RHOINV + PHI + PSI + TEMP
         ERRETM = EIGHT*( PHI-PSI ) + ERRETM + TWO*RHOINV + THREE*ABS( TEMP ) + ABS( TAU+ETA )*DW
!
         SWTCH = .FALSE.
         IF( ORGATI ) THEN
            IF( -W.GT.ABS( PREW ) / TEN ) SWTCH = .TRUE.
         ELSE
            IF( W.GT.ABS( PREW ) / TEN ) SWTCH = .TRUE.
         END IF
!
         TAU = TAU + ETA
!
!        Main loop to update the values of the array   DELTA
!
         ITER = NITER + 1
!
         DO 240 NITER = ITER, MAXIT
!
!           Test for convergence
!
            IF( ABS( W ).LE.EPS*ERRETM ) THEN
               IF( ORGATI ) THEN
                  DLAM = D( I ) + TAU
               ELSE
                  DLAM = D( IP1 ) + TAU
               END IF
               GO TO 250
            END IF
!
            IF( W.LE.ZERO ) THEN
               DLTLB = MAX( DLTLB, TAU )
            ELSE
               DLTUB = MIN( DLTUB, TAU )
            END IF
!
!           Calculate the new step
!
            IF( .NOT.SWTCH3 ) THEN
               IF( .NOT.SWTCH ) THEN
                  IF( ORGATI ) THEN
                     C = W - DELTA( IP1 )*DW - ( D( I )-D( IP1 ) )*( Z( I ) / DELTA( I ) )**2
                  ELSE
                     C = W - DELTA( I )*DW - ( D( IP1 )-D( I ) )* ( Z( IP1 ) / DELTA( IP1 ) )**2
                  END IF
               ELSE
                  TEMP = Z( II ) / DELTA( II )
                  IF( ORGATI ) THEN
                     DPSI = DPSI + TEMP*TEMP
                  ELSE
                     DPHI = DPHI + TEMP*TEMP
                  END IF
                  C = W - DELTA( I )*DPSI - DELTA( IP1 )*DPHI
               END IF
               A = ( DELTA( I )+DELTA( IP1 ) )*W - DELTA( I )*DELTA( IP1 )*DW
               B = DELTA( I )*DELTA( IP1 )*W
               IF( C.EQ.ZERO ) THEN
                  IF( A.EQ.ZERO ) THEN
                     IF( .NOT.SWTCH ) THEN
                        IF( ORGATI ) THEN
                           A = Z( I )*Z( I ) + DELTA( IP1 )* DELTA( IP1 )*( DPSI+DPHI )
                        ELSE
                           A = Z( IP1 )*Z( IP1 ) + DELTA( I )*DELTA( I )*( DPSI+DPHI )
                        END IF
                     ELSE
                        A = DELTA( I )*DELTA( I )*DPSI + DELTA( IP1 )*DELTA( IP1 )*DPHI
                     END IF
                  END IF
                  ETA = B / A
               ELSE IF( A.LE.ZERO ) THEN
                  ETA = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
               ELSE
                  ETA = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
               END IF
            ELSE
!
!              Interpolation using THREE most relevant poles
!
               TEMP = RHOINV + PSI + PHI
               IF( SWTCH ) THEN
                  C = TEMP - DELTA( IIM1 )*DPSI - DELTA( IIP1 )*DPHI
                  ZZ( 1 ) = DELTA( IIM1 )*DELTA( IIM1 )*DPSI
                  ZZ( 3 ) = DELTA( IIP1 )*DELTA( IIP1 )*DPHI
               ELSE
                  IF( ORGATI ) THEN
                     TEMP1 = Z( IIM1 ) / DELTA( IIM1 )
                     TEMP1 = TEMP1*TEMP1
                     C = TEMP - DELTA( IIP1 )*( DPSI+DPHI ) - ( D( IIM1 )-D( IIP1 ) )*TEMP1
                     ZZ( 1 ) = Z( IIM1 )*Z( IIM1 )
                     ZZ( 3 ) = DELTA( IIP1 )*DELTA( IIP1 )* ( ( DPSI-TEMP1 )+DPHI )
                  ELSE
                     TEMP1 = Z( IIP1 ) / DELTA( IIP1 )
                     TEMP1 = TEMP1*TEMP1
                     C = TEMP - DELTA( IIM1 )*( DPSI+DPHI ) - ( D( IIP1 )-D( IIM1 ) )*TEMP1
                     ZZ( 1 ) = DELTA( IIM1 )*DELTA( IIM1 )* ( DPSI+( DPHI-TEMP1 ) )
                     ZZ( 3 ) = Z( IIP1 )*Z( IIP1 )
                  END IF
               END IF
               CALL DLAED6( NITER, ORGATI, C, DELTA( IIM1 ), ZZ, W, ETA, INFO )
               IF( INFO.NE.0 ) GO TO 250
            END IF
!
!           Note, eta should be positive if w is negative, and
!           eta should be negative otherwise. However,
!           if for some reason caused by roundoff, eta*w > 0,
!           we simply use one Newton step instead. This way
!           will guarantee eta*w < 0.
!
            IF( W*ETA.GE.ZERO ) ETA = -W / DW
            TEMP = TAU + ETA
            IF( TEMP.GT.DLTUB .OR. TEMP.LT.DLTLB ) THEN
               IF( W.LT.ZERO ) THEN
                  ETA = ( DLTUB-TAU ) / TWO
               ELSE
                  ETA = ( DLTLB-TAU ) / TWO
               END IF
            END IF
!
            DO 210 J = 1, N
               DELTA( J ) = DELTA( J ) - ETA
  210       CONTINUE
!
            TAU = TAU + ETA
            PREW = W
!
!           Evaluate PSI and the derivative DPSI
!
            DPSI = ZERO
            PSI = ZERO
            ERRETM = ZERO
            DO 220 J = 1, IIM1
               TEMP = Z( J ) / DELTA( J )
               PSI = PSI + Z( J )*TEMP
               DPSI = DPSI + TEMP*TEMP
               ERRETM = ERRETM + PSI
  220       CONTINUE
            ERRETM = ABS( ERRETM )
!
!           Evaluate PHI and the derivative DPHI
!
            DPHI = ZERO
            PHI = ZERO
            DO 230 J = N, IIP1, -1
               TEMP = Z( J ) / DELTA( J )
               PHI = PHI + Z( J )*TEMP
               DPHI = DPHI + TEMP*TEMP
               ERRETM = ERRETM + PHI
  230       CONTINUE
!
            TEMP = Z( II ) / DELTA( II )
            DW = DPSI + DPHI + TEMP*TEMP
            TEMP = Z( II )*TEMP
            W = RHOINV + PHI + PSI + TEMP
            ERRETM = EIGHT*( PHI-PSI ) + ERRETM + TWO*RHOINV + THREE*ABS( TEMP ) + ABS( TAU )*DW
            IF( W*PREW.GT.ZERO .AND. ABS( W ).GT.ABS( PREW ) / TEN ) SWTCH = .NOT.SWTCH
!
  240    CONTINUE
!
!        Return with INFO = 1, NITER = MAXIT and not converged
!
         INFO = 1
         IF( ORGATI ) THEN
            DLAM = D( I ) + TAU
         ELSE
            DLAM = D( IP1 ) + TAU
         END IF
!
      END IF
!
  250 CONTINUE
!
      RETURN
!
!     End of DLAED4
!
      END
!> \brief \b DLAED8
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DLAED8 + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaed8.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaed8.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaed8.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAED8( ICOMPQ, K, N, QSIZ, D, Q, LDQ, INDXQ, RHO,
!                          CUTPNT, Z, DLAMDA, Q2, LDQ2, W, PERM, GIVPTR,
!                          GIVCOL, GIVNUM, INDXP, INDX, INFO )
! 
!       .. Scalar Arguments ..
!       integer(i4)  ::          CUTPNT, GIVPTR, ICOMPQ, INFO, K, LDQ, LDQ2, N,
!      $                   QSIZ
!       real(dp)    ::    RHO
!       ..
!       .. Array Arguments ..
!       integer(i4)  ::          GIVCOL( 2, * ), INDX( * ), INDXP( * ),
!      $                   INDXQ( * ), PERM( * )
!       real(dp)    ::    D( * ), DLAMDA( * ), GIVNUM( 2, * ),
!      $                   Q( LDQ, * ), Q2( LDQ2, * ), W( * ), Z( * )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAED8 merges the two sets of eigenvalues together into a single
!> sorted set.  Then it tries to deflate the size of the problem.
!> There are two ways in which deflation can occur:  when two or more
!> eigenvalues are close together or if there is a tiny element in the
!> Z vector.  For each such occurrence the order of the related secular
!> equation problem is reduced by one.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ICOMPQ
!> \verbatim
!>          ICOMPQ is INTEGER
!>          = 0:  Compute eigenvalues only.
!>          = 1:  Compute eigenvectors of original dense symmetric matrix
!>                also.  On entry, Q contains the orthogonal matrix used
!>                to reduce the original matrix to tridiagonal form.
!> \endverbatim
!>
!> \param[out] K
!> \verbatim
!>          K is INTEGER
!>         The number of non-deflated eigenvalues, and the order of the
!>         related secular equation.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         The dimension of the symmetric tridiagonal matrix.  N >= 0.
!> \endverbatim
!>
!> \param[in] QSIZ
!> \verbatim
!>          QSIZ is INTEGER
!>         The dimension of the orthogonal matrix used to reduce
!>         the full matrix to tridiagonal form.  QSIZ >= N if ICOMPQ = 1.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is real(dp)    ::  array, dimension (N)
!>         On entry, the eigenvalues of the two submatrices to be
!>         combined.  On exit, the trailing (N-K) updated eigenvalues
!>         (those which were deflated) sorted into increasing order.
!> \endverbatim
!>
!> \param[in,out] Q
!> \verbatim
!>          Q is real(dp)    ::  array, dimension (LDQ,N)
!>         If ICOMPQ = 0, Q is not referenced.  Otherwise,
!>         on entry, Q contains the eigenvectors of the partially solved
!>         system which has been previously updated in matrix
!>         multiplies with other partially solved eigensystems.
!>         On exit, Q contains the trailing (N-K) updated eigenvectors
!>         (those which were deflated) in its last N-K columns.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>         The leading dimension of the array Q.  LDQ >= max(1,N).
!> \endverbatim
!>
!> \param[in] INDXQ
!> \verbatim
!>          INDXQ is INTEGER array, dimension (N)
!>         The permutation which separately sorts the two sub-problems
!>         in D into ascending order.  Note that elements in the second
!>         half of this permutation must first have CUTPNT added to
!>         their values in order to be accurate.
!> \endverbatim
!>
!> \param[in,out] RHO
!> \verbatim
!>          RHO is real(dp)    :: 
!>         On entry, the off-diagonal element associated with the rank-1
!>         cut which originally split the two submatrices which are now
!>         being recombined.
!>         On exit, RHO has been modified to the value required by
!>         DLAED3.
!> \endverbatim
!>
!> \param[in] CUTPNT
!> \verbatim
!>          CUTPNT is INTEGER
!>         The location of the last eigenvalue in the leading
!>         sub-matrix.  min(1,N) <= CUTPNT <= N.
!> \endverbatim
!>
!> \param[in] Z
!> \verbatim
!>          Z is real(dp)    ::  array, dimension (N)
!>         On entry, Z contains the updating vector (the last row of
!>         the first sub-eigenvector matrix and the first row of the
!>         second sub-eigenvector matrix).
!>         On exit, the contents of Z are destroyed by the updating
!>         process.
!> \endverbatim
!>
!> \param[out] DLAMDA
!> \verbatim
!>          DLAMDA is real(dp)    ::  array, dimension (N)
!>         A copy of the first K eigenvalues which will be used by
!>         DLAED3 to form the secular equation.
!> \endverbatim
!>
!> \param[out] Q2
!> \verbatim
!>          Q2 is real(dp)    ::  array, dimension (LDQ2,N)
!>         If ICOMPQ = 0, Q2 is not referenced.  Otherwise,
!>         a copy of the first K eigenvectors which will be used by
!>         DLAED7 in a matrix multiply (DGEMM) to update the new
!>         eigenvectors.
!> \endverbatim
!>
!> \param[in] LDQ2
!> \verbatim
!>          LDQ2 is INTEGER
!>         The leading dimension of the array Q2.  LDQ2 >= max(1,N).
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is real(dp)    ::  array, dimension (N)
!>         The first k values of the final deflation-altered z-vector and
!>         will be passed to DLAED3.
!> \endverbatim
!>
!> \param[out] PERM
!> \verbatim
!>          PERM is INTEGER array, dimension (N)
!>         The permutations (from deflation and sorting) to be applied
!>         to each eigenblock.
!> \endverbatim
!>
!> \param[out] GIVPTR
!> \verbatim
!>          GIVPTR is INTEGER
!>         The number of Givens rotations which took place in this
!>         subproblem.
!> \endverbatim
!>
!> \param[out] GIVCOL
!> \verbatim
!>          GIVCOL is INTEGER array, dimension (2, N)
!>         Each pair of numbers indicates a pair of columns to take place
!>         in a Givens rotation.
!> \endverbatim
!>
!> \param[out] GIVNUM
!> \verbatim
!>          GIVNUM is real(dp)    ::  array, dimension (2, N)
!>         Each number indicates the S value to be used in the
!>         corresponding Givens rotation.
!> \endverbatim
!>
!> \param[out] INDXP
!> \verbatim
!>          INDXP is INTEGER array, dimension (N)
!>         The permutation used to place deflated values of D at the end
!>         of the array.  INDXP(1:K) points to the nondeflated D-values
!>         and INDXP(K+1:N) points to the deflated eigenvalues.
!> \endverbatim
!>
!> \param[out] INDX
!> \verbatim
!>          INDX is INTEGER array, dimension (N)
!>         The permutation used to sort the contents of D into ascending
!>         order.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup auxOTHERcomputational
!
!> \par Contributors:
!  ==================
!>
!> Jeff Rutter, Computer Science Division, University of California
!> at Berkeley, USA
!
!  =====================================================================
      SUBROUTINE DLAED8( ICOMPQ, K, N, QSIZ, D, Q, LDQ, INDXQ, RHO, &
                         CUTPNT, Z, DLAMDA, Q2, LDQ2, W, PERM, GIVPTR, &
                         GIVCOL, GIVNUM, INDXP, INDX, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes
      integer(i4)  ::          CUTPNT, GIVPTR, ICOMPQ, INFO, K, LDQ, LDQ2, N, QSIZ
      real(dp)    ::    RHO
!     ..
!     .. Array Arguments ..
      integer(i4)  ::          GIVCOL( 2, * ), INDX( * ), INDXP( * ), INDXQ( * ), PERM( * )
      real(dp)    ::    D( * ), DLAMDA( * ), GIVNUM( 2, * ), Q( LDQ, * ), Q2( LDQ2, * ), W( * ), Z( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      real(dp)    ::    MONE, ZERO, ONE, TWO, EIGHT
      PARAMETER          ( MONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, EIGHT = 8.0D0 )
!     ..
!     .. Local Scalars ..
!
      integer(i4)  ::          I, IMAX, J, JLAM, JMAX, JP, K2, N1, N1P1, N2
      real(dp)    ::    C, EPS, S, T, TAU, TOL
!     ..
!     .. External Functions ..
      integer(i4)  ::          IDAMAX
      real(dp)    ::    DLAMCH, DLAPY2
      EXTERNAL           IDAMAX, DLAMCH, DLAPY2
!     ..
!     .. External Subroutines ..
      EXTERNAL           DCOPY, DLACPY, DLAMRG, DROT, DSCAL, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
!
      IF( ICOMPQ.LT.0 .OR. ICOMPQ.GT.1 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( ICOMPQ.EQ.1 .AND. QSIZ.LT.N ) THEN
         INFO = -4
      ELSE IF( LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( CUTPNT.LT.MIN( 1, N ) .OR. CUTPNT.GT.N ) THEN
         INFO = -10
      ELSE IF( LDQ2.LT.MAX( 1, N ) ) THEN
         INFO = -14
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLAED8', -INFO )
         RETURN
      END IF
!
!     Need to initialize GIVPTR to O here in case of quick exit
!     to prevent an unspecified code behavior (usually sigfault) 
!     when IWORK array on entry to *stedc is not zeroed 
!     (or at least some IWORK entries which used in *laed7 for GIVPTR).
!
      GIVPTR = 0
!
!     Quick return if possible
!
      IF( N.EQ.0 ) RETURN
!
      N1 = CUTPNT
      N2 = N - N1
      N1P1 = N1 + 1
!
      IF( RHO.LT.ZERO ) THEN
         CALL DSCAL( N2, MONE, Z( N1P1 ), 1 )
      END IF
!
!     Normalize z so that norm(z) = 1
!
      T = ONE / SQRT( TWO )
      DO 10 J = 1, N
         INDX( J ) = J
   10 CONTINUE
      CALL DSCAL( N, T, Z, 1 )
      RHO = ABS( TWO*RHO )
!
!     Sort the eigenvalues into increasing order
!
      DO 20 I = CUTPNT + 1, N
         INDXQ( I ) = INDXQ( I ) + CUTPNT
   20 CONTINUE
      DO 30 I = 1, N
         DLAMDA( I ) = D( INDXQ( I ) )
         W( I ) = Z( INDXQ( I ) )
   30 CONTINUE
      I = 1
      J = CUTPNT + 1
      CALL DLAMRG( N1, N2, DLAMDA, 1, 1, INDX )
      DO 40 I = 1, N
         D( I ) = DLAMDA( INDX( I ) )
         Z( I ) = W( INDX( I ) )
   40 CONTINUE
!
!     Calculate the allowable deflation tolerence
!
      IMAX = IDAMAX( N, Z, 1 )
      JMAX = IDAMAX( N, D, 1 )
      EPS = DLAMCH( 'Epsilon' )
      TOL = EIGHT*EPS*ABS( D( JMAX ) )
!
!     If the rank-1 modifier is small enough, no more needs to be done
!     except to reorganize Q so that its columns correspond with the
!     elements in D.
!
      IF( RHO*ABS( Z( IMAX ) ).LE.TOL ) THEN
         K = 0
         IF( ICOMPQ.EQ.0 ) THEN
            DO 50 J = 1, N
               PERM( J ) = INDXQ( INDX( J ) )
   50       CONTINUE
         ELSE
            DO 60 J = 1, N
               PERM( J ) = INDXQ( INDX( J ) )
               CALL DCOPY( QSIZ, Q( 1, PERM( J ) ), 1, Q2( 1, J ), 1 )
   60       CONTINUE
            CALL DLACPY( 'A', QSIZ, N, Q2( 1, 1 ), LDQ2, Q( 1, 1 ), LDQ )
         END IF
         RETURN
      END IF
!
!     If there are multiple eigenvalues then the problem deflates.  Here
!     the number of equal eigenvalues are found.  As each equal
!     eigenvalue is found, an elementary reflector is computed to rotate
!     the corresponding eigensubspace so that the corresponding
!     components of Z are zero in this new basis.
!
      K = 0
      K2 = N + 1
      DO 70 J = 1, N
         IF( RHO*ABS( Z( J ) ).LE.TOL ) THEN
!
!           Deflate due to small z component.
!
            K2 = K2 - 1
            INDXP( K2 ) = J
            IF( J.EQ.N ) GO TO 110
         ELSE
            JLAM = J
            GO TO 80
         END IF
   70 CONTINUE
   80 CONTINUE
      J = J + 1
      IF( J.GT.N ) GO TO 100
      IF( RHO*ABS( Z( J ) ).LE.TOL ) THEN
!
!        Deflate due to small z component.
!
         K2 = K2 - 1
         INDXP( K2 ) = J
      ELSE
!
!        Check if eigenvalues are close enough to allow deflation.
!
         S = Z( JLAM )
         C = Z( J )
!
!        Find sqrt(a**2+b**2) without overflow or
!        destructive underflow.
!
         TAU = DLAPY2( C, S )
         T = D( J ) - D( JLAM )
         C = C / TAU
         S = -S / TAU
         IF( ABS( T*C*S ).LE.TOL ) THEN
!
!           Deflation is possible.
!
            Z( J ) = TAU
            Z( JLAM ) = ZERO
!
!           Record the appropriate Givens rotation
!
            GIVPTR = GIVPTR + 1
            GIVCOL( 1, GIVPTR ) = INDXQ( INDX( JLAM ) )
            GIVCOL( 2, GIVPTR ) = INDXQ( INDX( J ) )
            GIVNUM( 1, GIVPTR ) = C
            GIVNUM( 2, GIVPTR ) = S
            IF( ICOMPQ.EQ.1 ) THEN
               CALL DROT( QSIZ, Q( 1, INDXQ( INDX( JLAM ) ) ), 1, Q( 1, INDXQ( INDX( J ) ) ), 1, C, S )
            END IF
            T = D( JLAM )*C*C + D( J )*S*S
            D( J ) = D( JLAM )*S*S + D( J )*C*C
            D( JLAM ) = T
            K2 = K2 - 1
            I = 1
   90       CONTINUE
            IF( K2+I.LE.N ) THEN
               IF( D( JLAM ).LT.D( INDXP( K2+I ) ) ) THEN
                  INDXP( K2+I-1 ) = INDXP( K2+I )
                  INDXP( K2+I ) = JLAM
                  I = I + 1
                  GO TO 90
               ELSE
                  INDXP( K2+I-1 ) = JLAM
               END IF
            ELSE
               INDXP( K2+I-1 ) = JLAM
            END IF
            JLAM = J
         ELSE
            K = K + 1
            W( K ) = Z( JLAM )
            DLAMDA( K ) = D( JLAM )
            INDXP( K ) = JLAM
            JLAM = J
         END IF
      END IF
      GO TO 80
  100 CONTINUE
!
!     Record the last eigenvalue.
!
      K = K + 1
      W( K ) = Z( JLAM )
      DLAMDA( K ) = D( JLAM )
      INDXP( K ) = JLAM
!
  110 CONTINUE
!
!     Sort the eigenvalues and corresponding eigenvectors into DLAMDA
!     and Q2 respectively.  The eigenvalues/vectors which were not
!     deflated go into the first K slots of DLAMDA and Q2 respectively,
!     while those which were deflated go into the last N - K slots.
!
      IF( ICOMPQ.EQ.0 ) THEN
         DO 120 J = 1, N
            JP = INDXP( J )
            DLAMDA( J ) = D( JP )
            PERM( J ) = INDXQ( INDX( JP ) )
  120    CONTINUE
      ELSE
         DO 130 J = 1, N
            JP = INDXP( J )
            DLAMDA( J ) = D( JP )
            PERM( J ) = INDXQ( INDX( JP ) )
            CALL DCOPY( QSIZ, Q( 1, PERM( J ) ), 1, Q2( 1, J ), 1 )
  130    CONTINUE
      END IF
!
!     The deflated eigenvalues and their corresponding vectors go back
!     into the last N - K slots of D and Q respectively.
!
      IF( K.LT.N ) THEN
         IF( ICOMPQ.EQ.0 ) THEN
            CALL DCOPY( N-K, DLAMDA( K+1 ), 1, D( K+1 ), 1 )
         ELSE
            CALL DCOPY( N-K, DLAMDA( K+1 ), 1, D( K+1 ), 1 )
            CALL DLACPY( 'A', QSIZ, N-K, Q2( 1, K+1 ), LDQ2, Q( 1, K+1 ), LDQ )
         END IF
      END IF
!
      RETURN
!
!     End of DLAED8
!
      END
!> \brief \b DLAED5
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DLAED5 + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaed5.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaed5.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaed5.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAED5( I, D, Z, DELTA, RHO, DLAM )
! 
!       .. Scalar Arguments ..
!       integer(i4)  ::            I
!       real(dp)    ::   DLAM, RHO
!       ..
!       .. Array Arguments ..
!       real(dp)    ::   D( 2 ), DELTA( 2 ), Z( 2 )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> This subroutine computes the I-th eigenvalue of a symmetric rank-one
!> modification of a 2-by-2 diagonal matrix
!>
!>            diag( D )  +  RHO * Z * transpose(Z) .
!>
!> The diagonal elements in the array D are assumed to satisfy
!>
!>            D(i) < D(j)  for  i < j .
!>
!> We also assume RHO > 0 and that the Euclidean norm of the vector
!> Z is one.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] I
!> \verbatim
!>          I is integer(i4)  ::
!>         The index of the eigenvalue to be computed.  I = 1 or I = 2.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is real(dp)    :: array, dimension (2)
!>         The original eigenvalues.  We assume D(1) < D(2).
!> \endverbatim
!>
!> \param[in] Z
!> \verbatim
!>          Z is real(dp)    :: array, dimension (2)
!>         The components of the updating vector.
!> \endverbatim
!>
!> \param[out] DELTA
!> \verbatim
!>          DELTA is real(dp)    :: array, dimension (2)
!>         The vector DELTA contains the information necessary
!>         to construct the eigenvectors.
!> \endverbatim
!>
!> \param[in] RHO
!> \verbatim
!>          RHO is real(dp)    ::
!>         The scalar in the symmetric updating formula.
!> \endverbatim
!>
!> \param[out] DLAM
!> \verbatim
!>          DLAM is real(dp)    ::
!>         The computed lambda_I, the I-th updated eigenvalue.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup auxOTHERcomputational
!
!> \par Contributors:
!  ==================
!>
!>     Ren-Cang Li, Computer Science Division, University of California
!>     at Berkeley, USA
!>
!  =====================================================================
      SUBROUTINE DLAED5( I, D, Z, DELTA, RHO, DLAM )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes
      integer(i4)  ::            I
      real(dp)    ::   DLAM, RHO
!     ..
!     .. Array Arguments ..
      real(dp)    ::   D( 2 ), DELTA( 2 ), Z( 2 )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      real(dp)    ::   ZERO, ONE, TWO, FOUR
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, FOUR = 4.0D0 )
!     ..
!     .. Local Scalars ..
      real(dp)    ::   B, C, DEL, TAU, TEMP, W
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, SQRT
!     ..
!     .. Executable Statements ..
!
      DEL = D( 2 ) - D( 1 )
      IF( I.EQ.1 ) THEN
         W = ONE + TWO*RHO*( Z( 2 )*Z( 2 )-Z( 1 )*Z( 1 ) ) / DEL
         IF( W.GT.ZERO ) THEN
            B = DEL + RHO*( Z( 1 )*Z( 1 )+Z( 2 )*Z( 2 ) )
            C = RHO*Z( 1 )*Z( 1 )*DEL
!
!           B > ZERO, always
!
            TAU = TWO*C / ( B+SQRT( ABS( B*B-FOUR*C ) ) )
            DLAM = D( 1 ) + TAU
            DELTA( 1 ) = -Z( 1 ) / TAU
            DELTA( 2 ) = Z( 2 ) / ( DEL-TAU )
         ELSE
            B = -DEL + RHO*( Z( 1 )*Z( 1 )+Z( 2 )*Z( 2 ) )
            C = RHO*Z( 2 )*Z( 2 )*DEL
            IF( B.GT.ZERO ) THEN
               TAU = -TWO*C / ( B+SQRT( B*B+FOUR*C ) )
            ELSE
               TAU = ( B-SQRT( B*B+FOUR*C ) ) / TWO
            END IF
            DLAM = D( 2 ) + TAU
            DELTA( 1 ) = -Z( 1 ) / ( DEL+TAU )
            DELTA( 2 ) = -Z( 2 ) / TAU
         END IF
         TEMP = SQRT( DELTA( 1 )*DELTA( 1 )+DELTA( 2 )*DELTA( 2 ) )
         DELTA( 1 ) = DELTA( 1 ) / TEMP
         DELTA( 2 ) = DELTA( 2 ) / TEMP
      ELSE
!
!     Now I=2
!
         B = -DEL + RHO*( Z( 1 )*Z( 1 )+Z( 2 )*Z( 2 ) )
         C = RHO*Z( 2 )*Z( 2 )*DEL
         IF( B.GT.ZERO ) THEN
            TAU = ( B+SQRT( B*B+FOUR*C ) ) / TWO
         ELSE
            TAU = TWO*C / ( -B+SQRT( B*B+FOUR*C ) )
         END IF
         DLAM = D( 2 ) + TAU
         DELTA( 1 ) = -Z( 1 ) / ( DEL+TAU )
         DELTA( 2 ) = -Z( 2 ) / TAU
         TEMP = SQRT( DELTA( 1 )*DELTA( 1 )+DELTA( 2 )*DELTA( 2 ) )
         DELTA( 1 ) = DELTA( 1 ) / TEMP
         DELTA( 2 ) = DELTA( 2 ) / TEMP
      END IF
      RETURN
!
!     End OF DLAED5
!
      END
!> \brief \b DLAED6
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!> \htmlonly
!> Download DLAED6 + dependencies 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaed6.f"> 
!> [TGZ]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaed6.f"> 
!> [ZIP]</a> 
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaed6.f"> 
!> [TXT]</a>
!> \endhtmlonly 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAED6( KNITER, ORGATI, RHO, D, Z, FINIT, TAU, INFO )
! 
!       .. Scalar Arguments ..
!       LOGICAL            ORGATI
!       integer(i4)  ::            INFO, KNITER
!       real(dp)    ::   FINIT, RHO, TAU
!       ..
!       .. Array Arguments ..
!       real(dp)    ::   D( 3 ), Z( 3 )
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAED6 computes the positive or negative root (closest to the origin)
!> of
!>                  z(1)        z(2)        z(3)
!> f(x) =   rho + --------- + ---------- + ---------
!>                 d(1)-x      d(2)-x      d(3)-x
!>
!> It is assumed that
!>
!>       if ORGATI = .true. the root is between d(2) and d(3);
!>       otherwise it is between d(1) and d(2)
!>
!> This routine will be called by DLAED4 when necessary. In most cases,
!> the root sought is the smallest in magnitude, though it might not be
!> in some extremely rare situations.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] KNITER
!> \verbatim
!>          KNITER is integer(i4)  ::
!>               Refer to DLAED4 for its significance.
!> \endverbatim
!>
!> \param[in] ORGATI
!> \verbatim
!>          ORGATI is LOGICAL
!>               If ORGATI is true, the needed root is between d(2) and
!>               d(3); otherwise it is between d(1) and d(2).  See
!>               DLAED4 for further details.
!> \endverbatim
!>
!> \param[in] RHO
!> \verbatim
!>          RHO is real(dp)    ::
!>               Refer to the equation f(x) above.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is real(dp)    :: array, dimension (3)
!>               D satisfies d(1) < d(2) < d(3).
!> \endverbatim
!>
!> \param[in] Z
!> \verbatim
!>          Z is real(dp)    :: array, dimension (3)
!>               Each of the elements in z must be positive.
!> \endverbatim
!>
!> \param[in] FINIT
!> \verbatim
!>          FINIT is real(dp)    ::
!>               The value of f at 0. It is more accurate than the one
!>               evaluated inside this routine (if someone wants to do
!>               so).
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is real(dp)    ::
!>               The root of the equation f(x).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is integer(i4)  ::
!>               = 0: successful exit
!>               > 0: if INFO = 1, failure to converge
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup auxOTHERcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  10/02/03: This version has a few statements commented out for thread
!>  safety (machine parameters are computed on each entry). SJH.
!>
!>  05/10/06: Modified from a new version of Ren-Cang Li, use
!>     Gragg-Thornton-Warner cubic convergent scheme for better stability.
!> \endverbatim
!
!> \par Contributors:
!  ==================
!>
!>     Ren-Cang Li, Computer Science Division, University of California
!>     at Berkeley, USA
!>
!  =====================================================================
      SUBROUTINE DLAED6( KNITER, ORGATI, RHO, D, Z, FINIT, TAU, INFO )
!
!  -- LAPACK computational routine (version 3.4.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2011
!
!     .. Scalar Arguments ..
      use datatypes
      LOGICAL            ORGATI
      integer(i4)  ::            INFO, KNITER
      real(dp)    ::   FINIT, RHO, TAU
!     ..
!     .. Array Arguments ..
      real(dp)    ::   D( 3 ), Z( 3 )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      integer(i4)  ::            MAXIT
      PARAMETER          ( MAXIT = 40 )
      real(dp)    ::   ZERO, ONE, TWO, THREE, FOUR, EIGHT
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, THREE = 3.0D0, FOUR = 4.0D0, EIGHT = 8.0D0 )
!     ..
!     .. External Functions ..
      real(dp)    ::   DLAMCH
      EXTERNAL           DLAMCH
!     ..
!     .. Local Arrays ..
      real(dp)    ::   DSCALE( 3 ), ZSCALE( 3 )
!     ..
!     .. Local Scalars ..
      LOGICAL            SCALE
      integer(i4)  ::            I, ITER, NITER
      real(dp)    ::   A, B, BASE, C, DDF, DF, EPS, ERRETM, ETA, F, &
                         FC, SCLFAC, SCLINV, SMALL1, SMALL2, SMINV1, &
                         SMINV2, TEMP, TEMP1, TEMP2, TEMP3, TEMP4,  &
                         LBD, UBD
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, LOG, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
      INFO = 0
!
      IF( ORGATI ) THEN
         LBD = D(2)
         UBD = D(3)
      ELSE
         LBD = D(1)
         UBD = D(2)
      END IF
      IF( FINIT .LT. ZERO )THEN
         LBD = ZERO
      ELSE
         UBD = ZERO 
      END IF
!
      NITER = 1
      TAU = ZERO
      IF( KNITER.EQ.2 ) THEN
         IF( ORGATI ) THEN
            TEMP = ( D( 3 )-D( 2 ) ) / TWO
            C = RHO + Z( 1 ) / ( ( D( 1 )-D( 2 ) )-TEMP )
            A = C*( D( 2 )+D( 3 ) ) + Z( 2 ) + Z( 3 )
            B = C*D( 2 )*D( 3 ) + Z( 2 )*D( 3 ) + Z( 3 )*D( 2 )
         ELSE
            TEMP = ( D( 1 )-D( 2 ) ) / TWO
            C = RHO + Z( 3 ) / ( ( D( 3 )-D( 2 ) )-TEMP )
            A = C*( D( 1 )+D( 2 ) ) + Z( 1 ) + Z( 2 )
            B = C*D( 1 )*D( 2 ) + Z( 1 )*D( 2 ) + Z( 2 )*D( 1 )
         END IF
         TEMP = MAX( ABS( A ), ABS( B ), ABS( C ) )
         A = A / TEMP
         B = B / TEMP
         C = C / TEMP
         IF( C.EQ.ZERO ) THEN
            TAU = B / A
         ELSE IF( A.LE.ZERO ) THEN
            TAU = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
         ELSE
            TAU = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
         END IF
         IF( TAU .LT. LBD .OR. TAU .GT. UBD ) TAU = ( LBD+UBD )/TWO
         IF( D(1).EQ.TAU .OR. D(2).EQ.TAU .OR. D(3).EQ.TAU ) THEN
            TAU = ZERO
         ELSE
            TEMP = FINIT + TAU*Z(1)/( D(1)*( D( 1 )-TAU ) ) + &
                           TAU*Z(2)/( D(2)*( D( 2 )-TAU ) ) + &
                           TAU*Z(3)/( D(3)*( D( 3 )-TAU ) )
            IF( TEMP .LE. ZERO )THEN
               LBD = TAU
            ELSE
               UBD = TAU
            END IF
            IF( ABS( FINIT ).LE.ABS( TEMP ) ) TAU = ZERO
         END IF
      END IF
!
!     get machine parameters for possible scaling to avoid overflow
!
!     modified by Sven: parameters SMALL1, SMINV1, SMALL2,
!     SMINV2, EPS are not SAVEd anymore between one call to the
!     others but recomputed at each call
!
      EPS = DLAMCH( 'Epsilon' )
      BASE = DLAMCH( 'Base' )
      SMALL1 = BASE**( INT( LOG( DLAMCH( 'SafMin' ) ) / LOG( BASE ) / THREE ) )
      SMINV1 = ONE / SMALL1
      SMALL2 = SMALL1*SMALL1
      SMINV2 = SMINV1*SMINV1
!
!     Determine if scaling of inputs necessary to avoid overflow
!     when computing 1/TEMP**3
!
      IF( ORGATI ) THEN
         TEMP = MIN( ABS( D( 2 )-TAU ), ABS( D( 3 )-TAU ) )
      ELSE
         TEMP = MIN( ABS( D( 1 )-TAU ), ABS( D( 2 )-TAU ) )
      END IF
      SCALE = .FALSE.
      IF( TEMP.LE.SMALL1 ) THEN
         SCALE = .TRUE.
         IF( TEMP.LE.SMALL2 ) THEN
!
!        Scale up by power of radix nearest 1/SAFMIN**(2/3)
!
            SCLFAC = SMINV2
            SCLINV = SMALL2
         ELSE
!
!        Scale up by power of radix nearest 1/SAFMIN**(1/3)
!
            SCLFAC = SMINV1
            SCLINV = SMALL1
         END IF
!
!        Scaling up safe because D, Z, TAU scaled elsewhere to be O(1)
!
         DO 10 I = 1, 3
            DSCALE( I ) = D( I )*SCLFAC
            ZSCALE( I ) = Z( I )*SCLFAC
   10    CONTINUE
         TAU = TAU*SCLFAC
         LBD = LBD*SCLFAC
         UBD = UBD*SCLFAC
      ELSE
!
!        Copy D and Z to DSCALE and ZSCALE
!
         DO 20 I = 1, 3
            DSCALE( I ) = D( I )
            ZSCALE( I ) = Z( I )
   20    CONTINUE
      END IF
!
      FC = ZERO
      DF = ZERO
      DDF = ZERO
      DO 30 I = 1, 3
         TEMP = ONE / ( DSCALE( I )-TAU )
         TEMP1 = ZSCALE( I )*TEMP
         TEMP2 = TEMP1*TEMP
         TEMP3 = TEMP2*TEMP
         FC = FC + TEMP1 / DSCALE( I )
         DF = DF + TEMP2
         DDF = DDF + TEMP3
   30 CONTINUE
      F = FINIT + TAU*FC
!
      IF( ABS( F ).LE.ZERO ) GO TO 60
      IF( F .LE. ZERO )THEN
         LBD = TAU
      ELSE
         UBD = TAU
      END IF
!
!        Iteration begins -- Use Gragg-Thornton-Warner cubic convergent
!                            scheme
!
!     It is not hard to see that
!
!           1) Iterations will go up monotonically
!              if FINIT < 0;
!
!           2) Iterations will go down monotonically
!              if FINIT > 0.
!
      ITER = NITER + 1
!
      DO 50 NITER = ITER, MAXIT
!
         IF( ORGATI ) THEN
            TEMP1 = DSCALE( 2 ) - TAU
            TEMP2 = DSCALE( 3 ) - TAU
         ELSE
            TEMP1 = DSCALE( 1 ) - TAU
            TEMP2 = DSCALE( 2 ) - TAU
         END IF
         A = ( TEMP1+TEMP2 )*F - TEMP1*TEMP2*DF
         B = TEMP1*TEMP2*F
         C = F - ( TEMP1+TEMP2 )*DF + TEMP1*TEMP2*DDF
         TEMP = MAX( ABS( A ), ABS( B ), ABS( C ) )
         A = A / TEMP
         B = B / TEMP
         C = C / TEMP
         IF( C.EQ.ZERO ) THEN
            ETA = B / A
         ELSE IF( A.LE.ZERO ) THEN
            ETA = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
         ELSE
            ETA = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
         END IF
         IF( F*ETA.GE.ZERO ) THEN
            ETA = -F / DF
         END IF
!
         TAU = TAU + ETA
         IF( TAU .LT. LBD .OR. TAU .GT. UBD ) TAU = ( LBD + UBD )/TWO 
!
         FC = ZERO
         ERRETM = ZERO
         DF = ZERO
         DDF = ZERO
         DO 40 I = 1, 3
            TEMP = ONE / ( DSCALE( I )-TAU )
            TEMP1 = ZSCALE( I )*TEMP
            TEMP2 = TEMP1*TEMP
            TEMP3 = TEMP2*TEMP
            TEMP4 = TEMP1 / DSCALE( I )
            FC = FC + TEMP4
            ERRETM = ERRETM + ABS( TEMP4 )
            DF = DF + TEMP2
            DDF = DDF + TEMP3
   40    CONTINUE
         F = FINIT + TAU*FC
         ERRETM = EIGHT*( ABS( FINIT )+ABS( TAU )*ERRETM ) + ABS( TAU )*DF
         IF( ABS( F ).LE.EPS*ERRETM ) GO TO 60
         IF( F .LE. ZERO )THEN
            LBD = TAU
         ELSE
            UBD = TAU
         END IF
   50 CONTINUE
      INFO = 1
   60 CONTINUE
!
!     Undo scaling
!
      IF( SCALE ) TAU = TAU*SCLINV
      RETURN
!
!     End of DLAED6
!
      END
!
!  DIVIDE AND CONQUER RELATED ROUTINES BEFORE HERE
!
