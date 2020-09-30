        subroutine LoadMatrixTranspose(nd, n, mmax, m, A, f)
           ! Reads the n x m matrix A from the file named f, stored as its transpose
           implicit none
           integer nd, mmax, n, m, j
           real*8 A(nd, mmax)
           character*256 f
           open (2, file=f, status='old')
           m = 0
555        m = m + 1
           if (m .gt. mmax) stop 'DEATH ERROR: m>mmax in LoadVectorTranspose'
           read (2, *, end=666) (A(j, m), j=1, n)
           goto 555
666        close (2)
           m = m - 1
           print *, m, ' rows read from file ', f
           return
        end subroutine LoadMatrixTranspose
