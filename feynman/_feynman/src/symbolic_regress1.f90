        ! Max Tegmark 171119, 190128-31, 190506
        ! Loads templates.csv functions.dat and mystery.dat, returns winner.
        ! scp -P2222 symbolic_regress1.f euler@tor.mit.edu:FEYNMAN
        ! COMPILATION: a f 'f77 -O3 -o symbolic_regress1.x symbolic_regress1.f |& more'
        ! SAMPLE USAGE: call symbolic_regress1.x 10ops.txt arity2templates.txt mystery_constant.dat results.dat
        ! functions.dat contains a single line (say "0>+*-/") with the single-character symbols
        ! that will be used, drawn from this list:
        !
        ! Binary:
        ! +: add
        ! *: multiply
        ! -: subtract
        ! /: divide        (Put "D" instead of "/" in file, since f77 can't load backslash
        ! Unary:
        !  O: double    (x->2*x); note that this is the letter "O", not zero
        !  J: double+1  (x->2*x+1)
        !  >: increment (x -> x+1)
        !  <: decrement (x -> x-1)
        !  ~: negate          (x-> -x)
        !  \: invert    (x->1/x) (Put "I" instead of "\" in file, since f77 can't load backslash
        !  L: logaritm  (x-> ln(x)
        !  E: exponentiate (x->exp(x))
        !  S: sin:      (x->sin(x))
        !  C: cos:      (x->cos(x))
        !  A: abs:      (x->abs(x))
        !  N: arcsin    (x->arcsin(x))
        !  T: arctan    (x->arctan(x))
        !  R: sqrt        (x->sqrt(x))
        ! nonary:
        !  0
        !  1
        !  P: pi
        !  a, b, c, ...: input variables for function (need not be listed in functions.dat)

        subroutine symbolic_regress1() bind(C)
           implicit none
           character*256 opsfile, templatefile, mysteryfile, outfile, usedfuncs
           character*60 comline, functions, ops, formula
           integer arities(21), nvar, nvarmax, nmax, lnblnk
           parameter(nvarmax=20, nmax=10000000)
           real*8 newloss, minloss, maxloss, rmsloss, epsilon, DL, DL2, DL3
           real*8, allocatable :: xy(:,:)
           parameter(epsilon=0.00000001)
           data arities/2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0/
           data functions/"+*-/><~\OJLESCANTR01P"/
           integer nn(0:2)
           integer, allocatable :: ii(:), kk(:), radix(:)
           integer ndata, i, j, n
           integer*8 nformulas
           logical done
           character*60 func(0:2), template

           open (2, file='args.dat', status='old', err=666)
           read (2, *) opsfile, templatefile, mysteryfile, outfile
           close (2)

           nvar = 0
           write (*, '(1a24,i8)') 'Number of variables.....', nvar

           open (2, file=opsfile, status='old', err=668)
           read (2, *) usedfuncs
           close (2)
           nn(0) = 0
           nn(1) = 0
           nn(2) = 0
           do i = 1, lnblnk(usedfuncs)
              if (usedfuncs(i:i) .eq. 'D') usedfuncs(i:i) = '/'
              if (usedfuncs(i:i) .eq. 'I') usedfuncs(i:i) = '\'
              j = index(functions, usedfuncs(i:i))
              if (j .eq. 0) then
                 print *, 'DEATH ERROR: Unknown function requested: ', usedfuncs(i:i)
                 stop
              else
                 nn(arities(j)) = nn(arities(j)) + 1
                 func(arities(j)) (nn(arities(j)):nn(arities(j))) = functions(j:j)
              end if
           end do
           ! Add nonary ops to retrieve each of the input variables:
           do i = 1, nvar
              nn(0) = nn(0) + 1
              func(0) (nn(0):nn(0)) = char(96 + i)
           end do
           write (*, '(1a24,1a22)') 'Functions used..........', usedfuncs(1:lnblnk(usedfuncs))
           do i = 0, 2
              write (*, *) 'Arity ', i, ': ', func(i) (1:nn(i))
           end do

           write (*, '(1a24)') 'Loading mystery data....'
           allocate(xy(nvarmax + 1, nmax))
           allocate(ii(nmax), kk(nmax), radix(nmax))
           call LoadMatrixTranspose(nvarmax + 1, nvar + 1, nmax, ndata, xy, mysteryfile)
           write (*, '(1a24,i8)') 'Number of examples......', ndata

           print *, 'symbolic_regress1.f90: Searching for best fit...'
           nformulas = 0
           minloss = 1.e6
           template = ''
           ops = '===================='
           open (2, file=templatefile, status='old', err=670)
           open (3, file=outfile)
555        read (2, '(1a60)', end=665) template
           n = lnblnk(template)
           !print *,"template:",template(1:n),"#####"
           do i = 1, n
              ii(i) = ichar(template(i:i)) - 48
              radix(i) = nn(ii(i))
              kk(i) = 0
              !print *,'ASILOMAR ', i,ii(i),kk(i),radix(i)
           end do
           done = .false.
           do while ((minloss .gt. epsilon) .and. (.not. done))
              nformulas = nformulas + 1
              ! Analyze structure ii:
              do i = 1, n
                 ops(i:i) = func(ii(i)) (1 + kk(i):1 + kk(i))
                 !print *,'TEST ',i,ii(i), func(ii(i))
              end do
              !write(*,'(1f20.12,99i3)') minloss, (ii(i),i=1,n), (kk(i),i=1,n)
              !write(*,'(1a24)') ops(1:n)
              j = 1
              maxloss = 0.
              rmsloss = 0.
              do while ((maxloss .lt. minloss) .and. (j .le. ndata))
                 newloss = abs(xy(nvar + 1, j) - f(n, ii, ops, xy(1, j)))
                 rmsloss = rmsloss + newloss**2
            !!!!!print *,'newloss: ',j,newloss,xy(nvar,j),f(n,ii,ops,xy(1,j))
                 if (.not. ((newloss .ge. 0) .or. (newloss .le. 0))) newloss = 1.e30 ! This was a NaN :-)
                 if (maxloss .lt. newloss) maxloss = newloss
                 j = j + 1
              end do
              if (maxloss .lt. minloss) then ! We have a new best fit
                 minloss = maxloss
                 do j = j, ndata
                    newloss = xy(nvar + 1, j) - f(n, ii, ops, xy(1, j))
                    rmsloss = rmsloss + newloss**2
                 end do
                 rmsloss = sqrt(rmsloss/ndata)
                 DL = log(nformulas*max(1., rmsloss/epsilon))/log(2.)
                 DL2 = log(nformulas*max(1., rmsloss/1.e-15))/log(2.)
                 DL3 = (log(1.*nformulas) + sqrt(1.*ndata)*log(max(1., rmsloss/1.e-15)))/log(2.)
                 write (*, '(1f20.12,x,1a22,1i16,4f19.4)') minloss, ops(1:n), nformulas, rmsloss, DL, DL2, DL3
                 write (3, '(1f20.12,x,1a22,1i16,4f19.4)') minloss, ops(1:n), nformulas, rmsloss, DL, DL2, DL3
                 flush (3)
              end if
              call multiloop(n, radix, kk, done)
           end do
           goto 555
665        close (3)
           close (2)
           deallocate(xy)
           deallocate(ii, kk, radix)
           print *, 'All done: results in ', outfile
           return
666        stop 'DEATH ERROR: missing file args.dat'
668        print *, 'DEATH ERROR: missing file ', opsfile(1:lnblnk(opsfile))
           stop
670        print *, 'DEATH ERROR: missing file ', templatefile(1:lnblnk(templatefile))
           stop

        contains

           ! Binary:
           ! +: add
           ! *: multiply
           ! -: subtract
           ! /: divide        (Put "D" instead of "/" in file, since f77 can't load backslash
           ! Unary:
           !  >: increment (x -> x+1)
           !  <: decrement (x -> x-1)
           !  ~: negate          (x-> -x)
           !  \: invert    (x->1/x) (Put "I" instead of "\" in file, since f77 can't load backslash
           !  L: logaritm  (x-> ln(x)
           !  E: exponentiate (x->exp(x))
           !  S: sin:      (x->sin(x))
           !  C: cos:      (x->cos(x))
           !  A: abs:      (x->abs(x))
           !  N: arcsin    (x->arcsin(x))
           !  T: arctan    (x->arctan(x))
           !  R: sqrt        (x->sqrt(x))
           !  O: double    (x->2*x); note that this is the letter "O", not zero
           !  J: double+1  (x->2*x+1)
           ! nonary:
           !  0
           !  1
           !  P: pi
           real*8 function f(n, arities, ops, x) ! n=number of ops, x=arg vector
              implicit none
              integer nmax, n, i, j, arities(n), arity, lnblnk
              character*256 ops
              parameter(nmax=100)
              real*8 x(nmax), y, stack(nmax)
              character op
              !write(*,*) 'Evaluating function with ops = ',ops(1:n)
              !write(*,'(3f10.5,99i3)') (x(i),i=1,3), (arities(i),i=1,n)
              j = 0 ! Number of numbers on the stack
              do i = 1, n
                 arity = arities(i)
                 op = ops(i:i)
                 if (arity .eq. 0) then ! This is a nonary function
                    if (op .eq. "0") then
                       y = 0.
                    else if (op .eq. "1") then
                       y = 1.
                    else if (op .eq. "P") then
                       y = 4.*atan(1.) ! pi
                    else
                       y = x(ichar(op) - 96)
                    end if
                 else if (arity .eq. 1) then ! This is a unary function
                    if (op .eq. ">") then
                       y = stack(j) + 1
                    else if (op .eq. "<") then
                       y = stack(j) - 1
                    else if (op .eq. "~") then
                       y = -stack(j)
                    else if (op .eq. "\") then
                       y = 1./stack(j)
                    else if (op .eq. "L") then
                       y = log(stack(j))
                    else if (op .eq. "E") then
                       y = exp(stack(j))
                    else if (op .eq. "S") then
                       y = sin(stack(j))
                    else if (op .eq. "C") then
                       y = cos(stack(j))
                    else if (op .eq. "A") then
                       y = abs(stack(j))
                    else if (op .eq. "N") then
                       y = asin(stack(j))
                    else if (op .eq. "T") then
                       y = atan(stack(j))
                    else if (op .eq. "O") then
                       y = 2.*stack(j)
                    else if (op .eq. "J") then
                       y = 1 + 2.*stack(j)
                    else
                       y = sqrt(stack(j))
                    end if
                 else ! This is a binary function
                    if (op .eq. "+") then
                       y = stack(j - 1) + stack(j)
                    else if (op .eq. "-") then
                       y = stack(j - 1) - stack(j)
                    else if (op .eq. "*") then
                       y = stack(j - 1)*stack(j)
                    else
                       y = stack(j - 1)/stack(j)
                    end if
                 end if
                 j = j + 1 - arity
                 stack(j) = y
                 ! write(*,'(9f10.5)') (stack(k),k=1,j)
              end do
              if (j .ne. 1) stop 'DEATH ERROR: STACK UNBALANCED'
              f = stack(1)
              !write(*,'(9f10.5)') 666.,x(1),x(2),x(3),f
              return
           end function f

           subroutine multiloop(n, bases, i, done)
              ! Handles <n> nested loops with loop variables i(1),...i(n).
              ! Example: With n=3, bases=2, repeated calls starting with i=(000) will return
              ! 001, 010, 011, 100, 101, 110, 111, 000 (and done=.true. the last time).
              ! All it's doing is counting in mixed radix specified by the array <bases>.
              implicit none
              integer n, bases(n), i(n), k
              logical done
              done = .false.
              k = 1
555           i(k) = i(k) + 1
              if (i(k) .lt. bases(k)) return
              i(k) = 0
              k = k + 1
              if (k .le. n) goto 555
              done = .true.
              return
           end subroutine multiloop

           real*8 function limit(x)
              implicit none
              real*8 x, xmax
              parameter(xmax=666.)
              if (abs(x) .lt. xmax) then
                 limit = x
              else
                 limit = sign(xmax, x)
              end if
              return
           end function limit

           subroutine LoadMatrixTranspose(nd, n, mmax, m, A, f)
              ! Reads the n x m matrix A from the file named f, stored as its transpose
              implicit none
              integer nd, mmax, n, m, j
              real*8 A(nd, mmax)
              character*256 f
              open (2, file=f, status='old')
              m = 0
555           m = m + 1
              if (m .gt. mmax) stop 'DEATH ERROR: m>mmax in LoadVectorTranspose'
              read (2, *, end=666) (A(j, m), j=1, n)
              goto 555
666           close (2)
              m = m - 1
              print *, m, ' rows read from file ', f
              return
           end subroutine LoadMatrixTranspose

        end

