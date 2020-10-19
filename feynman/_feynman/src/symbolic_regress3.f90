        ! Max Tegmark 171119, 190128-31, 190218, 25
        ! Same as symbolic_regress3.f except that it fits for the symbolic formula plus an arbitrary constant.
        ! Loads templates.csv functions.dat and mystery.dat, returns winner.
        ! scp -P2222 symbolic_regress3.f euler@tor.mit.edu:FEYNMAN
        ! COMPILATION: a f 'f77 -O3 -o symbolic_regress3.x symbolic_regress3.f |& more'
        ! SAMPLE USAGE: call symbolic_regress3.x 6ops.txt arity2templates.txt mystery012.dat results.dat
        ! functions.dat contains a single line (say "0>+*-/") with the single-character symbols
        ! that will be used, drawn from this list:
        !
        ! Binary:
        ! +: add
        ! *: multiply
        ! -: subtract
        ! /: divide        (Put "D" instead of "/" in file, since f77 can't load backslash
        !
        ! Unary:
        !  >: increment (x -> x+1)
        !  <: decrement (x -> x-1)
        !  ~: negate          (x-> -x)
        !  \: invert    (x->1/x) (Put "I" instead of "\" in file, since f77 can't load backslash
        !  L: logaritm: (x-> ln(x)
        !  E: exponentiate (x->exp(x))
        !  S: sin:      (x->sin(x))
        !  C: cos:      (x->cos(x))
        !  A: abs:      (x->abs(x))
        !  N: arcsin:   (x->arcsin(x))
        !  T: arctan:   (x->arctan(x))
        !  R: sqrt        (x->sqrt(x))
        !
        ! nonary:
        !  0
        !  1
        !  a, b, c, ...: input variables for function (need not be listed in functions.dat)

        subroutine symbolic_regress3() bind(C)
           implicit none
           character*256 opsfile, templatefile, mysteryfile, outfile, usedfuncs
           character*60 comline, functions, ops, formula
           integer arities(21), nvar, nvarmax, nmax, lnblnk
           parameter(nvarmax=20, nmax=10000000)
           real*8 minloss, epsilon
           real*8, allocatable :: xy(:,:)
           real*8 ymin
           parameter(epsilon=0.00001)
           data arities/2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0/
           data functions/"+*-/><~\OJLESCANTR01P"/
           integer nn(0:2)
           integer, allocatable :: ii(:), radix(:)
           integer ndata, i, j, n, jmin
           integer*8 nformulas
           logical done
           character*60 func(0:2), template

           open (2, file='args.dat', status='old', err=666)
           read (2, *) opsfile, templatefile, mysteryfile, outfile
           close (2)

           comline = 'head -1 '//mysteryfile(1:lnblnk(mysteryfile))//' | wc > qaz.dat'
           if (system(comline) .ne. 0) stop 'DEATH ERROR counting columns'
           open (2, file='qaz.dat')
           read (2, *) i, nvar
           close (2)
           nvar = nvar - 1
           if (nvar .gt. nvarmax) stop 'DEATH ERROR: TOO MANY VARIABLES'
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
           allocate(ii(nmax), radix(nmax))
           call LoadMatrixTranspose(nvarmax + 1, nvar + 1, nmax, ndata, xy, mysteryfile)
           write (*, '(1a24,i8)') 'Number of examples......', ndata
           ! Find min(abs(y)) to use for offset estimation:
           jmin = 1
           ymin = abs(xy(1, nvar + 1))
           do j = 2, ndata
              if (ymin > abs(xy(nvar + 1, j))) then
                 ymin = abs(xy(nvar + 1, j))
                 jmin = j
              end if
           end do
           print *, 'Mystery data has largest magnitude ', ymin, ' at j=', jmin
           print *, 'symbolic_regress3.f90: Searching for best fit...'
           nformulas = 1
           minloss = 1.e6
           template = ''
           open (2, file=templatefile, status='old', err=670)
           open (3, file=outfile)
555        read (2, '(1a60)', end=665) template
           n = lnblnk(template)
           !print *,"template:",template(1:n),"#####"
           do i = 1, n
              ii(i) = ichar(template(i:i)) - 48
              radix(i) = nn(ii(i))
           end do
           done = .false.
           call multiloop(n, radix, loss_loop, loss_report, minloss, nformulas)
           goto 555
665        close (3)
           close (2)
           deallocate(xy)
           deallocate(ii, radix)
           print *, 'All done: results in ', outfile
           return
666        stop 'DEATH ERROR: missing file args.dat'
668        print *, 'DEATH ERROR: missing file ', opsfile(1:lnblnk(opsfile))
           stop
670        print *, 'DEATH ERROR: missing file ', templatefile(1:lnblnk(templatefile))
           stop

        contains

           subroutine loss_loop(kk, prefactor_out, minloss_out, rmsloss_out, ops_out)
              implicit none
              integer :: kk(*)
              real*8 :: prefactor_out, minloss_out, rmsloss_out
              character*60 :: ops, ops_out
              real*8 :: prefactor, newloss, maxloss, rmsloss, f3
              integer :: i, j

              ! Analyze structure ii:
              do i = 1, n
                 ops(i:i) = func(ii(i)) (1 + kk(i):1 + kk(i))
              end do

              prefactor = xy(nvar + 1, jmin) - f3(n, ii, ops, xy(1, jmin))
              maxloss = 0.
              rmsloss = 0.
              do j = 1, ndata
                 newloss = abs(xy(nvar + 1, j) - (prefactor + f3(n, ii, ops, xy(1, j))))
                 if (.not. ((newloss .ge. 0) .or. (newloss .le. 0))) newloss = 1.e30 ! This was a NaN :-)
                 if (newloss .ge. minloss) then
                    minloss_out = newloss
                 endif
                 rmsloss = rmsloss + newloss**2
                 maxloss = max(maxloss, newloss)
              end do
              if (maxloss .lt. minloss) then
                 ! We have a new best fit
                 prefactor_out = prefactor
                 minloss_out = maxloss
                 rmsloss_out = sqrt(rmsloss/ndata)
                 ops_out = ops
              end if

              ! TODO
              ! if (minloss .le. epsilon) then
              ! Stop multiloop
              ! endif

           end subroutine loss_loop

           subroutine loss_report(iformula, prefactor, minloss, rmsloss, ops)
              implicit none

              integer*8, value :: iformula
              real*8, value :: prefactor, minloss, rmsloss
              character*60 :: ops
              real*8 :: DL, DL2, DL3

              DL = log(iformula*max(1., minloss/epsilon))/log(2.)
              DL2 = log(iformula*max(1., minloss/1.e-15))/log(2.)
              DL3 = (log(1.*iformula) + sqrt(1.*ndata)*log(max(1., rmsloss/1.e-15)))/log(2.)
              write (*, '(2f20.12,x,1a22,1i16,4f19.4)') limit(minloss), limit(prefactor), ops(1:n), &
                 iformula, rmsloss, DL, DL2, DL3
              write (3, '(2f20.12,x,1a22,1i16,4f19.4)') limit(minloss), limit(prefactor), ops(1:n), &
                 iformula, rmsloss, DL, DL2, DL3
              flush (3)

           end subroutine loss_report

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

        end

