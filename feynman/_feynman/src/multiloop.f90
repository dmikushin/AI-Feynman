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
