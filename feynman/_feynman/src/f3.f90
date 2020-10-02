        real*8 function f3(n, arities, ops, x) ! n=number of ops, x=arg vector
           implicit none
           integer nmax, n, i, j, arities(n), arity, lnblnk
           character*60 ops
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
        end
