!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File: ct20eval.f90
! Auxiliary subroutines for returning constraint and derivative values
!
! 29 Dec 2008: First version for v4 of snctrl.
! 09 Feb 2010: v5.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module ct20eval
  use precision
  use ct15usr,   only : odeconS, odeconA, odeconD, algconS, algconA, algconD

  implicit none

  private
  public  :: sncEvalTR, & ! contains funConTR, derConTR
             sncEvalHS    ! contains funConHS, derConHS

contains

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine sncEvalTR( Status, mode, nY, nU, nP, nC, nPhs, ndPtr, step, s, &
                        ytype, ctype, neJ, neG, &
                        x, nnJac, fCon, nnCon, gCon, neJac, &
                        usrode, usralg, &
                        cu, lencu, iu, leniu, ru, lenru, &
                        cw, lencw, iw, leniw, rw, lenrw )

    integer(ip), intent(in)    :: Status, mode, nY, nU, nP, nC, nPhs, &
                                  ndPtr(nPhs+1), neJ(nP), neG(nP), &
                                  ytype(nY,nPhs), ctype(nC,nPhs), &
                                  s, nnJac, nnCon, neJac, &
                                  leniu, lencu, lenru, leniw, lencw, lenrw
    real(rp),    intent(in)    :: step(s), x(nnJac)

    integer(ip), intent(inout) :: iu(leniu), iw(leniw)
    real(rp),    intent(inout) :: ru(lenru), rw(lenrw)
    character(8),intent(inout) :: cu(lencu), cw(lencw)

    real(rp),    intent(out)   :: fCon(nnCon), gCon(neJac)
    external     :: usrode, usralg

    !===========================================================================
    ! Evaluates TR-discretized constraints and derivatives.
    !
    ! 02 Jan 2009: First version sncEvalB.
    ! 19 Jan 2009: Only does the TR version. Renamed sncEvalTR.
    ! 28 Jan 2009: Added dense version.
    ! 09 Feb 2010: v5.
    !===========================================================================

    if ( iw(iIntf) == 0 ) then   ! S version

       ! Function values
       if ( mode == 0 .or. mode == 2 ) then
          call funconTR ( Status, nY, nU, nP, nC, nPhs, ndPtr, step, s, &
                          ytype, ctype, &
                          odeconS, algconS, usrode, usralg, &
                          x, nnJac, fCon, nnCon, &
                          cu, lencu, iu, leniu, ru, lenru, &
                          cw, lencw, iw, leniw, rw, lenrw )
       end if

       ! Derivatives
       if ( mode == 1 .or. mode == 2 ) then
          call derConTR ( Status, nY, nU, nP, nC, nPhs, ndPtr, step, s, &
                          ytype, ctype, &
                          odeconS, algconS, usrode, usralg, &
                          x, nnJac, gCon, neJac, &
                          cu, lencu, iu, leniu, ru, lenru, &
                          cw, lencw, iw, leniw, rw, lenrw )
       end if

    else if ( iw(iIntf) == 1 ) then   ! D version

       ! Function values
       if ( mode == 0 .or. mode == 2 ) then
          call funconTR ( Status, nY, nU, nP, nC, nPhs, ndPtr, step, s, &
                          ytype, ctype, &
                          odeconD, algconD, usrode, usralg, &
                          x, nnJac, fCon, nnCon, &
                          cu, lencu, iu, leniu, ru, lenru, &
                          cw, lencw, iw, leniw, rw, lenrw )
       end if

       ! Derivatives
       if ( mode == 1 .or. mode == 2 ) then
          call derConTR ( Status, nY, nU, nP, nC, nPhs, ndPtr, step, s, &
                          ytype, ctype, &
                          odeconD, algconD, usrode, usralg, &
                          x, nnJac, gCon, neJac, &
                          cu, lencu, iu, leniu, ru, lenru, &
                          cw, lencw, iw, leniw, rw, lenrw )
       end if

    else if ( iw(iIntf) == 2 ) then   ! A version

       ! Function values
       if ( mode == 0 .or. mode == 2 ) then
          call funconTR ( Status, nY, nU, nP, nC, nPhs, ndPtr, step, s, &
                          ytype, ctype, &
                          odeconA, algconA, usrode, usralg, &
                          x, nnJac, fCon, nnCon, &
                          cu, lencu, iu, leniu, ru, lenru, &
                          cw, lencw, iw, leniw, rw, lenrw )
       end if

       ! Derivatives
       if ( mode == 1 .or. mode == 2 ) then
          call derConTR ( Status, nY, nU, nP, nC, nPhs, ndPtr, step, s, &
                          ytype, ctype, &
                          odeconA, algconA, usrode, usralg, &
                          x, nnJac, gCon, neJac, &
                          cu, lencu, iu, leniu, ru, lenru, &
                          cw, lencw, iw, leniw, rw, lenrw )
       end if
    end if

  contains

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine funConTR ( Status, nY, nU, nP, nC, nPhs, ndPtr, step, s, &
                          ytype, ctype, &
                          odecon, algcon, usrode, usralg, &
                          x, nnJac, fCon, nnCon, &
                          cu, lencu, iu, leniu, ru, lenru, &
                          cw, lencw, iw, leniw, rw, lenrw )
      integer(ip), intent(in)    :: Status, nY, nU, nP, nC, nPhs, s, nnJac, &
                                    ndPtr(nPhs+1), &
                                    ytype(nY,nPhs), ctype(nC,nPhs), nnCon, &
                                    leniu, lencu, lenru, leniw, lencw, lenrw
      real(rp),    intent(in)    :: step(s), x(nnJac)

      integer(ip), intent(inout) :: iu(leniu), iw(leniw)
      real(rp),    intent(inout) :: ru(lenru), rw(lenrw)
      character(8),intent(inout) :: cu(lencu), cw(lencw)

      real(rp),    intent(out)   :: fCon(nnCon)
      external     :: odecon, algcon, usrode, usralg

      !=========================================================================
      ! This routine gets the user-defined state equations and then
      ! computes the TR-discretized constraints.
      ! Does the ODE constraints first and if necessary, the algebraic
      ! constraints.  Then the results are combined.
      !
      ! 29 Dec 2008: First version of funConTR for snctrl v4.
      ! 24 Jan 2009: Reorganized a bit for sparse structures and phases.
      ! 09 Feb 2010: v5.
      !=========================================================================
      integer(ip) :: nVar, nJ, nG, pind, nNodes, ndS, ndE, indS, indE, &
                     frwind, crwind, ind, rowind, dFlag(nY), i, j, p
      real(rp)    :: ftCon(nnCon), hstep
      integer(ip), allocatable :: Jrow(:), Grow(:), col(:)
      real(rp),    allocatable :: F(:,:), C(:,:), Jval(:,:), Gval(:,:), ctCon(:)


      nVar   = nY+nU
      pind   = nnJac - nP + 1
      fCon   = 0

      allocate ( col(nY+nU+nP+1) )

      !-------------------------------------------------------------------------
      ! Assign ODE function values
      !-------------------------------------------------------------------------
      ftCon  = 0
      frwind = 0
      ind    = 0
      indS   = 1

      do p = 1, nPhs

         nJ  = neJ(p)
         ndS = ndPtr(p)
         ndE = ndPtr(p+1)
         nNodes = nde - ndS + 1
         indE   = indS-1 + nNodes*(nY+nU)

         ! Get user-defined state equations
         allocate( F(nY,ndS:ndE) )
         allocate( Jrow(nJ), Jval(nJ,ndS:ndE) )
         call odecon( Status, usrode, p, nPhs, nY, nU, nP, nNodes, dFlag, &
                      F, Jrow, Jval, col, nJ, x(indS:indE), x(pind), &
                      1_ip, 0_ip, &
                      cu, lencu, iu, leniu, ru, lenru, &
                      cw, lencw, iw, leniw, rw, lenrw )
         indS = indE + 1

         ! Compute the TR-discretized constraints
         do j = ndS, ndE-1
            hstep = step(j)
            do i = 1, nY
               if ( ytype(i,p) == 0 ) then
                  frwind = frwind + 1
                  ftCon(frwind) = -x(i+nVar+ind) + x(i+ind) + &
                       (hstep/2.0)*( F(i,j) + F(i,j+1) )
               end if
            end do
            ind = ind + nVar
         end do

         ind = ind + nVar ! for phase continuity vars
         deallocate ( F, Jrow, Jval )
      end do

      if ( nC == 0 ) then
         fCon = ftCon

         deallocate ( col )
         return
      end if

      !-------------------------------------------------------------------
      ! Assign algebraic constraint function values
      !-------------------------------------------------------------------
      allocate ( ctCon(nnCon) )

      ctCon  = 0
      crwind = 0
      ind    = 0
      indS   = 1

      do p = 1, nPhs

         nG  = neG(p)
         ndS = ndPtr(p)
         ndE = ndPtr(p+1)
         nNodes = nde - ndS + 1
         indE   = indS-1 + nNodes*(nY+nU)

         ! Get user-defined state equations
         allocate( C(nC,ndS:ndE) )
         allocate( Grow(nG), Gval(nG,ndS:ndE) )
         call algcon( Status, usralg, p, nPhs, nC, nY, nU, nP, nNodes, &
                      C, Grow, Gval, col, nG, x(indS:indE), x(pind), &
                      1_ip, 0_ip, &
                      cu, lencu, iu, leniu, ru, lenru, &
                      cw, lencw, iw, leniw, rw, lenrw )

         indS = indE + 1

         do j = ndS, ndE
            do i = 1, nC
               if ( ctype(i,p) == 0 ) then
                  crwind = crwind + 1
                  ctCon(crwind) = C(i,j)
               end if
            end do
         end do

         deallocate ( C, Grow, Gval )
      end do

      deallocate ( col )

      !-------------------------------------------------------------------
      ! Combine algebraic and ode constraints
      !-------------------------------------------------------------------
      rowind = 0
      crwind = 0
      frwind = 0

      do p = 1, nPhs
         do j = ndPtr(p), ndPtr(p+1)-1

            do i = 1, nC
               if ( ctype(i,p) == 0 ) then
                  rowind = rowind + 1
                  crwind = crwind + 1
                  fCon(rowind) = ctCon(crwind)
               end if
            end do

            do i = 1, nY
               if ( ytype(i,p) == 0 ) then
                  rowind = rowind + 1
                  frwind = frwind + 1
                  fCon(rowind) = ftCon(frwind)
               end if
            end do
         end do

         j = ndPtr(p+1)
         do i = 1, nC
            if ( ctype(i,p) == 0 ) then
               rowind = rowind + 1
               crwind = crwind + 1
               fCon(rowind) = ctCon(crwind)
            end if
         end do
      end do

      deallocate ( ctCon )

    end subroutine funConTR

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine derConTR ( Status, nY, nU, nP, nC, nPhs, ndPtr, step, s, &
                          ytype, ctype, &
                          odecon, algcon, usrode, usralg, x, nnJac, &
                          gCon, neJac, &
                          cu, lencu, iu, leniu, ru, lenru, &
                          cw, lencw, iw, leniw, rw, lenrw )
      integer(ip), intent(in)    :: Status, nY, nU, nP, nC, nPhs, s, nnJac, &
                                    ytype(nY,nPhs), ctype(nC,nPhs), &
                                    ndPtr(nPhs+1), neJac, &
                                    leniu, lencu, lenru, leniw, lencw, lenrw
      real(rp),    intent(in)    :: step(s), x(nnJac)

      integer(ip), intent(inout) :: iu(leniu), iw(leniw)
      real(rp),    intent(inout) :: ru(lenru), rw(lenrw)
      character(8),intent(inout) :: cu(lencu), cw(lencw)

      real(rp),    intent(out) :: gCon(neJac)
      external     :: odecon, algcon, usrode, usralg

      !=========================================================================
      ! Returns derivative values for the TR-discretized constraints.
      ! Does the ODE constraints first and if necessary, the algebraic
      ! constraints.  Then the results are combined.
      !
      ! 29 Dec 2008: First version of derConTR for snctrl v4.
      ! 24 Jan 2009: Reorganized a bit for sparse structures and phases.
      ! 09 Feb 2010: v5.
      !=========================================================================
      integer(ip) :: nVar, pind, count, nJ, nG, nNodes, ndS, ndE, indS, indE, &
                     dFlag(nY), i, j, k, p, q
      real(rp)    :: hstep, phstep
      integer(ip), allocatable :: Jrow(:), Grow(:), Jcol(:), Gcol(:)
      real(rp),    allocatable :: F(:,:), C(:,:), Jval(:,:), Gval(:,:)


      nVar  = nY+nU
      pind = nnJac - nP + 1
      gCon  = 0

      if ( nC == 0 ) then
         allocate ( Jcol(nY+nU+nP+1) )

         !----------------------------------------------------------------
         ! Assign ODE derivatives
         !----------------------------------------------------------------
         count = 0
         indS  = 1

         do p = 1, nPhs

            nJ  = neJ(p)
            ndS = ndPtr(p)
            ndE = ndPtr(p+1)
            nNodes = ndE - ndS + 1
            indE   = indS-1 + nNodes*(nY+nU)

            ! Get user-defined Jacobian elements
            allocate( F(nY,ndS:ndE) )
            allocate( Jrow(nJ), Jval(nJ,ndS:ndE) )
            call odecon( Status, usrode, p, nPhs, nY, nU, nP, nNodes, dFlag, &
                         F, Jrow, Jval, Jcol, nJ, x(indS:indE), x(pind), &
                         0_ip, 1_ip, &
                         cu, lencu, iu, leniu, ru, lenru, &
                         cw, lencw, iw, leniw, rw, lenrw )
            indS = indE+1

            ! At first node of the phase:
            j     = ndS
            hstep = step(j)

            do k = 1, nY+nU
               do q = Jcol(k), Jcol(k+1)-1
                  i = Jrow(q)

                  if ( ytype(i,p) == 0 ) then
                     count = count + 1
                     gCon(count) = (hstep/2.0)*Jval(q,j)
                     if ( i == k ) &
                          gCon(count) = gCon(count) + 1.0
                  end if
               end do

               if ( k <= nY ) then
                  if ( dFlag(k) == 0 .and. ytype(k,p) == 0 ) then
                     count = count + 1
                     gCon(count) = 1.0
                  end if
               end if

            end do

            ! At inner nodes:
            do j = ndS+1, ndE-1
               phstep = hstep
               hstep  = step(j)

               do k = 1, nY+nU
                  do q = Jcol(k), Jcol(k+1)-1
                     i = Jrow(q)

                     if ( ytype(i,p) == 0 ) then
                        count = count + 1
                        gCon(count) = (phstep/2.0)*Jval(q,j)
                        if ( k == i) &
                             gCon(count) = gCon(count) - 1.0

                        count = count + 1
                        gCon(count) = (hstep/2.0)*Jval(q,j)
                        if ( k == i) &
                             gCon(count) = gCon(count) + 1.0
                     end if
                  end do

                  if ( k <= nY ) then
                     if ( dFlag(k) == 0 .and. ytype(k,p) == 0 ) then
                        count = count + 1
                        gCon(count) = -1.0

                        count = count + 1
                        gCon(count) = 1.0
                     end if
                  end if
               end do
            end do

            ! At final node of phase:
            j = ndE
            do k = 1, nY+nU
               do q = Jcol(k), Jcol(k+1)-1
                  i = Jrow(q)

                  if ( ytype(i,p) == 0 ) then
                     count = count + 1
                     gCon(count) = (hstep/2.0)*Jval(q,j)
                     if ( k == i) &
                          gCon(count) = gCon(count) - 1.0
                  end if
               end do

               if ( k <= nY ) then
                  if ( dFlag(k) == 0 .and. ytype(k,p) == 0 ) then
                     count = count + 1
                     gCon(count) = -1.0
                  end if
               end if
            end do

            deallocate ( F, Jrow, Jval )
         end do

         ! Parameters
         indS = 1

         do k = 1, nP
            do p = 1, nPhs

               nJ  = neJ(p)
               ndS = ndPtr(p)
               ndE = ndPtr(p+1)
               nNodes = ndE - ndS + 1
               indE   = indS-1 + nNodes*(nY+nU)

               ! Get user-defined Jacobian elements
               allocate( F(nY,ndS:ndE) )
               allocate( Jrow(nJ), Jval(nJ,ndS:ndE) )
               call odecon( Status, usrode, p, nPhs, nY, nU, nP, nNodes, dFlag, &
                            F, Jrow, Jval, Jcol, nJ, x(indS:indE), x(pind), &
                            0_ip, 1_ip, &
                            cu, lencu, iu, leniu, ru, lenru, &
                            cw, lencw, iw, leniw, rw, lenrw )
               indS = indE+1

               do j = ndS, ndE-1
                  hstep = step(j)

                  do q = Jcol(k+nVar), Jcol(k+nVar+1)-1
                     i = Jrow(q)
                     if ( ytype(i,p) == 0 ) then
                        count = count + 1
                        gCon(count) = (hstep/2.0)*( Jval(q,j) + Jval(q,j+1) )
                     end if
                  end do
               end do

               deallocate ( F, Jrow, Jval )
            end do
         end do

         deallocate ( Jcol )
         return
      end if

      !-------------------------------------------------------------------
      ! Assign ODE derivatives and algebraic constraint derivatives
      !-------------------------------------------------------------------
      count = 0
      indS  = 1

      allocate ( Jcol(nY+nU+nP+1), Gcol(nY+nU+nP+1) )
      do p = 1, nPhs

         nJ  = neJ(p)
         nG  = neG(p)
         ndS = ndPtr(p)
         ndE = ndPtr(p+1)
         nNodes = ndE - ndS + 1
         indE   = indS-1 + nNodes*(nY+nU)

         ! Get user-defined Jacobian elements
         allocate( F(nY,ndS:ndE) )
         allocate( Jrow(nJ), Jval(nJ,ndS:ndE) )
         call odecon( Status, usrode, p, nPhs, nY, nU, nP, nNodes, dFlag, &
                      F, Jrow, Jval, Jcol, nJ, x(indS:indE), x(pind), &
                      0_ip, 1_ip, &
                      cu, lencu, iu, leniu, ru, lenru, &
                      cw, lencw, iw, leniw, rw, lenrw )


         ! Get user-defined state equations
         allocate( C(nC,ndS:ndE) )
         allocate( Grow(nG), Gval(nG,ndS:ndE) )
         call algcon( Status, usralg, p, nPhs, nC, nY, nU, nP, nNodes, &
                      C, Grow, Gval, Gcol, nG, x(indS:indE), x(pind), &
                      0_ip, 1_ip, &
                      cu, lencu, iu, leniu, ru, lenru, &
                      cw, lencw, iw, leniw, rw, lenrw )
         indS = indE + 1

         ! At first node of the phase:
         j     = ndS
         hstep = step(j)

         do k = 1, nY+nU
            do q = Gcol(k), Gcol(k+1)-1
               i = Grow(q)
               if ( ctype(i,p) == 0 ) then
                  count = count + 1
                  gCon(count) = Gval(q,j)
               end if
            end do

            do q = Jcol(k), Jcol(k+1)-1
               i = Jrow(q)

               if ( ytype(i,p) == 0 ) then
                  count = count + 1
                  gCon(count) = (hstep/2d+0)*Jval(q,j)
                  if ( i == k ) &
                       gCon(count) = gCon(count) + 1.0
               end if
            end do

            if ( k <= nY ) then
               if ( dFlag(k) == 0 .and. ytype(k,p) == 0 ) then
                  count = count + 1
                  gCon(count) = 1.0
               end if
            end if
         end do

         ! At inner nodes:
         do j = ndS+1, ndE-1
            phstep = hstep
            hstep  = step(j)

            do k = 1, nY+nU
               do q = Gcol(k), Gcol(k+1)-1
                  i = Grow(q)
                  if ( ctype(i,p) == 0 ) then
                     count = count + 1
                     gCon(count) = Gval(q,j)
                  end if
               end do

               do q = Jcol(k), Jcol(k+1)-1
                  i = Jrow(q)

                  if ( ytype(i,p) == 0 ) then
                     count = count + 1
                     gCon(count) = (phstep/2.0)*Jval(q,j)
                     if ( k == i) &
                          gCon(count) = gCon(count) - 1.0

                     count = count + 1
                     gCon(count) = (hstep/2.0d+0)*Jval(q,j)
                     if ( k == i) &
                          gCon(count) = gCon(count) + 1.0

                  end if
               end do

               if ( k <= nY ) then
                  if ( dFlag(k) == 0 .and. ytype(k,p) == 0 ) then
                     count = count + 1
                     gCon(count) = -1.0

                     count = count + 1
                     gCon(count) = 1.0
                  end if
               end if
            end do
         end do

         ! At final node of phase:
         j = ndE
         do k = 1, nY+nU
            do q = Gcol(k), Gcol(k+1)-1
               i = Grow(q)
               if ( ctype(i,p) == 0 ) then
                  count = count + 1
                  gCon(count) = Gval(q,j)
               end if
            end do

            do q = Jcol(k), Jcol(k+1)-1
               i = Jrow(q)

               if ( ytype(i,p) == 0 ) then
                  count = count + 1
                  gCon(count) = (hstep/2.0)*Jval(q,j)
                  if ( k == i) &
                       gCon(count) = gCon(count) - 1.0
               end if
            end do

            if ( k <= nY ) then
               if ( dFlag(k) == 0 .and. ytype(k,p) == 0 ) then
                  count = count + 1
                  gCon(count) = -1.0
               end if
            end if
         end do

         deallocate ( F, Jrow, Jval )
         deallocate ( C, Grow, Gval )
      end do

      ! Parameters
      indS = 1

      do k = 1, nP
         do p = 1, nPhs

            nJ  = neJ(p)
            nG  = neG(p)
            ndS = ndPtr(p)
            ndE = ndPtr(p+1)
            nNodes = ndE - ndS + 1
            indE   = indS-1 + nNodes*(nY+nU)

            ! Get user-defined Jacobian elements
            allocate( F(nY,ndS:ndE) )
            allocate( Jrow(nJ), Jval(nJ,ndS:ndE) )
            call odecon( Status, usrode, p, nPhs, nY, nU, nP, nNodes, dFlag, &
                         F, Jrow, Jval, Jcol, nJ, x(indS:indE), x(pind), &
                         0_ip, 1_ip, &
                         cu, lencu, iu, leniu, ru, lenru, &
                         cw, lencw, iw, leniw, rw, lenrw )

            ! Get user-defined state equations
            allocate( C(nC,ndS:ndE) )
            allocate( Grow(nG), Gval(nG,ndS:ndE) )
            call algcon( Status, usralg, p, nPhs, nC, nY, nU, nP, nNodes, &
                         C, Grow, Gval, Gcol, nG, x(indS:indE), x(pind), &
                         0_ip, 1_ip, &
                         cu, lencu, iu, leniu, ru, lenru, &
                         cw, lencw, iw, leniw, rw, lenrw )

            indS = indE+1

            ! Alg constraints
            do j = ndS, ndE
               do q = Gcol(k+nVar), Gcol(k+nVar+1)-1
                  i = Grow(q)
                  if ( ctype(i,p) == 0 ) then
                     count = count + 1
                     gCon(count) = Gval(q,j)
                  end if
               end do
            end do

            ! ODE constraints
            do j = ndS, ndE-1
               hstep = step(j)
               do q = Jcol(k+nVar), Jcol(k+nVar+1)-1
                  i = Jrow(q)
                  if ( ytype(i,p) == 0 ) then
                     count = count + 1
                     gCon(count) = (hstep/2.0)*( Jval(q,j) + Jval(q,j+1) )
                  end if
               end do
            end do

            deallocate ( F, Jrow, Jval )
            deallocate ( C, Grow, Gval )
         end do
      end do

      deallocate ( Jcol, Gcol )

    end subroutine derConTR

  end subroutine sncEvalTR

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine sncEvalHS( Status, mode, nY, nU, nP, nC, nPhs, ndPtr, step, s, &
                        ytype, ctype, neJ, neG, &
                        x, nnJac, fCon, nnCon, gCon, neJac, &
                        usrode, usralg, &
                        cu, lencu, iu, leniu, ru, lenru, &
                        cw, lencw, iw, leniw, rw, lenrw )
    integer(ip),  intent(in)   :: Status, mode, nY, nU, nP, nC, nPhs, &
                                  ndPtr(nPhs+1), neJ(nP), neG(nP), &
                                  s, nnJac, nnCon, neJac, &
                                  ytype(nY,nPhs), ctype(nC,nPhs), &
                                  leniu, lencu, lenru, leniw, lencw, lenrw
    real(rp),    intent(in)    :: step(s), x(nnJac)

    integer(ip), intent(inout) :: iu(leniu), iw(leniw)
    real(rp),    intent(inout) :: ru(lenru), rw(lenrw)
    character(8),intent(inout) :: cu(lencu), cw(lencw)

    real(rp),    intent(out)   :: fCon(nnCon), gCon(neJac)
    external     :: usrode, usralg

    !===========================================================================
    ! Evaluates HS-discretized constraints and derivatives.
    !
    ! 02 Jan 2009: First version sncEvalB.
    ! 19 Jan 2009: Only does HS.  Renamed sncEvalHS.
    ! 28 Jan 2009: Added dense version.
    ! 09 Feb 2010: v5.
    !===========================================================================

    if ( iw(iIntf) == 0 ) then   ! S version

       ! Function values
       if ( mode == 0 .or. mode == 2 ) then
          call funconHS( Status, nY, nU, nP, nC, nPhs, ndPtr, step, s, &
                         ytype, ctype, &
                         odeconS, algconS, usrode, usralg, &
                         x, nnJac, fCon, nnCon, &
                         cu, lencu, iu, leniu, ru, lenru, &
                         cw, lencw, iw, leniw, rw, lenrw )
       end if

       ! Derivatives
       if ( mode == 1 .or. mode == 2 ) then
          call derConHS( Status, nY, nU, nP, nC, nPhs, ndPtr, step, s, &
                         ytype, ctype, &
                         odeconS, algconS, usrode, usralg, &
                         x, nnJac, gCon, neJac, &
                         cu, lencu, iu, leniu, ru, lenru, &
                         cw, lencw, iw, leniw, rw, lenrw )
       end if

    else if ( iw(iIntf) == 1 ) then   ! D version

       ! Function values
       if ( mode == 0 .or. mode == 2 ) then
          call funconHS( Status, nY, nU, nP, nC, nPhs, ndPtr, step, s, &
                         ytype, ctype, &
                         odeconD, algconD, usrode, usralg, &
                         x, nnJac, fCon, nnCon, &
                         cu, lencu, iu, leniu, ru, lenru, &
                         cw, lencw, iw, leniw, rw, lenrw )
       end if

       ! Derivatives
       if ( mode == 1 .or. mode == 2 ) then
          call derConHS( Status, nY, nU, nP, nC, nPhs, ndPtr, step, s, &
                         ytype, ctype, &
                         odeconD, algconD, usrode, usralg, &
                         x, nnJac, gCon, neJac, &
                         cu, lencu, iu, leniu, ru, lenru, &
                         cw, lencw, iw, leniw, rw, lenrw )
       end if

    else if ( iw(iIntf) == 2 ) then   ! A version

       ! Function values
       if ( mode == 0 .or. mode == 2 ) then
          call funconHS( Status, nY, nU, nP, nC, nPhs, ndPtr, step, s, &
                         ytype, ctype, &
                         odeconA, algconA, usrode, usralg, &
                         x, nnJac, fCon, nnCon, &
                         cu, lencu, iu, leniu, ru, lenru, &
                         cw, lencw, iw, leniw, rw, lenrw )
       end if

       ! Derivatives
       if ( mode == 1 .or. mode == 2 ) then
          call derConHS( Status, nY, nU, nP, nC, nPhs, ndPtr, step, s, &
                         ytype, ctype, &
                         odeconA, algconA, usrode, usralg, &
                         x, nnJac, gCon, neJac, &
                         cu, lencu, iu, leniu, ru, lenru, &
                         cw, lencw, iw, leniw, rw, lenrw )
       end if
    end if

  contains

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine funConHS( Status, nY, nU, nP, nC, nPhs, ndPtr, step, s, &
                         ytype, ctype, &
                         odecon, algcon, usrode, usralg, &
                         x, nnJac, fCon, nnCon, &
                         cu, lencu, iu, leniu, ru, lenru, &
                         cw, lencw, iw, leniw, rw, lenrw )

      integer(ip), intent(in)    :: Status, nY, nU, nP, nC, nPhs, s, nnJac, &
                                    ndPtr(nPhs+1), nnCon, &
                                    ytype(nY,nPhs), ctype(nC,nPhs), &
                                    leniu, lencu, lenru, leniw, lencw, lenrw
      real(rp),    intent(in)    :: step(s), x(nnJac)

      integer(ip), intent(inout) :: iu(leniu), iw(leniw)
      real(rp),    intent(inout) :: ru(lenru), rw(lenrw)
      character(8),intent(inout) :: cu(lencu), cw(lencw)

      real(rp),    intent(out)   :: fCon(nnCon)
      external     :: odecon, algcon, usrode, usralg

      !=========================================================================
      ! Returns functions values for HS-discretized constraints.
      ! Does the ODE constraints first and if necessary, the algebraic
      ! constraints.  Then the results are combined.
      !
      ! 29 Dec 2008: First version of funConHS for snctrl v4.
      ! 09 Feb 2010: v5.
      !=========================================================================
      integer(ip) :: nVar, nJ, nG, pind, nNodes, ndS, ndE, indS, indE, &
                     frwind, crwind, ind, rowind, dFlag(nY), i, j, p
      real(rp)    :: ftCon(nnCon), hstep
      integer(ip), allocatable :: Jrow(:), Grow(:), col(:)
      real(rp),    allocatable :: F(:,:), C(:,:), Jval(:,:), Gval(:,:), ctCon(:)

      nVar   = nY+nU
      pind  = nnJac - nP + 1
      fCon   = 0

      allocate ( col(nY+nU+nP+1) )
      !-------------------------------------------------------------------------
      ! Assign ODE function values
      !-------------------------------------------------------------------------
      ftCon  = 0
      frwind = 0
      ind    = 0
      indS   = 1

      do p = 1, nPhs

         nJ  = neJ(p)
         ndS = ndPtr(p)
         ndE = ndPtr(p+1)
         nNodes = nde - ndS + 1
         indE   = indS-1 + nNodes*(nY+nU)

         ! Get user-defined state equations
         allocate ( F(nY,ndS:ndE) )
         allocate ( Jrow(nJ), Jval(nJ,ndS:ndE) )
         call odecon( Status, usrode, p, nPhs, nY, nU, nP, nNodes, dFlag, &
                      F, Jrow, Jval, col, nJ, x(indS:indE), x(pind), &
                      1_ip, 0_ip, &
                      cu, lencu, iu, leniu, ru, lenru, &
                      cw, lencw, iw, leniw, rw, lenrw )

         indS = indE + 1

         do j = ndS, ndE-2, 2
            hstep = step((j+1)/2)

            do i = 1, nY     ! ode constraint 1
               if ( ytype(i,p) == 0 ) then
                  frwind = frwind + 1
                  ftCon(frwind) = -x(i+2*nVar+ind) + x(i+ind) + &
                       (hstep/6.0d+0)*( F(i,j) + 4*F(i,j+1) + F(i,j+2) )
               end if
            end do

            do i = 1, nY     ! ode constraint 2
               if ( ytype(i,p) == 0 ) then
                  frwind = frwind + 1
                  ftCon(frwind) = -x(i+nVar+ind) + &
                       0.5d+0*( x(i+ind) + x(i+ind+2*nVar) ) + &
                       (hstep/8.0d+0)*( F(i,j) - F(i,j+2) )
               end if
            end do

            ind = ind + 2*nVar
         end do

         ind = ind + nVar ! for phase continuity vars
         deallocate ( F, Jrow, Jval )
      end do

      if ( nC == 0 ) then
         fCon(1:nnCon) = ftCon(1:nnCon)
         deallocate ( col )
         return
      end if

      !-------------------------------------------------------------------------
      ! Assign algebraic constraint function values
      !-------------------------------------------------------------------------
      allocate ( ctCon(nnCon) )

      ctCon  = 0
      crwind = 0
      ind    = 0
      indS   = 1

      do p = 1, nPhs

         nG  = neG(p)
         ndS = ndPtr(p)
         ndE = ndPtr(p+1)
         nNodes = nde - ndS + 1
         indE   = indS-1 + nNodes*(nY+nU)

         ! Get user-defined state equations
         allocate ( C(nC,ndS:ndE) )
         allocate ( Grow(nG), Gval(nG,ndS:ndE) )
         call algcon( Status, usralg, p, nPhs, nC, nY, nU, nP, nNodes, &
                      C, Grow, Gval, col, nG, x(indS:indE), x(pind), &
                      1_ip, 0_ip, &
                      cu, lencu, iu, leniu, ru, lenru, &
                      cw, lencw, iw, leniw, rw, lenrw )

         indS = indE + 1

         do j = ndS, ndE, 2
            do i = 1, nC
               if ( ctype(i,p) == 0 ) then
                  crwind = crwind + 1
                  ctCon(crwind) = C(i,j)
               end if
            end do
         end do

         deallocate ( C, Grow, Gval )
      end do

      deallocate ( col )

      !-------------------------------------------------------------------------
      ! Combine algebraic and ode constraints
      !-------------------------------------------------------------------------
      rowind = 0
      crwind = 0
      frwind = 0

      do p = 1, nPhs
         do j = ndPtr(p), ndPtr(p+1)-2, 2

            do i = 1, nC
               if ( ctype(i,p) == 0 ) then
                  rowind = rowind + 1
                  crwind = crwind + 1
                  fCon(rowind) = ctCon(crwind)
               end if
            end do

            do i = 1, nY
               if ( ytype(i,p) == 0 ) then
                  rowind = rowind + 1
                  frwind = frwind + 1
                  fCon(rowind) = ftCon(frwind)
               end if
            end do

            do i = 1, nY
               if ( ytype(i,p) == 0 ) then
                  rowind = rowind + 1
                  frwind = frwind + 1
                  fCon(rowind) = ftCon(frwind)
               end if
            end do
         end do

         j = ndPtr(p+1)
         do i = 1, nC
            if ( ctype(i,p) == 0 ) then
               rowind = rowind + 1
               crwind = crwind + 1
               fCon(rowind) = ctCon(crwind)
            end if
         end do
      end do

      deallocate ( ctCon )

    end subroutine funConHS

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine derConHS( Status, nY, nU, nP, nC, nPhs, ndPtr, step, s, &
                         ytype, ctype, &
                         odecon, algcon, usrode, usralg, &
                         x, nnJac, gCon, neJac, &
                         cu, lencu, iu, leniu, ru, lenru, &
                         cw, lencw, iw, leniw, rw, lenrw )
      integer(ip), intent(in)    :: Status, nY, nU, nP, nC, nPhs, s, nnJac, &
                                    ndPtr(nPhs+1), neJac, &
                                    ytype(nY,nPhs), ctype(nC,nPhs), &
                                    leniu, lencu, lenru, leniw, lencw, lenrw
      real(rp),    intent(in)    :: step(s), x(nnJac)

      integer(ip), intent(inout) :: iu(leniu), iw(leniw)
      real(rp),    intent(inout) :: ru(lenru), rw(lenrw)
      character(8),intent(inout) :: cu(lencu), cw(lencw)

      real(rp),    intent(out) :: gCon(neJac)
      external     :: odecon, algcon, usrode, usralg

      !=========================================================================
      ! Returns derivative values for the HS-discretized constraints.
      ! Does the ODE constraints first and if necessary, the algebraic
      ! constraints.  Then the results are combined.
      !
      ! 29 Dec 2008: First version of derConHS for snctrl v4.
      ! 27 Jan 2009: Reorganized for sparse structures and phases.
      ! 09 Feb 2010: v5.
      !=========================================================================
      integer(ip) :: nVar, pind, count, nJ, nG, nNodes, ndS, ndE, &
                     indS, indE, dFlag(nY), i, j, k, p, q
      real(rp)    :: hstep, phstep
      integer(ip), allocatable :: Jrow(:), Jcol(:), Grow(:), Gcol(:)
      real(rp),    allocatable :: F(:,:), C(:,:), Jval(:,:), Gval(:,:)


      nVar = nY+nU
      pind = nnJac - nP + 1
      gCon = 0

      if ( nC == 0 ) then
         !----------------------------------------------------------------------
         ! Assign ODE derivatives
         !----------------------------------------------------------------------
         count = 0
         indS  = 1

         allocate ( Jcol(nY+nU+nP+1) )

         do p = 1, nPhs
            nJ  = neJ(p)
            ndS = ndPtr(p)
            ndE = ndPtr(p+1)
            nNodes = ndE - ndS + 1
            indE   = indS-1 + nNodes*(nY+nU)

            ! Get user-defined Jacobian elements
            allocate ( F(nY,ndS:ndE) )
            allocate ( Jrow(nJ), Jval(nJ,ndS:ndE) )
            call odecon( Status, usrode, p, nPhs, nY, nU, nP, nNodes, dFlag, &
                         F, Jrow, Jval, Jcol, nJ, x(indS:indE), x(pind), &
                         0_ip, 1_ip, &
                         cu, lencu, iu, leniu, ru, lenru, &
                         cw, lencw, iw, leniw, rw, lenrw )

            indS = indE+1

            ! At first node of the phase:
            j     = ndS
            hstep = step( (j+1)/2 )

            do k = 1, nY+nU
               do q = Jcol(k), Jcol(k+1)-1
                  i = Jrow(q)

                  if ( ytype(i,p) == 0 ) then
                     count = count + 1
                     gCon(count) = (hstep/6.0)*Jval(q,j)
                     if ( k == i ) &
                          gCon(count) = gCon(count) + 1.0

                     count = count + 1
                     gCon(count) = (hstep/8.0)*Jval(q,j)
                     if ( k == i ) &
                          gCon(count) = gCon(count) + 0.5
                  end if
               end do

               if ( k <= nY ) then
                  if ( dFlag(k) == 0 .and. ytype(k,p) == 0 ) then
                     count = count + 1
                     gCon(count) = 1.0

                     count = count + 1
                     gCon(count) = 0.5
                  end if
               end if

            end do

            ! Midpoint
            do k = 1, nY+nU
               do q = Jcol(k), Jcol(k+1)-1
                  i = Jrow(q)

                  if ( ytype(i,p) == 0 ) then
                     count = count + 1
                     gCon(count) = (2*hstep/3.0)*Jval(q,j+1)

                     count = count + 1
                     gCon(count) = 0.0
                     if ( k == i ) &
                          gCon(count) = -1.0
                  end if
               end do

               if ( k <= nY ) then
                  if ( dFlag(k) == 0 .and. ytype(k,p) == 0 ) then
                     count = count + 1
                     gCon(count) = -1.0
                  end if
               end if
            end do

            ! At inner nodes:
            do j = ndS+2, ndE-2, 2
               phstep = hstep
               hstep  = step((j+1)/2)

               do k = 1, nY+nU
                  do q = Jcol(k), Jcol(k+1)-1
                     i = Jrow(q)

                     if ( ytype(i,p) == 0 ) then
                        count = count + 1
                        gCon(count) = (phstep/6.0)*Jval(q,j)
                        if ( k == i ) &
                             gCon(count) = gCon(count) - 1.0

                        count = count + 1
                        gCon(count) = -(phstep/8.0)*Jval(q,j)
                        if ( k == i ) &
                             gCon(count) = gCon(count) + 0.5


                        count = count + 1
                        gCon(count) = (hstep/6.0)*Jval(q,j)
                        if ( k == i ) &
                             gCon(count) = gCon(count) + 1.0

                        count = count + 1
                        gCon(count) = (hstep/8.0)*Jval(q,j)
                        if ( k == i ) &
                             gCon(count) = gCon(count) + 0.5
                     end if
                  end do

                  if ( k <= nY ) then
                     if ( dFlag(k) == 0 .and. ytype(k,p) == 0 ) then
                        count = count + 1
                        gCon(count) = -1.0

                        count = count + 1
                        gCon(count) = 0.5

                        count = count + 1
                        gCon(count) = 1.0

                        count = count + 1
                        gCon(count) = 0.5
                     end if
                  end if
               end do

               ! Midpoint
               do k = 1, nY+nU
                  do q = Jcol(k), Jcol(k+1)-1
                     i = Jrow(q)

                     if ( ytype(i,p) == 0 ) then
                        count = count + 1
                        gCon(count) = (2*hstep/3.0)*Jval(q,j+1)

                        count = count + 1
                        gCon(count) = 0.0
                        if ( k == i) &
                             gCon(count) = -1.0
                     end if
                  end do

                  if ( k <= nY ) then
                     if ( dFlag(k) == 0 .and. ytype(k,p) == 0 ) then
                        count = count + 1
                        gCon(count) = -1.0
                     end if
                  end if
               end do
            end do

            ! At final node of phase:
            j = ndE
            phstep = hstep

            do k = 1, nY+nU
               do q = Jcol(k), Jcol(k+1)-1
                  i = Jrow(q)

                  if ( ytype(i,p) == 0 ) then
                     count = count + 1
                     gCon(count) = (phstep/6.0)*Jval(q,j)
                     if ( k == i ) &
                          gCon(count) = gCon(count) - 1.0

                     count = count + 1
                     gCon(count) = -(phstep/8.0)*Jval(q,j)
                     if ( k == i ) &
                          gCon(count) = gCon(count) + 0.5
                  end if
               end do

               if ( k <= nY ) then
                  if ( dFlag(k) == 0 .and. ytype(k,p) == 0 ) then
                     count = count + 1
                     gCon(count) = -1.0

                     count = count + 1
                     gCon(count) = 0.5
                  end if
               end if
            end do

            deallocate ( F, Jrow, Jval )
         end do

         ! Parameters
         indS = 1
         do k = 1, nP
            do p = 1, nPhs

               nJ  = neJ(p)
               ndS = ndPtr(p)
               ndE = ndPtr(p+1)
               nNodes = ndE - ndS + 1
               indE   = indS-1 + nNodes*(nY+nU)

               ! Get user-defined Jacobian elements
               allocate ( F(nY,ndS:ndE) )
               allocate ( Jrow(nJ), Jval(nJ,ndS:ndE) )
               call odecon( Status, usrode, p, nPhs, nY, nU, nP, nNodes, dFlag, &
                            F, Jrow, Jval, Jcol, nJ, x(indS:indE), x(pind), &
                            0_ip, 1_ip, &
                            cu, lencu, iu, leniu, ru, lenru, &
                            cw, lencw, iw, leniw, rw, lenrw )

               indS = indE+1

               do j = ndS, ndE-2, 2
                  hstep = step( (j+1)/2 )

                  ! ODE Constraint 1 and 2
                  do q = Jcol(k+nVar), Jcol(k+nVar+1)-1
                     i = Jrow(q)

                     if ( ytype(i,p) == 0 ) then
                        count = count + 1
                        gCon(count) = (hstep/6.0)*( Jval(q,j) + &
                             4*Jval(q,j+1) + Jval(q,j+2) )
                        count = count + 1
                        gCon(count) = (hstep/8.0)*( Jval(q,j) - Jval(q,j+2) )
                     end if
                  end do
               end do

               deallocate ( F, Jrow, Jval )
            end do
         end do

         deallocate ( Jcol )
         return
      end if


      !-------------------------------------------------------------------------
      ! Assign ODE derivatives and algebraic constraint derivatives
      !-------------------------------------------------------------------------
      count = 0
      indS  = 1

      allocate ( Jcol(nY+nU+nP+1), Gcol(nY+nU+nP+1) )

      do p = 1, nPhs

         nJ  = neJ(p)
         nG  = neG(p)
         ndS = ndPtr(p)
         ndE = ndPtr(p+1)
         nNodes = ndE - ndS + 1
         indE   = indS-1 + nNodes*(nY+nU)

         ! Get user-defined Jacobian elements
         allocate ( F(nY,ndS:ndE) )
         allocate ( Jrow(nJ), Jval(nJ,ndS:ndE) )
         call odecon( Status, usrode, p, nPhs, nY, nU, nP, nNodes, dFlag, &
                      F, Jrow, Jval, Jcol, nJ, x(indS:indE), x(pind), &
                      0_ip, 1_ip, &
                      cu, lencu, iu, leniu, ru, lenru, &
                      cw, lencw, iw, leniw, rw, lenrw )

         ! Get user-defined state equations
         allocate ( C(nC,ndS:ndE) )
         allocate ( Grow(nG), Gval(nG,ndS:ndE) )
         call algcon( Status, usralg, p, nPhs, nC, nY, nU, nP, nNodes, &
                      C, Grow, Gval, Gcol, nG, x(indS:indE), x(pind), &
                      0_ip, 1_ip, &
                      cu, lencu, iu, leniu, ru, lenru, &
                      cw, lencw, iw, leniw, rw, lenrw )

         indS = indE+1

         ! At first node of the phase:
         j     = ndS
         hstep = step( (j+1)/2 )

         do k = 1, nY+nU
            do q = Gcol(k), Gcol(k+1)-1
               i = Grow(q)

               if ( ctype(i,p) == 0 ) then
                  count = count + 1
                  gCon(count) = Gval(q,j)
               end if
            end do

            do q = Jcol(k), Jcol(k+1)-1
               i = Jrow(q)

               if ( ytype(i,p) == 0 ) then
                  count = count + 1
                  gCon(count) = (hstep/6.0)*Jval(q,j)
                  if ( k == i ) &
                       gCon(count) = gCon(count) + 1.0

                  count = count + 1
                  gCon(count) = (hstep/8.0)*Jval(q,j)
                  if ( k == i ) &
                       gCon(count) = gCon(count) + 0.5
               end if
            end do

            if ( k <= nY ) then
               if ( dFlag(k) == 0 .and. ytype(k,p) == 0 ) then
                  count = count + 1
                  gCon(count) = 1.0

                  count = count + 1
                  gCon(count) = 0.5
               end if
            end if
         end do

         ! Midpoint
         do k = 1, nY+nU
            do q = Jcol(k), Jcol(k+1)-1
               i = Jrow(q)

               if ( ytype(i,p) == 0 ) then
                  count = count + 1
                  gCon(count) = (2*hstep/3.0)*Jval(q,j+1)

                  count = count + 1
                  gCon(count) = 0.0
                  if ( k == i ) &
                       gCon(count) = -1.0
               end if
            end do

            if ( k <= nY ) then
               if ( dFlag(k) == 0 .and. ytype(k,p) == 0 ) then
                  count = count + 1
                  gCon(count) = -1.0
               end if
            end if
         end do

         ! At inner nodes:
         do j = ndS+2, ndE-2, 2
            phstep = hstep
            hstep  = step((j+1)/2)

            do k = 1, nY+nU
               do q = Gcol(k), Gcol(k+1)-1
                  i = Grow(q)

                  if ( ctype(i,p) == 0 ) then
                     count = count + 1
                     gCon(count) = Gval(q,j)
                  end if
               end do

               do q = Jcol(k), Jcol(k+1)-1
                  i = Jrow(q)

                  if ( ytype(i,p) == 0 ) then
                     count = count + 1
                     gCon(count) = (phstep/6.0)*Jval(q,j)
                     if ( k == i ) &
                          gCon(count) = gCon(count) - 1.0

                     count = count + 1
                     gCon(count) = -(phstep/8.0)*Jval(q,j)
                     if ( k == i ) &
                          gCon(count) = gCon(count) + 0.5


                     count = count + 1
                     gCon(count) = (hstep/6.0)*Jval(q,j)
                     if ( k == i ) &
                          gCon(count) = gCon(count) + 1.0

                     count = count + 1
                     gCon(count) = (hstep/8.0)*Jval(q,j)
                     if ( k == i ) &
                          gCon(count) = gCon(count) + 0.5
                  end if
               end do

               if ( k <= nY ) then
                  if ( dFlag(k) == 0 .and. ytype(k,p) == 0 ) then
                     count = count + 1
                     gCon(count) = -1.0

                     count = count + 1
                     gCon(count) = 0.5


                     count = count + 1
                     gCon(count) = 1.0

                     count = count + 1
                     gCon(count) = 0.5
                  end if
               end if
            end do

            ! Midpoint
            do k = 1, nY+nU
               do q = Jcol(k), Jcol(k+1)-1
                  i = Jrow(q)

                  if ( ytype(i,p) == 0 ) then
                     count = count + 1
                     gCon(count) = (2*hstep/3.0)*Jval(q,j+1)

                     count = count + 1
                     gCon(count) = 0.0
                     if ( k == i ) &
                          gCon(count) = -1.0
                  end if
               end do

               if ( k <= nY ) then
                  if ( dFlag(k) == 0 .and. ytype(k,p) == 0 ) then
                     count = count + 1
                     gCon(count) = -1.0
                  end if
               end if
            end do
         end do

         ! At final node of phase:
         j = ndE
         phstep = hstep

         do k = 1, nY+nU
            do q = Gcol(k), Gcol(k+1)-1
               i = Grow(q)

               if ( ctype(i,p) == 0 ) then
                  count = count + 1
                  gCon(count) = Gval(q,j)
               end if
            end do

            do q = Jcol(k), Jcol(k+1)-1
               i = Jrow(q)

               if ( ytype(i,p) == 0 ) then
                  count = count + 1
                  gCon(count) = (phstep/6.0)*Jval(q,j)
                  if ( k == i ) &
                       gCon(count) = gCon(count) - 1.0

                  count = count + 1
                  gCon(count) = -(phstep/8.0)*Jval(q,j)
                  if ( k == i ) &
                       gCon(count) = gCon(count) + 0.5
               end if
            end do

            if ( k <= nY ) then
               if ( dFlag(k) == 0 .and. ytype(k,p) == 0 ) then
                  count = count + 1
                  gCon(count) = -1.0

                  count = count + 1
                  gCon(count) = 0.5
               end if
            end if
         end do

         deallocate ( F, Jrow, Jval )
         deallocate ( C, Grow, Gval )
      end do

      ! Parameters
      indS = 1
      do k = 1, nP
         do p = 1, nPhs

            nJ  = neJ(p)
            nG  = neG(p)
            ndS = ndPtr(p)
            ndE = ndPtr(p+1)
            nNodes = ndE - ndS + 1
            indE   = indS-1 + nNodes*(nY+nU)

            ! Get user-defined Jacobian elements
            allocate ( F(nY,ndS:ndE) )
            allocate ( Jrow(nJ), Jval(nJ,ndS:ndE) )
            call odecon( Status, usrode, p, nPhs, nY, nU, nP, nNodes, dFlag, &
                         F, Jrow, Jval, Jcol, nJ, x(indS:indE), x(pind), &
                         0_ip, 1_ip, &
                         cu, lencu, iu, leniu, ru, lenru, &
                         cw, lencw, iw, leniw, rw, lenrw )

            ! Get user-defined state equations
            allocate ( C(nC,ndS:ndE) )
            allocate ( Grow(nG), Gval(nG,ndS:ndE) )
            call algcon( Status, usralg, p, nPhs, nC, nY, nU, nP, nNodes, &
                         C, Grow, Gval, Gcol, nG, x(indS:indE), x(pind), &
                         0_ip, 1_ip, &
                         cu, lencu, iu, leniu, ru, lenru, &
                         cw, lencw, iw, leniw, rw, lenrw )

            indS = indE+1

            ! Alg constraints
            do j = ndS, ndE, 2
               do q = Gcol(k+nVar), Gcol(k+nVar+1)-1
                  i = Grow(q)

                  if ( ctype(i,p) == 0 ) then
                     count = count + 1
                     gCon(count) = Gval(q,j)
                  end if
               end do
            end do

            ! ODE Constraint 1 and 2
            do j = ndS, ndE-2, 2
               hstep = step( (j+1)/2 )
               do q = Jcol(k+nVar), Jcol(k+nVar+1)-1
                  i = Jrow(q)

                  if ( ytype(i,p) == 0 ) then
                     count = count + 1
                     gCon(count) = (hstep/6.0)*( Jval(q,j) + &
                          4.0*Jval(q,j+1) + Jval(q,j+2) )

                     count = count + 1
                     gCon(count) = (hstep/8.0)*( Jval(q,j) - Jval(q,j+2) )
                  end if
               end do
            end do

            deallocate ( F, Jrow, Jval )
            deallocate ( C, Grow, Gval )
         end do
      end do

      deallocate ( Jcol, Gcol )

    end subroutine derConHS

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  end subroutine sncEvalHS

end module ct20eval
