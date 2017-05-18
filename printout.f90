!------------------------------------------------------------------
!---                    Subroutine printout                     ---
!------------------------------------------------------------------

!    iprint = 1 ---> Last printout before exit.
!    iprint = 2 ---> Print due a change of Paflov parameter.
!    iprint = 3 ---> Security copy.
!
!  If we use negative values the output files will be written
!  in binary format.


subroutine printout(iprint,namefile,namefile1,             &
                    psi,elem,nx,ny,nz,hx,hy,hz,limp,       &
                    xmax,ymax,zmax,ximp,yimp,zimp,psix,    &
                    paflv,iter)

implicit none
logical                        :: limp
integer   (kind=4), intent(in) :: iter,iprint
integer   (kind=4), intent(in) :: nx,ny,nz
complex   (kind=8), intent(in) :: psi(nx,ny,nz)
real      (kind=8), intent(in) :: psix(nx,ny,nz)
real      (kind=8), intent(in) :: xmax,ymax,zmax,hx,hy,hz
real      (kind=8), intent(in) :: ximp,yimp,zimp,paflv
character (len=3)              :: elem
character (len=60)             :: namefile,namefile1


!-----------------------------------------------
! If iprint < 0 then print binary files...   ---
!-----------------------------------------------
!
if(iprint.lt.0) then
   open(10,file=namefile,form='UNFORMATTED')
   write(10) xmax,ymax,zmax,hx,hy,hz,nx,ny,nz,limp,ximp,yimp,zimp
   write(10) 0.0d0, 0.0d0, 0.0d0
   write(10) psi
   close(10)
   if(limp) then
      open(11,file=namefile1,form='UNFORMATTED')
      write(11) xmax,ymax,zmax,hx,hy,hz,nx,ny,nz,limp,ximp,yimp,zimp
      write(11) psix
      close(11)
   endif
else
!------------------------------------------------------------
! If iprint > 0 then print 'normal files'                 ---
!               the first lines of the file are comments  ---
!               telling you the kind of backup file       ---
!------------------------------------------------------------
!
  open(10,file=namefile,form='FORMATTED')
  if(limp) then
    write(10,1010) iter,paflv,elem
  else
    write(10,1020) iter,paflv
  end if
  select case(iprint)
    case(1)  !........................................ last printout
      write(10,1030)
    case(2)  !............ print due of a change of Paflov parameter
      write(10,1040)
    case(3)  !........................... print due to a backup copy
      write(10,1050)
  end select
!...(Print the density of helium)
  write(10,*) xmax,ymax,zmax,hx,hy,hz,nx,ny,nz,limp,ximp,yimp,zimp
  write(10,*) 0.0d0, 0.0d0, 0.0d0
  write(10,*) psi
  close(10)

!... If (limp.eq.true) then print the wave function
  if(limp) then
     open(10,file=namefile1,form='FORMATTED')
     write(10,1010) iter,paflv,elem
     select case(iprint)
       case(1)  !....................................... last printout
         write(10,1030)
       case(2)  ! print due of a change of Paflov parameter
         write(10,1040)
       case(3)  ! print due to a backup copy
         write(10,1050)
     end select
     write(10,*) xmax,ymax,zmax,hx,hy,hz,nx,ny,nz,limp,ximp,yimp,zimp
     write(10,*) psix
     close(10)
  end if
end if

return

1010 format('#  Density after ',I17,' iterations.',/,   &
            '#  Paflov factor ',0P,F6.3,/,             &
            '#  Impurity      ',A)
1020 format('#  Density after ',I17,' iterations.',/,   &
            '#  Paflov factor ',0P,F6.3)
1030 format('#  Last printout of the run.')
1040 format('#  Change of Paflov parameter.')
1050 format('#  Backup copy.')
end
