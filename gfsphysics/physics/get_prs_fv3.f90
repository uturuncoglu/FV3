module gfs_fv3_needs

   use machine,  only: kind_phys
   use physcons, only: con_fvirt
   use GFS_typedefs, only: GFS_grid_type

!--- public declarations
   public get_prs_fv3, get_phi_fv3

!--- local variables
   real(kind=kind_phys), parameter :: zero = 0.0_kind_phys
   real(kind=kind_phys), parameter :: half = 0.5_kind_phys

contains

   subroutine get_prs_fv3(ix, levs, ntrac, phii, prsi, tgrs, qgrs, del, del_gz)
     integer, intent(in) :: ix, levs, ntrac
     real(kind=kind_phys), dimension(ix,levs+1),     intent(in)    :: phii
     real(kind=kind_phys), dimension(ix,levs+1),     intent(in)    :: prsi
     real(kind=kind_phys), dimension(ix,levs),       intent(in)    :: tgrs
     real(kind=kind_phys), dimension(ix,levs,ntrac), intent(in)    :: qgrs
     real(kind=kind_phys), dimension(ix,levs),       intent(inout) :: del
     real(kind=kind_phys), dimension(ix,levs+1),     intent(inout) :: del_gz

! SJL: Adjust the geopotential height hydrostatically in a way consistent with FV3 discretization
! del_gz is a temp array recording the old info before (t,q) are adjusted
     do k=1,levs
       do i=1,ix
            del(i,k) = prsi(i,k) - prsi(i,k+1)
         del_gz(i,k) = (phii(i,k+1) - phii(i,k)) /                    &
                        (tgrs(i,k)*(1.+con_fvirt*max(zero,qgrs(i,k,1))))
       enddo
     enddo

   end subroutine get_prs_fv3


   subroutine get_phi_fv3(ix, levs, ntrac, gt0, gq0, del_gz, phii, phil,  &
                          Grid, testp)

     integer, intent(in) :: ix, levs, ntrac
     real(kind=kind_phys), dimension(ix,levs),       intent(in)    :: gt0
     real(kind=kind_phys), dimension(ix,levs,ntrac), intent(in)    :: gq0
     real(kind=kind_phys), dimension(ix,levs+1),     intent(inout) :: del_gz
     real(kind=kind_phys), dimension(ix,levs+1),     intent(inout) :: phii
     real(kind=kind_phys), dimension(ix,levs),       intent(inout) :: phil

     type(GFS_grid_type),                            intent(in)    :: Grid
     real(kind=kind_phys),intent(IN) :: testp(2)	! lat/lon in radians

! SJL: Adjust the heighz hydrostatically in a way consistent with FV3 discretization
     do i=1,ix
        phii(i,1) = zero
     enddo
     do k=1,levs
       do i=1,ix
         del_gz(i,k) = del_gz(i,k)*gt0(i,k) *                          &
     &                 (1.+con_fvirt*max(zero,gq0(i,k,1)))
         phii(i,k+1) = phii(i,k) + del_gz(i,k)
         phil(i,k)   = half*(phii(i,k) + phii(i,k+1))

         if (min(phii(i,k),phii(i,k+1)).lt.-1.e11 .or.			&
          abs(grid%xlat(i)-testp(1))+abs(grid%xlon(i)-testp(2))<1.e-4)	&
           print 99,'inside get_phi_fv3  i,k,lat,lon=',i,k,		&
           grid%xlat(i),grid%xlon(i),					&
           'del_gz',del_gz(i,k),'gt0',gt0(i,k),'gq0',gq0(i,k,1),	&
           'phiik',phii(i,k),'phiik+1',phii(i,k+1)
 99      format (a,2i4,2f9.5/4(a9,'=',es10.3))

       enddo
     enddo

!!   do i=1,ix
!!     print 97,'get_phi_fv3',i,minval(phii(i,:)),maxval(phii(i,:))
!!97   format ('inside ',a12,': min/maxval(phii), i=',i3,2es11.3)
!!     if (minval(phii(i,:)).lt.-1.e11) print '(5(i5,es11.3))',		&
!!       (k,phii(i,k),k=1,levs+1)
!!   end do

   end subroutine get_phi_fv3

    logical function NaN(arg)
! --- check whether 'arg' is NaN
    implicit none
    real,intent(IN) :: arg
    character       :: string*9
    write (string,'(es9.1)') arg
    NaN = index(string,'NaN') > 0
    return
    end function NaN
end module gfs_fv3_needs
