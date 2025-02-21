module GFS_restart

  use machine,          only: kind_phys
  use GFS_typedefs,     only: GFS_control_type,  GFS_statein_type,  &
                              GFS_stateout_type, GFS_sfcprop_type,  &
                              GFS_coupling_type, GFS_grid_type,     &
                              GFS_tbd_type,      GFS_cldprop_type,  &
                              GFS_radtend_type,  GFS_diag_type,     &
                              GFS_init_type   

  type var_subtype
    real(kind=kind_phys), pointer :: var2p(:)   => null()  !< 2D data saved in packed format [dim(ix)]
    real(kind=kind_phys), pointer :: var3p(:,:) => null()  !< 3D data saved in packed format [dim(ix,levs)]
  end type var_subtype

  type GFS_restart_type
    integer           :: num2d                    !< current number of registered 2D restart variables
    integer           :: num3d                    !< current number of registered 3D restart variables
    character(len=32), allocatable :: name2d(:)   !< variable name as it will appear in the restart file
    character(len=32), allocatable :: name3d(:)   !< variable name as it will appear in the restart file
    type(var_subtype), allocatable :: data(:,:)   !< holds pointers to data in packed format (allocated to (nblks,max(2d/3dfields))
  end type GFS_restart_type

  public GFS_restart_type, GFS_restart_populate

  CONTAINS
!*******************************************************************************************

!---------------------
! GFS_restart_populate
!---------------------
  subroutine GFS_restart_populate (Restart, Model, Statein, Stateout, Sfcprop, &
                                   Coupling, Grid, Tbd, Cldprop, Radtend, IntDiag, Init_parm)
!----------------------------------------------------------------------------------------!
!   RESTART_METADATA                                                                         !
!     Restart%num2d          [int*4  ]  number of 2D variables to output             !
!     Restart%num3d          [int*4  ]  number of 3D variables to output             !
!     Restart%name2d         [char=32]  variable name in restart file                !
!     Restart%name3d         [char=32]  variable name in restart file                !
!     Restart%fld2d(:,:,:)   [real*8 ]  pointer to 2D data (im,nblks,MAX_RSTRT)      !
!     Restart%fld3d(:,:,:,:) [real*8 ]  pointer to 3D data (im,levs,nblks,MAX_RSTRT) !
!----------------------------------------------------------------------------------------!
    type(GFS_restart_type),     intent(inout) :: Restart
    type(GFS_control_type),     intent(in)    :: Model
    type(GFS_statein_type),     intent(in)    :: Statein(:)
    type(GFS_stateout_type),    intent(in)    :: Stateout(:)
    type(GFS_sfcprop_type),     intent(in)    :: Sfcprop(:)
    type(GFS_coupling_type),    intent(in)    :: Coupling(:)
    type(GFS_grid_type),        intent(in)    :: Grid(:)
    type(GFS_tbd_type),         intent(in)    :: Tbd(:)
    type(GFS_cldprop_type),     intent(in)    :: Cldprop(:)
    type(GFS_radtend_type),     intent(in)    :: Radtend(:)
    type(GFS_diag_type),        intent(in)    :: IntDiag(:)
    type(GFS_init_type),        intent(in)    :: Init_parm

    !--- local variables
    integer :: nblks, num, nb, max_rstrt, offset
    character(len=2) :: c2 = ''
    
    nblks = size(Init_parm%blksz)
    max_rstrt = size(Restart%name2d)

    Restart%num2d = 3 + Model%ntot2d + Model%nctp
    Restart%num3d = Model%ntot3d

    allocate (Restart%name2d(Restart%num2d))
    allocate (Restart%name3d(Restart%num3d))
    allocate (Restart%data(nblks,max(Restart%num2d,Restart%num3d)))

    Restart%name2d(:) = ' '
    Restart%name3d(:) = ' '

    !--- Cldprop variables
    Restart%name2d(1) = 'cv'
    Restart%name2d(2) = 'cvt'
    Restart%name2d(3) = 'cvb'
    do nb = 1,nblks
      Restart%data(nb,1)%var2p => Cldprop(nb)%cv(:)
      Restart%data(nb,2)%var2p => Cldprop(nb)%cvt(:)
      Restart%data(nb,3)%var2p => Cldprop(nb)%cvb(:)
    enddo

    !--- phy_f2d variables
    offset = 3
    do num = 1,Model%ntot2d
       !--- set the variable name
      write(c2,'(i2.2)') num
      Restart%name2d(num+offset) = 'phy_f2d_'//c2
      do nb = 1,nblks
        Restart%data(nb,num+offset)%var2p => Tbd(nb)%phy_f2d(:,num)
      enddo
    enddo

    !--- phy_fctd variables
    offset = offset + Model%ntot2d
    do num = 1, Model%nctp
       !--- set the variable name
      write(c2,'(i2.2)') num
      Restart%name2d(num+offset) = 'phy_fctd_'//c2
      do nb = 1,nblks
        Restart%data(nb,num+offset)%var2p => Tbd(nb)%phy_fctd(:,num)
      enddo
    enddo

    !--- phy_f3d variables
    do num = 1,Model%ntot3d
       !--- set the variable name
      write(c2,'(i2.2)') num
      Restart%name3d(num) = 'phy_f3d_'//c2
      do nb = 1,nblks
        Restart%data(nb,num)%var3p => Tbd(nb)%phy_f3d(:,:,num)
      enddo
    enddo

  end subroutine GFS_restart_populate

end module GFS_restart
