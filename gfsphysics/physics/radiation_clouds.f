!>  \file radiation_clouds.f
!!  This file contains routines to compute cloud related quantities
!!  for radiation computations.
!              module_radiation_clouds description             !!!!!
!  ==========================================================  !!!!!
!                                                                      !
!    the 'radiation_clouds.f' contains:                                !
!                                                                      !
!       'module_radiation_clouds' ---  compute cloud related quantities!
!                for radiation computations                            !
!                                                                      !
!    the following are the externally accessable subroutines:          !
!                                                                      !
!       'cld_init'           --- initialization routine                !
!          inputs:                                                     !
!           ( si, NLAY, imp_physics,  me )                                           !
!          outputs:                                                    !
!           ( none )                                                   !
!                                                                      !
!       'progcld1'           --- zhao/moorthi prognostic cloud scheme  !
!          inputs:                                                     !
!           (plyr,plvl,tlyr,tvly,qlyr,qstl,rhly,clw,                   !
!            xlat,xlon,slmsk,dz,delp,                                  !
!            IX, NLAY, NLP1,                                           !
!            uni_cld, lmfshal, lmfdeep2, cldcov,                       !
!            effrl,effri,effrr,effrs,effr_in,                          !
!          outputs:                                                    !
!            clouds,clds,mtop,mbot,de_lgth)                            !
!                                                                      !
!       'progcld2'           --- ferrier prognostic cloud microphysics !
!          inputs:                                                     !
!           (plyr,plvl,tlyr,tvly,qlyr,qstl,rhly,clw,                   !
!            xlat,xlon,slmsk,dz,delp, f_ice,f_rain,r_rime,flgmin,      !
!            IX, NLAY, NLP1, lmfshal, lmfdeep2,                        !
!          outputs:                                                    !
!            clouds,clds,mtop,mbot,de_lgth)                            !
!                                                                      !
!       'progcld3'           --- zhao/moorthi prognostic cloud + pdfcld!
!          inputs:                                                     !
!           (plyr,plvl,tlyr,tvly,qlyr,qstl,rhly,clw, cnvw,cnvc,        !
!            xlat,xlon,slmsk, dz, delp,                                !
!            ix, nlay, nlp1,                                           !
!            deltaq,sup,kdt,me,                                        !
!          outputs:                                                    !
!            clouds,clds,mtop,mbot,de_lgth)                            !
!                                                                      !
!       'progcld4'           --- gfdl-lin cloud microphysics           !
!          inputs:                                                     !
!           (plyr,plvl,tlyr,tvly,qlyr,qstl,rhly,clw,cnvw,cnvc,         !
!            xlat,xlon,slmsk, dz, delp,                                !
!            ix, nlay, nlp1,                                           !
!          outputs:                                                    !
!            clouds,clds,mtop,mbot,de_lgth)                            !
!                                                                      !
!       'progcld4o'          --- inactive                              !
!                                                                      !
!       'progcld5'           --- thompson/wsm6 cloud microphysics      !
!          inputs:                                                     !
!           (plyr,plvl,tlyr,qlyr,qstl,rhly,clw,                        !
!            xlat,xlon,slmsk, dz, delp,                                !
!            ntrac,ntcw,ntiw,ntrw,ntsw,ntgl,                           !
!            ix, nlay, nlp1,                                           !
!            uni_cld, lmfshal, lmfdeep2, cldcov,                       !
!            re_cloud,re_ice,re_snow,                                  !
!          outputs:                                                    !
!            clouds,clds,mtop,mbot,de_lgth)                            !
!                                                                      !
!       'progclduni'           --- for unified clouds with MG microphys!
!          inputs:                                                     !
!           (plyr,plvl,tlyr,tvly,ccnd,ncnd,                            !
!            xlat,xlon,slmsk,dz,delp, IX, NLAY, NLP1, cldtot,          !
!            effrl,effri,effrr,effrs,effr_in,                          !
!          outputs:                                                    !
!            clouds,clds,mtop,mbot,de_lgth)                            !
!                                                                      !
!    internal accessable only subroutines:                             !
!       'gethml'             --- get diagnostic hi, mid, low clouds    !
!                                                                      !
!                                                                      !
!    cloud array description:                                          !
!          clouds(:,:,1)  -  layer total cloud fraction                !
!          clouds(:,:,2)  -  layer cloud liq water path                !
!          clouds(:,:,3)  -  mean effective radius for liquid cloud    !
!          clouds(:,:,4)  -  layer cloud ice water path                !
!          clouds(:,:,5)  -  mean effective radius for ice cloud       !
!          clouds(:,:,6)  -  layer rain drop water path                !
!          clouds(:,:,7)  -  mean effective radius for rain drop       !
!   **     clouds(:,:,8)  -  layer snow flake water path               !
!          clouds(:,:,9)  -  mean effective radius for snow flake      !
!   ** fu's scheme need to be normalized by snow density (g/m**3/1.0e6)!
!                                                                      !
!    external modules referenced:                                      !
!                                                                      !
!       'module physparam'           in 'physparam.f'                  !
!       'module physcons'            in 'physcons.f'                   !
!       'module module_microphysics' in 'module_bfmicrophysics.f'      !
!                                                                      !
!                                                                      !
! program history log:                                                 !
!      nov 1992,   y.h., k.a.c, a.k. - cloud parameterization          !
!        'cldjms' patterned after slingo and slingo's work (jgr,       !
!        1992), stratiform clouds are allowed in any layer except      !
!        the surface and upper stratosphere. the relative humidity     !
!        criterion may cery in different model layers.                 !
!      mar 1993,   kenneth campana   - created original crhtab tables  !
!        for critical rh look up references.
!      feb 1994,   kenneth campana   - modified to use only one table  !
!        for all forecast hours.                                       !
!      oct 1995,   kenneth campana   - tuned cloud rh curves           !
!        rh-cld relation from tables created using mitchell-hahn       !
!        tuning technique on airforce rtneph observations.             !
!      nov 1995,   kenneth campana   - the bl relationships used       !
!        below llyr, except in marine stratus regions.                 !
!      apr 1996,   kenneth campana   - save bl cld amt in cld(,5)      !
!      aug 1997,   kenneth campana   - smooth out last bunch of bins   !
!        of the rh look up tables                                      !
!      dec 1998,   s. moorthi        - added prognostic cloud method   !
!      apr 2003,   yu-tai hou        - rewritten in frotran 90         !
!        modulized form 'module_rad_clouds' from combining the original!
!        subroutines 'cldjms', 'cldprp', and 'gcljms'. and seperated   !
!        prognostic and diagnostic methods into two packages.          !
!      --- 2003,   s. moorthi        - adapted b. ferrier's prognostic !
!        cloud scheme to ncep gfs model as additional option.          !
!      apr 2004,   yu-tai hou        - separated calculation of the    !
!        averaged h,m,l,bl cloud amounts from each of the cld schemes  !
!        to become an shared individule subprogram 'gethml'.           !
!      may 2004,   yu-tai hou        - rewritten ferrier's scheme as a !
!        separated program 'progcld2' in the cloud module.             !
!      apr 2005,   yu-tai hou        - modified cloud array and module !
!        structures.                                                   !
!      dec 2008,   yu-tai hou        - changed low-cld calculation,    !
!        now cantains clds from sfc layer and upward to the low/mid    !
!        boundary (include bl-cld). h,m,l clds domain boundaries are   !
!        adjusted for better agreement with observations.              !
!      jan 2011,   yu-tai hou        - changed virtual temperature     !
!        as input variable instead of originally computed inside the   !
!        two prognostic cld schemes 'progcld1' and 'progcld2'.         !
!      aug 2012,   yu-tai hou        - modified subroutine cld_init    !
!        to pass all fixed control variables at the start. and set     !
!        their correponding internal module variables to be used by    !
!        module subroutines.                                           !
!      nov 2012,   yu-tai hou        - modified control parameters     !
!        thru module 'physparam'.                                      !
!      apr 2013,   h.lin/y.hou       - corrected error in calculating  !
!        llyr for different vertical indexing directions.              !
!      jul 2013, r. sun/h. pan       - modified to use pdf cloud and   !
!        convective cloud cover and water for radiation                !
!                                                                      !
!      jul 2014 s. moorthi - merging with gfs version                  !
!      feb 2017 a. cheng   - add odepth output, effective radius input !
!      Jan 2018 S Moorthi  - update to include physics from ipdv4      !
!      jun 2018 h-m lin/y-t hou   - removed the legacy subroutine      !
!        'diagcld1' for diagnostic cloud scheme, added new cloud       !
!        overlapping method of de-correlation length, and optimized    !
!        the code structure.                                           !
!                                                                      !
!!!!!  ==========================================================  !!!!!
!!!!!                       end descriptions                       !!!!!
!!!!!  ==========================================================  !!!!!

!> \ingroup rad
!! \defgroup module_radiation_clouds module_radiation_clouds
!! @{
!> This module computes cloud related quantities for radiation
!! computations.
!!
!! Knowledge of cloud properties and their vertical structure is
!! important for meteorological studies due to their impact on both the
!! Earth's radiation budget and adiabatic heating within the atmosphere.
!! Cloud properties in the US National Oceanic and Atmospheric
!! Administration National Centers for Environmental Prediction Global
!! Forecast System (GFS) include (i) cloud liquid/ice water path; (ii)
!! the fraction of clouds; (iii) effective radius of water/ice droplet:
!!\version NCEP-Radiation_clouds    v5.1  Nov 2012
!!
!! This module has three externally accessible subroutines:
!!  - cld_init()           --- initialization routine
!!  - progcld1()           --- zhao/moorthi prognostic cloud scheme 
!!  - progcld2()           --- ferrier prognostic cloud microphysics
!!  - progcld3()           --- zhao/moorthi prognostic cloud + pdfcld
!!  - progcld4()           --- gfdl-lin cloud microphysics
!!  - progcld5()           --- thompson/wsm6 cloud microphysics
!!  - progclduni()         --- unified clouds with MG microphys
!!
!!  and one internally accessable only subroutines:                             
!!  - gethml()             --- get diagnostic hi, mid, low,total,BL clouds    
!!
!> \section gen_al  General Algorithm
!! @{
!! -# Cloud Liquid/Ice Water Path (LWP,IWP)
!!\n We define the fraction of liquid and ice cloud as:
!!\n Fraction of ice cloud (F): \f$F=(273.16K-T)/20\f$
!!\n LWP = total cloud condensate path X (1-F)
!!\n IWP = total clod condensate path X F 
!!
!! -# GFS Cloud Fraction
!! \n The cloud fraction in a given grid box of the GFS model is
!! computed using the parameterization scheme of Xu and Randall(1996)
!! \cite xu_and_randall_1996 :
!!  \f[
!!  \sigma =RH^{k_{1}}\left[1-exp\left(-\frac{k_{2}q_{l}}{\left[\left(1-RH\right)q_{s}\right]^{k_{3}}}\right)\right]
!!  \f]
!!  Where \f$RH\f$ is relative humidity, \f$q_{l}\f$ is the cloud
!!  condensate, \f$q_{s}\f$ is saturation specific humidity,
!! \f$k_{1}(=0.25)\f$, \f$k_{2}(=100)\f$, \f$k_{3}(=0.49)\f$ are the
!! empirical parameters. The cloud condensate is partitioned into
!!  cloud water and ice in radiation based on temperature. Cloud drop
!! effective radius ranges 5-10 microns over land depending on
!!  temperature. Ice crystal radius is function of ice water content
!! (Heymsfield and McFarquhar (1996) \cite heymsfield_and_mcfarquhar_1996).
!!  Maximum-randomly cloud overlapping is used in both long-wave
!! radiation and short-wave radiation. Convective clouds are not
!! considered in radiation.
!!\n
!! -# The parameterization of effective radius of water/ice droplet
!!    (\f$r_{e}\f$)
!>\n Cloud properties in the GFS model are derived from atmospheric
!! condition and cloud condensate amount(NCEP Office Note 441). Cloud
!! condensate information can be provided by one of the choices of
!! different cloud microphysics built in the model. The legacy version
!! of diagnostic cloud scheme in the early GFS has been discontinued.
!! \n For the parameterization of effective radius,\f$r_{ew}\f$, of
!! water droplet, we fix \f$r_{ew}\f$ to a value of \f$10\mu m\f$ over
!! the oceans. Over the land, \f$\f$ is defined as:
!!\f[
!! r_{ew} = 5+5\times F
!!\f]
!! Thus, the effective radius of cloud water droplets will reach to a
!! minimum values of \f$5\mu m\f$ when F=0, and to a maximum value of 
!! \f$10\mu m\f$ when the ice fraction is increasing.
!! \n For ice clouds, following Heymsfield and McFarquhar (1996)
!!  \cite heymsfield_and_mcfarquhar_1996,
!! we have made the effective ice droplet radius to be an empirical
!! function of ice water concentration (IWC) and environmental temperature as:
!! \f[
!! r_{ei}=\begin{cases}(1250/9.917)IWC^{0.109} & T <-50^0C \\(1250/9.337)IWC^{0.080} & -50^0C \leq T<-40^0C\\(1250/9.208)IWC^{0.055} & -40^0C\leq T <-30^0C\\(1250/9.387)IWC^{0.031} & -30^0C \leq T\end{cases}
!! \f]
!! where IWC and IWP satisfy:
!! \f[
!! IWP_{\triangle Z}=\int_{\triangle Z} IWCdZ
!! \f]
!! @}
!========================================!
      module module_radiation_clouds     !
!........................................!
!
      use physparam,           only : icldflg, iovrsw, iovrlw,          &
     &                                lcrick, lcnorm, lnoprec,          &
     &                                ivflip, kind_phys, kind_io4
      use physcons,            only : con_fvirt, con_ttp, con_rocp,     &
     &                                con_t0c, con_pi, con_g, con_rd,   &
     &                                con_thgni
      use module_microphysics, only : rsipath2
      use module_iounitdef,    only : NICLTUN
!
      implicit   none
!
      private

!  ---  version tag and last revision date
      character(40), parameter ::                                       &
     &   VTAGCLD='NCEP-Radiation_clouds    v5.1  Nov 2012 '
!    &   VTAGCLD='NCEP-Radiation_clouds    v5.0  Aug 2012 '

!  ---  set constant parameters
      real (kind=kind_phys), parameter :: gfac=1.0e5/con_g              &
     &,                                   gord=con_g/con_rd
!> number of fields in cloud array
      integer, parameter, public :: NF_CLDS = 9
!> number of cloud vertical domains
      integer, parameter, public :: NK_CLDS = 3  

!> pressure limits of cloud domain interfaces (low,mid,high) in mb (0.1kPa)
      real (kind=kind_phys), save :: ptopc(NK_CLDS+1,2)

!org  data ptopc / 1050., 642., 350., 0.0,  1050., 750., 500., 0.0 /
      data ptopc / 1050., 650., 400., 0.0,  1050., 750., 500., 0.0 /

!     real (kind=kind_phys), parameter :: climit = 0.01
      real (kind=kind_phys), parameter :: climit = 0.001, climit2=0.05
      real (kind=kind_phys), parameter :: ovcst  = 1.0 - 1.0e-8

!> default liq radius to 10 micron
      real (kind=kind_phys), parameter :: reliq_def = 10.0    
!> default ice radius to 50 micron
      real (kind=kind_phys), parameter :: reice_def = 50.0   
!> default rain radius to 1000 micron
      real (kind=kind_phys), parameter :: rrain_def = 1000.0 
!> default snow radius to 250 micron
      real (kind=kind_phys), parameter :: rsnow_def = 250.0 

!> default cld single scat albedo
      real (kind=kind_phys), parameter :: cldssa_def = 0.99  
!> default cld asymmetry factor
      real (kind=kind_phys), parameter :: cldasy_def = 0.84 

!> upper limit of boundary layer clouds
      integer  :: llyr   = 2          
!> maximum-random cloud overlapping method
      integer  :: iovr   = 1           

      public progcld1, progcld2, progcld3, progcld4, progclduni,        &
     &                 cld_init, progcld5, progcld4o


! =================
      contains
! =================


!> This subroutine is an initialization program for cloud-radiation
!! calculations and sets up boundary layer cloud top.
!!\param si              model vertical sigma layer interface
!!\param NLAY            vertical layer number
!!\param me              print control flag
!!\section gen_cld_init General Algorithm
!> @{
!-----------------------------------
      subroutine cld_init                                               &
     &     ( si, NLAY, imp_physics, me ) !  ---  inputs
!  ---  outputs:
!          ( none )

!  ===================================================================  !
!                                                                       !
! abstract: cld_init is an initialization program for cloud-radiation   !
!   calculations. it sets up boundary layer cloud top.                  !
!                                                                       !
!                                                                       !
! inputs:                                                               !
!   si (L+1)        : model vertical sigma layer interface              !
!   NLAY            : vertical layer number                             !
!   imp_physics     : MP identifier                                     !
!   me              : print control flag                                !
!                                                                       !
!  outputs: (none)                                                      !
!           to module variables                                         !
!                                                                       !
!  external module variables: (in physparam)                            !
!   icldflg         : cloud optical property scheme control flag        !
!                     =0: abort! diagnostic cloud method discontinued   !
!                     =1: model use prognostic cloud method             !
!   imp_physics         : cloud microphysics scheme control flag        !
!                     =99: zhao/carr/sundqvist microphysics cloud       !
!                     =98: zhao/carr/sundqvist microphysics cloud+pdfcld!
!                     =11: GFDL microphysics cloud                      !
!                     =8: Thompson microphysics                         !
!                     =6: WSM6 microphysics                             !
!                     =10: MG microphysics                              !
!   iovrsw/iovrlw   : sw/lw control flag for cloud overlapping scheme   !
!                     =0: random overlapping clouds                     !
!                     =1: max/ran overlapping clouds                    !
!                     =2: maximum overlap clouds       (mcica only)     !
!                     =3: decorrelation-length overlap (mcica only)     !
!   ivflip          : control flag for direction of vertical index      !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!  usage:       call cld_init                                           !
!                                                                       !
!  subroutines called:    rhtable                                       !
!                                                                       !
!  ===================================================================  !
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: NLAY, me, imp_physics

      real (kind=kind_phys), intent(in) :: si(:)

!  ---  outputs: (none)

!  ---  locals:
      integer :: k, kl, ier

!
!===> ...  begin here
!
!  ---  set up module variables

      iovr    = max( iovrsw, iovrlw )    !cld ovlp used for diag HML cld output

      if (me == 0) print *, VTAGCLD      !print out version tag

      if ( icldflg == 0 ) then
        print *,' - Diagnostic Cloud Method has been discontinued'
        stop

      else
        if (me == 0) then
          print *,' - Using Prognostic Cloud Method'
          if (imp_physics == 99) then
            print *,'   --- Zhao/Carr/Sundqvist microphysics'
          elseif (imp_physics == 98) then
            print *,'   --- zhao/carr/sundqvist + pdf cloud'
          elseif (imp_physics == 11) then
            print *,'   --- GFDL Lin cloud microphysics'
          elseif (imp_physics == 8) then
            print *,'   --- Thompson cloud microphysics'
          elseif (imp_physics == 6) then
            print *,'   --- WSM6 cloud microphysics'
          elseif (imp_physics == 10) then
            print *,'   --- MG cloud microphysics'
          else
            print *,'  !!! ERROR in cloud microphysc specification!!!', &
     &              '  imp_physics (NP3D) =',imp_physics
            stop
          endif
        endif
      endif

!> -# Compute the top of BL cld (llyr), which is the topmost non 
!!    cld(low) layer for stratiform (at or above lowest 0.1 of the
!!     atmosphere).

      if ( ivflip == 0 ) then    ! data from toa to sfc
        lab_do_k0 : do k = NLAY, 2, -1
          kl = k
          if (si(k) < 0.9e0) exit lab_do_k0
        enddo  lab_do_k0

        llyr = kl
      else                      ! data from sfc to top
        lab_do_k1 : do k = 2, NLAY
          kl = k
          if (si(k) < 0.9e0) exit lab_do_k1
        enddo  lab_do_k1

        llyr = kl - 1
      endif                     ! end_if_ivflip

!
      return
!...................................
      end subroutine cld_init
!-----------------------------------
!> @}

!> This subroutine computes cloud related quantities using 
!! zhao/moorthi's prognostic cloud microphysics scheme.
!!\param plyr    (IX,NLAY), model layer mean pressure in mb (100Pa)
!!\param plvl    (IX,NLP1), model level pressure in mb (100Pa)
!!\param tlyr    (IX,NLAY), model layer mean temperature in K
!!\param tvly    (IX,NLAY), model layer virtual temperature in K
!!\param qlyr    (IX,NLAY), layer specific humidity in gm/gm
!!\param qstl    (IX,NLAY), layer saturate humidity in gm/gm
!!\param rhly    (IX,NLAY), layer relative humidity \f$ (=qlyr/qstl) \f$
!!\param clw     (IX,NLAY), layer cloud condensate amount
!!\param xlat    (IX), grid latitude in radians, default to pi/2 -> 
!!               -pi/2 range, otherwise see in-line comment
!!\param xlon    (IX), grid longitude in radians  (not used)
!!\param slmsk   (IX), sea/land mask array (sea:0,land:1,sea-ice:2)
!!\param dz      (IX,NLAY), layer thickness (km)
!!\param delp    (IX,NLAY), model layer pressure thickness in mb (100Pa)
!!\param IX           horizontal dimention
!!\param NLAY,NLP1    vertical layer/level dimensions
!!\param clouds      (IX,NLAY,NF_CLDS), cloud profiles
!!\n                 (:,:,1) - layer total cloud fraction
!!\n                 (:,:,2) - layer cloud liq water path \f$(g/m^2)\f$
!!\n                 (:,:,3) - mean eff radius for liq cloud (micron)
!!\n                 (:,:,4) - layer cloud ice water path \f$(g/m^2)\f$
!!\n                 (:,:,5) - mean eff radius for ice cloud (micron)
!!\n                 (:,:,6) - layer rain drop water path (not assigned)
!!\n                 (:,:,7) - mean eff radius for rain drop (micron)
!!\n                 (:,:,8) - layer snow flake water path (not assigned)
!!\n                 (:,:,9) - mean eff radius for snow flake (micron)
!!\n  *** fu's scheme need to be normalized by snow density \f$ (g/m^3/1.0e6)\f$
!!\param clds       (IX,5), fraction of clouds for low, mid, hi, tot, bl
!!\param mtop       (IX,3), vertical indices for low, mid, hi cloud tops
!!\param mbot       (IX,3), vertical indices for low, mid, hi cloud bases
!!\param de_lgth    (IX),   clouds decorrelation length (km)
!>\section gen_progcld1 General Algorithm
!> @{
!-----------------------------------
      subroutine progcld1                                               &
     &     ( plyr,plvl,tlyr,tvly,qlyr,qstl,rhly,clw,                    &    !  ---  inputs:
     &       xlat,xlon,slmsk,dz,delp, IX, NLAY, NLP1,                   &
     &       uni_cld, lmfshal, lmfdeep2, cldcov,                        &
     &       effrl,effri,effrr,effrs,effr_in,                           &
     &       clouds,clds,mtop,mbot,de_lgth                              &    !  ---  outputs:
     &      )

! =================   subprogram documentation block   ================ !
!                                                                       !
! subprogram:    progcld1    computes cloud related quantities using    !
!   zhao/moorthi's prognostic cloud microphysics scheme.                !
!                                                                       !
! abstract:  this program computes cloud fractions from cloud           !
!   condensates, calculates liquid/ice cloud droplet effective radius,  !
!   and computes the low, mid, high, total and boundary layer cloud     !
!   fractions and the vertical indices of low, mid, and high cloud      !
!   top and base.  the three vertical cloud domains are set up in the   !
!   initial subroutine "cld_init".                                      !
!                                                                       !
! usage:         call progcld1                                          !
!                                                                       !
! subprograms called:   gethml                                          !
!                                                                       !
! attributes:                                                           !
!   language:   fortran 90                                              !
!   machine:    ibm-sp, sgi                                             !
!                                                                       !
!                                                                       !
!  ====================  definition of variables  ====================  !
!                                                                       !
! input variables:                                                      !
!   plyr  (IX,NLAY) : model layer mean pressure in mb (100Pa)           !
!   plvl  (IX,NLP1) : model level pressure in mb (100Pa)                !
!   tlyr  (IX,NLAY) : model layer mean temperature in k                 !
!   tvly  (IX,NLAY) : model layer virtual temperature in k              !
!   qlyr  (IX,NLAY) : layer specific humidity in gm/gm                  !
!   qstl  (IX,NLAY) : layer saturate humidity in gm/gm                  !
!   rhly  (IX,NLAY) : layer relative humidity (=qlyr/qstl)              !
!   clw   (IX,NLAY) : layer cloud condensate amount                     !
!   xlat  (IX)      : grid latitude in radians, default to pi/2 -> -pi/2!
!                     range, otherwise see in-line comment              !
!   xlon  (IX)      : grid longitude in radians  (not used)             !
!   slmsk (IX)      : sea/land mask array (sea:0,land:1,sea-ice:2)      !
!   dz    (ix,nlay) : layer thickness (km)                              !
!   delp  (ix,nlay) : model layer pressure thickness in mb (100Pa)      !
!   IX              : horizontal dimention                              !
!   NLAY,NLP1       : vertical layer/level dimensions                   !
!   uni_cld         : logical - true for cloud fraction from shoc       !
!   lmfshal         : logical - true for mass flux shallow convection   !
!   lmfdeep2        : logical - true for mass flux deep convection      !
!   cldcov          : layer cloud fraction (used when uni_cld=.true.    !
!                                                                       !
! output variables:                                                     !
!   clouds(IX,NLAY,NF_CLDS) : cloud profiles                            !
!      clouds(:,:,1) - layer total cloud fraction                       !
!      clouds(:,:,2) - layer cloud liq water path         (g/m**2)      !
!      clouds(:,:,3) - mean eff radius for liq cloud      (micron)      !
!      clouds(:,:,4) - layer cloud ice water path         (g/m**2)      !
!      clouds(:,:,5) - mean eff radius for ice cloud      (micron)      !
!      clouds(:,:,6) - layer rain drop water path         not assigned  !
!      clouds(:,:,7) - mean eff radius for rain drop      (micron)      !
!  *** clouds(:,:,8) - layer snow flake water path        not assigned  !
!      clouds(:,:,9) - mean eff radius for snow flake     (micron)      !
!  *** fu's scheme need to be normalized by snow density (g/m**3/1.0e6) !
!   clds  (IX,5)    : fraction of clouds for low, mid, hi, tot, bl      !
!   mtop  (IX,3)    : vertical indices for low, mid, hi cloud tops      !
!   mbot  (IX,3)    : vertical indices for low, mid, hi cloud bases     !
!   de_lgth(ix)     : clouds decorrelation length (km)                  !
!                                                                       !
! module variables:                                                     !
!   ivflip          : control flag of vertical index direction          !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!   lmfshal         : mass-flux shallow conv scheme flag                !
!   lmfdeep2        : scale-aware mass-flux deep conv scheme flag       !
!   lcrick          : control flag for eliminating CRICK                !
!                     =t: apply layer smoothing to eliminate CRICK      !
!                     =f: do not apply layer smoothing                  !
!   lcnorm          : control flag for in-cld condensate                !
!                     =t: normalize cloud condensate                    !
!                     =f: not normalize cloud condensate                !
!                                                                       !
!  ====================    end of description    =====================  !
!
      implicit none

!  ---  inputs
      integer,  intent(in) :: IX, NLAY, NLP1

      logical, intent(in)  :: uni_cld, lmfshal, lmfdeep2, effr_in

      real (kind=kind_phys), dimension(:,:), intent(in) :: plvl, plyr,  &
     &       tlyr,  tvly,  qlyr,  qstl, rhly, clw, cldcov, delp, dz,    &
     &       effrl, effri, effrr, effrs

      real (kind=kind_phys), dimension(:),   intent(in) :: xlat, xlon,  &
     &       slmsk

!  ---  outputs
      real (kind=kind_phys), dimension(:,:,:), intent(out) :: clouds

      real (kind=kind_phys), dimension(:,:),   intent(out) :: clds
      real (kind=kind_phys), dimension(:),     intent(out) :: de_lgth

      integer,               dimension(:,:),   intent(out) :: mtop,mbot

!  ---  local variables:
      real (kind=kind_phys), dimension(IX,NLAY) :: cldtot, cldcnv,      &
     &       cwp, cip, crp, csp, rew, rei, res, rer, tem2d, clwf

      real (kind=kind_phys) :: ptop1(IX,NK_CLDS+1), rxlat(ix)

      real (kind=kind_phys) :: clwmin, clwm, clwt, onemrh, value,       &
     &       tem1, tem2, tem3

      integer :: i, k, id, nf

!  ---  constant values
!     real (kind=kind_phys), parameter :: xrc3 = 200.
      real (kind=kind_phys), parameter :: xrc3 = 100.

!
!===> ... begin here
!
      do nf=1,nf_clds
        do k=1,nlay
          do i=1,ix
            clouds(i,k,nf) = 0.0
          enddo
        enddo
      enddo
!     clouds(:,:,:) = 0.0

      if(effr_in) then
        do k = 1, NLAY
          do i = 1, IX
            cldtot(i,k) = 0.0
            cldcnv(i,k) = 0.0
            cwp   (i,k) = 0.0
            cip   (i,k) = 0.0
            crp   (i,k) = 0.0
            csp   (i,k) = 0.0
            rew   (i,k) = effrl (i,k)
            rei   (i,k) = effri (i,k)
            rer   (i,k) = effrr (i,k)
            res   (i,k) = effrs (i,k)
            tem2d (i,k) = min(1.0, max(0.0,(con_ttp-tlyr(i,k))*0.05))
            clwf(i,k)   = 0.0
          enddo
        enddo
      else
        do k = 1, NLAY
          do i = 1, IX
            cldtot(i,k) = 0.0
            cldcnv(i,k) = 0.0
            cwp   (i,k) = 0.0
            cip   (i,k) = 0.0
            crp   (i,k) = 0.0
            csp   (i,k) = 0.0
            rew   (i,k) = reliq_def            ! default liq radius to 10 micron
            rei   (i,k) = reice_def            ! default ice radius to 50 micron
            rer   (i,k) = rrain_def            ! default rain radius to 1000 micron
            res   (i,k) = rsnow_def            ! default snow radius to 250 micron
            tem2d (i,k) = min(1.0, max(0.0, (con_ttp-tlyr(i,k))*0.05))
            clwf(i,k)   = 0.0
          enddo
        enddo
      endif
!
      if ( lcrick ) then
        do i = 1, IX
          clwf(i,1)    = 0.75*clw(i,1)    + 0.25*clw(i,2)
          clwf(i,nlay) = 0.75*clw(i,nlay) + 0.25*clw(i,nlay-1)
        enddo
        do k = 2, NLAY-1
          do i = 1, IX
            clwf(i,K) = 0.25*clw(i,k-1) + 0.5*clw(i,k) + 0.25*clw(i,k+1)
          enddo
        enddo
      else
        do k = 1, NLAY
          do i = 1, IX
            clwf(i,k) = clw(i,k)
          enddo
        enddo
      endif

!> -# Find top pressure for each cloud domain for given latitude.
!     ptopc(k,i): top presure of each cld domain (k=1-4 are sfc,L,m,h;
!  ---  i=1,2 are low-lat (<45 degree) and pole regions)

      do i =1, IX
        rxlat(i) = abs( xlat(i) / con_pi )      ! if xlat in pi/2 -> -pi/2 range
!       rxlat(i) = abs(0.5 - xlat(i)/con_pi)    ! if xlat in 0 -> pi range
      enddo

      do id = 1, 4
        tem1 = ptopc(id,2) - ptopc(id,1)

        do i =1, IX
          ptop1(i,id) = ptopc(id,1) + tem1*max( 0.0, 4.0*rxlat(i)-1.0 )
        enddo
      enddo

!> -# Compute cloud liquid/ice condensate path in \f$ g/m^2 \f$ .

        do k = 1, NLAY
          do i = 1, IX
            clwt     = max(0.0, clwf(i,k)) * gfac * delp(i,k)
            cip(i,k) = clwt * tem2d(i,k)
            cwp(i,k) = clwt - cip(i,k)
          enddo
        enddo

!> -# Compute effective liquid cloud droplet radius over land.

      if(.not. effr_in) then
        do i = 1, IX
          if (nint(slmsk(i)) == 1) then
            do k = 1, NLAY
              rew(i,k) = 5.0 + 5.0 * tem2d(i,k)
            enddo
          endif
        enddo
      endif

      if (uni_cld) then     ! use unified sgs clouds generated outside
        do k = 1, NLAY
          do i = 1, IX
            cldtot(i,k) = cldcov(i,k)
          enddo
        enddo

      else

!> -# Calculate layer cloud fraction.

        clwmin = 0.0
        if (.not. lmfshal) then
          do k = 1, NLAY
          do i = 1, IX
            clwt = 1.0e-6 * (plyr(i,k)*0.001)
!           clwt = 2.0e-6 * (plyr(i,k)*0.001)

            if (clwf(i,k) > clwt) then

              onemrh= max( 1.e-10, 1.0-rhly(i,k) )
              clwm  = clwmin / max( 0.01, plyr(i,k)*0.001 )

              tem1  = min(max(sqrt(sqrt(onemrh*qstl(i,k))),0.0001),1.0)
              tem1  = 2000.0 / tem1

!             tem1  = 1000.0 / tem1

              value = max( min( tem1*(clwf(i,k)-clwm), 50.0 ), 0.0 )
              tem2  = sqrt( sqrt(rhly(i,k)) )

              cldtot(i,k) = max( tem2*(1.0-exp(-value)), 0.0 )
            endif
          enddo
          enddo
        else
          do k = 1, NLAY
          do i = 1, IX
            clwt = 1.0e-6 * (plyr(i,k)*0.001)
!           clwt = 2.0e-6 * (plyr(i,k)*0.001)

            if (clwf(i,k) > clwt) then
              onemrh= max( 1.e-10, 1.0-rhly(i,k) )
              clwm  = clwmin / max( 0.01, plyr(i,k)*0.001 )
!
              tem1  = min(max((onemrh*qstl(i,k))**0.49,0.0001),1.0)  !jhan
              if (lmfdeep2) then
                tem1  = xrc3 / tem1
              else
                tem1  = 100.0 / tem1
              endif
!
              value = max( min( tem1*(clwf(i,k)-clwm), 50.0 ), 0.0 )
              tem2  = sqrt( sqrt(rhly(i,k)) )

              cldtot(i,k) = max( tem2*(1.0-exp(-value)), 0.0 )
            endif
          enddo
          enddo
        endif

      endif                                ! if (uni_cld) then

      do k = 1, NLAY
        do i = 1, IX
          if (cldtot(i,k) < climit) then
            cldtot(i,k) = 0.0
            cwp(i,k)    = 0.0
            cip(i,k)    = 0.0
            crp(i,k)    = 0.0
            csp(i,k)    = 0.0
          endif
        enddo
      enddo

      if ( lcnorm ) then
        do k = 1, NLAY
          do i = 1, IX
            if (cldtot(i,k) >= climit) then
              tem1 = 1.0 / max(climit2, cldtot(i,k))
              cwp(i,k) = cwp(i,k) * tem1
              cip(i,k) = cip(i,k) * tem1
              crp(i,k) = crp(i,k) * tem1
              csp(i,k) = csp(i,k) * tem1
            endif
          enddo
        enddo
      endif

!> -# Compute effective ice cloud droplet radius following Heymsfield 
!!    and McFarquhar (1996) \cite heymsfield_and_mcfarquhar_1996.

      if(.not.effr_in) then
        do k = 1, NLAY
          do i = 1, IX
            tem2 = tlyr(i,k) - con_ttp
  
            if (cip(i,k) > 0.0) then
              tem3 = gord * cip(i,k) * plyr(i,k) / (delp(i,k)*tvly(i,k))

              if (tem2 < -50.0) then
                rei(i,k) = (1250.0/9.917) * tem3 ** 0.109
              elseif (tem2 < -40.0) then
                rei(i,k) = (1250.0/9.337) * tem3 ** 0.08
              elseif (tem2 < -30.0) then
                rei(i,k) = (1250.0/9.208) * tem3 ** 0.055
              else
                rei(i,k) = (1250.0/9.387) * tem3 ** 0.031
              endif
!             rei(i,k)   = max(20.0, min(rei(i,k), 300.0))
!             rei(i,k)   = max(10.0, min(rei(i,k), 100.0))
              rei(i,k)   = max(10.0, min(rei(i,k), 150.0))
!             rei(i,k)   = max(5.0,  min(rei(i,k), 130.0))
            endif
          enddo
        enddo
      endif

!
      do k = 1, NLAY
        do i = 1, IX
          clouds(i,k,1) = cldtot(i,k)
          clouds(i,k,2) = cwp(i,k)
          clouds(i,k,3) = rew(i,k)
          clouds(i,k,4) = cip(i,k)
          clouds(i,k,5) = rei(i,k)
!         clouds(i,k,6) = 0.0
          clouds(i,k,7) = rer(i,k)
!         clouds(i,k,8) = 0.0
          clouds(i,k,9) = res(i,k)
        enddo
      enddo

!  --- ...  estimate clouds decorrelation length in km
!           this is only a tentative test, need to consider change later

      if ( iovr == 3 ) then
        do i = 1, ix
          de_lgth(i) = max( 0.6, 2.78-4.6*rxlat(i) )
        enddo
      endif

!> -# Call gethml() to compute low,mid,high,total, and boundary layer
!!    cloud fractions and clouds top/bottom layer indices for low, mid,
!!    and high clouds.
!  ---  compute low, mid, high, total, and boundary layer cloud fractions
!       and clouds top/bottom layer indices for low, mid, and high clouds.
!       The three cloud domain boundaries are defined by ptopc.  The cloud
!       overlapping method is defined by control flag 'iovr', which may
!       be different for lw and sw radiation programs.

      call gethml                                                       &
!  ---  inputs:
     &     ( plyr, ptop1, cldtot, cldcnv, dz, de_lgth,                  &
     &       IX,NLAY,                                                   &
!  ---  outputs:
     &       clds, mtop, mbot                                           &
     &     )


!
      return
!...................................
      end subroutine progcld1
!-----------------------------------
!> @}

!> This subroutine computes cloud related quantities using ferrier's
!! prognostic cloud microphysics scheme.
!!\param plyr    (IX,NLAY), model layer mean pressure in mb (100Pa)
!!\param plvl    (IX,NLP1), model level pressure in mb (100Pa)
!!\param tlyr    (IX,NLAY), model layer mean temperature in K
!!\param tvly    (IX,NLAY), model layer virtual temperature in K
!!\param qlyr    (IX,NLAY), layer specific humidity in gm/gm
!!\param qstl    (IX,NLAY), layer saturate humidity in gm/gm
!!\param rhly    (IX,NLAY), layer relative humidity (=qlyr/qstl)
!!\param clw     (IX,NLAY), layer cloud condensate amount
!!\param f_ice   (IX,NLAY), fraction of layer cloud ice  (ferrier micro-phys)
!!\param f_rain  (IX,NLAY), fraction of layer rain water (ferrier micro-phys)
!!\param r_rime  (IX,NLAY), mass ratio of total ice to unrimed ice (>=1)
!!\param flgmin  (IX), minimum large ice fraction
!!\param xlat    (IX), grid latitude in radians, default to pi/2 -> 
!!               -pi/2 range, otherwise see in-line comment
!!\param xlon    (IX), grid longitude in radians  (not used)
!!\param slmsk   (IX), sea/land mask array (sea:0,land:1,sea-ice:2)
!!\param dz      (IX,NLAY), layer thickness (km)
!!\param delp    (IX,NLAY), model layer pressure thickness in mb (100Pa)
!!\param IX      horizontal dimention
!!\param NLAY,NLP1    vertical layer/level dimensions
!!\param clouds      (IX,NLAY,NF_CLDS), cloud profiles
!!\n                 (:,:,1) - layer total cloud fraction
!!\n                 (:,:,2) - layer cloud liq water path  \f$(g/m^2)\f$
!!\n                 (:,:,3) - mean eff radius for liq cloud (micron)
!!\n                 (:,:,4) - layer cloud ice water path  \f$(g/m^2)\f$
!!\n                 (:,:,5) - mean eff radius for ice cloud (micron)
!!\n                 (:,:,6) - layer rain drop water path  \f$(g/m^2)\f$
!!\n                 (:,:,7) - mean eff radius for rain drop (micron)
!!\n                 (:,:,8) - layer snow flake water path \f$(g/m^2)\f$
!!\n                 (:,:,9) - mean eff radius for snow flake (micron)
!!\n  *** fu's scheme need to be normalized by snow density (g/m**3/1.0e6)
!!\param clds  (IX,5), fraction of clouds for low, mid, hi, tot, bl
!!\param mtop  (IX,3), vertical indices for low, mid, hi cloud tops
!!\param mbot  (IX,3), vertical indices for low, mid, hi cloud bases
!!\param de_lgth   (IX),   clouds decorrelation length (km)
!>\section gen_progcld2 General Algorithm
!> @{
!-----------------------------------
      subroutine progcld2                                               &
     &     ( plyr,plvl,tlyr,tvly,qlyr,qstl,rhly,clw,                    &    !  ---  inputs:
     &       xlat,xlon,slmsk,dz,delp, f_ice,f_rain,r_rime,flgmin,       &
     &       IX, NLAY, NLP1, lmfshal, lmfdeep2,                         &
     &       clouds,clds,mtop,mbot,de_lgth                              &    !  ---  outputs:
     &      )

! =================   subprogram documentation block   ================ !
!                                                                       !
! subprogram:    progcld2    computes cloud related quantities using    !
!   ferrier's prognostic cloud microphysics scheme.                     !
!                                                                       !
! abstract:  this program computes cloud fractions from cloud           !
!   condensates, calculates liquid/ice cloud droplet effective radius,  !
!   and computes the low, mid, high, total and boundary layer cloud     !
!   fractions and the vertical indices of low, mid, and high cloud      !
!   top and base.  the three vertical cloud domains are set up in the   !
!   initial subroutine "cld_init".                                      !
!                                                                       !
! usage:         call progcld2                                          !
!                                                                       !
! subprograms called:   gethml                                          !
!                                                                       !
! attributes:                                                           !
!   language:   fortran 90                                              !
!   machine:    ibm-sp, sgi                                             !
!                                                                       !
!                                                                       !
!  ====================  definition of variables  ====================  !
!                                                                       !
! input variables:                                                      !
!   plyr  (IX,NLAY) : model layer mean pressure in mb (100Pa)           !
!   plvl  (IX,NLP1) : model level pressure in mb (100Pa)                !
!   tlyr  (IX,NLAY) : model layer mean temperature in k                 !
!   tvly  (IX,NLAY) : model layer virtual temperature in k              !
!   qlyr  (IX,NLAY) : layer specific humidity in gm/gm                  !
!   qstl  (IX,NLAY) : layer saturate humidity in gm/gm                  !
!   rhly  (IX,NLAY) : layer relative humidity (=qlyr/qstl)              !
!   clw   (IX,NLAY) : layer cloud condensate amount                     !
!   f_ice (IX,NLAY) : fraction of layer cloud ice  (ferrier micro-phys) !
!   f_rain(IX,NLAY) : fraction of layer rain water (ferrier micro-phys) !
!   r_rime(IX,NLAY) : mass ratio of total ice to unrimed ice (>=1)      !
!   flgmin(IX)      : minimim large ice fraction                        !
!   xlat  (IX)      : grid latitude in radians, default to pi/2 -> -pi/2!
!                     range, otherwise see in-line comment              !
!   xlon  (IX)      : grid longitude in radians  (not used)             !
!   slmsk (IX)      : sea/land mask array (sea:0,land:1,sea-ice:2)      !
!   dz    (ix,nlay) : layer thickness (km)                              !
!   delp  (ix,nlay) : model layer pressure thickness in mb (100Pa)      !
!   IX              : horizontal dimention                              !
!   NLAY,NLP1       : vertical layer/level dimensions                   !
!                                                                       !
! output variables:                                                     !
!   clouds(IX,NLAY,NF_CLDS) : cloud profiles                            !
!      clouds(:,:,1) - layer total cloud fraction                       !
!      clouds(:,:,2) - layer cloud liq water path         (g/m**2)      !
!      clouds(:,:,3) - mean eff radius for liq cloud      (micron)      !
!      clouds(:,:,4) - layer cloud ice water path         (g/m**2)      !
!      clouds(:,:,5) - mean eff radius for ice cloud      (micron)      !
!      clouds(:,:,6) - layer rain drop water path         (g/m**2)      !
!      clouds(:,:,7) - mean eff radius for rain drop      (micron)      !
!  *** clouds(:,:,8) - layer snow flake water path        (g/m**2)      !
!      clouds(:,:,9) - mean eff radius for snow flake     (micron)      !
!  *** fu's scheme need to be normalized by snow density (g/m**3/1.0e6) !
!   clds  (IX,5)    : fraction of clouds for low, mid, hi, tot, bl      !
!   mtop  (IX,3)    : vertical indices for low, mid, hi cloud tops      !
!   mbot  (IX,3)    : vertical indices for low, mid, hi cloud bases     !
!   de_lgth(ix)     : clouds decorrelation length (km)                  !
!                                                                       !
! external module variables:                                            !
!   ivflip          : control flag of vertical index direction          !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!   lmfshal         : mass-flux shallow conv scheme flag                !
!   lmfdeep2        : scale-aware mass-flux deep conv scheme flag       !
!   lcrick          : control flag for eliminating CRICK                !
!                     =t: apply layer smoothing to eliminate CRICK      !
!                     =f: do not apply layer smoothing                  !
!   lcnorm          : control flag for in-cld condensate                !
!                     =t: normalize cloud condensate                    !
!                     =f: not normalize cloud condensate                !
!   lnoprec         : precip effect in radiation flag (ferrier scheme)  !
!                     =t: snow/rain has no impact on radiation          !
!                     =f: snow/rain has impact on radiation             !
!                                                                       !
!  ====================    end of description    =====================  !
!
      implicit none

!  ---  constants

!  ---  inputs
      integer,  intent(in) :: IX, NLAY, NLP1

      logical, intent(in)  :: lmfshal, lmfdeep2

      real (kind=kind_phys), dimension(:,:), intent(in) :: plvl, plyr,  &
     &       tlyr, tvly, qlyr, qstl, rhly, clw, f_ice, f_rain, r_rime,  &
     &       dz, delp

      real (kind=kind_phys), dimension(:),   intent(in) :: xlat, xlon,  &
     &       slmsk
      real (kind=kind_phys), dimension(:), intent(in) :: flgmin

!  ---  outputs
      real (kind=kind_phys), dimension(:,:,:), intent(out) :: clouds

      real (kind=kind_phys), dimension(:,:),   intent(out) :: clds
      real (kind=kind_phys), dimension(:),     intent(out) :: de_lgth

      integer,               dimension(:,:),   intent(out) :: mtop,mbot

!  ---  local variables:
      real (kind=kind_phys), dimension(IX,NLAY) :: cldtot, cldcnv,      &
     &       cwp, cip, crp, csp, rew, rei, res, rer, tem2d, clw2,       &
     &       qcwat, qcice, qrain, fcice, frain, rrime, rsden, clwf

      real (kind=kind_phys) :: ptop1(IX,NK_CLDS+1), rxlat(ix)

      real (kind=kind_phys) :: clwmin, clwm, clwt, onemrh, value,       &
     &       tem1, tem2, tem3

      integer :: i, k, id

!  ---  constant values
!     real (kind=kind_phys), parameter :: xrc3 = 200.
      real (kind=kind_phys), parameter :: xrc3 = 100.

!
!===> ... begin here
!
!     clouds(:,:,:) = 0.0

      do k = 1, NLAY
        do i = 1, IX
          cldtot(i,k) = 0.0
          cldcnv(i,k) = 0.0
          cwp   (i,k) = 0.0
          cip   (i,k) = 0.0
          crp   (i,k) = 0.0
          csp   (i,k) = 0.0
          rew   (i,k) = reliq_def            ! default liq radius to 10 micron
          rei   (i,k) = reice_def            ! default ice radius to 50 micron
          rer   (i,k) = rrain_def            ! default rain radius to 1000 micron
          res   (i,k) = rsnow_def            ! default snow radius to 250 micron
          fcice (i,k) = max(0.0, min(1.0, f_ice(i,k)))
          frain (i,k) = max(0.0, min(1.0, f_rain(i,k)))
          rrime (i,k) = max(1.0, r_rime(i,k))
          tem2d (i,k) = tlyr(i,k) - con_t0c
        enddo
      enddo
!
      if ( lcrick ) then
        do i = 1, IX
          clwf(i,1)    = 0.75*clw(i,1)    + 0.25*clw(i,2)
          clwf(i,nlay) = 0.75*clw(i,nlay) + 0.25*clw(i,nlay-1)
        enddo
        do k = 2, NLAY-1
          do i = 1, IX
            clwf(i,K) = 0.25*clw(i,k-1) + 0.5*clw(i,k) + 0.25*clw(i,k+1)
          enddo
        enddo
      else
        do k = 1, NLAY
          do i = 1, IX
            clwf(i,k) = clw(i,k)
          enddo
        enddo
      endif

!> -# Find top pressure (ptopc) for each cloud domain for given latitude.
!    - ptopc(k,i): top pressure of each cld domain (k=1-4 are sfc,l,m,
!!     h; i=1,2 are low-lat (<45 degree) and pole regions)

      do i =1, IX
        rxlat(i) = abs( xlat(i) / con_pi )      ! if xlat in pi/2 -> -pi/2 range
!       rxlat(i) = abs(0.5 - xlat(i)/con_pi)    ! if xlat in 0 -> pi range
      enddo

      do id = 1, 4
        tem1 = ptopc(id,2) - ptopc(id,1)

        do i =1, IX
          ptop1(i,id) = ptopc(id,1) + tem1*max( 0.0, 4.0*rxlat(i)-1.0 )
        enddo
      enddo

!> -# Seperate cloud condensate into liquid, ice, and rain types, and
!! save the liquid+ice condensate in array clw2 for later calculation
!!  of cloud fraction.

      do k = 1, NLAY
        do i = 1, IX
          if (tem2d(i,k) > -40.0) then
            qcice(i,k) = clwf(i,k) * fcice(i,k)
            tem1       = clwf(i,k) - qcice(i,k)
            qrain(i,k) = tem1 * frain(i,k)
            qcwat(i,k) = tem1 - qrain(i,k)
            clw2 (i,k) = qcwat(i,k) + qcice(i,k)
          else
            qcice(i,k) = clwf(i,k)
            qrain(i,k) = 0.0
            qcwat(i,k) = 0.0
            clw2 (i,k) = clwf(i,k)
          endif
        enddo
      enddo

!> -# Call module_microphysics::rsipath2(), in Ferrier's scheme, to
!! compute layer's cloud liquid, ice, rain, and snow water condensate
!! path and the partical effective radius for liquid droplet, rain drop, 
!! and snow flake.
      call  rsipath2                                                    &
!  ---  inputs:
     &     ( plyr, plvl, tlyr, qlyr, qcwat, qcice, qrain, rrime,        &
     &       IX, NLAY, ivflip, flgmin,                                  &
!  ---  outputs:
     &       cwp, cip, crp, csp, rew, rer, res, rsden                   &
     &     )


        do k = 1, NLAY
          do i = 1, IX
            tem2d(i,k) = (con_g * plyr(i,k))                            &
     &                 / (con_rd* delp(i,k))
          enddo
        enddo

!> -# Calculate layer cloud fraction.

        clwmin = 0.0e-6
        if (.not. lmfshal) then
          do k = 1, NLAY
          do i = 1, IX
!           clwt = 1.0e-7 * (plyr(i,k)*0.001)
!           clwt = 1.0e-6 * (plyr(i,k)*0.001)
            clwt = 2.0e-6 * (plyr(i,k)*0.001)
!           clwt = 5.0e-6 * (plyr(i,k)*0.001)
!           clwt = 5.0e-6

            if (clw2(i,k) > clwt) then
              onemrh= max( 1.e-10, 1.0-rhly(i,k) )
              clwm  = clwmin / max( 0.01, plyr(i,k)*0.001 )

!             tem1  = min(max(sqrt(onemrh*qstl(i,k)),0.0001),1.0)
!             tem1  = 100.0 / tem1

              tem1  = min(max(sqrt(sqrt(onemrh*qstl(i,k))),0.0001),1.0)
              tem1  = 2000.0 / tem1
!             tem1  = 2400.0 / tem1
!cnt          tem1  = 2500.0 / tem1
!             tem1  = min(max(sqrt(onemrh*qstl(i,k)),0.0001),1.0)
!             tem1  = 2000.0 / tem1
!             tem1  = 1000.0 / tem1
!             tem1  = 100.0 / tem1

              value = max( min( tem1*(clw2(i,k)-clwm), 50.0 ), 0.0 )
              tem2  = sqrt( sqrt(rhly(i,k)) )

              cldtot(i,k) = max( tem2*(1.0-exp(-value)), 0.0 )
            endif
          enddo
          enddo
        else
          do k = 1, NLAY
          do i = 1, IX
!           clwt = 1.0e-6 * (plyr(i,k)*0.001)
            clwt = 2.0e-6 * (plyr(i,k)*0.001)

            if (clw2(i,k) > clwt) then
              onemrh= max( 1.e-10, 1.0-rhly(i,k) )
              clwm  = clwmin / max( 0.01, plyr(i,k)*0.001 )
!
              tem1  = min(max((onemrh*qstl(i,k))**0.49,0.0001),1.0)   !jhan
              if (lmfdeep2) then
                tem1  = xrc3 / tem1
              else
                tem1  = 100.0 / tem1
              endif
!
              value = max( min( tem1*(clw2(i,k)-clwm), 50.0 ), 0.0 )
              tem2  = sqrt( sqrt(rhly(i,k)) )

              cldtot(i,k) = max( tem2*(1.0-exp(-value)), 0.0 )
            endif
          enddo
          enddo
        endif

      do k = 1, NLAY
        do i = 1, IX
          if (cldtot(i,k) < climit) then
            cldtot(i,k) = 0.0
            cwp(i,k)    = 0.0
            cip(i,k)    = 0.0
            crp(i,k)    = 0.0
            csp(i,k)    = 0.0
          endif
        enddo
      enddo

!     When lnoprec = .true. snow/rain has no impact on radiation
      if ( lnoprec ) then
        do k = 1, NLAY
          do i = 1, IX
            crp(i,k) = 0.0
            csp(i,k) = 0.0
          enddo
        enddo
      endif
!
      if ( lcnorm ) then
        do k = 1, NLAY
          do i = 1, IX
            if (cldtot(i,k) >= climit) then
              tem1 = 1.0 / max(climit2, cldtot(i,k))
              cwp(i,k) = cwp(i,k) * tem1
              cip(i,k) = cip(i,k) * tem1
              crp(i,k) = crp(i,k) * tem1
              csp(i,k) = csp(i,k) * tem1
            endif
          enddo
        enddo
      endif

!> -# Calculate effective ice cloud droplet radius.

      do k = 1, NLAY
        do i = 1, IX
          tem1 = tlyr(i,k) - con_ttp
          tem2 = cip(i,k)

          if (tem2 > 0.0) then
            tem3 = tem2d(i,k) * tem2 / tvly(i,k)

            if (tem1 < -50.0) then
              rei(i,k) = (1250.0/9.917) * tem3 ** 0.109
            elseif (tem1 < -40.0) then
              rei(i,k) = (1250.0/9.337) * tem3 ** 0.08
            elseif (tem1 < -30.0) then
              rei(i,k) = (1250.0/9.208) * tem3 ** 0.055
            else
              rei(i,k) = (1250.0/9.387) * tem3 ** 0.031
            endif

!           if (lprnt .and. k == l) print *,' reiL=',rei(i,k),' icec=', &
!    &        icec,' cip=',cip(i,k),' tem=',tem,' delt=',delt

            rei(i,k)   = max(10.0, min(rei(i,k), 300.0))
!           rei(i,k)   = max(20.0, min(rei(i,k), 300.0))
!!!!        rei(i,k)   = max(30.0, min(rei(i,k), 300.0))
!           rei(i,k)   = max(50.0, min(rei(i,k), 300.0))
!           rei(i,k)   = max(100.0, min(rei(i,k), 300.0))
          endif
        enddo
      enddo
!
      do k = 1, NLAY
        do i = 1, IX
          clouds(i,k,1) = cldtot(i,k)
          clouds(i,k,2) = cwp(i,k)
          clouds(i,k,3) = rew(i,k)
          clouds(i,k,4) = cip(i,k)
          clouds(i,k,5) = rei(i,k)
          clouds(i,k,6) = crp(i,k)
          clouds(i,k,7) = rer(i,k)
!         clouds(i,k,8) = csp(i,k)               !ncar scheme
          clouds(i,k,8) = csp(i,k) * rsden(i,k)  !fu's scheme
          clouds(i,k,9) = res(i,k)
        enddo
      enddo

!  --- ...  estimate clouds decorrelation length in km
!           this is only a tentative test, need to consider change later

      if ( iovr == 3 ) then
        do i = 1, ix
          de_lgth(i) = max( 0.6, 2.78-4.6*rxlat(i) )
        enddo
      endif

!> -# Call gethml(), to compute low, mid, high, total, and boundary
!! layer cloud fractions and clouds top/bottom layer indices for low,
!! mid, and high clouds.
!       The three cloud domain boundaries are defined by ptopc.  The cloud
!       overlapping method is defined by control flag 'iovr', which may
!       be different for lw and sw radiation programs.

      call gethml                                                       &
!  ---  inputs:
     &     ( plyr, ptop1, cldtot, cldcnv, dz, de_lgth,                  &
     &       IX,NLAY,                                                   &
!  ---  outputs:
     &       clds, mtop, mbot                                           &
     &     )


!
      return
!...................................
      end subroutine progcld2
!-----------------------------------
!> @}

!> This subroutine computes cloud related quantities using 
!! zhao/moorthi's prognostic cloud microphysics scheme + pdfcld.
!!\param plyr       (ix,nlay), model layer mean pressure in mb (100pa)
!!\param plvl       (ix,nlp1), model level pressure in mb (100pa)
!!\param tlyr       (ix,nlay), model layer mean temperature in K
!!\param tvly       (ix,nlay), model layer virtual temperature in K
!!\param qlyr       (ix,nlay), layer specific humidity in gm/gm
!!\param qstl       (ix,nlay), layer saturate humidity in gm/gm
!!\param rhly       (ix,nlay), layer relative humidity (=qlyr/qstl)
!!\param clw        (ix,nlay), layer cloud condensate amount
!!\param xlat       (ix), grid latitude in radians, default to pi/2 ->
!!                  -pi/2 range, otherwise see in-line comment
!!\param xlon       (ix), grid longitude in radians  (not used)
!!\param slmsk      (ix), sea/land mask array (sea:0,land:1,sea-ice:2)
!!\param dz         (IX,NLAY), layer thickness (km)
!!\param delp       (IX,NLAY), model layer pressure thickness in mb (100Pa)
!!\param ix         horizontal dimention
!!\param nlay,nlp1  vertical layer/level dimensions
!!\param deltaq     (ix,nlay), half total water distribution width
!!\param sup        supersaturation
!!\param kdt           
!!\param me         print control flag
!!\param clouds    (ix,nlay,nf_clds), cloud profiles
!!\n               (:,:,1) - layer total cloud fraction
!!\n               (:,:,2) - layer cloud liq water path (g/m**2)
!!\n               (:,:,3) - mean eff radius for liq cloud (micron)
!!\n               (:,:,4) - layer cloud ice water path (g/m**2)
!!\n               (:,:,5) - mean eff radius for ice cloud (micron)
!!\n               (:,:,6) - layer rain drop water path         not assigned
!!\n               (:,:,7) - mean eff radius for rain drop (micron)
!!\n               (:,:,8) - layer snow flake water path        not assigned
!!\n               (:,:,9) - mean eff radius for snow flake(micron)
!!\param clds      (ix,5), fraction of clouds for low, mid, hi, tot, bl
!!\param mtop      (ix,3), vertical indices for low, mid, hi cloud tops
!!\param mbot      (ix,3), vertical indices for low, mid, hi cloud bases
!!\param de_lgth   (ix),   clouds decorrelation length (km)
!!\section gen_progcld3 General Algorithm
!> @{
!-----------------------------------
      subroutine progcld3                                               &
     &     ( plyr,plvl,tlyr,tvly,qlyr,qstl,rhly,clw,cnvw,cnvc,          &    !  ---  inputs:
     &       xlat,xlon,slmsk, dz, delp,                                 &
     &       ix, nlay, nlp1,                                            &
     &       deltaq,sup,kdt,me,                                         &
     &       clouds,clds,mtop,mbot,de_lgth                              &    !  ---  outputs:
     &      )

! =================   subprogram documentation block   ================ !
!                                                                       !
! subprogram:    progcld3    computes cloud related quantities using    !
!   zhao/moorthi's prognostic cloud microphysics scheme.                !
!                                                                       !
! abstract:  this program computes cloud fractions from cloud           !
!   condensates, calculates liquid/ice cloud droplet effective radius,  !
!   and computes the low, mid, high, total and boundary layer cloud     !
!   fractions and the vertical indices of low, mid, and high cloud      !
!   top and base.  the three vertical cloud domains are set up in the   !
!   initial subroutine "cld_init".                                      !
!                                                                       !
! usage:         call progcld3                                          !
!                                                                       !
! subprograms called:   gethml                                          !
!                                                                       !
! attributes:                                                           !
!   language:   fortran 90                                              !
!   machine:    ibm-sp, sgi                                             !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
! input variables:                                                      !
!   plyr  (ix,nlay) : model layer mean pressure in mb (100pa)           !
!   plvl  (ix,nlp1) : model level pressure in mb (100pa)                !
!   tlyr  (ix,nlay) : model layer mean temperature in k                 !
!   tvly  (ix,nlay) : model layer virtual temperature in k              !
!   qlyr  (ix,nlay) : layer specific humidity in gm/gm                  !
!   qstl  (ix,nlay) : layer saturate humidity in gm/gm                  !
!   rhly  (ix,nlay) : layer relative humidity (=qlyr/qstl)              !
!   clw   (ix,nlay) : layer cloud condensate amount                     !
!   xlat  (ix)      : grid latitude in radians, default to pi/2 -> -pi/2!
!                     range, otherwise see in-line comment              !
!   xlon  (ix)      : grid longitude in radians  (not used)             !
!   slmsk (ix)      : sea/land mask array (sea:0,land:1,sea-ice:2)      !
!   dz    (ix,nlay) : layer thickness (km)                              !
!   delp  (ix,nlay) : model layer pressure thickness in mb (100Pa)      !
!   ix              : horizontal dimention                              !
!   nlay,nlp1       : vertical layer/level dimensions                   !
!   cnvw  (ix,nlay) : layer convective cloud condensate                 !
!   cnvc  (ix,nlay) : layer convective cloud cover                      !
!   deltaq(ix,nlay) : half total water distribution width               !
!   sup             : supersaturation                                   !

!                                                                       !
! output variables:                                                     !
!   clouds(ix,nlay,nf_clds) : cloud profiles                            !
!      clouds(:,:,1) - layer total cloud fraction                       !
!      clouds(:,:,2) - layer cloud liq water path         (g/m**2)      !
!      clouds(:,:,3) - mean eff radius for liq cloud      (micron)      !
!      clouds(:,:,4) - layer cloud ice water path         (g/m**2)      !
!      clouds(:,:,5) - mean eff radius for ice cloud      (micron)      !
!      clouds(:,:,6) - layer rain drop water path         not assigned  !
!      clouds(:,:,7) - mean eff radius for rain drop      (micron)      !
!  *** clouds(:,:,8) - layer snow flake water path        not assigned  !
!      clouds(:,:,9) - mean eff radius for snow flake     (micron)      !
!  *** fu's scheme need to be normalized by snow density (g/m**3/1.0e6) !
!   clds  (ix,5)    : fraction of clouds for low, mid, hi, tot, bl      !
!   mtop  (ix,3)    : vertical indices for low, mid, hi cloud tops      !
!   mbot  (ix,3)    : vertical indices for low, mid, hi cloud bases     !
!   de_lgth(ix)     : clouds decorrelation length (km)                  !
!                                                                       !
! module variables:                                                     !
!   ivflip          : control flag of vertical index direction          !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!   lcrick          : control flag for eliminating crick                !
!                     =t: apply layer smoothing to eliminate crick      !
!                     =f: do not apply layer smoothing                  !
!   lcnorm          : control flag for in-cld condensate                !
!                     =t: normalize cloud condensate                    !
!                     =f: not normalize cloud condensate                !
!                                                                       !
!  ====================    end of description    =====================  !
!
      implicit none

!  ---  inputs
      integer,  intent(in) :: ix, nlay, nlp1,kdt

      real (kind=kind_phys), dimension(:,:), intent(in) :: plvl, plyr,    &
     &       tlyr, tvly, qlyr, qstl, rhly, clw, dz, delp
!     &       tlyr, tvly, qlyr, qstl, rhly, clw, cnvw, cnvc
!      real (kind=kind_phys), dimension(:,:), intent(in) :: deltaq
      real (kind=kind_phys), dimension(:,:) :: deltaq, cnvw, cnvc
      real (kind=kind_phys) qtmp,qsc,rhs
      real (kind=kind_phys), intent(in) :: sup
      real (kind=kind_phys), parameter :: epsq = 1.0e-12

      real (kind=kind_phys), dimension(:),   intent(in) :: xlat, xlon,    &
     &       slmsk
      integer :: me

!  ---  outputs
      real (kind=kind_phys), dimension(:,:,:), intent(out) :: clouds

      real (kind=kind_phys), dimension(:,:),   intent(out) :: clds
      real (kind=kind_phys), dimension(:),     intent(out) :: de_lgth

      integer,               dimension(:,:),   intent(out) :: mtop,mbot

!  ---  local variables:
      real (kind=kind_phys), dimension(ix,nlay) :: cldtot, cldcnv,      &
     &       cwp, cip, crp, csp, rew, rei, res, rer, tem2d, clwf

      real (kind=kind_phys) :: ptop1(ix,nk_clds+1), rxlat(ix)

      real (kind=kind_phys) :: clwmin, clwm, clwt, onemrh, value,       &
     &       tem1, tem2, tem3

      integer :: i, k, id, nf

!
!===> ... begin here
!
      do nf=1,nf_clds
        do k=1,nlay
          do i=1,ix
            clouds(i,k,nf) = 0.0
          enddo
        enddo
      enddo
!     clouds(:,:,:) = 0.0

      do k = 1, nlay
        do i = 1, ix
          cldtot(i,k) = 0.0
          cldcnv(i,k) = 0.0
          cwp   (i,k) = 0.0
          cip   (i,k) = 0.0
          crp   (i,k) = 0.0
          csp   (i,k) = 0.0
          rew   (i,k) = reliq_def            ! default liq radius to 10 micron
          rei   (i,k) = reice_def            ! default ice radius to 50 micron
          rer   (i,k) = rrain_def            ! default rain radius to 1000 micron
          res   (i,k) = rsnow_def            ! default snow radius to 250 micron
          tem2d (i,k) = min( 1.0, max( 0.0, (con_ttp-tlyr(i,k))*0.05 ) )
          clwf(i,k)   = 0.0
        enddo
      enddo
!
      if ( lcrick ) then
        do i = 1, ix
          clwf(i,1)    = 0.75*clw(i,1)    + 0.25*clw(i,2)
          clwf(i,nlay) = 0.75*clw(i,nlay) + 0.25*clw(i,nlay-1)
        enddo
        do k = 2, nlay-1
          do i = 1, ix
            clwf(i,k) = 0.25*clw(i,k-1) + 0.5*clw(i,k) + 0.25*clw(i,k+1)
          enddo
        enddo
      else
        do k = 1, nlay
          do i = 1, ix
            clwf(i,k) = clw(i,k)
          enddo
        enddo
      endif

      if(kdt==1) then
        do k = 1, nlay
          do i = 1, ix
            deltaq(i,k) = (1.-0.95)*qstl(i,k)
          enddo
        enddo
      endif

!> -# Find top pressure (ptopc) for each cloud domain for given latitude.
!     ptopc(k,i): top presure of each cld domain (k=1-4 are sfc,l,m,h;
!  ---  i=1,2 are low-lat (<45 degree) and pole regions)

      do i =1, IX
        rxlat(i) = abs( xlat(i) / con_pi )      ! if xlat in pi/2 -> -pi/2 range
!       rxlat(i) = abs(0.5 - xlat(i)/con_pi)    ! if xlat in 0 -> pi range
      enddo

      do id = 1, 4
        tem1 = ptopc(id,2) - ptopc(id,1)

        do i =1, ix
          ptop1(i,id) = ptopc(id,1) + tem1*max( 0.0, 4.0*rxlat(i)-1.0 )
        enddo
      enddo

!> -# Calculate liquid/ice condensate path in \f$ g/m^2 \f$

        do k = 1, nlay
          do i = 1, ix
            clwt = max(0.0,(clwf(i,k)+cnvw(i,k))) * gfac * delp(i,k)
            cip(i,k) = clwt * tem2d(i,k)
            cwp(i,k) = clwt - cip(i,k)
          enddo
        enddo

!> -# Calculate effective liquid cloud droplet radius over land.

      do i = 1, ix
        if (nint(slmsk(i)) == 1) then
          do k = 1, nlay
            rew(i,k) = 5.0 + 5.0 * tem2d(i,k)
          enddo
        endif
      enddo

!> -# Calculate layer cloud fraction.

          do k = 1, nlay
          do i = 1, ix
            tem1 = tlyr(i,k) - 273.16
            if(tem1 < con_thgni) then  ! for pure ice, has to be consistent with gscond
              qsc = sup * qstl(i,k)
              rhs = sup
            else
              qsc = qstl(i,k)
              rhs = 1.0
            endif
           if(rhly(i,k) >= rhs) then
             cldtot(i,k) = 1.0
           else
             qtmp = qlyr(i,k) + clwf(i,k) - qsc
             if(deltaq(i,k) > epsq) then
!              if(qtmp <= -deltaq(i,k) .or. cwmik < epsq) then
               if(qtmp <= -deltaq(i,k)) then
                 cldtot(i,k) = 0.0
               elseif(qtmp >= deltaq(i,k)) then
                 cldtot(i,k) = 1.0
               else
                 cldtot(i,k) = 0.5*qtmp/deltaq(i,k) + 0.5
                 cldtot(i,k) = max(cldtot(i,k),0.0)
                 cldtot(i,k) = min(cldtot(i,k),1.0)
               endif
             else
               if(qtmp > 0.) then
                 cldtot(i,k) = 1.0
               else
                 cldtot(i,k) = 0.0
               endif
             endif
           endif
           cldtot(i,k) = cnvc(i,k) + (1-cnvc(i,k))*cldtot(i,k)
           cldtot(i,k) = max(cldtot(i,k),0.)
           cldtot(i,k) = min(cldtot(i,k),1.)

          enddo
          enddo

      do k = 1, nlay
        do i = 1, ix
          if (cldtot(i,k) < climit) then
            cldtot(i,k) = 0.0
            cwp(i,k)    = 0.0
            cip(i,k)    = 0.0
            crp(i,k)    = 0.0
            csp(i,k)    = 0.0
          endif
        enddo
      enddo

      if ( lcnorm ) then
        do k = 1, nlay
          do i = 1, ix
            if (cldtot(i,k) >= climit) then
              tem1 = 1.0 / max(climit2, cldtot(i,k))
              cwp(i,k) = cwp(i,k) * tem1
              cip(i,k) = cip(i,k) * tem1
              crp(i,k) = crp(i,k) * tem1
              csp(i,k) = csp(i,k) * tem1
            endif
          enddo
        enddo
      endif

!> -# Calculate effective ice cloud droplet radius following Heymsfield
!!    and McFarquhar (1996) \cite heymsfield_and_mcfarquhar_1996.

      do k = 1, nlay
        do i = 1, ix
          tem2 = tlyr(i,k) - con_ttp

          if (cip(i,k) > 0.0) then
!           tem3 = gord * cip(i,k) * (plyr(i,k)/delp(i,k)) / tvly(i,k)
            tem3 = gord * cip(i,k) * plyr(i,k) / (delp(i,k)*tvly(i,k))

            if (tem2 < -50.0) then
              rei(i,k) = (1250.0/9.917) * tem3 ** 0.109
            elseif (tem2 < -40.0) then
              rei(i,k) = (1250.0/9.337) * tem3 ** 0.08
            elseif (tem2 < -30.0) then
              rei(i,k) = (1250.0/9.208) * tem3 ** 0.055
            else
              rei(i,k) = (1250.0/9.387) * tem3 ** 0.031
            endif
!           rei(i,k)   = max(20.0, min(rei(i,k), 300.0))
!           rei(i,k)   = max(10.0, min(rei(i,k), 100.0))
            rei(i,k)   = max(10.0, min(rei(i,k), 150.0))
!           rei(i,k)   = max(5.0, min(rei(i,k), 130.0))
          endif
        enddo
      enddo

!
      do k = 1, nlay
        do i = 1, ix
          clouds(i,k,1) = cldtot(i,k)
          clouds(i,k,2) = cwp(i,k)
          clouds(i,k,3) = rew(i,k)
          clouds(i,k,4) = cip(i,k)
          clouds(i,k,5) = rei(i,k)
!         clouds(i,k,6) = 0.0
          clouds(i,k,7) = rer(i,k)
!         clouds(i,k,8) = 0.0
          clouds(i,k,9) = res(i,k)
        enddo
      enddo

!  --- ...  estimate clouds decorrelation length in km
!           this is only a tentative test, need to consider change later

      if ( iovr == 3 ) then
        do i = 1, ix
          de_lgth(i) = max( 0.6, 2.78-4.6*rxlat(i) )
        enddo
      endif

!> -# Call gethml() to compute low,mid,high,total, and boundary layer
!! cloud fractions and clouds top/bottom layer indices for low, mid,
!! and high clouds.
!       the three cloud domain boundaries are defined by ptopc.  the cloud
!       overlapping method is defined by control flag 'iovr', which may
!       be different for lw and sw radiation programs.


      call gethml                                                       &
!  ---  inputs:
     &     ( plyr, ptop1, cldtot, cldcnv, dz, de_lgth,                  &
     &       ix,nlay,                                                   &
!  ---  outputs:
     &       clds, mtop, mbot                                           &
     &     )


!
      return
!...................................
      end subroutine progcld3
!-----------------------------------


!-----------------------------------
      subroutine progcld4                                               &
!...................................

!  ---  inputs:
     &     ( plyr,plvl,tlyr,tvly,qlyr,qstl,rhly,clw,cnvw,cnvc,          &
     &       xlat,xlon,slmsk,cldtot, dz, delp,                          &
     &       IX, NLAY, NLP1,                                            &
!  ---  outputs:
     &       clouds,clds,mtop,mbot,de_lgth                              &
     &      )

! =================   subprogram documentation block   ================ !
!                                                                       !
! subprogram:    progcld4    computes cloud related quantities using    !
!   GFDL Lin MP prognostic cloud microphysics scheme.                   !
!                                                                       !
! abstract:  this program computes cloud fractions from cloud           !
!   condensates, calculates liquid/ice cloud droplet effective radius,  !
!   and computes the low, mid, high, total and boundary layer cloud     !
!   fractions and the vertical indices of low, mid, and high cloud      !
!   top and base.  the three vertical cloud domains are set up in the   !
!   initial subroutine "cld_init".                                      !
!                                                                       !
! usage:         call progcld4                                          !
!                                                                       !
! subprograms called:   gethml                                          !
!                                                                       !
! attributes:                                                           !
!   language:   fortran 90                                              !
!   machine:    ibm-sp, sgi                                             !
!                                                                       !
!                                                                       !
!  ====================  definition of variables  ====================  !
!                                                                       !
! input variables:                                                      !
!   plyr  (IX,NLAY) : model layer mean pressure in mb (100Pa)           !
!   plvl  (IX,NLP1) : model level pressure in mb (100Pa)                !
!   tlyr  (IX,NLAY) : model layer mean temperature in k                 !
!   tvly  (IX,NLAY) : model layer virtual temperature in k              !
!   qlyr  (IX,NLAY) : layer specific humidity in gm/gm                  !
!   qstl  (IX,NLAY) : layer saturate humidity in gm/gm                  !
!   rhly  (IX,NLAY) : layer relative humidity (=qlyr/qstl)              !
!   clw   (IX,NLAY) : layer cloud condensate amount                     !
!   cnvw  (IX,NLAY) : layer convective cloud condensate                 !
!   cnvc  (IX,NLAY) : layer convective cloud cover                      !
!   xlat  (IX)      : grid latitude in radians, default to pi/2 -> -pi/2!
!                     range, otherwise see in-line comment              !
!   xlon  (IX)      : grid longitude in radians  (not used)             !
!   slmsk (IX)      : sea/land mask array (sea:0,land:1,sea-ice:2)      !
!   dz    (ix,nlay) : layer thickness (km)                              !
!   delp  (ix,nlay) : model layer pressure thickness in mb (100Pa)      !
!   IX              : horizontal dimention                              !
!   NLAY,NLP1       : vertical layer/level dimensions                   !
!                                                                       !
! output variables:                                                     !
!   clouds(IX,NLAY,NF_CLDS) : cloud profiles                            !
!      clouds(:,:,1) - layer total cloud fraction                       !
!      clouds(:,:,2) - layer cloud liq water path         (g/m**2)      !
!      clouds(:,:,3) - mean eff radius for liq cloud      (micron)      !
!      clouds(:,:,4) - layer cloud ice water path         (g/m**2)      !
!      clouds(:,:,5) - mean eff radius for ice cloud      (micron)      !
!      clouds(:,:,6) - layer rain drop water path         not assigned  !
!      clouds(:,:,7) - mean eff radius for rain drop      (micron)      !
!  *** clouds(:,:,8) - layer snow flake water path        not assigned  !
!      clouds(:,:,9) - mean eff radius for snow flake     (micron)      !
!  *** fu's scheme need to be normalized by snow density (g/m**3/1.0e6) !
!   clds  (IX,5)    : fraction of clouds for low, mid, hi, tot, bl      !
!   mtop  (IX,3)    : vertical indices for low, mid, hi cloud tops      !
!   mbot  (IX,3)    : vertical indices for low, mid, hi cloud bases     !
!   de_lgth(ix)     : clouds decorrelation length (km)                  !
!                                                                       !
! module variables:                                                     !
!   ivflip          : control flag of vertical index direction          !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!   lsashal         : control flag for shallow convection               !
!   lcrick          : control flag for eliminating CRICK                !
!                     =t: apply layer smoothing to eliminate CRICK      !
!                     =f: do not apply layer smoothing                  !
!   lcnorm          : control flag for in-cld condensate                !
!                     =t: normalize cloud condensate                    !
!                     =f: not normalize cloud condensate                !
!                                                                       !
!  ====================    end of description    =====================  !
!
      implicit none

!  ---  inputs
      integer,  intent(in) :: IX, NLAY, NLP1

      real (kind=kind_phys), dimension(:,:), intent(in) :: plvl, plyr,  &
     &       tlyr, tvly, qlyr, qstl, rhly, clw, cldtot, cnvw, cnvc,     &
     &       delp, dz

      real (kind=kind_phys), dimension(:),   intent(in) :: xlat, xlon,  &
     &       slmsk

!  ---  outputs
      real (kind=kind_phys), dimension(:,:,:), intent(out) :: clouds

      real (kind=kind_phys), dimension(:,:),   intent(out) :: clds
      real (kind=kind_phys), dimension(:),     intent(out) :: de_lgth

      integer,               dimension(:,:),   intent(out) :: mtop,mbot

!  ---  local variables:
      real (kind=kind_phys), dimension(IX,NLAY) :: cldcnv,              &
     &       cwp, cip, crp, csp, rew, rei, res, rer, tem2d, clwf

      real (kind=kind_phys) :: ptop1(IX,NK_CLDS+1), rxlat(ix)

      real (kind=kind_phys) :: clwmin, clwm, clwt, onemrh, value,       &
     &       tem1, tem2, tem3

      integer :: i, k, id, nf

!
!===> ... begin here
!
      do nf=1,nf_clds
        do k=1,nlay
          do i=1,ix
            clouds(i,k,nf) = 0.0
          enddo
        enddo
      enddo
!     clouds(:,:,:) = 0.0

      do k = 1, NLAY
        do i = 1, IX
          cldcnv(i,k) = 0.0
          cwp   (i,k) = 0.0
          cip   (i,k) = 0.0
          crp   (i,k) = 0.0
          csp   (i,k) = 0.0
          rew   (i,k) = reliq_def            ! default liq radius to 10 micron
          rei   (i,k) = reice_def            ! default ice radius to 50 micron
          rer   (i,k) = rrain_def            ! default rain radius to 1000 micron
          res   (i,k) = rsnow_def            ! default snow radius to 250 micron
          tem2d (i,k) = min( 1.0, max( 0.0, (con_ttp-tlyr(i,k))*0.05 ) )
          clwf(i,k)   = 0.0
        enddo
      enddo
!
      if ( lcrick ) then
        do i = 1, IX
          clwf(i,1)    = 0.75*clw(i,1)    + 0.25*clw(i,2)
          clwf(i,nlay) = 0.75*clw(i,nlay) + 0.25*clw(i,nlay-1)
        enddo
        do k = 2, NLAY-1
          do i = 1, IX
            clwf(i,K) = 0.25*clw(i,k-1) + 0.5*clw(i,k) + 0.25*clw(i,k+1)
          enddo
        enddo
      else
        do k = 1, NLAY
          do i = 1, IX
            clwf(i,k) = clw(i,k)
          enddo
        enddo
      endif

!  ---  find top pressure for each cloud domain for given latitude
!       ptopc(k,i): top presure of each cld domain (k=1-4 are sfc,L,m,h;
!  ---  i=1,2 are low-lat (<45 degree) and pole regions)

      do i =1, IX
        rxlat(i) = abs( xlat(i) / con_pi )      ! if xlat in pi/2 -> -pi/2 range
!       rxlat(i) = abs(0.5 - xlat(i)/con_pi)    ! if xlat in 0 -> pi range
      enddo

      do id = 1, 4
        tem1 = ptopc(id,2) - ptopc(id,1)

        do i =1, IX
          ptop1(i,id) = ptopc(id,1) + tem1*max( 0.0, 4.0*rxlat(i)-1.0 )
        enddo
      enddo

!  ---  compute liquid/ice condensate path in g/m**2

        do k = 1, NLAY
          do i = 1, IX
            clwt     = max(0.0,(clwf(i,k)+cnvw(i,k))) * gfac * delp(i,k)
            cip(i,k) = clwt * tem2d(i,k)
            cwp(i,k) = clwt - cip(i,k)
          enddo
        enddo

!  ---  effective liquid cloud droplet radius over land

      do i = 1, IX
        if (nint(slmsk(i)) == 1) then
          do k = 1, NLAY
            rew(i,k) = 5.0 + 5.0 * tem2d(i,k)
          enddo
        endif
      enddo

      do k = 1, NLAY
        do i = 1, IX
          if (cldtot(i,k) < climit) then
            cwp(i,k)    = 0.0
            cip(i,k)    = 0.0
            crp(i,k)    = 0.0
            csp(i,k)    = 0.0
          endif
        enddo
      enddo

      if ( lcnorm ) then
        do k = 1, NLAY
          do i = 1, IX
            if (cldtot(i,k) >= climit) then
              tem1 = 1.0 / max(climit2, cldtot(i,k))
              cwp(i,k) = cwp(i,k) * tem1
              cip(i,k) = cip(i,k) * tem1
              crp(i,k) = crp(i,k) * tem1
              csp(i,k) = csp(i,k) * tem1
            endif
          enddo
        enddo
      endif

!  ---  effective ice cloud droplet radius

      do k = 1, NLAY
        do i = 1, IX
          tem2 = tlyr(i,k) - con_ttp

          if (cip(i,k) > 0.0) then
            tem3 = gord * cip(i,k) * plyr(i,k) / (delp(i,k)*tvly(i,k))

            if (tem2 < -50.0) then
              rei(i,k) = (1250.0/9.917) * tem3 ** 0.109
            elseif (tem2 < -40.0) then
              rei(i,k) = (1250.0/9.337) * tem3 ** 0.08
            elseif (tem2 < -30.0) then
              rei(i,k) = (1250.0/9.208) * tem3 ** 0.055
            else
              rei(i,k) = (1250.0/9.387) * tem3 ** 0.031
            endif
!           rei(i,k)   = max(20.0, min(rei(i,k), 300.0))
!           rei(i,k)   = max(10.0, min(rei(i,k), 100.0))
            rei(i,k)   = max(10.0, min(rei(i,k), 150.0))
!           rei(i,k)   = max(5.0,  min(rei(i,k), 130.0))
          endif
        enddo
      enddo

!
      do k = 1, NLAY
        do i = 1, IX
          clouds(i,k,1) = cldtot(i,k)
          clouds(i,k,2) = cwp(i,k)
          clouds(i,k,3) = rew(i,k)
          clouds(i,k,4) = cip(i,k)
          clouds(i,k,5) = rei(i,k)
!         clouds(i,k,6) = 0.0
          clouds(i,k,7) = rer(i,k)
!         clouds(i,k,8) = 0.0
          clouds(i,k,9) = res(i,k)
        enddo
      enddo

!  --- ...  estimate clouds decorrelation length in km
!           this is only a tentative test, need to consider change later

      if ( iovr == 3 ) then
        do i = 1, ix
          de_lgth(i) = max( 0.6, 2.78-4.6*rxlat(i) )
        enddo
      endif

!  ---  compute low, mid, high, total, and boundary layer cloud fractions
!       and clouds top/bottom layer indices for low, mid, and high clouds.
!       The three cloud domain boundaries are defined by ptopc.  The cloud
!       overlapping method is defined by control flag 'iovr', which may
!       be different for lw and sw radiation programs.

      call gethml                                                       &
!  ---  inputs:
     &     ( plyr, ptop1, cldtot, cldcnv, dz, de_lgth,                  &
     &       IX,NLAY,                                                   &
!  ---  outputs:
     &       clds, mtop, mbot                                           &
     &     )


!
      return
!...................................
      end subroutine progcld4
!-----------------------------------

!-----------------------------------
      subroutine progcld4o                                              &
!...................................

!  ---  inputs:
     &     ( plyr,plvl,tlyr,tvly,qlyr,qstl,rhly,clw,                    &
     &       xlat,xlon,slmsk, dz, delp,                                 &
     &       ntrac,ntcw,ntiw,ntrw,ntsw,ntgl,ntclamt,                    &
     &       IX, NLAY, NLP1,                                            &
!  ---  outputs:
     &       clouds,clds,mtop,mbot,de_lgth                              &
     &      )

! =================   subprogram documentation block   ================ !
!                                                                       !
! subprogram:    progcld4o   computes cloud related quantities using    !
!   GFDL Lin MP prognostic cloud microphysics scheme. Moist species     !
!   from MP are fed into the corresponding arrays for calcuation of     !
!                                                                       !
! abstract:  this program computes cloud fractions from cloud           !
!   condensates, calculates liquid/ice cloud droplet effective radius,  !
!   and computes the low, mid, high, total and boundary layer cloud     !
!   fractions and the vertical indices of low, mid, and high cloud      !
!   top and base.  the three vertical cloud domains are set up in the   !
!   initial subroutine "cld_init".                                      !
!                                                                       !
! usage:         call progcld4o                                         !
!                                                                       !
! subprograms called:   gethml                                          !
!                                                                       !
! attributes:                                                           !
!   language:   fortran 90                                              !
!   machine:    ibm-sp, sgi                                             !
!                                                                       !
!                                                                       !
!  ====================  definition of variables  ====================  !
!                                                                       !
! input variables:                                                      !
!   plyr  (IX,NLAY) : model layer mean pressure in mb (100Pa)           !
!   plvl  (IX,NLP1) : model level pressure in mb (100Pa)                !
!   tlyr  (IX,NLAY) : model layer mean temperature in k                 !
!   tvly  (IX,NLAY) : model layer virtual temperature in k              !
!   qlyr  (IX,NLAY) : layer specific humidity in gm/gm                  !
!   qstl  (IX,NLAY) : layer saturate humidity in gm/gm                  !
!   rhly  (IX,NLAY) : layer relative humidity (=qlyr/qstl)              !
!   clw   (IX,NLAY,NTRAC) : layer cloud condensate amount               !
!   xlat  (IX)      : grid latitude in radians, default to pi/2 -> -pi/2!
!                     range, otherwise see in-line comment              !
!   xlon  (IX)      : grid longitude in radians  (not used)             !
!   slmsk (IX)      : sea/land mask array (sea:0,land:1,sea-ice:2)      !
!   dz    (ix,nlay) : layer thickness (km)                              !
!   delp  (ix,nlay) : model layer pressure thickness in mb (100Pa)      !
!   IX              : horizontal dimention                              !
!   NLAY,NLP1       : vertical layer/level dimensions                   !
!                                                                       !
! output variables:                                                     !
!   clouds(IX,NLAY,NF_CLDS) : cloud profiles                            !
!      clouds(:,:,1) - layer total cloud fraction                       !
!      clouds(:,:,2) - layer cloud liq water path         (g/m**2)      !
!      clouds(:,:,3) - mean eff radius for liq cloud      (micron)      !
!      clouds(:,:,4) - layer cloud ice water path         (g/m**2)      !
!      clouds(:,:,5) - mean eff radius for ice cloud      (micron)      !
!      clouds(:,:,6) - layer rain drop water path         not assigned  !
!      clouds(:,:,7) - mean eff radius for rain drop      (micron)      !
!  *** clouds(:,:,8) - layer snow flake water path        not assigned  !
!      clouds(:,:,9) - mean eff radius for snow flake     (micron)      !
!  *** fu's scheme need to be normalized by snow density (g/m**3/1.0e6) !
!   clds  (IX,5)    : fraction of clouds for low, mid, hi, tot, bl      !
!   mtop  (IX,3)    : vertical indices for low, mid, hi cloud tops      !
!   mbot  (IX,3)    : vertical indices for low, mid, hi cloud bases     !
!   de_lgth(ix)     : clouds decorrelation length (km)                  !
!                                                                       !
! module variables:                                                     !
!   ivflip          : control flag of vertical index direction          !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!   lsashal         : control flag for shallow convection               !
!   lcrick          : control flag for eliminating CRICK                !
!                     =t: apply layer smoothing to eliminate CRICK      !
!                     =f: do not apply layer smoothing                  !
!   lcnorm          : control flag for in-cld condensate                !
!                     =t: normalize cloud condensate                    !
!                     =f: not normalize cloud condensate                !
!                                                                       !
!  ====================    end of description    =====================  !
!
      implicit none

!  ---  inputs
      integer,  intent(in) :: IX, NLAY, NLP1
      integer,  intent(in) :: ntrac, ntcw, ntiw, ntrw, ntsw, ntgl,      &
     &		 		ntclamt

      real (kind=kind_phys), dimension(:,:), intent(in) :: plvl, plyr,  &
     &       tlyr, tvly, qlyr, qstl, rhly, delp, dz


      real (kind=kind_phys), dimension(:,:,:), intent(in) :: clw
      real (kind=kind_phys), dimension(:),   intent(in) :: xlat, xlon,  &
     &       slmsk

!  ---  outputs
      real (kind=kind_phys), dimension(:,:,:), intent(out) :: clouds

      real (kind=kind_phys), dimension(:,:),   intent(out) :: clds
      real (kind=kind_phys), dimension(:),     intent(out) :: de_lgth

      integer,               dimension(:,:),   intent(out) :: mtop,mbot

!  ---  local variables:
      real (kind=kind_phys), dimension(IX,NLAY) :: cldcnv,              &
     &       cwp, cip, crp, csp, rew, rei, res, rer, tem2d

      real (kind=kind_phys) :: ptop1(IX,NK_CLDS+1), rxlat(ix)

      real (kind=kind_phys) :: clwmin, clwm, clwt, onemrh, value,       &
     &       tem1, tem2, tem3
      real (kind=kind_phys), dimension(IX,NLAY) :: cldtot

      integer :: i, k, id, nf

!
!===> ... begin here
!
      do nf=1,nf_clds
        do k=1,nlay
          do i=1,ix
            clouds(i,k,nf) = 0.0
          enddo
        enddo
      enddo
!     clouds(:,:,:) = 0.0

      do k = 1, NLAY
        do i = 1, IX
          cldcnv(i,k) = 0.0
          cwp   (i,k) = 0.0
          cip   (i,k) = 0.0
          crp   (i,k) = 0.0
          csp   (i,k) = 0.0
          rew   (i,k) = reliq_def            ! default liq radius to 10 micron
          rei   (i,k) = reice_def            ! default ice radius to 50 micron
          rer   (i,k) = rrain_def            ! default rain radius to 1000 micron
          res   (i,k) = rsnow_def            ! default snow radius to 250 micron
          tem2d (i,k) = min( 1.0, max( 0.0, (con_ttp-tlyr(i,k))*0.05 ) )
          cldtot(i,k) = clw(i,k,ntclamt)
        enddo
      enddo

!  ---  find top pressure for each cloud domain for given latitude
!       ptopc(k,i): top presure of each cld domain (k=1-4 are sfc,L,m,h;
!  ---  i=1,2 are low-lat (<45 degree) and pole regions)

      do i =1, IX
        rxlat(i) = abs( xlat(i) / con_pi )      ! if xlat in pi/2 -> -pi/2 range
!       rxlat(i) = abs(0.5 - xlat(i)/con_pi)    ! if xlat in 0 -> pi range
      enddo

      do id = 1, 4
        tem1 = ptopc(id,2) - ptopc(id,1)

        do i =1, IX
          ptop1(i,id) = ptopc(id,1) + tem1*max( 0.0, 4.0*rxlat(i)-1.0 )
        enddo
      enddo

!  ---  compute liquid/ice condensate path in g/m**2

        do k = 1, NLAY
          do i = 1, IX
            cwp(i,k) = max(0.0, clw(i,k,ntcw) * gfac * delp(i,k))
            cip(i,k) = max(0.0, clw(i,k,ntiw) * gfac * delp(i,k))
            crp(i,k) = max(0.0, clw(i,k,ntrw) * gfac * delp(i,k))
            csp(i,k) = max(0.0, (clw(i,k,ntsw)+clw(i,k,ntgl)) *         &
     &                  gfac * delp(i,k))
          enddo
        enddo

!  ---  effective liquid cloud droplet radius over land

      do i = 1, IX
        if (nint(slmsk(i)) == 1) then
          do k = 1, NLAY
            rew(i,k) = 5.0 + 5.0 * tem2d(i,k)
          enddo
        endif
      enddo

      do k = 1, NLAY
        do i = 1, IX
          if (cldtot(i,k) < climit) then
            cwp(i,k)    = 0.0
            cip(i,k)    = 0.0
            crp(i,k)    = 0.0
            csp(i,k)    = 0.0
          endif
        enddo
      enddo

      if ( lcnorm ) then
        do k = 1, NLAY
          do i = 1, IX
            if (cldtot(i,k) >= climit) then
              tem1 = 1.0 / max(climit2, cldtot(i,k))
              cwp(i,k) = cwp(i,k) * tem1
              cip(i,k) = cip(i,k) * tem1
              crp(i,k) = crp(i,k) * tem1
              csp(i,k) = csp(i,k) * tem1
            endif
          enddo
        enddo
      endif

!  ---  effective ice cloud droplet radius

      do k = 1, NLAY
        do i = 1, IX
          tem2 = tlyr(i,k) - con_ttp

          if (cip(i,k) > 0.0) then
            tem3 = gord * cip(i,k) * plyr(i,k) / (delp(i,k)*tvly(i,k))

            if (tem2 < -50.0) then
              rei(i,k) = (1250.0/9.917) * tem3 ** 0.109
            elseif (tem2 < -40.0) then
              rei(i,k) = (1250.0/9.337) * tem3 ** 0.08
            elseif (tem2 < -30.0) then
              rei(i,k) = (1250.0/9.208) * tem3 ** 0.055
            else
              rei(i,k) = (1250.0/9.387) * tem3 ** 0.031
            endif
!           rei(i,k)   = max(20.0, min(rei(i,k), 300.0))
!           rei(i,k)   = max(10.0, min(rei(i,k), 100.0))
            rei(i,k)   = max(10.0, min(rei(i,k), 150.0))
!           rei(i,k)   = max(5.0,  min(rei(i,k), 130.0))
          endif
        enddo
      enddo

!
      do k = 1, NLAY
        do i = 1, IX
          clouds(i,k,1) = cldtot(i,k)
          clouds(i,k,2) = cwp(i,k)
          clouds(i,k,3) = rew(i,k)
          clouds(i,k,4) = cip(i,k)
          clouds(i,k,5) = rei(i,k)
          clouds(i,k,6) = crp(i,k) 
          clouds(i,k,7) = rer(i,k)
          clouds(i,k,8) = csp(i,k)
          clouds(i,k,9) = rei(i,k)
        enddo
      enddo

!  --- ...  estimate clouds decorrelation length in km
!           this is only a tentative test, need to consider change later

      if ( iovr == 3 ) then
        do i = 1, ix
          de_lgth(i) = max( 0.6, 2.78-4.6*rxlat(i) )
        enddo
      endif

!  ---  compute low, mid, high, total, and boundary layer cloud fractions
!       and clouds top/bottom layer indices for low, mid, and high clouds.
!       The three cloud domain boundaries are defined by ptopc.  The cloud
!       overlapping method is defined by control flag 'iovr', which may
!       be different for lw and sw radiation programs.

      call gethml                                                       &
!  ---  inputs:
     &     ( plyr, ptop1, cldtot, cldcnv, dz, de_lgth,                  &
     &       IX,NLAY,                                                   &
!  ---  outputs:
     &       clds, mtop, mbot                                           &
     &     )


!
      return
!...................................
      end subroutine progcld4o
!-----------------------------------

!-----------------------------------

      subroutine progcld5                                               &
     &     ( plyr,plvl,tlyr,qlyr,qstl,rhly,clw,                         &    !  ---  inputs:
     &       xlat,xlon,slmsk,dz,delp,                                   & 
     &       ntrac,ntcw,ntiw,ntrw,ntsw,ntgl,                            &            
     &       IX, NLAY, NLP1,                                            &
     &       uni_cld, lmfshal, lmfdeep2, cldcov,                        &    
     &       re_cloud,re_ice,re_snow,                                   & 
     &       clouds,clds,mtop,mbot,de_lgth                              &    !  ---  outputs:
     &      )

! =================   subprogram documentation block   ================ !
!                                                                       !
! subprogram:    progcld5    computes cloud related quantities using    !
!   Thompson/WSM6 cloud microphysics scheme.                !
!                                                                       !
! abstract:  this program computes cloud fractions from cloud           !
!   condensates,                                                        !
!   and computes the low, mid, high, total and boundary layer cloud     !
!   fractions and the vertical indices of low, mid, and high cloud      !
!   top and base.  the three vertical cloud domains are set up in the   !
!   initial subroutine "cld_init".                                      !
!                                                                       !
! usage:         call progcld5                                          !
!                                                                       !
! subprograms called:   gethml                                          !
!                                                                       !
! attributes:                                                           !
!   language:   fortran 90                                              !
!   machine:    ibm-sp, sgi                                             !
!                                                                       !
!                                                                       !
!  ====================  definition of variables  ====================  !
!                                                                       !
! input variables:                                                      !
!   plyr  (IX,NLAY) : model layer mean pressure in mb (100Pa)           !
!   plvl  (IX,NLP1) : model level pressure in mb (100Pa)                !
!   tlyr  (IX,NLAY) : model layer mean temperature in k                 !
!   tvly  (IX,NLAY) : model layer virtual temperature in k              !
!   qlyr  (IX,NLAY) : layer specific humidity in gm/gm                  !
!   qstl  (IX,NLAY) : layer saturate humidity in gm/gm                  !
!   rhly  (IX,NLAY) : layer relative humidity (=qlyr/qstl)              !
!   clw   (IX,NLAY,ntrac) : layer cloud condensate amount               !
!   xlat  (IX)      : grid latitude in radians, default to pi/2 -> -pi/2!
!                     range, otherwise see in-line comment              !
!   xlon  (IX)      : grid longitude in radians  (not used)             !
!   slmsk (IX)      : sea/land mask array (sea:0,land:1,sea-ice:2)      !
!   dz    (ix,nlay) : layer thickness (km)                              !
!   delp  (ix,nlay) : model layer pressure thickness in mb (100Pa)      !
!   IX              : horizontal dimention                              !
!   NLAY,NLP1       : vertical layer/level dimensions                   !
!   uni_cld         : logical - true for cloud fraction from shoc       !
!   lmfshal         : logical - true for mass flux shallow convection   !
!   lmfdeep2        : logical - true for mass flux deep convection      !
!   cldcov          : layer cloud fraction (used when uni_cld=.true.    !
!                                                                       !
! output variables:                                                     !
!   clouds(IX,NLAY,NF_CLDS) : cloud profiles                            !
!      clouds(:,:,1) - layer total cloud fraction                       !
!      clouds(:,:,2) - layer cloud liq water path         (g/m**2)      !
!      clouds(:,:,3) - mean eff radius for liq cloud      (micron)      !
!      clouds(:,:,4) - layer cloud ice water path         (g/m**2)      !
!      clouds(:,:,5) - mean eff radius for ice cloud      (micron)      !
!      clouds(:,:,6) - layer rain drop water path         not assigned  !
!      clouds(:,:,7) - mean eff radius for rain drop      (micron)      !
!  *** clouds(:,:,8) - layer snow flake water path        not assigned  !
!      clouds(:,:,9) - mean eff radius for snow flake     (micron)      !
!  *** fu's scheme need to be normalized by snow density (g/m**3/1.0e6) !
!   clds  (IX,5)    : fraction of clouds for low, mid, hi, tot, bl      !
!   mtop  (IX,3)    : vertical indices for low, mid, hi cloud tops      !
!   mbot  (IX,3)    : vertical indices for low, mid, hi cloud bases     !
!   de_lgth(ix)     : clouds decorrelation length (km)                  !
!                                                                       !
! module variables:                                                     !
!   ivflip          : control flag of vertical index direction          !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!   lmfshal         : mass-flux shallow conv scheme flag                !
!   lmfdeep2        : scale-aware mass-flux deep conv scheme flag       !
!   lcrick          : control flag for eliminating CRICK                !
!                     =t: apply layer smoothing to eliminate CRICK      !
!                     =f: do not apply layer smoothing                  !
!   lcnorm          : control flag for in-cld condensate                !
!                     =t: normalize cloud condensate                    !
!                     =f: not normalize cloud condensate                !
!                                                                       !
!  ====================    end of description    =====================  !
!
      implicit none

!  ---  inputs
      integer,  intent(in) :: IX, NLAY, NLP1
      integer,  intent(in) :: ntrac, ntcw, ntiw, ntrw, ntsw, ntgl

      logical, intent(in)  :: uni_cld, lmfshal, lmfdeep2

      real (kind=kind_phys), dimension(:,:), intent(in) :: plvl, plyr,  &
     &       tlyr, qlyr, qstl, rhly, cldcov, delp, dz,                  &
     &       re_cloud, re_ice, re_snow 

      real (kind=kind_phys), dimension(:,:,:), intent(in) :: clw

      real (kind=kind_phys), dimension(:),   intent(in) :: xlat, xlon,  &
     &       slmsk

!  ---  outputs
      real (kind=kind_phys), dimension(:,:,:), intent(out) :: clouds

      real (kind=kind_phys), dimension(:,:),   intent(out) :: clds
      real (kind=kind_phys), dimension(:),     intent(out) :: de_lgth

      integer,               dimension(:,:),   intent(out) :: mtop,mbot

!  ---  local variables:
      real (kind=kind_phys), dimension(IX,NLAY) :: cldtot, cldcnv,      &
     &       cwp, cip, crp, csp, rew, rei, res, rer, tem2d, clwf

      real (kind=kind_phys) :: ptop1(IX,NK_CLDS+1), rxlat(ix)

      real (kind=kind_phys) :: clwmin, clwm, clwt, onemrh, value,       &
     &       tem1, tem2, tem3

      integer :: i, k, id, nf

!  ---  constant values
!     real (kind=kind_phys), parameter :: xrc3 = 200.
      real (kind=kind_phys), parameter :: xrc3 = 100.

!
!===> ... begin here
!
      do nf=1,nf_clds
        do k=1,nlay
          do i=1,ix
            clouds(i,k,nf) = 0.0
          enddo
        enddo
      enddo
!     clouds(:,:,:) = 0.0

      do k = 1, NLAY
        do i = 1, IX
          cldtot(i,k) = 0.0
          cldcnv(i,k) = 0.0
          cwp   (i,k) = 0.0
          cip   (i,k) = 0.0
          crp   (i,k) = 0.0
          csp   (i,k) = 0.0
          rew   (i,k) = re_cloud(i,k)
          rei   (i,k) = re_ice(i,k) 
          rer   (i,k) = rrain_def            ! default rain radius to 1000 micron
          res   (i,k) = re_snow(i,K) 
!         tem2d (i,k) = min( 1.0, max( 0.0, (con_ttp-tlyr(i,k))*0.05 ) )
          clwf(i,k)   = 0.0
        enddo
      enddo
!
!      if ( lcrick ) then
!        do i = 1, IX
!          clwf(i,1)    = 0.75*clw(i,1)    + 0.25*clw(i,2)
!          clwf(i,nlay) = 0.75*clw(i,nlay) + 0.25*clw(i,nlay-1)
!        enddo
!        do k = 2, NLAY-1
!          do i = 1, IX
!            clwf(i,K) = 0.25*clw(i,k-1) + 0.5*clw(i,k) + 0.25*clw(i,k+1)
!          enddo
!        enddo
!      else
!        do k = 1, NLAY
!          do i = 1, IX
!            clwf(i,k) = clw(i,k)
!          enddo
!        enddo
!      endif

        do k = 1, NLAY
          do i = 1, IX
            clwf(i,k) = clw(i,k,ntcw) +  clw(i,k,ntiw) + clw(i,k,ntsw)
          enddo
        enddo
!> -# Find top pressure for each cloud domain for given latitude.
!     ptopc(k,i): top presure of each cld domain (k=1-4 are sfc,L,m,h;
!  ---  i=1,2 are low-lat (<45 degree) and pole regions)

      do i =1, IX
        rxlat(i) = abs( xlat(i) / con_pi )      ! if xlat in pi/2 -> -pi/2 range
!       rxlat(i) = abs(0.5 - xlat(i)/con_pi)    ! if xlat in 0 -> pi range
      enddo

      do id = 1, 4
        tem1 = ptopc(id,2) - ptopc(id,1)

        do i =1, IX
          ptop1(i,id) = ptopc(id,1) + tem1*max( 0.0, 4.0*rxlat(i)-1.0 )
        enddo
      enddo

!> -# Compute cloud liquid/ice condensate path in \f$ g/m^2 \f$ .

        do k = 1, NLAY
          do i = 1, IX
            cwp(i,k) = max(0.0, clw(i,k,ntcw) * gfac * delp(i,k))
            cip(i,k) = max(0.0, clw(i,k,ntiw) * gfac * delp(i,k))
            crp(i,k) = max(0.0, clw(i,k,ntrw) * gfac * delp(i,k))
            csp(i,k) = max(0.0, (clw(i,k,ntsw)+clw(i,k,ntgl)) *         &
     &                  gfac * delp(i,k))
          enddo
        enddo

      if (uni_cld) then     ! use unified sgs clouds generated outside
        do k = 1, NLAY
          do i = 1, IX
            cldtot(i,k) = cldcov(i,k)
          enddo
        enddo

      else

!> -# Calculate layer cloud fraction.

        clwmin = 0.0
        if (.not. lmfshal) then
          do k = 1, NLAY
          do i = 1, IX
            clwt = 1.0e-6 * (plyr(i,k)*0.001)
!           clwt = 2.0e-6 * (plyr(i,k)*0.001)

            if (clwf(i,k) > clwt) then

              onemrh= max( 1.e-10, 1.0-rhly(i,k) )
              clwm  = clwmin / max( 0.01, plyr(i,k)*0.001 )

              tem1  = min(max(sqrt(sqrt(onemrh*qstl(i,k))),0.0001),1.0)
              tem1  = 2000.0 / tem1

!             tem1  = 1000.0 / tem1

              value = max( min( tem1*(clwf(i,k)-clwm), 50.0 ), 0.0 )
              tem2  = sqrt( sqrt(rhly(i,k)) )

              cldtot(i,k) = max( tem2*(1.0-exp(-value)), 0.0 )
            endif
          enddo
          enddo
        else
          do k = 1, NLAY
          do i = 1, IX
            clwt = 1.0e-6 * (plyr(i,k)*0.001)
!           clwt = 2.0e-6 * (plyr(i,k)*0.001)

            if (clwf(i,k) > clwt) then
              onemrh= max( 1.e-10, 1.0-rhly(i,k) )
              clwm  = clwmin / max( 0.01, plyr(i,k)*0.001 )
!
              tem1  = min(max((onemrh*qstl(i,k))**0.49,0.0001),1.0)  !jhan
              if (lmfdeep2) then
                tem1  = xrc3 / tem1
              else
                tem1  = 100.0 / tem1
              endif
!
              value = max( min( tem1*(clwf(i,k)-clwm), 50.0 ), 0.0 )
              tem2  = sqrt( sqrt(rhly(i,k)) )

              cldtot(i,k) = max( tem2*(1.0-exp(-value)), 0.0 )
            endif
          enddo
          enddo
        endif

      endif                                ! if (uni_cld) then

      do k = 1, NLAY
        do i = 1, IX
          if (cldtot(i,k) < climit) then
            cldtot(i,k) = 0.0
            cwp(i,k)    = 0.0
            cip(i,k)    = 0.0
            crp(i,k)    = 0.0
            csp(i,k)    = 0.0
          endif
        enddo
      enddo

      if ( lcnorm ) then
        do k = 1, NLAY
          do i = 1, IX
            if (cldtot(i,k) >= climit) then
              tem1 = 1.0 / max(climit2, cldtot(i,k))
              cwp(i,k) = cwp(i,k) * tem1
              cip(i,k) = cip(i,k) * tem1
              crp(i,k) = crp(i,k) * tem1
              csp(i,k) = csp(i,k) * tem1
            endif
          enddo
        enddo
      endif

!
      do k = 1, NLAY
        do i = 1, IX
          clouds(i,k,1) = cldtot(i,k)
          clouds(i,k,2) = cwp(i,k)
          clouds(i,k,3) = rew(i,k)
          clouds(i,k,4) = cip(i,k)
          clouds(i,k,5) = rei(i,k)
          clouds(i,k,6) = crp(i,k)  ! added for Thompson 
          clouds(i,k,7) = rer(i,k)
          clouds(i,k,8) = csp(i,k)  ! added for Thompson 
          clouds(i,k,9) = rei(i,k)
        enddo
      enddo

!  --- ...  estimate clouds decorrelation length in km
!           this is only a tentative test, need to consider change later

      if ( iovr == 3 ) then
        do i = 1, ix
          de_lgth(i) = max( 0.6, 2.78-4.6*rxlat(i) )
        enddo
      endif

!> -# Call gethml() to compute low,mid,high,total, and boundary layer
!!    cloud fractions and clouds top/bottom layer indices for low, mid,
!!    and high clouds.
!  ---  compute low, mid, high, total, and boundary layer cloud fractions
!       and clouds top/bottom layer indices for low, mid, and high clouds.
!       The three cloud domain boundaries are defined by ptopc.  The cloud
!       overlapping method is defined by control flag 'iovr', which may
!       be different for lw and sw radiation programs.

      call gethml                                                       &
!  ---  inputs:
     &     ( plyr, ptop1, cldtot, cldcnv, dz, de_lgth,                  &
     &       IX,NLAY,                                                   &
!  ---  outputs:
     &       clds, mtop, mbot                                           &
     &     )


!
      return
!...................................
      end subroutine progcld5
!...................................

!> @}

!> This subroutine computes cloud related quantities using 
!! zhao/moorthi's prognostic cloud microphysics scheme.
!!\param plyr    (IX,NLAY), model layer mean pressure in mb (100Pa)
!!\param plvl    (IX,NLP1), model level pressure in mb (100Pa)
!!\param tlyr    (IX,NLAY), model layer mean temperature in K
!!\param tvly    (IX,NLAY), model layer virtual temperature in K
!!\param qlyr    (IX,NLAY), layer specific humidity in gm/gm
!!\param clw     (IX,NLAY), layer cloud liquid water amount
!!\param ciw     (IX,NLAY), layer cloud ice water amount
!!\param xlat    (IX), grid latitude in radians, default to pi/2 -> 
!!               -pi/2 range, otherwise see in-line comment
!!\param xlon    (IX), grid longitude in radians  (not used)
!!\param slmsk   (IX), sea/land mask array (sea:0,land:1,sea-ice:2)
!!\param dz      (IX,NLAY), layer thickness (km)
!!\param delp    (IX,NLAY), model layer pressure thickness in mb (100Pa)
!!\param IX           horizontal dimention
!!\param NLAY,NLP1    vertical layer/level dimensions
!!\param clouds      (IX,NLAY,NF_CLDS), cloud profiles
!!\n                 (:,:,1) - layer total cloud fraction
!!\n                 (:,:,2) - layer cloud liq water path \f$(g/m^2)\f$
!!\n                 (:,:,3) - mean eff radius for liq cloud (micron)
!!\n                 (:,:,4) - layer cloud ice water path \f$(g/m^2)\f$
!!\n                 (:,:,5) - mean eff radius for ice cloud (micron)
!!\n                 (:,:,6) - layer rain drop water path (not assigned)
!!\n                 (:,:,7) - mean eff radius for rain drop (micron)
!!\n                 (:,:,8) - layer snow flake water path (not assigned)
!!\n                 (:,:,9) - mean eff radius for snow flake (micron)
!!\n  *** fu's scheme need to be normalized by snow density \f$ (g/m^3/1.0e6)\f$
!!\param clds       (IX,5), fraction of clouds for low, mid, hi, tot, bl
!!\param mtop       (IX,3), vertical indices for low, mid, hi cloud tops
!!\param mbot       (IX,3), vertical indices for low, mid, hi cloud bases
!!\param de_lgth    (IX),   clouds decorrelation length (km)
!>\section gen_progclduni General Algorithm
!> @{
!-----------------------------------
      subroutine progclduni                                             &
     &     ( plyr,plvl,tlyr,tvly,ccnd,ncnd,                             &    !  ---  inputs:
     &       xlat,xlon,slmsk,dz,delp, IX, NLAY, NLP1, cldtot,           &
     &       effrl,effri,effrr,effrs,effr_in,                           &
     &       clouds,clds,mtop,mbot,de_lgth                              &    !  ---  outputs:
     &      )

! =================   subprogram documentation block   ================ !
!                                                                       !
! subprogram:    progclduni    computes cloud related quantities using    !
!   for unified cloud microphysics scheme.                !
!                                                                       !
! abstract:  this program computes cloud fractions from cloud           !
!   condensates, calculates liquid/ice cloud droplet effective radius,  !
!   and computes the low, mid, high, total and boundary layer cloud     !
!   fractions and the vertical indices of low, mid, and high cloud      !
!   top and base.  the three vertical cloud domains are set up in the   !
!   initial subroutine "cld_init".                                      !
!                                                                       !
! usage:         call progclduni                                          !
!                                                                       !
! subprograms called:   gethml                                          !
!                                                                       !
! attributes:                                                           !
!   language:   fortran 90                                              !
!   machine:    ibm-sp, sgi                                             !
!                                                                       !
!                                                                       !
!  ====================  definition of variables  ====================  !
!                                                                       !
! input variables:                                                      !
!   plyr  (IX,NLAY) : model layer mean pressure in mb (100Pa)           !
!   plvl  (IX,NLP1) : model level pressure in mb (100Pa)                !
!   tlyr  (IX,NLAY) : model layer mean temperature in k                 !
!   tvly  (IX,NLAY) : model layer virtual temperature in k              !
!   ccnd  (IX,NLAY) : layer cloud condensate amount                     !
!   ncnd            : number of layer cloud condensate types            !
!   xlat  (IX)      : grid latitude in radians, default to pi/2 -> -pi/2!
!                     range, otherwise see in-line comment              !
!   xlon  (IX)      : grid longitude in radians  (not used)             !
!   slmsk (IX)      : sea/land mask array (sea:0,land:1,sea-ice:2)      !
!   IX              : horizontal dimention                              !
!   NLAY,NLP1       : vertical layer/level dimensions                   !
!   cldtot          : unified cloud fracrion from moist physics         !
!   effrl (ix,nlay) : effective radius for liquid water                 !
!   effri (ix,nlay) : effective radius for ice water                    !
!   effrr (ix,nlay) : effective radius for rain water                   !
!   effrs (ix,nlay) : effective radius for snow water                   !
!   effr_in         : logical - if .true. use input effective radii     !
!   dz    (ix,nlay) : layer thickness (km)                              !
!   delp  (ix,nlay) : model layer pressure thickness in mb (100Pa)      !
!                                                                       !
! output variables:                                                     !
!   clouds(IX,NLAY,NF_CLDS) : cloud profiles                            !
!      clouds(:,:,1) - layer total cloud fraction                       !
!      clouds(:,:,2) - layer cloud liq water path         (g/m**2)      !
!      clouds(:,:,3) - mean eff radius for liq cloud      (micron)      !
!      clouds(:,:,4) - layer cloud ice water path         (g/m**2)      !
!      clouds(:,:,5) - mean eff radius for ice cloud      (micron)      !
!      clouds(:,:,6) - layer rain drop water path         not assigned  !
!      clouds(:,:,7) - mean eff radius for rain drop      (micron)      !
!  *** clouds(:,:,8) - layer snow flake water path        not assigned  !
!      clouds(:,:,9) - mean eff radius for snow flake     (micron)      !
!  *** fu's scheme need to be normalized by snow density (g/m**3/1.0e6) !
!   clds  (IX,5)    : fraction of clouds for low, mid, hi, tot, bl      !
!   mtop  (IX,3)    : vertical indices for low, mid, hi cloud tops      !
!   mbot  (IX,3)    : vertical indices for low, mid, hi cloud bases     !
!   de_lgth(ix)     : clouds decorrelation length (km)                  !
!                                                                       !
! module variables:                                                     !
!   ivflip          : control flag of vertical index direction          !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!   lmfshal         : mass-flux shallow conv scheme flag                !
!   lmfdeep2        : scale-aware mass-flux deep conv scheme flag       !
!   lcrick          : control flag for eliminating CRICK                !
!                     =t: apply layer smoothing to eliminate CRICK      !
!                     =f: do not apply layer smoothing                  !
!   lcnorm          : control flag for in-cld condensate                !
!                     =t: normalize cloud condensate                    !
!                     =f: not normalize cloud condensate                !
!                                                                       !
!  ====================    end of description    =====================  !
!
      implicit none

!  ---  inputs
      integer,  intent(in) :: IX, NLAY, NLP1, ncnd
      logical,  intent(in) :: effr_in

      real (kind=kind_phys), dimension(:,:,:), intent(in) :: ccnd
      real (kind=kind_phys), dimension(:,:),   intent(in) :: plvl, plyr,&
     &       tlyr, tvly, cldtot, effrl, effri, effrr, effrs, dz, delp

      real (kind=kind_phys), dimension(:),   intent(in) :: xlat, xlon,  &
     &       slmsk

!  ---  outputs
      real (kind=kind_phys), dimension(:,:,:), intent(out) :: clouds

      real (kind=kind_phys), dimension(:,:),   intent(out) :: clds

      real (kind=kind_phys), dimension(:),     intent(out) :: de_lgth

      integer,               dimension(:,:),   intent(out) :: mtop,mbot

!  ---  local variables:
      real (kind=kind_phys), dimension(IX,NLAY) :: cldcnv, cwp, cip,    &
     &       crp, csp, rew, rei, res, rer, tem2d
      real (kind=kind_phys), dimension(IX,NLAY,ncnd) :: cndf

      real (kind=kind_phys) :: ptop1(IX,NK_CLDS+1), rxlat(ix)

      real (kind=kind_phys) :: tem1, tem2, tem3

      integer :: i, k, id, nf, n

!
!===> ... begin here
!
      do nf=1,nf_clds
        do k=1,nlay
          do i=1,ix
            clouds(i,k,nf) = 0.0
          enddo
        enddo
      enddo
!     clouds(:,:,:) = 0.0

      if (effr_in) then
        do k = 1, NLAY
          do i = 1, IX
            cldcnv(i,k) = 0.0
            cwp   (i,k) = 0.0
            cip   (i,k) = 0.0
            crp   (i,k) = 0.0
            csp   (i,k) = 0.0
            rew   (i,k) = effrl (i,k)
            rei   (i,k) = max(10.0, min(150.0,effri (i,k)))
            rer   (i,k) = effrr (i,k)
            res   (i,k) = effrs (i,k)
            tem2d (i,k) = min( 1.0, max( 0.0,(con_ttp-tlyr(i,k))*0.05))
          enddo
        enddo
      else
        do k = 1, NLAY
          do i = 1, IX
            cldcnv(i,k) = 0.0
            cwp   (i,k) = 0.0
            cip   (i,k) = 0.0
            crp   (i,k) = 0.0
            csp   (i,k) = 0.0
            rew   (i,k) = reliq_def            ! default liq radius to 10 micron
            rei   (i,k) = reice_def            ! default ice radius to 50 micron
            rer   (i,k) = rrain_def            ! default rain radius to 1000 micron
            res   (i,k) = rsnow_def            ! default snow radius to 250 micron
            tem2d (i,k) = min(1.0, max(0.0, (con_ttp-tlyr(i,k))*0.05))
          enddo
        enddo
      endif
!
      do n=1,ncnd
        do k = 1, NLAY
          do i = 1, IX
            cndf(i,k,n) = ccnd(i,k,n)
          enddo
        enddo
      enddo
      if ( lcrick ) then
        do n=1,ncnd
          do i = 1, IX
            cndf(i,1,n)    = 0.75*ccnd(i,1,n)    + 0.25*ccnd(i,2,n)
            cndf(i,nlay,n) = 0.75*ccnd(i,nlay,n) + 0.25*ccnd(i,nlay-1,n)
          enddo
          do k = 2, NLAY-1
            do i = 1, IX
              cndf(i,K,n) = 0.25 * (ccnd(i,k-1,n) + ccnd(i,k+1,n))      &
     &                    + 0.5  *  ccnd(i,k,n)
            enddo
          enddo
        enddo
      endif

!> -# Find top pressure for each cloud domain for given latitude.
!     ptopc(k,i): top presure of each cld domain (k=1-4 are sfc,L,m,h;
!  ---  i=1,2 are low-lat (<45 degree) and pole regions)

      do i =1, IX
        rxlat(i) = abs( xlat(i) / con_pi )     ! if xlat in pi/2 -> -pi/2 range
!       rxlat(i) = abs(0.5 - xlat(i)/con_pi)   ! if xlat in 0 -> pi range
      enddo

      do id = 1, 4
        tem1 = ptopc(id,2) - ptopc(id,1)

        do i =1, IX
          ptop1(i,id) = ptopc(id,1) + tem1*max( 0.0, 4.0*rxlat(i)-1.0 )
        enddo
      enddo

!> -# Compute cloud liquid/ice condensate path in \f$ g/m^2 \f$ .

        if (ncnd == 2) then
          do k = 1, NLAY
            do i = 1, IX
              tem1      = gfac * delp(i,k)
              cwp(i,k)  = cndf(i,k,1) * tem1
              cip(i,k)  = cndf(i,k,2) * tem1
            enddo
          enddo
        elseif (ncnd == 4 .or. ncnd == 5) then
          do k = 1, NLAY
            do i = 1, IX
              tem1      = gfac * delp(i,k)
              cwp(i,k)  = cndf(i,k,1) * tem1
              cip(i,k)  = cndf(i,k,2) * tem1
              crp(i,k)  = cndf(i,k,3) * tem1
              csp(i,k)  = cndf(i,k,4) * tem1
            enddo
          enddo
        endif

!> -# Compute effective liquid cloud droplet radius over land.

      if(.not. effr_in) then
        do i = 1, IX
          if (nint(slmsk(i)) == 1) then
            do k = 1, NLAY
              rew(i,k) = 5.0 + 5.0 * tem2d(i,k)
            enddo
          endif
        enddo
      endif

      do k = 1, NLAY
        do i = 1, IX
          if (cldtot(i,k) < climit) then
            cwp(i,k)    = 0.0
            cip(i,k)    = 0.0
            crp(i,k)    = 0.0
            csp(i,k)    = 0.0
          endif
        enddo
      enddo

      if ( lcnorm ) then
        do k = 1, NLAY
          do i = 1, IX
            if (cldtot(i,k) >= climit) then
              tem1 = 1.0 / max(climit2, cldtot(i,k))
              cwp(i,k) = cwp(i,k) * tem1
              cip(i,k) = cip(i,k) * tem1
              crp(i,k) = crp(i,k) * tem1
              csp(i,k) = csp(i,k) * tem1
            endif
          enddo
        enddo
      endif

!> -# Compute effective ice cloud droplet radius following Heymsfield 
!!    and McFarquhar (1996) \cite heymsfield_and_mcfarquhar_1996.

      if(.not. effr_in) then
        do k = 1, NLAY
          do i = 1, IX
            tem2 = tlyr(i,k) - con_ttp

            if (cip(i,k) > 0.0) then
              tem3 = gord * cip(i,k) * plyr(i,k) / (delp(i,k)*tvly(i,k))

              if (tem2 < -50.0) then
                rei(i,k) = (1250.0/9.917) * tem3 ** 0.109
              elseif (tem2 < -40.0) then
                rei(i,k) = (1250.0/9.337) * tem3 ** 0.08
              elseif (tem2 < -30.0) then
                rei(i,k) = (1250.0/9.208) * tem3 ** 0.055
              else
                rei(i,k) = (1250.0/9.387) * tem3 ** 0.031
              endif
!             rei(i,k)   = max(20.0, min(rei(i,k), 300.0))
!             rei(i,k)   = max(10.0, min(rei(i,k), 100.0))
              rei(i,k)   = max(10.0, min(rei(i,k), 150.0))
!             rei(i,k)   = max(5.0,  min(rei(i,k), 130.0))
            endif
          enddo
        enddo
      endif

!
      do k = 1, NLAY
        do i = 1, IX
          clouds(i,k,1) = cldtot(i,k)
          clouds(i,k,2) = cwp(i,k)
          clouds(i,k,3) = rew(i,k)
          clouds(i,k,4) = cip(i,k)
          clouds(i,k,5) = rei(i,k)
          clouds(i,k,6) = crp(i,k)
          clouds(i,k,7) = rer(i,k)
          clouds(i,k,8) = csp(i,k)
          clouds(i,k,9) = res(i,k)
        enddo
      enddo

!  --- ...  estimate clouds decorrelation length in km
!           this is only a tentative test, need to consider change later

      if ( iovr == 3 ) then
        do i = 1, ix
          de_lgth(i) = max( 0.6, 2.78-4.6*rxlat(i) )
        enddo
      endif

!> -# Call gethml() to compute low,mid,high,total, and boundary layer
!!    cloud fractions and clouds top/bottom layer indices for low, mid,
!!    and high clouds.
!  ---  compute low, mid, high, total, and boundary layer cloud fractions
!       and clouds top/bottom layer indices for low, mid, and high clouds.
!       The three cloud domain boundaries are defined by ptopc.  The cloud
!       overlapping method is defined by control flag 'iovr', which may
!       be different for lw and sw radiation programs.

      call gethml                                                       &
!  ---  inputs:
     &     ( plyr, ptop1, cldtot, cldcnv, dz, de_lgth,                  &
     &       IX,NLAY,                                                   &
!  ---  outputs:
     &       clds, mtop, mbot                                           &
     &     )


!
      return
!...................................
      end subroutine progclduni
!-----------------------------------
!> @}


!> This subroutine computes high, mid, low, total, and boundary cloud
!! fractions and cloud top/bottom layer indices for model diagnostic
!! output. The three cloud domain boundaries are defined by ptopc. The
!! cloud overlapping method is defined by control flag 'iovr', which is
!! also used by LW and SW radiation programs.                          
!> \param plyr    (IX,NLAY), model layer mean pressure in mb (100Pa)
!> \param ptop1   (IX,4), pressure limits of cloud domain interfaces
!!                    (sfc,low,mid,high) in mb (100Pa)
!> \param cldtot  (IX,NLAY), total or stratiform cloud profile in fraction
!> \param cldcnv  (IX,NLAY), convective cloud (for diagnostic scheme only)
!> \param dz      (IX,NLAY), layer thickness (km)
!> \param de_lgth (IX),  clouds decorrelation length (km)
!> \param IX      horizontal dimension
!> \param NLAY    vertical layer dimensions
!> \param clds   (IX,5), fraction of clouds for low, mid, hi, tot, bl 
!> \param mtop   (IX,3),vertical indices for low, mid, hi cloud tops 
!> \param mbot   (IX,3),vertical indices for low, mid, hi cloud bases
!!
!>\section detail Detailed Algorithm
!! @{
!-----------------------------------                                    !
      subroutine gethml                                                 &
     &     ( plyr, ptop1, cldtot, cldcnv, dz, de_lgth,                  &       !  ---  inputs:
     &       IX, NLAY,                                                  &
     &       clds, mtop, mbot                                           &       !  ---  outputs:
     &     )

!  ===================================================================  !
!                                                                       !
! abstract: compute high, mid, low, total, and boundary cloud fractions !
!   and cloud top/bottom layer indices for model diagnostic output.     !
!   the three cloud domain boundaries are defined by ptopc.  the cloud  !
!   overlapping method is defined by control flag 'iovr', which is also !
!   used by lw and sw radiation programs.                               !
!                                                                       !
! usage:         call gethml                                            !
!                                                                       !
! subprograms called:  none                                             !
!                                                                       !
! attributes:                                                           !
!   language:   fortran 90                                              !
!   machine:    ibm-sp, sgi                                             !
!                                                                       !
!                                                                       !
!  ====================  definition of variables  ====================  !
!                                                                       !
! input variables:                                                      !
!   plyr  (IX,NLAY) : model layer mean pressure in mb (100Pa)           !
!   ptop1 (IX,4)    : pressure limits of cloud domain interfaces        !
!                     (sfc,low,mid,high) in mb (100Pa)                  !
!   cldtot(IX,NLAY) : total or straiform cloud profile in fraction      !
!   cldcnv(IX,NLAY) : convective cloud (for diagnostic scheme only)     !
!   dz    (ix,nlay) : layer thickness (km)                              !
!   de_lgth(ix)     : clouds vertical de-correlation length (km)        !
!   IX              : horizontal dimention                              !
!   NLAY            : vertical layer dimensions                         !
!                                                                       !
! output variables:                                                     !
!   clds  (IX,5)    : fraction of clouds for low, mid, hi, tot, bl      !
!   mtop  (IX,3)    : vertical indices for low, mid, hi cloud tops      !
!   mbot  (IX,3)    : vertical indices for low, mid, hi cloud bases     !
!                                                                       !
! external module variables:  (in physparam)                            !
!   ivflip          : control flag of vertical index direction          !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!                                                                       !
! internal module variables:                                            !
!   iovr            : control flag for cloud overlap                    !
!                     =0 random overlapping clouds                      !
!                     =1 max/ran overlapping clouds                     !
!                     =2 maximum overlapping  ( for mcica only )        !
!                     =3 decorr-length ovlp   ( for mcica only )        !
!                                                                       !
!  ====================    end of description    =====================  !
!
      implicit none!

!  ---  inputs:
      integer, intent(in) :: IX, NLAY

      real (kind=kind_phys), dimension(:,:), intent(in) :: plyr, ptop1, &
     &       cldtot, cldcnv, dz
      real (kind=kind_phys), dimension(:),   intent(in) :: de_lgth

!  ---  outputs
      real (kind=kind_phys), dimension(:,:), intent(out) :: clds

      integer,               dimension(:,:), intent(out) :: mtop, mbot

!  ---  local variables:
      real (kind=kind_phys) :: cl1(IX), cl2(IX), dz1(ix)
      real (kind=kind_phys) :: pcur, pnxt, ccur, cnxt, alfa

      integer, dimension(IX):: idom, kbt1, kth1, kbt2, kth2
      integer :: i, k, id, id1, kstr, kend, kinc

!
!===> ... begin here
!
      clds(:,:) = 0.0

      do i = 1, IX
        cl1(i) = 1.0
        cl2(i) = 1.0
      enddo

!  ---  total and bl clouds, where cl1, cl2 are fractions of clear-sky view
!       layer processed from surface and up

!> -# Calculate total and BL cloud fractions (maximum-random cloud
!!    overlapping is operational).

      if ( ivflip == 0 ) then                   ! input data from toa to sfc
        kstr = NLAY
        kend = 1
        kinc = -1
      else                                      ! input data from sfc to toa
        kstr = 1
        kend = NLAY
        kinc = 1
      endif                                     ! end_if_ivflip

      if ( iovr == 0 ) then                     ! random overlap

        do k = kstr, kend, kinc
          do i = 1, IX
            ccur = min( ovcst, max( cldtot(i,k), cldcnv(i,k) ))
            if (ccur >= climit) cl1(i) = cl1(i) * (1.0 - ccur)
          enddo

          if (k == llyr) then
            do i = 1, IX
              clds(i,5) = 1.0 - cl1(i)          ! save bl cloud
            enddo
          endif
        enddo

        do i = 1, IX
          clds(i,4) = 1.0 - cl1(i)              ! save total cloud
        enddo

      elseif ( iovr == 1 ) then                 ! max/ran overlap

        do k = kstr, kend, kinc
          do i = 1, IX
            ccur = min( ovcst, max( cldtot(i,k), cldcnv(i,k) ))
            if (ccur >= climit) then             ! cloudy layer
              cl2(i) = min( cl2(i), (1.0 - ccur) )
            else                                ! clear layer
              cl1(i) = cl1(i) * cl2(i)
              cl2(i) = 1.0
            endif
          enddo

          if (k == llyr) then
            do i = 1, IX
              clds(i,5) = 1.0 - cl1(i) * cl2(i) ! save bl cloud
            enddo
          endif
        enddo

        do i = 1, IX
          clds(i,4) = 1.0 - cl1(i) * cl2(i)     ! save total cloud
        enddo

      elseif ( iovr == 2 ) then                 ! maximum overlap all levels

        cl1(:) = 0.0

        do k = kstr, kend, kinc
          do i = 1, IX
            ccur = min( ovcst,  max( cldtot(i,k), cldcnv(i,k) ))
            if (ccur >= climit) cl1(i) = max( cl1(i), ccur )
          enddo

          if (k == llyr) then
            do i = 1, IX
              clds(i,5) = cl1(i)    ! save bl cloud
            enddo
          endif
        enddo

        do i = 1, IX
          clds(i,4) = cl1(i)        ! save total cloud
        enddo

      elseif ( iovr == 3 ) then                 ! random if clear-layer divided,
                                                ! otherwise de-corrlength method
        do i = 1, ix
          dz1(i) = - dz(i,kstr)
        enddo

        do k = kstr, kend, kinc
          do i = 1, ix
            ccur = min( ovcst, max( cldtot(i,k), cldcnv(i,k) ))
            if (ccur >= climit) then                           ! cloudy layer
              alfa = exp( -0.5*((dz1(i)+dz(i,k)))/de_lgth(i) )
              dz1(i) = dz(i,k)
              cl2(i) =    alfa      * min(cl2(i), (1.0 - ccur))         & ! maximum part
     &               + (1.0 - alfa) * (cl2(i) * (1.0 - ccur))             ! random part
            else                                               ! clear layer
              cl1(i) = cl1(i) * cl2(i)
              cl2(i) = 1.0
              if (k /= kend) dz1(i) = -dz(i,k+kinc)
            endif
          enddo

          if (k == llyr) then
            do i = 1, ix
              clds(i,5) = 1.0 - cl1(i) * cl2(i) ! save bl cloud
            enddo
          endif
        enddo

        do i = 1, ix
          clds(i,4) = 1.0 - cl1(i) * cl2(i)     ! save total cloud
        enddo

      endif                                     ! end_if_iovr

!  ---  high, mid, low clouds, where cl1, cl2 are cloud fractions
!       layer processed from one layer below llyr and up
!  ---  change! layer processed from surface to top, so low clouds will
!       contains both bl and low clouds.

!> -# Calculte high, mid, low cloud fractions and vertical indices of
!!    cloud tops/bases.
      if ( ivflip == 0 ) then                   ! input data from toa to sfc

        do i = 1, IX
          cl1 (i) = 0.0
          cl2 (i) = 0.0
          kbt1(i) = NLAY
          kbt2(i) = NLAY
          kth1(i) = 0
          kth2(i) = 0
          idom(i) = 1
          mbot(i,1) = NLAY
          mtop(i,1) = NLAY
          mbot(i,2) = NLAY - 1
          mtop(i,2) = NLAY - 1
          mbot(i,3) = NLAY - 1
          mtop(i,3) = NLAY - 1
        enddo

!org    do k = llyr-1, 1, -1
        do k = NLAY, 1, -1
          do i = 1, IX
            id = idom(i)
            id1= id + 1

            pcur = plyr(i,k)
            ccur = min( ovcst, max( cldtot(i,k), cldcnv(i,k) ))

            if (k > 1) then
              pnxt = plyr(i,k-1)
              cnxt = min( ovcst, max( cldtot(i,k-1), cldcnv(i,k-1) ))
            else
              pnxt = -1.0
              cnxt = 0.0
            endif

            if (pcur < ptop1(i,id1)) then
              id = id + 1
              id1= id1 + 1
              idom(i) = id
            endif

            if (ccur >= climit) then
              if (kth2(i) == 0) kbt2(i) = k
              kth2(i) = kth2(i) + 1

              if ( iovr == 0 ) then
                cl2(i) = cl2(i) + ccur - cl2(i)*ccur
              else
                cl2(i) = max( cl2(i), ccur )
              endif

              if (cnxt < climit .or. pnxt < ptop1(i,id1)) then
                kbt1(i) = nint( (cl1(i)*kbt1(i) + cl2(i)*kbt2(i) )      &
     &                  / (cl1(i) + cl2(i)) )
                kth1(i) = nint( (cl1(i)*kth1(i) + cl2(i)*kth2(i) )      &
     &                  / (cl1(i) + cl2(i)) )
                cl1 (i) = cl1(i) + cl2(i) - cl1(i)*cl2(i)

                kbt2(i) = k - 1
                kth2(i) = 0
                cl2 (i) = 0.0
              endif   ! end_if_cnxt_or_pnxt
            endif     ! end_if_ccur

            if (pnxt < ptop1(i,id1)) then
              clds(i,id) = cl1(i)
              mtop(i,id) = min( kbt1(i), kbt1(i)-kth1(i)+1 )
              mbot(i,id) = kbt1(i)

              cl1 (i) = 0.0
              kbt1(i) = k - 1
              kth1(i) = 0

              if (id1 <= NK_CLDS) then
                mbot(i,id1) = kbt1(i)
                mtop(i,id1) = kbt1(i)
              endif
            endif     ! end_if_pnxt

          enddo       ! end_do_i_loop
        enddo         ! end_do_k_loop

      else                                      ! input data from sfc to toa

        do i = 1, IX
          cl1 (i) = 0.0
          cl2 (i) = 0.0
          kbt1(i) = 1
          kbt2(i) = 1
          kth1(i) = 0
          kth2(i) = 0
          idom(i) = 1
          mbot(i,1) = 1
          mtop(i,1) = 1
          mbot(i,2) = 2
          mtop(i,2) = 2
          mbot(i,3) = 2
          mtop(i,3) = 2
        enddo

!org    do k = llyr+1, NLAY
        do k = 1, NLAY
          do i = 1, IX
            id = idom(i)
            id1= id + 1

            pcur = plyr(i,k)
            ccur = min( ovcst, max( cldtot(i,k), cldcnv(i,k) ))

            if (k < NLAY) then
              pnxt = plyr(i,k+1)
              cnxt = min( ovcst, max( cldtot(i,k+1), cldcnv(i,k+1) ))
            else
              pnxt = -1.0
              cnxt = 0.0
            endif

            if (pcur < ptop1(i,id1)) then
              id = id + 1
              id1= id1 + 1
              idom(i) = id
            endif

            if (ccur >= climit) then
              if (kth2(i) == 0) kbt2(i) = k
              kth2(i) = kth2(i) + 1

              if ( iovr == 0 ) then
                cl2(i) = cl2(i) + ccur - cl2(i)*ccur
              else
                cl2(i) = max( cl2(i), ccur )
              endif

              if (cnxt < climit .or. pnxt < ptop1(i,id1)) then
                kbt1(i) = nint( (cl1(i)*kbt1(i) + cl2(i)*kbt2(i))       &
     &                  / (cl1(i) + cl2(i)) )
                kth1(i) = nint( (cl1(i)*kth1(i) + cl2(i)*kth2(i))       &
     &                  / (cl1(i) + cl2(i)) )
                cl1 (i) = cl1(i) + cl2(i) - cl1(i)*cl2(i)

                kbt2(i) = k + 1
                kth2(i) = 0
                cl2 (i) = 0.0
              endif     ! end_if_cnxt_or_pnxt
            endif       ! end_if_ccur

            if (pnxt < ptop1(i,id1)) then
              clds(i,id) = cl1(i)
              mtop(i,id) = max( kbt1(i), kbt1(i)+kth1(i)-1 )
              mbot(i,id) = kbt1(i)

              cl1 (i) = 0.0
              kbt1(i) = min(k+1, nlay)
              kth1(i) = 0

              if (id1 <= NK_CLDS) then
                mbot(i,id1) = kbt1(i)
                mtop(i,id1) = kbt1(i)
              endif
            endif     ! end_if_pnxt

          enddo       ! end_do_i_loop
        enddo         ! end_do_k_loop

      endif                                     ! end_if_ivflip

!
      return
!...................................
      end subroutine gethml
!-----------------------------------
!! @}

!
!........................................!
      end module module_radiation_clouds !
!========================================!
!> @}
