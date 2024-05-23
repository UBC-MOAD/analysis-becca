!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! - Ariane - (May - 2007)
!! 
!! bruno.blanke@univ-brest.fr and nicolas.grima@univ-brest.fr
!! 
!! This software is a computer program whose purpose is 
!! the computation of 3D streamlines in a given velocity field 
!! (as the output of an Ocean General Circulation Model) and 
!! subsequent water masses analyses.
!! 
!! This software is governed by the CeCILL license under French law and
!! abiding by the rules of distribution of free software.  You can  use, 
!! modify and/ or redistribute the software under the terms of the CeCILL
!! license as circulated by CEA, CNRS and INRIA at the following URL
!! "http://www.cecill.info". 
!! 
!! As a counterpart to the access to the source code and  rights to copy,
!! modify and redistribute granted by the license, users are provided only
!! with a limited warranty  and the software's author,  the holder of the
!! economic rights,  and the successive licensors  have only  limited
!! liability. 
!! 
!! In this respect, the user's attention is drawn to the risks associated
!! with loading,  using,  modifying and/or developing or reproducing the
!! software by the user in light of its specific status of free software,
!! that may mean  that it is complicated to manipulate,  and  that  also
!! therefore means  that it is reserved for developers  and  experienced
!! professionals having in-depth computer knowledge. Users are therefore
!! encouraged to load and test the software's suitability as regards their
!! requirements in conditions enabling the security of their systems and/or 
!! data to be ensured and,  more generally, to use and operate it in the 
!! same conditions as regards security. 
!! 
!! The fact that you are presently reading this means that you have had
!! knowledge of the CeCILL license and that you accept its terms.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!****f* ariane/trajec
!! NAME
!!   trajec (trajec.f90 - Main program - Fortran 90)
!!
!! USAGE
!!   
!!
!! FUNCTION
!! --- Particle trajectories from the analytical computation of 3D
!!     streamlines in a 3D velocity fiels.
!!
!! --- Qualitative or quantitative diagnostics
!!
!! --- Computations done over periods equal to the available sampling
!!     for the velocity field, and over which the velocity is assumed
!!     constant
!!
!! --- The program works with the OPA tensorial formalism
!!
!! --- input: simplified meshmask, 3D TRANSPORT field (velocity x
!!     section), parameters (on various files)
!!
!! --- output: individual trajectories or quantitative outputs
!!
!! --- no-vectorized version
!!
!! --- some comments (yes !) in English
!!
!! --- copyright: Bruno Blanke (Bruno.Blanke@univ-brest.fr)
!! first version: Summer 1992
!! this version: March 2007
!! http://www.ifremer.fr/lpo/blanke/ARIANE
!! 
!! AUTHOR
!!   * Origin  : Bruno Blanke  (1992)
!!   * F77toF90: Nicolas Grima (April-May 2005)
!! 
!! CREATION DATE
!!   * April-May 2005
!!
!! HISTORY
!!   Date (dd/mm/yyyy/) - Modification(s)
!!
!! RESULT
!!   
!!
!! EXAMPLES
!!   
!!
!! NOTES
!!
!! (see mod_lun)
!! current input files:
!! -------------------

!! current output files:
!! ---------------------
!!
!! TODO
!!   
!!
!! PORTABILITY
!!         Machine-OS    - Fortran90/95 compiler
!!   * i686-pc-linux-gnu -         ifort
!!
!! SEE ALSO
!!
!! SOURCE
!!=========================================================================
SUBROUTINE trajec_seq()

  USE mod_precision
  USE mod_cst
  USE mod_namelist
  USE mod_input_grid
  USE mod_seq
  USE mod_input_data
  USE mod_output_data
  USE mod_trajec_subs
  USE mod_init_particules
  USE mod_posin
  USE mod_flags
  USE mod_stats
  USE mod_stati
  USE mod_txt
  USE mod_lun
  USE mod_fx
  USE mod_fy
  USE mod_fz
  USE mod_zinter
  USE mod_sigma
  USE mod_quant
  USE mod_reducmem
  USE mod_orca
  USE mod_save_netcdf

  !--------------!
  ! DECLARATIONS !
  !--------------!
  IMPLICIT NONE 

  INTEGER(kind=iprec) :: &
       nsect2      

  INTEGER(kind=iprec) :: &
       CASE            , & !
       i, j, n         , & !
       is, js, ks      , & !
       nis, njs, nks   , & !
       i1s, j1s, k1s   , & !
       i2s, j2s, k2s   , & !
       ii, jj, kk      , & !
       ndone           , & !
       iter            , & !
       lref     = iOne , & !
       kk0             , & !
       ifl             , & !
                                !!       cnt      = 1    , & ! counter
       iter_min        , & !
       iter_max        , & ! 
       imin_tmp        , & !
       imax_tmp        , & !
       indice_traj        , & !! NG 29 april 2010
       flag_write_position    !! NG 29 april 2010

  INTEGER(kind=iprec), DIMENSION(:), ALLOCATABLE ::  &
       igo   , & !
       icolim, & !
       inext , & !
       i0    , & !
       j0    , & !
       k0    , & !
       l0    , & !
       i1    , & !
       i2    , & !
       j1    , & !
       j2    , & !
       k1    , & !
       k2    , & !
       iflag , & !
       iperfl, & !
       jperfl, & !
       ind_traj, & !
       nfin

  REAL(kind=rprec) :: &
       subtime =   -1._rprec  , & !
       tfic    =    rZero  , & !
       tlap    =    rZero  , & !
       tfin    =    rZero  , & !
       trtot1  =    rZero  , & !
       prof    =    rZero  , & !
       tol     =    2e-3    , & !
       x0   , & !
       y0   , & !
       z0   , & !
       xs   , & !
       ys   , & !
       zs   , & !
       ft   , & !
       xt   , & !
       yt   , & !
       zt   , & !
       st   , & !
       r_lmt, & !
       dumm, & !
       timestep, & !
       icount

  REAL(kind=rprec), DIMENSION(:), ALLOCATABLE :: & !
       age  , & !
       agenew, & !
       truage, & !
       gi   , & !
       gj   , & !
       gk   , & !
       gl   , & !
       hi   , & !
       hj   , & !
       hk   , & !
       hl   , & !
       tnext, & !
       tswap, & !
       zbuoy, & !
       ti   , & !
       si   , & !
       ri   , & !
       u    , & !
       v    , & !
       w    , & !
       du   , & !
       dv   , & !
       dw   , & !
       tx   , & !
       ty   , & !
       tz   , & !
       tf   , & !
       sf   , & !
       rf   , & !
       tfac , & !
       init_depth , & ! !! NG: 15_09_2008
       final_depth, & !! NG: 15_09_2008
       dtime

  REAL(kind = qprec), DIMENSION(:), ALLOCATABLE :: ttt

  REAL(kind=rprec), DIMENSION(:,:,:), ALLOCATABLE :: & !
       array_qual_traj
  REAL(kind=rprec), DIMENSION(:,:), ALLOCATABLE :: array_qual_traj_2D

  LOGICAL :: open_door=.TRUE., first_pass=.TRUE., first_call=.TRUE.
  LOGICAL :: Logic_flag
  LOGICAL :: horizontal_eddy, vert_NS_eddy, vert_EW_eddy

  !-------------------------------------------!
  !- DYNAMIC ALLOCATIONS AND INITIALIZATIONS -!
  !-------------------------------------------!
  CALL sub_posin_alloc(nmax)

  ALLOCATE(igo(nmax))   ; igo(:)   = iZero
  CALL sub_memory(size(igo),'i','igo','trajec_seq')
  IF (key_2D) THEN
    ALLOCATE(icolim(nmax))
    icolim(:)   = iZero
    CALL sub_memory(size(icolim),'i','icolim','trajec_seq')
  ENDIF
  ALLOCATE(i0(nmax))    ; i0(:)    = iZero
  CALL sub_memory(size(i0),'i','i0','trajec_seq')
  ALLOCATE(j0(nmax))    ; j0(:)    = iZero
  CALL sub_memory(size(j0),'i','j0','trajec_seq')
  ALLOCATE(k0(nmax))    ; k0(:)    = iZero
  CALL sub_memory(size(k0),'i','k0','trajec_seq')
  ALLOCATE(l0(nmax))    ; l0(:)    = iZero
  CALL sub_memory(size(l0),'i','l0','trajec_seq')
  ALLOCATE(i1(nmax))    ; i1(:)    = iZero
  CALL sub_memory(size(i1),'i','i1','trajec_seq')
  ALLOCATE(i2(nmax))    ; i2(:)    = iZero
  CALL sub_memory(size(i2),'i','i2','trajec_seq')
  ALLOCATE(j1(nmax))    ; j1(:)    = iZero
  CALL sub_memory(size(j1),'i','j1','trajec_seq')
  ALLOCATE(j2(nmax))    ; j2(:)    = iZero
  CALL sub_memory(size(j2),'i','j2','trajec_seq')
  ALLOCATE(k1(nmax))    ; k1(:)    = iZero
  CALL sub_memory(size(k1),'i','k1','trajec_seq')
  ALLOCATE(k2(nmax))    ; k2(:)    = iZero
  CALL sub_memory(size(k2),'i','k2','trajec_seq')
  ALLOCATE(iflag(nmax)) ; iflag(:) = iZero
  CALL sub_memory(size(iflag),'i','iflag','trajec_seq')
  ALLOCATE(iperfl(nmax)); iperfl(:)= iZero
  CALL sub_memory(size(iperfl),'i','iperfl','trajec_seq')
  ALLOCATE(jperfl(nmax)); jperfl(:)= iZero
  CALL sub_memory(size(jperfl),'i','jperfl','trajec_seq')
  ALLOCATE(nfin(nmax))  ; nfin(:)  = -iOne
  CALL sub_memory(size(nfin),'i','nfin','trajec_seq')

  ALLOCATE(age(nmax))   ; age(:)   = rZero
  CALL sub_memory(size(age),'r','age','trajec_seq')
  ALLOCATE(agenew(nmax)); agenew(:)= rZero
  CALL sub_memory(size(agenew),'r','agenew','trajec_seq')
  ALLOCATE(truage(nmax)); truage(:)= rZero
  CALL sub_memory(size(truage),'r','truage','trajec_seq')
  ALLOCATE(dtime(nmax)); dtime(:)= 0._rprec
  CALL sub_memory(size(dtime),'r','dtime','trajec_seq')

  tage(:) = rZero ! define is mod_posin.f90

  ALLOCATE(gi(nmax))   ; gi(:)   = rZero
  CALL sub_memory(size(gi),'r','gi','trajec_seq')
  ALLOCATE(gj(nmax))   ; gj(:)   = rZero
  CALL sub_memory(size(gj),'r','gj','trajec_seq')
  ALLOCATE(gk(nmax))   ; gk(:)   = rZero
  CALL sub_memory(size(gk),'r','gk','trajec_seq')
  ALLOCATE(gl(nmax))   ; gl(:)   = rZero
  CALL sub_memory(size(gl),'r','gl','trajec_seq')
  ALLOCATE(hi(nmax))   ; hi(:)   = rZero
  CALL sub_memory(size(hi),'r','hi','trajec_seq')
  ALLOCATE(hj(nmax))   ; hj(:)   = rZero
  CALL sub_memory(size(hj),'r','hj','trajec_seq')
  ALLOCATE(hk(nmax))   ; hk(:)   = rZero
  CALL sub_memory(size(hk),'r','hk','trajec_seq')
  ALLOCATE(hl(nmax))   ; hl(:)   = rZero
  CALL sub_memory(size(hl),'r','hl','trajec_seq')
  ALLOCATE(inext(nmax)); inext(:)= iZero
  CALL sub_memory(size(inext),'i','inext','trajec_seq')
  ALLOCATE(tnext(nmax)); tnext(:)= rZero
  CALL sub_memory(size(tnext),'r','tnext','trajec_seq')
  ALLOCATE(tswap(nmax)); tswap(:)= rZero
  CALL sub_memory(size(tswap),'r','tswap','trajec_seq')
  ALLOCATE(zbuoy(nmax)); zbuoy(:)= 1.e9_rprec
  CALL sub_memory(size(zbuoy),'r','zbuoy','trajec_seq')
  ALLOCATE(u(nmax))    ; u(:)    = rZero
  CALL sub_memory(size(u),'r','u','trajec_seq')
  ALLOCATE(v(nmax))    ; v(:)    = rZero
  CALL sub_memory(size(v),'r','v','trajec_seq')
  ALLOCATE(w(nmax))    ; w(:)    = rZero
  CALL sub_memory(size(w),'r','w','trajec_seq')
  ALLOCATE(du(nmax))   ; du(:)   = rZero
  CALL sub_memory(size(du),'r','du','trajec_seq')
  ALLOCATE(dv(nmax))   ; dv(:)   = rZero
  CALL sub_memory(size(dv),'r','dv','trajec_seq')
  ALLOCATE(dw(nmax))   ; dw(:)   = rZero
  CALL sub_memory(size(dw),'r','dw','trajec_seq')
  ALLOCATE(tx(nmax))   ; tx(:)   = rZero
  CALL sub_memory(size(tx),'r','tx','trajec_seq')
  ALLOCATE(ty(nmax))   ; ty(:)   = rZero
  CALL sub_memory(size(ty),'r','ty','trajec_seq')
  ALLOCATE(tz(nmax))   ; tz(:)   = rZero
  CALL sub_memory(size(tz),'r','tz','trajec_seq')
  ALLOCATE(tf(nmax))   ; tf(:)   = rZero
  CALL sub_memory(size(tf),'r','tf','trajec_seq')
  ALLOCATE(sf(nmax))   ; sf(:)   = rZero
  CALL sub_memory(size(sf),'r','sf','trajec_seq')
  ALLOCATE(rf(nmax))   ; rf(:)   = rZero
  CALL sub_memory(size(rf),'r','rf','trajec_seq')
  ALLOCATE(tfac(nmax)) ; tfac(:) = rZero
  CALL sub_memory(size(tfac),'r','tfac','trajec_seq')
  ALLOCATE(ttt(nmax))  ; ttt(:)  = 0._qprec
  CALL sub_memory(size(ttt),'r','ttt','trajec_seq')

  !! NG: 15_09_2008
  ALLOCATE(init_depth(nmax))  ; init_depth(:)  = rZero
  CALL sub_memory(size(init_depth),'r','init_depth','trajec_seq')
  ALLOCATE(final_depth(nmax)) ; final_depth(:) = rZero
  CALL sub_memory(size(final_depth),'r','final_depth','trajec_seq')

  IF (key_alltracers) THEN
    ALLOCATE(ti(nmax)); ti(:) = rZero
    CALL sub_memory(size(ti),'r','ti','trajec_seq')
    ALLOCATE(si(nmax)); si(:) = rZero
    CALL sub_memory(size(si),'r','si','trajec_seq')
    ALLOCATE(ri(nmax)); ri(:) = rZero
    CALL sub_memory(size(ri),'r','ri','trajec_seq')
  ENDIF

  !======================!
  !- Read region limits -!
  !======================!----------------------------!
  ! This is done in qualitative and quantitative mode !
  !---------------------------------------------------!
  CALL sub_reducmem_read_reg_lim()

  !==========================!
  ! Open ASCII output files -! (if key_ascii_outputs=True)
  !==========================!
  CALL sub_output_data_open_files()

  !-------------------!
  !- READ INPUT MESH -!
  !-------------------!---------------------------------------------------!
  !-- Allocate and read coordinates and scale factors 
  !-- xx_tt, xx_uu, xx_vv, xx_ff
  !-- yy_tt, yy_uu, yy_vv, yy_ff
  !-- zz_tt, zz_ww
  !-- e2u
  !-- e1v
  !-- e1t, e2t, e3t
  !-- tmask
  !-----------------------------------------------------------------------!
  CALL sub_input_grid()

  !==========================
  !- QUALITATIVE VARIABLES -!
  !=================================================================
  ! sampling time (in seconds) for the available transport field
  !=================================================================
  tfic = tunit   * REAL(ntfic, kind=rprec)

  !=================================================================
  ! sampling time (in seconds) between two successive output positions
  !=================================================================
  tlap = ABS(delta_t) * REAL(frequency, kind=rprec)

  !=================================================================
  ! maximum integration time (im seconds) for individual trajectories
  !=================================================================
  tfin = tlap    * REAL(ABS(nb_output), kind=rprec)

  !==============================================
  ! mask is written (for graphical trajectories)
  !==============================================
  IF ((TRIM(mode) == 'qualitative').AND.(mask)) THEN
    DO j = 1, jmt
      DO i = 1, imt

        !!NG: 20 may 2010 IF (tmask(i,j,1,1) == rZero) THEN
        IF (tmask(i,j,1,1) <= 0.5_rprec) THEN

          IF (key_alltracers) THEN
            WRITE(lun_traj,7055)0,xx_tt(i,j,1,1),yy_tt(i,j,1,1),0.,0.,0.,0.,0.
          ELSE
            WRITE(lun_traj,7055)0,xx_tt(i,j,1,1),yy_tt(i,j,1,1),0.,0.
          ENDIF

        ENDIF

      END DO
    END DO
  ENDIF

  !!NG: BEGIN SEQUENTIAL DATA ACCESS
  !=============================================!
  !- Allocate and Initialize input data arrays -!
  !=============================================!
  CALL sub_input_data_alloc_init_seq()
  !!NG: END SEQUENTIAL DATA ACCESS

  !=================================================================
  ! we initialize the *new* particle position 
  !=================================================================
  ! BBL __________________________________________________________________________
  ! BBL BBL: Initialisation a zero des prochaines positions a calculer sur la
  ! BBL        grille spatiale (I, J, K) et temporelle (L)
  ! BBL
  hi(:) = rZero
  hj(:) = rZero
  hk(:) = rZero
  hl(:) = rZero

  !=======================================!
  !- Read or Compute particule positions -!
  !=======================================!
  CALL sub_init_particules_positions()

  !=======================!
  !- NetCDF data outputs -!
  !=======================!
  !---------------------------!
  !- Open Netcdf output file -!
  !---------------------------!
  CALL sub_save_netcdf_data_init()

  !!NG  !------------------------------------------------!
  !!NG  !- Save Initial Positions in netcdf output file -!
  !!NG  !------------------------------------------------!
  !!NG  CALL sub_save_netcdf_init_pos()
  !!NG
  !!NG  !- Deallocate init_temp, init_salt, init_dens (mod_posin) -!
  !!NG  CALL sub_posin_dealloc_init_tracers()

  IF (TRIM(mode) == 'quantitative') THEN
    CALL sub_quant_alloc(imt, jmt, kmt)
    CALL sub_quant_init() !  initialize a few arrays

    !-------------------------------------------------------------------------
    ! we do not yet know whether a final criterion (criter1) will be satisfied
    !-------------------------------------------------------------------------
    icrit1 = 0

    secname(nsect+1) = 'Criter1'

  ELSEIF (TRIM(mode) == 'qualitative') THEN

    !! NG 29 april 2010 
    IF (small_virtual_memory) THEN

      ALLOCATE(array_qual_traj_2D(ntraj, nstat))
      CALL sub_memory(size(array_qual_traj_2D),'r','array_qual_traj_2D','trajec_seq')

      array_qual_traj_2D(:,:)=mask_value

    ELSE

      ALLOCATE(ind_traj(ntraj))
      CALL sub_memory(size(ind_traj),'i','ind_traj','trajec_seq')

      ALLOCATE(array_qual_traj(ntraj,nb_output+1,nstat))
      CALL sub_memory(size(array_qual_traj),'r','array_qual_traj','trajec_seq')

      array_qual_traj(:,:,:)=mask_value

    END IF

  ELSE
    STOP
  ENDIF

  !=============================================================!
  != First loop on the number of trajectories; initializations =!
  !=============================================================!
  DO n = 1, ntraj

    !====================================================================
    ! time index must be within [0.5, lmt+.5[; otherwise time translation
    !NG time index must be within ]0.5, lmt+.5]; otherwise time translation
    !====================================================================
    !NG !!! To work "lmt" must be greater than zero !!! (else infinit loop)
    !NG !!! "lmt" must be tested before to use it.  !!!
    !NG !!! If "lmt" is tested, GOTO could be removed. !!!
    ! BBL __________________________________________________________________________
    ! BBL BBL: Recadrage du parametre temporel FL dans la fenetre [0.5, lmt+0.5[
    ! BBL
100 IF (tfl(n) < 0.5_rprec) THEN
      tfl(n) = tfl(n) + REAL(lmt,kind=rprec)
      GOTO 100
    ENDIF ! (fl(n) < 0.5)

    IF (tfl(n) >= (REAL(lmt,kind=rprec) + 0.5_rprec)) THEN
      tfl(n) = tfl(n) - REAL(lmt,kind=rprec)
      GOTO 100
    ENDIF

    !====================================================================
    ! the program is able to deal with constant-depth particles
    !  (qualitative experiments only) [QUALI]
    !====================================================================
    ! BBL __________________________________________________________________________
    ! BBL BBL: La variable ZBUOY, initialisee par defaut a 1.e9, n'est active que
    ! BBL        dans le cas de calculs en mode QUALITATIVE, et seulement alors pour
    ! BBL        les particules dont on veut calculer des trajectoires a profondeur
    ! BBL        CONSTANTE
    ! BBL       Sa valeur est alors egale a l'immersion initiale de chaque particule
    ! BBL        concernee.
    ! BBL       NOTE: dans le cas d'un maillage variable dans le temps (evolution
    ! BBL        d'une surface libre par exemple), le calcul de l'immersion de
    ! BBL        reference n'est pas parfait car le maillage vertical utilise dans
    ! BBL        la fonction FZ est un maillage statique (ZETA = 0) [la dynamique du
    ! BBL        modele, incluant ZETA, n'a en effet pas encore ete lue a ce stade
    ! BBL        du code]
    ! BBL
    zbuoy(n)=1.e9_rprec

    IF (tfk(n) <= rZero) THEN

      IF (TRIM(mode) == 'quantitative') THEN 
        WRITE(lun_error,*)''
        WRITE(lun_error,*)'---    ERROR    ---'
        WRITE(lun_error,*)'Isobaric particles allowed in QUALITATIVE experiments ONLY '
        WRITE(lun_error,*)'Particule number n      = ', n
        WRITE(lun_error,*)'                 tfk(n) = ', tfk(n)
        WRITE(lun_error,*)'--- STOP ARIANE ---'
        STOP
      ENDIF

      IF (key_roms) THEN
        If (first_pass) THEN
          CALL sub_input_approx_zz_ww_roms()
          first_pass=.FALSE.
        ENDIF
        zbuoy(n)= fz(gi=tfi(n), gj=tfj(n), gk=tfk(n), il=1)
      ELSEIF (key_mars) THEN
        If (first_pass) THEN
          CALL sub_input_approx_zz_ww_mars()
          first_pass=.FALSE.
        ENDIF
        zbuoy(n)= fz(gi=tfi(n), gj=tfj(n), gk=tfk(n), il=1)
      ELSEIF (key_symphonie) THEN
        WRITE(lun_error,*)'       Isobaric particles'
        WRITE(lun_error,*)' NOT AVAILABLE FOR SYMPHONIE - STOP'
        STOP
      ELSE
        zbuoy(n)= fz(gi=tfi(n), gj=tfj(n), gk=tfk(n))
      ENDIF

      tfk(n)   = -tfk(n)

      !!NG: to debug___ write(*,*) zbuoy(n)

    ENDIF

    !====================================================================
    ! time (tfl(n)) and position (VELOCITY grid: fi, tfj(n), tfk(n))
    ! indices used to  detect the relevant embedding temperature
    !  gridcell (i, j, k, l).
    !====================================================================
    ! BBL __________________________________________________________________________
    ! BBL BBL: Definition de la maille TEMPERATURE (I0, J0, K0, L0) contenant la
    ! BBL        position initiale de chaque particule
    ! BBL
    i0(n) =  INT(tfi(n)) + 1
    j0(n) =  INT(tfj(n)) + 1
    k0(n) =  INT(tfk(n))
    l0(n) = NINT(tfl(n))


    !===========================!
    !- Periodicity in OPA-ORCA -!
    !===========================!
    IF (.NOT.(key_roms.OR.key_mars.OR.key_symphonie)) THEN

      CALL sub_orca_north_pole(i0(n), j0(n))

      CALL sub_orca_east_west_periodic(i0(n))

    ENDIF

    !===================================================================
    ! Compute diagnostics for initial tracers values
    !===================================================================
    ! BBL __________________________________________________________________________
    ! BBL BBL: La maille TEMPERATURE contenant la position initiale ne saurait
    ! BBL        etre une maille TERRE
    ! BBL
    CALL sub_reducmem_shift_or_not_ind(i0(n),j0(n),k0(n),is,js,ks)

    IF (tmask(is,js,ks,1) <= .5_rprec) THEN

      WRITE(lun_error,7000)'   false start', &
           n,i0(n),j0(n),k0(n),l0(n),tfi(n),tfj(n),tfk(n),tfl(n)

      STOP
    ENDIF ! (tmask(i0(n),j0(n),k0(n),1)

    !------------------------------------------------------------------!
    !- fi, fj, fk and fl refer to the initial position                -!
    !- gi, gj, gk and gl refer to the current position                -!
    !- BEWARE: we use a reverse vertical axis for easier computations -!
    !------------------------------------------------------------------!
    ! BBL __________________________________________________________________________
    ! BBL BBL: Positionnement initial de chaque particule a partir des valeurs
    ! BBL        FI, FJ, FK et FL lues sur fichier ou calculees par POSINI
    ! BBL       NOTE: les calculs faits dans le coeur du code utilisent des indices
    ! BBL        verticaux NEGATIFS, d'ou le changement de signe opere sur gk
    gi(n) =  tfi(n)
    gj(n) =  tfj(n)
    gk(n) = -tfk(n)
    gl(n) =  tfl(n)

    !----------------------!
    !- age initialization -!
    !----------------------!
    ! BBL __________________________________________________________________________
    ! BBL BBL: La variable AGE definit l'age (en secondes) de chaque particule
    ! BBL        depuis le debut de son integration temporelle
    ! BBL
    age(n) = tage(n)

    ! BBL __________________________________________________________________________
    ! BBL BBL: La variable IGO est utilisee pour referencer les etats de chaque
    ! BBL        particule, et en particulier les conditions suivantes:
    ! BBL         - pas encore partie
    ! BBL         - en cours d'integration
    ! BBL         - en attente d'une valeur actualisee du champ de transport
    !             - arrivee a destination
    igo(n) = iZero
    if (key_2D) icolim(n) = iZero

    !
    ! next output instant (for positions in qualitative experiments)
    !
    ! BBL __________________________________________________________________________
    ! BBL BBL: En mode QUALITATIVE la variable TNEXT definit le prochain instant
    ! BBL        (en secondes) ou la sortie d'une position sur le fichier TRAJ.QL
    ! BBL        sera necessaire
    ! BBL       NOTE: ce temps sera utilise par comparaison a l'AGE de chaque
    ! BBL        particule; dans la mesure ou on echantillonne les trajectoires de
    ! BBL        toutes les particules avec la meme frequence, la variable TNEXT est
    ! BBL        initialisee a la meme valeur pour toutes les particules
    ! BBL
    inext(n) = iOne
    tnext(n) = tlap

    !------------------------------------------------------------------
    !- next time swap instant (between two samples of velocity field) -
    !------------------------------------------------------------------
    ! BBL __________________________________________________________________________
    ! BBL BBL: La variable TSWAP definit le prochain instant (en secondes) ou une
    ! BBL        actualisation du champ de transport sera necessaire
    ! BBL       NOTE: ce temps sera utilise par comparaison avec l'AGE de chaque
    ! BBL        particule; dans la mesure ou les particules ne sont pas toutes
    ! BBL        initialisees au meme moment sur l'axe temporel, la valeur de TSWAP
    ! BBL        differe selon la particule consideree; son calcul depend aussi de
    ! BBL        l'orientation temporelle chosie (BACKWARD ou FORWARD)
    ! BBL
    IF (TRIM(forback) == 'forward') THEN
      tswap(n)=tfic*ABS( tfl(n)-(l0(n)+.5_rprec) ) + tage(n)
    ELSE
      tswap(n)=tfic*ABS( tfl(n)-(l0(n)-.5_rprec) ) + tage(n)
    ENDIF

    !!BBL - march 2007
    !!NG: option to follow an experience interrupted. Age is stored
    !!NG: in the final position file.
    !!NG: Dont USE this option in backward mode !!!
    IF (key_read_age) THEN
      dumm = tswap(n) / tfic
      tswap(n) = NINT(tswap(n)-NINT(dumm)*tfic)+NINT(dumm)*tfic
    ENDIF
    !
    ! end initializations
    !
  ENDDO ! n=1,ntraj

  !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
  !!===============================================================================!!
  !!_______________________________________________________________________________!!

  !-------------------------------------------!
  !- boucle sur les pas de temps disponibles -!
  !-------------------------------------------!
  ! BBL __________________________________________________________________________
  ! BBL BBL: Pour des calculs SEQUENTIELS en mode CLIMATOLOGIQUE, il est utile de
  ! BBL       definir le nombre maximal de cycles autorises sur l'axe temporel, de
  ! BBL       maniere a pouvoir filtrer les particules *trop lentes*
  ! BBL
  !! TODO: maxcyl should be a namelist parameter...
  !!NG: namelist => maxcycl0 = 300
  !!NG: namelist => maxcycl  = maxcycl0
  ! BBL __________________________________________________________________________
  ! BBL BBL: Dans le cas LMT=1, tout le champ de transport est charge en memoire,
  ! BBL       et la restriction precedente peut etre relachee
  ! BBL
  !!NG : IF (lmt == 1) maxcycl=12599999
  ! BBL __________________________________________________________________________
  ! BBL BBL: Boucle sur la duree definie par l'archivage utilise pour le champ de
  ! BBL      vitesse
  ! BBL


  !!========================================================================!!
  !!============================== DO LOOP: CYCLES =========================!!
  !!========================================================================!!
  !- Comments -!
  WRITE(lun_standard,*)''
  WRITE(lun_standard,*)'=================================='
  WRITE(lun_standard,*)'= ENTER >>>>> MAIN LOOP ON CYCLE ='
  WRITE(lun_standard,*)'=================================='

  !! TODO: infini do loop and exit.
  DO iter = 1, maxcycles

    if (iter==1) THEN
      first_call=.TRUE.
    ELSE
      first_call=.FALSE.
    ENDIF

    !- Initialize some variables need to the sequential mode -!
    iter_min=1
    iter_max=lmt
    IF (iter == 1) THEN
      IF ( (lmin /= 1).AND.(TRIM(forback) == 'forward')) THEN
        CALL sub_seq_init(lmin)
        iter_min=lmin
      ELSEIF (TRIM(forback) == 'backward') THEN
        CALL sub_seq_init(lmax)
        iter_max=lmax
      ELSE
        CALL sub_seq_init()
      ENDIF
    ELSE
      CALL sub_seq_init()
    ENDIF

    ! BBL __________________________________________________________________________
    ! BBL BBL: On utlise NDONE pour comptabiliser le nombre de particules ayant
    ! BBL       atteint leur etat final; le calcul n'est fait que pour informer
    ! BBL       l'utilisateur sur les progres effectifs realises par les calculs,
    ! BBL       cycle temporel apres cycle temporel
    ! BBL      Le calcul est fait a chaque iteration ITER si LMT > 1, et seulement
    ! BBL       au debut de la premiere iteration si LMT=1 (dans ce dernier cas en
    ! BBL       effet tout le champ de transport peut etre charge en memoire; les
    ! BBL       calculs sont donc beaucoup plus rapides, et il n'est pas necessaire
    ! BBL       d'informer l'utilisateur...)
    ! BBL
    IF ( ((iter == 1).AND.(lmt == 1)).OR. &
         ((iter /= 1) .AND.(lmt > 1)) ) THEN

      IF (TRIM(mode) == 'quantitative') THEN
        ndone = 0 ! initialize to 0 the number of trajectorie(s) done

        DO  n = 1, ntraj

          ! BBL ___________________________________________________________________
          ! BBL -- BBL: IGO=3 ou -3 signifie que l'integration de la particule est
          ! BBL          terminee
          ! BBL
          IF (ABS(igo(n)) >= 3) ndone = ndone + 1

        ENDDO

        WRITE(lun_standard,*)
        WRITE(lun_standard,*)'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
        WRITE(lun_standard,*)'HOW MANY PARTICLES REACH A SECTION: ', &
             INT(100._rprec*REAL(ndone,kind=rprec)/REAL(ntraj,kind=rprec)), '%'
        WRITE(lun_standard,*)'    - number of particules         : ', ntraj
        WRITE(lun_standard,*)'    - done                         : ', ndone
        WRITE(lun_standard,*)'    - rest                         : ', ntraj - ndone
        WRITE(lun_standard,*)'    - cycle                        : ', iter-1
        WRITE(lun_standard,*)'<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'

      ELSE

        WRITE(lun_standard,*)
        WRITE(lun_standard,*)'>>>>>>>>>>>>>>>>>>>>>>>>'
        WRITE(lun_standard,*)'>>  Cycle : ', iter-1
        WRITE(lun_standard,*)'<<<<<<<<<<<<<<<<<<<<<<<<'

      ENDIF

      id_comments = .FALSE.

    ELSE

      WRITE(lun_standard,*)
      WRITE(lun_standard,*)'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      WRITE(lun_standard,*)' PARTICULE TRAJECTORIES will be computed:'
      WRITE(lun_standard,*)'    - number of trajectories         : ', ntraj
      WRITE(lun_standard,*)'<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'

      !- To add comments when input data are reading -!
      !- comments only at the first reading -!
      id_comments = .TRUE.

    ENDIF ! ( ((iter == 1).AND.(lmt == 1)).OR.((iter /= 1) .AND.(lmt > 1)) )

    ! BBL __________________________________________________________________________
    ! BBL BBL: LREF designe l'indice actif (entre 1 et LMT) sur l'axe temporel,
    ! BBL        c'est-a-dire l'indice temporel du champ de transport a charger en
    ! BBL        memoire en debut de cycle
    ! BBL       LREF vaut donc 1 en mode FORWARD et LMT en mode BACKWARD
    ! BBL

    lref = iter_min

    IF (TRIM(forback) == 'backward') lref = iter_max

    !!    cnt = iter_min                                    ! counter
    !!    IF (TRIM(forback) == 'backward') cnt = iter_max ! Backward mode

    ! BBL __________________________________________________________________________
    ! BBL BBL: Le chargement effectif en memoire s'effectue obligatoirement a la
    ! BBL        premiere iteration (il faut faire au moins une lecture!), et a
    ! BBL        chaque modification de l'indice temporel si jamais LMT > 1
    ! BBL

9998 CONTINUE

    IF ( (iter == 1).OR.(lmt > 1) ) THEN

      !!--------------------------------!!
      !!- READ INPUT DATA SEQUENTIALLY -!!
      !!--------------------------------!!

      !- Comments -!
      WRITE(lun_standard,*)'---------------------------------'
      WRITE(lun_standard,*)'- READ INPUT DATA :', lref, '-'

      IF (key_roms) THEN
        CALL sub_input_data_seq_roms( &
             ncids(:)               , &
             varids(:)              , &
             new_file(:)            , &
             ind_file(:)            , &
             ind_time(:)            , &
             ind_time_size(:)       , &
             sdimsorders(:,:)        )

      ELSEIF (key_symphonie) THEN
        CALL sub_input_data_seq_symphonie( &
             ncids(:)                    , &
             varids(:)                   , &
             new_file(:)                 , &
             ind_file(:)                 , &
             ind_time(:)                 , &
             ind_time_size(:)            , &
             sdimsorders(:,:)            )

      ELSEIF (key_mars) THEN
        CALL sub_input_data_seq_mars( &
             ncids(:)               , &
             varids(:)              , &
             new_file(:)            , &
             ind_file(:)            , &
             ind_time(:)            , &
             ind_time_size(:)       , &
             sdimsorders(:,:)       , &
             lref)

      ELSEIF (key_B2C_grid) THEN

        CALL sub_input_data_B2C_gridz_seq( &
             ncids(:)                    , &
             varids(:)                   , &
             new_file(:)                 , &
             ind_file(:)                 , &
             ind_time(:)                 , &
             ind_time_size(:)            , &
             sdimsorders(:,:)            , &
             first_call                  , &
             lref                        )

        !!NG: 16/04/2009 WRITE(lun_error,*)''
        !!NG: 16/04/2009 WRITE(lun_error,*)'STOP === key_B2C_grid in sequential mode is not available === STOP'
        !!NG: 16/04/2009 STOP

      ELSE ! OPA-NEMO
        CALL sub_input_data_seq_opa( &
             ncids(:)              , &
             varids(:)             , &
             new_file(:)           , &
             ind_file(:)           , &
             ind_time(:)           , &
             ind_time_size(:)      , &
             sdimsorders(:,:)       )
      ENDIF

      IF (TRIM(forback) == 'backward') THEN
        !!        cnt = cnt - 1
        delta_t = - delta_t
        uu(:,:,:,:) = -uu(:,:,:,:)
        vv(:,:,:,:) = -vv(:,:,:,:)
        ww(:,:,:,:) = -ww(:,:,:,:)
        !!      ELSE
        !!        cnt = cnt + 1
      ENDIF

      IF (open_door) THEN
        !------------------------------------------------!
        !- Save Initial Positions in netcdf output file -!
        !------------------------------------------------!
        CALL sub_save_netcdf_init_pos()

        !- Deallocate init_temp, init_salt, init_dens (mod_posin) -!
        CALL sub_posin_dealloc_init_tracers()

        open_door=.FALSE.
      ENDIF


    ENDIF ! ((iter == 1).OR.(lmt > 1))

    !!----------------------------------------------------------------------------!!
    !!---1111111111111111111---- DO LOOP : ON TRAJECTORIES ----111111111111111----!!
    !!----------------------------------------------------------------------------!!
    ! marquage des particules
    ! code utilise:
    !  igo=-3: integration de la particule terminee
    !  igo=-2: attente de la disponibilite d'un nouveau champ de transport
    !  igo=0: pas encore partie
    !  igo=1: active
    !  igo=2: interception temporelle (swap) venant d'etre realisee
    !  igo=3: interception spatiale/criter1 venant d'etre realisee
    !

    DO n = 1, ntraj
      ! BBL __________________________________________________________________________
      ! BBL BBL: Si l'integration de la particule s'est terminee lors de
      ! BBL        l'utilisation du champ de transport precedent (IGO=3), on la place
      ! BBL        definitivement dans l'etat IGO=-3
      ! BBL
      IF (ABS(igo(n)) == 3) igo(n) = -3

      ! BBL __________________________________________________________________________
      ! BBL BBL: Si l'integration de la particule avait ete stoppee du fait d'une
      ! BBL        mise a jour necesaire du champ de transport (IGO=+/-2), on la remet
      ! BBL        en activite (IGO=1)
      ! BBL
      IF (ABS(igo(n)) == 2) igo(n) = 1

      IF ((igo(n) == 0).AND.(l0(n) == lref)) THEN

        ! BBL __________________________________________________________________________
        ! BBL BBL: Si la particule n'est pas encore partie (IGO=0) et si l'instant
        ! BBL        LREF qui vient d'etre charge en en memoire l'autorisee, elle est
        ! BBL        mise en activite (IGO=1)
        ! BBL
        igo(n) = 1
        !
        ! diagnostic for initial tracers values
        !
        ! BBL __________________________________________________________________________
        ! BBL - BBL: On calcule la valeur des traceurs sur la position initiale

        CALL sub_reducmem_shift_or_not_ind(i0(n),j0(n),k0(n),is,js,ks)

        IF (key_alltracers) THEN

          IF (key_nointerpolstats) THEN
            ti(n)=tt(is,js,ks,1)
            si(n)=ss(is,js,ks,1)
            ri(n)=rr(is,js,ks,1)
          ELSE
            ti(n)=zinter(tt,gi(n),gj(n),gk(n),1._rprec)
            si(n)=zinter(ss,gi(n),gj(n),gk(n),1._rprec)

            IF (key_approximatesigma) THEN
              ri(n)=zinter(rr,gi(n),gj(n),gk(n),1._rprec)

            ELSE
              ri(n)=sigma(zsigma,si(n),ti(n))

            ENDIF ! key_approximatesigma

          ENDIF ! key_nointerpolstats

        ENDIF ! key_alltracers

        !-------------------------------!
        !-QUANTITATIVE INITIALIZATIONS -!
        !-------------------------------!
        IF (TRIM(mode) == 'quantitative') THEN

          ! BBL __________________________________________________________________________
          ! BBL -- BBL: En mode QUANTITATIVE, on met a jour si necessaire la valeur des
          ! BBL          champs de transport 2D au voisinage de la section initiale
          ! BBL

          IF (key_roms.OR.key_mars.OR.key_symphonie) THEN
            prof=fz(gi=tfi(n), gj=tfj(n), gk=-tfk(n), il=1)
          ELSE
            prof=fz(gi=tfi(n), gj=tfj(n), gk=-tfk(n))
          ENDIF

          !! NG: 15_09_2008
          !! NG: Here we stored the depth of the particle at the intial position.
          !! NG: Because with ROMS or symphonie the zz_ww in mod_fz is different at
          !! NG: time step.

          init_depth(n)=prof


          CALL sub_reducmem_shift_or_not_ind( NINT(tfi(n)), NINT(tfj(n)), NINT(tfk(n)), &
               nis, njs, nks)

          IF ((ABS(ANINT(tfi(n))-tfi(n)) <= 1.e-5_rprec).AND.(uu(nis, js, ks,1) > rZero)) THEN
            IF (key_alltracers) THEN
              CALL sub_quant_initside_u(NINT(tfi(n)), j0(n), k0(n), 1, &
                   prof, ttr(n), ti(n), si(n), ri(n))
            ELSE
              CALL sub_quant_initside_u(NINT(tfi(n)), j0(n), k0(n), 1, &
                   prof, ttr(n))
            ENDIF
          ENDIF

          IF ((ABS(ANINT(tfj(n))-tfj(n)) <= 1.e-5_rprec).AND.(vv(is,njs,ks,1) > rZero)) THEN
            IF (key_alltracers) THEN
              CALL sub_quant_initside_v(i0(n), NINT(tfj(n)), k0(n), 1, &
                   prof, ttr(n), ti(n), si(n), ri(n))
            ELSE
              CALL sub_quant_initside_v(i0(n), NINT(tfj(n)), k0(n), 1, &
                   prof, ttr(n))
            ENDIF
          ENDIF

          IF ((ABS(ANINT(tfj(n))-tfj(n)) <= 1.e-5_rprec).AND.(ww(is,js,nks,1) < rZero)) THEN
            CALL sub_quant_initside_w(i0(n), j0(n), NINT(tfk(n)), 1, ttr(n))
          ENDIF

        ENDIF

        !---------------!
        !- QUALITATIVE -!
        !---------------!
        IF (TRIM(mode) == 'qualitative') THEN

          !------------------------------------------------!
          !- the initial position is written on "traj.ql" -!
          !------------------------------------------------!
          ! BBL ________________________________________________________________________
          ! BBL -- BBL: Calcul de la position geographique initiale
          ! BBL         La position verticale n'est pas calculee dans le cas de flotteurs
          ! BBL          isobares, mais elle est prise egale a la valeur de reference
          ! BBL          calculee a l'initialisation
          ! BBL
          x0 = fx(gi(n),gj(n))
          y0 = fy(gi(n),gj(n))

          IF (zbuoy(n) > 1.e8_rprec) THEN
            IF (key_roms.OR.key_mars.OR.key_symphonie) THEN
              z0= fz(gi=gi(n), gj=gj(n), gk=gk(n), il=1)
            ELSE
              z0=fz(gi(n),gj(n),gk(n))
            ENDIF
          ELSE
            z0=zbuoy(n)
          ENDIF

          !! NG 29 april 2010
          IF (small_virtual_memory) THEN

            indice_traj = 1_iprec

            array_qual_traj_2D(n,1) = x0
            array_qual_traj_2D(n,2) = y0
            array_qual_traj_2D(n,3) = z0
            array_qual_traj_2D(n,4) = rZero

            IF (key_alltracers) THEN

              array_qual_traj_2D(n,5) = ti(n)
              array_qual_traj_2D(n,6) = si(n)
              array_qual_traj_2D(n,7) = ri(n)

              IF (key_iU_jV_kW) THEN

                array_qual_traj_2D(n,8)  = gi(n)
                array_qual_traj_2D(n,9)  = gj(n)
                array_qual_traj_2D(n,10) = gk(n)
              END IF

            ELSE
              IF (key_iU_jV_kW) THEN

                array_qual_traj_2D(n,5) = gi(n)
                array_qual_traj_2D(n,6) = gj(n)
                array_qual_traj_2D(n,7) = gk(n)
              END IF

            ENDIF ! key_alltracers

          ELSE

            ind_traj(n) = 1

            array_qual_traj(n,ind_traj(n),1) = x0
            array_qual_traj(n,ind_traj(n),2) = y0
            array_qual_traj(n,ind_traj(n),3) = z0
            array_qual_traj(n,ind_traj(n),4) = rZero

            IF (key_alltracers) THEN

              array_qual_traj(n,ind_traj(n),5) = ti(n)
              array_qual_traj(n,ind_traj(n),6) = si(n)
              array_qual_traj(n,ind_traj(n),7) = ri(n)

              IF (key_iU_jV_kW) THEN

                array_qual_traj(n,ind_traj(n),8)  = gi(n)
                array_qual_traj(n,ind_traj(n),9)  = gj(n)
                array_qual_traj(n,ind_traj(n),10) = gk(n)

              END IF

            ELSE

              IF (key_iU_jV_kW) THEN

                array_qual_traj(n,ind_traj(n),5) = gi(n)
                array_qual_traj(n,ind_traj(n),6) = gj(n)
                array_qual_traj(n,ind_traj(n),7) = gk(n)

              END IF

            ENDIF ! key_alltracers

          END IF
          !! NG 29 april 2010

        ENDIF

      ENDIF

    ENDDO

    !! NG 29 april 2010
    !! Write output trajectories: The initial positions
    if (small_virtual_memory) THEN

      IF (key_alltracers.AND.key_iU_jV_kW) THEN

        CALL sub_save_netcdf_trajectories_generic(             &
             indice_traj                                     , &
             array_qual_traj_2D(:,1), array_qual_traj_2D(:,2), &
             array_qual_traj_2D(:,3), array_qual_traj_2D(:,4), &
             array_qual_traj_2D(:,5), array_qual_traj_2D(:,6), &
             array_qual_traj_2D(:,7), array_qual_traj_2D(:,8), &
             array_qual_traj_2D(:,9), array_qual_traj_2D(:,10) )

      ELSEIF (key_alltracers.AND.(.NOT.key_iU_jV_kW)) THEN

        CALL sub_save_netcdf_trajectories_generic(             &
             indice_traj                                     , &
             array_qual_traj_2D(:,1), array_qual_traj_2D(:,2), &
             array_qual_traj_2D(:,3), array_qual_traj_2D(:,4), &
             traj_temp=array_qual_traj_2D(:,5), &
             traj_salt=array_qual_traj_2D(:,6), &
             traj_dens=array_qual_traj_2D(:,7))

      ELSEIF ((.NOT.key_alltracers).AND.key_iU_jV_kW) THEN
        CALL sub_save_netcdf_trajectories_generic(             &
             indice_traj                                     , &
             array_qual_traj_2D(:,1), array_qual_traj_2D(:,2), &
             array_qual_traj_2D(:,3), array_qual_traj_2D(:,4), &
             traj_iU=array_qual_traj_2D(:,5), &
             traj_jV=array_qual_traj_2D(:,6), &
             traj_kW=array_qual_traj_2D(:,7))  
      ELSE
        CALL sub_save_netcdf_trajectories_generic(             &
             indice_traj                                     , &
             array_qual_traj_2D(:,1), array_qual_traj_2D(:,2), &
             array_qual_traj_2D(:,3), array_qual_traj_2D(:,4)  )
      ENDIF

      array_qual_traj_2D(:,:) = mask_value

      indice_traj = indice_traj + 1_iprec

    ELSE

    ENDIF
    !! NG 29 april 2010

     icount = 0
     DO n = 1, ntraj
       dtime(n) = age(n)
     ENDDO
2222 CONTINUE

     icount = icount + 1
     if (icount == 49) then
        DO n = 1, ntraj
           IF (igo(n) == 1) then
              IF (age(n) - dtime(n) > 10) then  ! we are moving just slowly
                 dtime(n) = age(n)
              ELSE
                 write (lun_standard, *) "Some Kind of Eddy"
                 write (lun_standard, *) gi(n), gj(n), gk(n)
                 ii = int(gi(n) + tol)
                 jj = int(gj(n) + tol)
                 kk = -int(gk(n) - tol)
                 !ii = nint(gi(n)) ! SEA 2024/05
                 !jj = nint(gj(n))
                 !kk = -nint(gk(n))
                 write (lun_standard, *) ii, jj, kk, ' ii, jj, kk'

                 ! now get their equivalents on the reduced grid RB 2022/07/20
                 CALL sub_reducmem_shift_or_not_ind(ii,jj,kk,is,js,ks)
                 write (lun_standard, *) is, js, ks, ' reduced grid'
                 write (lun_standard, *) tol, ' tol'

                 ! SEA 2024/05
                 write (lun_standard, *) 'Horizontal Eddy test', uu(is, js, ks, 1)*uu(is, js+1, ks, 1)
                 write (lun_standard, *) vv(is, js, ks, 1)*vv(is+1, js, ks, 1)
                 write (lun_standard, *) 'Vertical NS Eddy test', uu(is, js+1, ks-1, 1)*uu(is, js+1, ks, 1)
                 write (lun_standard, *) ww(is, js+1, ks, 1)*ww(is+1, js+1, ks, 1)
                 write (lun_standard, *) 'Vertical EW Eddy test', vv(is+1, js, ks-1, 1)*vv(is+1, js, ks, 1)
                 write (lun_standard, *) ww(is+1, js, ks, 1)*ww(is+1, js+1, ks, 1)
                 horizontal_eddy = (uu(is, js, ks, 1)*uu(is, js+1, ks, 1) < 0.) .AND. &
                      (vv(is, js, ks, 1)*vv(is+1, js, ks, 1) < 0.)
                 vert_NS_eddy = (uu(is, js+1, ks-1, 1)*uu(is, js+1, ks, 1) < 0.) .AND. &
                      (ww(is, js+1, ks, 1)*ww(is+1, js+1, ks, 1) < 0.)
                 vert_EW_eddy = (vv(is+1, js, ks-1, 1)*vv(is+1, js, ks, 1) < 0.) .AND. &
                      (ww(is+1, js, ks, 1)*ww(is+1, js+1, ks, 1) < 0.)

                 IF ( .not. horizontal_eddy .and. .not. vert_NS_eddy .and. .not. vert_EW_eddy) THEN
                    ! Eddy problem
                    CASE = 4
                 ELSEIF ( horizontal_eddy .NEQV. (vert_NS_eddy .NEQV. vert_EW_eddy) ) THEN
                    IF (horizontal_eddy) THEN
                       CASE = 1
                    ELSEIF (vert_NS_eddy) THEN
                       CASE = 2
                    ELSE
                       CASE = 3
                    ENDIF
                 ELSE
                    xs = abs(gi(n) - ii )
                    ys = abs(gj(n) - jj )
                    zs = abs(-gk(n) + kk)
                    IF ((xs > ys) .AND. (xs > zs)) THEN
                       CASE = 3
                    ELSEIF ((ys > xs) .AND. (ys > zs)) THEN
                       CASE = 2 
                    ELSE
                       CASE = 1
                    ENDIF
                 ENDIF

                 IF (CASE == 1) THEN
                    write (lun_standard, *) 'Horizontal Eddy', n
                    gi(n) = gi(n) + timestep * 0.5 * (uu(is, js, ks, 1) + uu(is, js+1, ks, 1)) / tfac(n)
                    gj(n) = gj(n) + timestep * 0.5 * (vv(is, js, ks, 1) + vv(is+1, js, ks, 1)) / tfac(n)
                    w(n) = ww(is, js, ks, 1)   - &
                         (gk(n) + REAL(k0(n), kind=rprec))   * (ww(is, js, ks+1, 1) - ww(is, js, ks, 1))
                    gk(n) = gk(n) - 10 * w(n) / tfac(n)
                    age(n) = age(n) + timestep
                    IF (gi(n) > int(gi(n))) then
                       i0(n) = int(gi(n)) + 1
                    ELSE
                       i0(n) = int(gi(n))
                    ENDIF
                    IF (gj(n) > int(gj(n))) then
                       j0(n) = int(gj(n)) + 1
                    ELSE
                       j0(n) = int(gj(n))
                    ENDIF
                 ELSEIF (CASE == 2) THEN
                    write (lun_standard, *) 'Vertical Eddy, NS axis', n
                    write (lun_standard, *) gi(n), gj(n), gk(n), ks
                    write (lun_standard, *) uu(is, js, ks, 1), uu(is, js, ks+1, 1), &
                         ww(is, js, ks, 1), ww(is+1, js, ks, 1)
                    gi(n) = gi(n) + 10 * 0.5 * (uu(is, js, ks, 1) + uu(is, js, ks+1, 1)) / tfac(n)
                    gk(n) = gk(n) - 10 * 0.5 * (ww(is, js, ks, 1) + ww(is+1, js, ks, 1)) / tfac(n)
                    v(n) = vv(is, js-1, ks, 1) + &
                         (gj(n) - REAL(j0(n)-1, kind=rprec)) * (vv(is, js, ks, 1) - vv(is, js-1, ks, 1))
                    gj(n) = gj(n) + 10 * v(n) / tfac(n)
                    age(n) = age(n) + 10
                    IF (gi(n) > int(gi(n))) then
                       i0(n) = int(gi(n)) + 1
                    ELSE
                       i0(n) = int(gi(n))
                    ENDIF
                    IF (gk(n) > int(gk(n))) then
                       k0(n) = - int(gk(n)) - 1
                    ELSE
                       k0(n) = - int(gk(n))
                    ENDIF
                 ELSEIF (CASE == 3) THEN
                    write (lun_standard, *) 'Vertical Eddy, EW axis', n
                    write (lun_standard, *) gi(n), gj(n), gk(n), kk
                    write (lun_standard, *) vv(is, js, ks, 1), vv(is, js, ks+1, 1), &
                         ww(is, js, ks, 1), ww(is, js+1, ks, 1)
                    gj(n) = gj(n) + 10 * 0.5 * (vv(is, js, ks, 1) + vv(is, js, ks+1, 1)) / tfac(n)
                    gk(n) = gk(n) - 10 * 0.5 * (ww(is, js, ks, 1) + ww(is, js+1, ks, 1)) / tfac(n)
                    u(n) = uu(is-1, js, ks, 1) + &
                         (gi(n) - REAL(i0(n)-1, kind=rprec)) * (uu(is, js, ks, 1) - uu(is, js-1, ks, 1))
                    gi(n) = gi(n) + 10 * u(n) / tfac(n)
                    age(n) = age(n) + 10
                    IF (gj(n) > int(gj(n))) then
                       j0(n) = int(gj(n)) + 1
                    ELSE
                       j0(n) = int(gj(n))
                    ENDIF
                    IF (gk(n) > int(gk(n))) then
                       k0(n) = - int(gk(n)) - 1
                    ELSE
                       k0(n) = - int(gk(n))
                    ENDIF
                 ELSEIF ( abs(gi(n) - nint(gi(n))) < tol .and. abs(gj(n) - nint(gj(n))) < tol .and. &
                      abs(gk(n) - nint(gk(n))) < tol ) then
                    write (lun_standard, *) 'Corner Problem', n
                    gi(n) = gi(n) + 10 * 0.25 * (uu(is, js, ks, 1) + uu(is, js+1, ks, 1) + &
                         uu(is, js, ks+1, 1) + uu(is, js+1, ks+1, 1)) / tfac(n)
                    gj(n) = gj(n) + 10 * 0.25 * (vv(is, js, ks, 1) + vv(is+1, js, ks, 1) + &
                         vv(is, js, ks+1, 1) + vv(is+1, js, ks+1, 1)) / tfac(n)
                    gk(n) = gk(n) - 10 * 0.25 * (ww(is, js, ks, 1) + ww(is, js+1, ks, 1) + &
                         ww(is+1, js, ks, 1) + ww(is+1, js+1, ks, 1)) / tfac(n)
                    age(n) = age(n) + 10
                    write (lun_standard, *) 'Move to ', gi(n), gj(n), gk(n)
                    IF (gi(n) > int(gi(n))) then
                       i0(n) = int(gi(n)) + 1
                    ELSE
                       i0(n) = int(gi(n))
                    ENDIF
                    IF (gj(n) > int(gj(n))) then
                       j0(n) = int(gj(n)) + 1
                    ELSE
                       j0(n) = int(gj(n))
                    ENDIF
                    IF (gk(n) > int(gk(n))) then
                       k0(n) = - int(gk(n)) - 1
                    ELSE
                       k0(n) = - int(gk(n))
                    ENDIF
                 ELSE
                    write (lun_standard, *) 'Eddy Problem', n, age(n), dtime(n)
                    write (lun_standard, *) is, js, ks, ' i, j, k (reduced grid)'

                    write (lun_standard, *) uu(is, js, ks, 1), uu(is, js+1, ks, 1), uu(is, js-1, ks, 1), ' uu, uu(j+1), uu(j-1)'
                    write (lun_standard, *) uu(is-1, js, ks, 1), uu(is-1, js+1, ks, 1), ' uu(i-1), uu(i-1,j+1)'
                    write (lun_standard, *) uu(is, js, ks, 1), uu(is, js, ks+1, 1), uu(is, js, ks-1, 1), ' uu, uu(k+1), uu(k-1)'
                    write (lun_standard, *)
                    write (lun_standard, *) vv(is, js, ks, 1), vv(is+1, js, ks, 1), vv(is-1, js, ks, 1), ' vv, vv(i+1), vv(i-1)'
                    write (lun_standard, *) vv(is, js, ks, 1), vv(is, js, ks+1, 1), vv(is, js, ks-1, 1), ' vv, vv(k+1), vv(k-1)'
                    write (lun_standard, *)
                    write (lun_standard, *) ww(is, js, ks, 1), ww(is+1, js, ks, 1), ww(is-1, js, ks, 1), ' ww, ww(i+1) ww(i-1)'
                    write (lun_standard, *) ww(is, js, ks, 1), ww(is, js+1, ks, 1), ww(is, js-1, ks, 1), ' ww, ww(j+1) ww(j-1)'
                    stop
                 ENDIF
              ENDIF
            ENDIF
        ENDDO
        icount = 1
     ENDIF
    !!----------------------------------------------------------------------------!!
    !!---2222222222222222222---- DO LOOP : ON TRAJECTORIES ----222222222222222----!!
    !!----------------------------------------------------------------------------!!

    
    DO  n = 1, ntraj

      IF (igo(n) == 1) THEN
        ! BBL __________________________________________________________________________
        ! BBL - BBL: On ne fait les calculs suivants que pour les particules *actives*
        ! BBL         (etat IGO=1)
        !
        ! linear interpolation for the local velocity
        !
        ! BBL __________________________________________________________________________
        ! BBL - BBL: Le champ de transport (UU, VV, WW) est interpole lineairement,
        ! BBL         composante par composante, sur la position courante (GI, GJ, GK)
        ! BBL
        CALL sub_reducmem_shift_or_not_ind(i0(n),j0(n),k0(n),is,js,ks)

        u(n) = uu(is-1,js,ks,1) + &
             (gi(n) - REAL(i0(n)-1,kind=rprec)) * (uu(is,js,ks,1)   - uu(is-1,js,ks,1))

        v(n) = vv(is,js-1,ks,1) + &
             (gj(n) - REAL(j0(n)-1,kind=rprec)) * (vv(is,js,ks,1)   - vv(is,js-1,ks,1))

        w(n) = ww(is,js,ks,1)   - &
             (gk(n) + REAL(k0(n),kind=rprec))   * (ww(is,js,ks+1,1) - ww(is,js,ks,1))
        !
        ! entry and exit indices (depending on the sign of the velocity)
        !
        ! BBL __________________________________________________________________________
        ! BBL - BBL: Pour chaque composante du transport, on calcule les indices des
        ! BBL         facettes *ENTREE* et *SORTIE*, a partir du signe de la composante
        ! BBL         du transport interpolee localement
        ! BBL
        i1(n) = i0(n) - 1
        i2(n) = i0(n)
        IF (u(n) < rZero) THEN
          i2(n) = i0(n) - 1
          i1(n) = i0(n)
        ENDIF

        j1(n) = j0(n) - 1
        j2(n) = j0(n)
        IF (v(n) < rZero) THEN
          j2(n) = j0(n) - 1
          j1(n) = j0(n)
        ENDIF

        k1(n) = k0(n) + 1
        k2(n) = k0(n)
        IF (w(n) < rZero) THEN
          k2(n) = k0(n) + 1
          k1(n) = k0(n)
        ENDIF

        !!NG      ENDIF
        !!NG    ENDDO

        !!----------------------------------------------------------------------------!!
        !!---3333333333333333333---- DO LOOP : ON TRAJECTORIES ---3333333333333333----!!
        !!----------------------------------------------------------------------------!!

        !!NG    DO  n=1,ntraj
        !!NG      IF (igo(n) == 1) THEN

        ! BBL __________________________________________________________________________
        ! BBL - BBL: On ne fait les calculs suivants que pour les particules *actives*
        ! BBL         (etat IGO=1)
        !
        ! times to get to each edge of the (i, j, k) grid cell:
        ! - infinite, if Uexit * Ulocal < 0
        ! - linear relationship, if Uexit = Ulocal
        ! - logarithmic formulation, in all other cases
        !  BEWARE: use of log(u)-log(v) instead of log (u/v),
        !   better for numerical results
        !
        ! BBL __________________________________________________________________________
        ! BBL - BBL: DU, DV et DW sont des tableaux intermediaires
        ! BBL
        CALL sub_reducmem_shift_or_not_ind(i0(n),j0(n),k0(n),is ,js ,ks )
        CALL sub_reducmem_shift_or_not_ind(i1(n),j1(n),k1(n),i1s,j1s,k1s)
        CALL sub_reducmem_shift_or_not_ind(i2(n),j2(n),k2(n),i2s,j2s,k2s)

        du(n) = ( uu(i2s,js ,ks ,1) - uu(i1s,js ,ks ,1) ) * REAL(i2(n)-i1(n),kind=rprec)
        dv(n) = ( vv(is ,j2s,ks ,1) - vv(is ,j1s,ks ,1) ) * REAL(j2(n)-j1(n),kind=rprec)
        dw(n) = ( ww(is ,js ,k2s,1) - ww(is ,js ,k1s,1) ) * REAL(k1(n)-k2(n),kind=rprec)

        ! BBL __________________________________________________________________________
        ! BBL - BBL: Les temps ARIANE (TX, TY, TZ, TTT) sont normalises, et le facteur
        ! BBL         de normalisation pour obtenir des temps *physiques* fait
        ! BBL         intervenir le volume de la maille TEMPERATURE contenant chaque
        ! BBL         particule
        ! BBL
        !!NG        tfac(n)=e3t(is, js, ks,1)*e1t(is, js,1,1)*e2t(is, js,1,1)

        IF (key_roms.OR.key_mars.OR.key_symphonie) THEN

          tfac(n) = e3t(is,js,ks,1) * e1t(is,js,1,1) * e2t(is,js,1,1)

        ELSE
          !- OPA case => E3T is not dependent of time -!
          IF (key_partialsteps) THEN
            tfac(n) = e3t(is,js,ks,1) * e1t(is,js,1,1) * e2t(is,js,1,1)
          ELSE
            tfac(n) = e3t(1,1,ks,1) * e1t(is,js,1,1) * e2t(is,js,1,1)
          ENDIF

        ENDIF

        !
        ! BBL __________________________________________________________________________
        ! BBL - BBL: Si le signe de la composante zonale du transport est different du
        ! BBL         signe de la composante zonale du transport sur la facette de
        ! BBL         sortie, le temps d'acces vers cette facette est infini
        ! BBL
        ! BBL __________________________________________________________________________
        ! BBL - BBL: Si la composante zonale du transport est la meme sur les facettes
        ! BBL         d'entree et de sortie, le temps d'acces est defini par une simple
        ! BBL         relation lineaire
        ! BBL
        ! BBL __________________________________________________________________________
        ! BBL - BBL: Temps d'acces (vers la facette de sortie) dans le cas general
        ! BBL

        IF (( u(n) * uu(i2s,js,ks,1)) <= rZero) THEN
          tx(n) = 1.e35_rprec
          !!NG: 20 may 2010 ELSEIF (du(n) == rZero) THEN
        ELSEIF (abs(du(n)/u(n)) <= 1.e-11_rprec) THEN
          tx(n) = ( REAL(i2(n),  kind=rprec) - gi(n) ) / u(n)
        ELSE
          tx(n) = ( LOG(ABS(uu(i2s,js,ks,1))) - LOG(ABS(u(n))) ) / du(n)
        ENDIF

        !
        ! BBL __________________________________________________________________________
        ! BBL - BBL: Les calculs sont equivalents pour la composante meridienne du
        ! BBL         transport
        ! BBL

        IF ((v(n) * vv(is,j2s,ks,1)) <= rZero) THEN
          ty(n) = 1.e35_rprec
          !!NG: 20 may 2010 ELSEIF (dv(n) == rZero) THEN
        ELSEIF (abs(dv(n)/v(n)) <= 1.e-11_rprec) THEN
          ty(n) = ( REAL(j2(n), kind=rprec) - gj(n) ) / v(n)
        ELSE
          ty(n) = ( LOG(ABS(vv(is,j2s,ks,1))) - LOG(ABS(v(n))) ) / dv(n)
        ENDIF

        ! BBL __________________________________________________________________________
        ! BBL - BBL: Les calculs sont equivalents pour la composante verticale du
        ! BBL         transport, avec un traitement specifique des flotteurs isobares
        ! BBL         (pour lesquels les temps d'acces vers les facettes de sortie sont
        ! BBL         forces a une valeur *infinie*)
        ! BBL
        !!NG: pour BBL
        IF ((key_2dquant).OR.(zbuoy(n) <= 1.e8_rprec).OR.((w(n)*ww(is,js,k2s,1)) <= rZero)) THEN
          tz(n) = 1.e35_rprec
          !!NG: 20 may 2010 ELSEIF (dw(n) == rZero) THEN
        ELSEIF (abs(dw(n)/w(n)) <= 1.e-11_rprec) THEN
          tz(n) = -( REAL(k2(n), kind=rprec) + gk(n) ) / w(n)
        ELSE
          tz(n) = ( LOG(ABS(ww(is,js,k2s,1))) - LOG(ABS(w(n))) ) / dw(n)
        ENDIF

        !-----------------------------------------------------------------------------
        ! the true exit time is given by the shortest times (of tx, ty, tz)
        ! entry and exit indices are then updated
        !-----------------------------------------------------------------------------

        ! RB printing the exit time of the particle i'm having trouble with at each timestep - 2024/03
        IF ((n) == 1053126) THEN
          write (lun_standard, *) gi(n), gj(n), gk(n), 'gi, gj, gk - problem particle position'
          write (lun_standard, *) is, js, ks, 'is, js, ks'
          write (lun_standard, *) i1s, j1s, k1s, 'i1s, j1s, k1s'
          write (lun_standard, *) i2s, j2s, k2s, 'i2s, j2s, k2s'
          !ii = int(gi(n) + tol)
          !jj = int(gj(n) + tol)
          !kk = -int(gk(n) + tol)
          !write (lun_standard, *) ii, jj, kk, ' ii, jj, kk'
          write (lun_standard, *) tx(n), ty(n), tz(n), 'tx, ty, tz'
          write (lun_standard, *) uu(is, js, ks, 1), ww(is, js, ks, 1), ' uu(i,j,k), ww(i,j,k)'
          write (lun_standard, *) uu(is-1, js, ks, 1), ww(is-1, js, ks, 1), ' uu(i-1,j,k), ww(i-1,j,k)'
          write (lun_standard, *) uu(is+1, js, ks, 1), ww(is+1, js, ks, 1), ' uu(i+1,j,k), ww(i+1,j,k)'
          write (lun_standard, *) uu(is, js, ks-1, 1), ww(is, js, ks-1, 1), ' uu(i,j,k-1), ww(i,j,k-1)'
          write (lun_standard, *) uu(is, js, ks+1, 1), ww(is, js, ks+1, 1), ' uu(i,j,k+1), ww(i,j,k+1)'

        ENDIF

        ! BBL __________________________________________________________________________
        ! BBL - BBL: Si une particule sort de la maille temperature qui la contient,
        ! BBL         elle le fait en un temps TTT, egal au minimum des 3 temps d'acces
        ! BBL         precedents (TX, TY, TZ)
        ! BBL
        ttt(n)=MIN(tx(n),ty(n),tz(n))

        IF (n == 40) THEN !BNB 
          write (lun_standard, *) tx(n), ty(n), tz(n), ttt(n), ' tx(n), ty(n), tz(n), ttt(n)'
          write (lun_standard, *) v(n), dv(n), vv(is ,j2s,ks ,1), vv(is ,j1s,ks ,1), ' v(n), dv(n), vv(j2s), vv(j1s)'
          write (lun_standard, *) w(n), dw(n), ww(is ,js ,k2s,1), ww(is ,js ,k1s,1), ' w(n), dw(n), ww(k2s), ww(k1s)'
        ENDIF


        !!NG: pour BBL 
        IF (.NOT.key_2dquant) THEN
          !
          ! BBL __________________________________________________________________________
          ! BBL - BBL: Dans certains cas, les 3 temps d'acces (TX, TY, TZ) sont infinis...
          ! BBL
          IF (ttt(n) >= 1.e34_qprec) THEN
            ! BBL __________________________________________________________________________
            ! BBL - BBL: ... cela peut notamment se produire pour des particules isobares
            ! BBL         (pour lesquelles on a force TZ a etre infini); dans ce cas, c'est
            ! BBL         *normal* et on continue les calculs
            ! BBL
            !! NG: 12/05/2009 IF (zbuoy(n) <= 1.e8_rprec) GOTO 204
            ! BBL __________________________________________________________________________
            ! BBL - BBL: ... si la particule n'est pas isobare, c'est qu'elle se trouve
            ! BBL         piegee dans une maille ou les 6 composantes du transport sont
            ! BBL         *entrantes*; cela peut se produire dans le cas d'un maillage
            ! BBL         vertical qui *respire* (avec la montee et la descente de la
            ! BBL         surface libre);
            ! BBL        La valeur de TTT est ramenee a une valeur *finie*, choisie de telle
            ! BBL         facon qu'elle excede la duree documentee par le champ de
            ! BBL         transport present en memoire
            ! BBL        NOTE: cette precaution est en fait inutile, et on pourrait tout
            ! BBL         aussi bien continuer avec la valeur *infinie* de TTT, comme pour
            ! BBL         les particules isobares; a la place de la ligne de code suivante,
            ! BBL         il faudrait en fait introduire un test du type:
            ! BBL          - si le maillage vertical peut evoluer dans le temps, une valeur
            ! BBL           de TTT *infinie*  n'est pas incongrue, et on continue les
            ! BBL           calculs sans modifier TTT
            ! BBL          - si le maillage vertical est statique, TTT *infini* est une
            ! BBL           aberration, et il faut stopper les calculs avec un message
            ! BBL           d'erreur 'negative time'
            ! BBL
            !! NG: 12/05/2009 ttt(n)=2._rprec*(tswap(n)-age(n))/tfac(n)
            IF (zbuoy(n) > 1.e8_rprec)  THEN
              ttt(n)=2._qprec*(tswap(n)-age(n))/tfac(n)
            ENDIF

          ENDIF

        ENDIF
        !
        ! BBL __________________________________________________________________________
        ! BBL - BBL: Dans de tres rares cas, TTT prend une valeur negative!
        ! BBL
        IF (ttt(n) < 0._qprec) THEN
          ! BBL __________________________________________________________________________
          ! BBL - BBL: Si TTT est *significativement* negatif (> 1 seconde), il y a un
          ! BBL         probleme et on stoppe l'integration avec un message d'erreur
          ! BBL
          IF (ABS(tfac(n)*ttt(n)) > 1._rprec) THEN
            WRITE(lun_error,7000)' negative time', &
                 n,i0(n),j0(n),k0(n),l0(n),tfi(n),tfj(n),tfk(n),tfl(n)
            STOP
          ENDIF
          ! BBL __________________________________________________________________________
          ! BBL - BBL: Sinon, il s'agit d'un probleme d'arrondi machine, et on force TTT
          ! BBL         a etre egal a 0
          ! BBL
          ttt(n)=0._qprec
        ENDIF

        IF (key_2D) THEN
          ! BBL_DEBUG_BBL START ontologie phenomenologique de la particule colimacon
          !   avec (pour SLAM) neutralisation additionelle du test sur zbuoy (car on tourne avec un "W" actif fictif)
          !         IF (zbuoy(n) <= 1.e8) THEN
          IF ( ((abs(gi(n)-int(gi(n))).le.1e-3).and.((tfac(n)*ty(n)).le.864)).or. &
               ((abs(gj(n)-int(gj(n))).le.1e-3).and.((tfac(n)*tx(n)).le.864)) ) THEN
            icolim(n)=icolim(n)+1
            IF (icolim(n).gt.99) THEN
              ttt(n)=(tswap(n)-age(n))/tfac(n)
              igo(n)= iTwo
              WRITE(lun_error,7000)'     colimacon', &
                   n,i0(n),j0(n),k0(n),l0(n),tfi(n),tfj(n),tfk(n),tfl(n)
              icolim(n)= iZero
            ENDIF
          ELSE
            icolim(n)= iZero
          ENDIF
          !         ENDIF
          ! BBL_DEBUG_BBL END ontologie phenomenologique de la particule colimacon
        ENDIF

        ! BBL __________________________________________________________________________
        ! BBL - BBL: On ne fait les calculs suivants que pour les particules *actives*
        ! BBL         (etat IGO=1)
        !
        ! do we need to update the velocity field WITHIN the gridcell ?
        !
        ! update needed
        !
        IF (key_2D) THEN
          logic_flag = ((age(n)+ttt(n)*tfac(n)) > tswap(n))
        ELSE
          logic_flag = ((lmt > 1).AND.((age(n)+ttt(n)*tfac(n)) > tswap(n)))
        ENDIF
        IF (logic_flag) THEN
          ! BBL __________________________________________________________________________
          ! BBL - BBL: Si LMT > 1, le temps d'acces vers la facette de sortie NE SAURAIT
          ! BBL         exceder la duree documentee par le champ de transport present en
          ! BBL         memoire
          ! BBL        Si tel est le cas, on seuille TTT au temps que le champ present en
          ! BBL         memoire permet encore de documenter, et la particule passe dans
          ! BBL         l'etat IGO=2: une mise a jour du champ de transport s'avere
          ! BBL         necessaire avant de la faire *repartir*
          ! BBL        NOTE: pour comparer les temps ARIANE aux temps physiques (AGE,
          ! BBL         TSWAP, TFIN, ...), l'utilisation du facteur de normalisation est
          ! BBL         necessaire
          ! BBL
          ttt(n)=(tswap(n)-age(n))/tfac(n)

          !! NG: 20 may 2010
          !! NG: tmp = age(n)+ttt(n)*tfac(n)

          !! NG: 20 may 2010
          !! NG: IF (tmp /= tswap(n)) THEN
          !! NG:   write(lun_error,*)'trajec_seq.f90: 1 - precision could be a problem !!!'
          !! NG:    write(lun_error,*)'          n:', n
          !! NG:    write(lun_error,'(a,F30.15)')'   tswap(n):', tswap(n)
          !! NG:    write(lun_error,'(a,F30.15)')'        tmp:', tmp
          !! NG:    write(lun_error,'(a,F30.15)')'     ttt(n):',ttt(n)
          !! NG:    write(lun_error,'(a,F30.15)')'     age(n):',age(n)
          !! NG:    write(lun_error,'(a,F30.15)')'    tfac(n):',tfac(n)
          !! NG: END IF

          igo(n) = 2
        ENDIF
        !
        ! end of the integration
        !
        IF ((TRIM(mode) == 'qualitative').AND.((age(n)+ttt(n)*tfac(n)) > tfin)) THEN
          ! BBL __________________________________________________________________________
          ! BBL - BBL: Dans un cas QUALITATIVE, la valeur retenue pour TTT ne doit pas
          ! BBL         conduire a un age de la particule superieur au temps maximal
          ! BBL         d'integration demande
          ! BBL        Si tel est le cas, on seuille TTT au temps maximal admissible, et
          ! BBL         la particule passe dans l'etat IGO=3 (integration terminee)
          ! BBL

          ttt(n)=(tfin-age(n))/tfac(n)

          !!NG : 14 march 2011 tmp = age(n) + ttt(n) * tfac(n)

          !! NG : 20 may 2010
          !! NG: IF (tmp /= tfin) THEN
          !! NG:    write(lun_error,*)'trajec_seq.f90: 2 - precision could be a problem !!!'
          !! NG:    write(lun_error,*)'      n:', n
          !! NG:    write(lun_error,'(a,F30.15)')'      tfin:', tfin
          !! NG:    write(lun_error,'(a,F30.15)')'       tmp:', tmp
          !! NG:    write(lun_error,'(a,F30.15)')'    ttt(n):', ttt(n)
          !! NG:    write(lun_error,'(a,F30.15)')'    age(n):', age(n)
          !! NG:    write(lun_error,'(a,F30.15)')'   tfac(n):', tfac(n)
          !! NG: END IF
          igo(n) = 3
        ENDIF

        !====================================================================
        ! case 1: we do not reach the side of the cell (for a given direction)
        ! diagnostic of the final position
        ! - linear formulation, if Uexit = Ulocal
        ! - exponential formulation, in other cases
        !====================================================================
        IF (tx(n) > ttt(n)) THEN

          ! BBL __________________________________________________________________________
          ! BBL - BBL: On traite le cas ou le temps de sortie dans la direction zonale est
          ! BBL          plus grand que le temps TTT retenu pour les calculs
          ! BBL        On calcule donc une nouvelle position dans la maille TEMPERATURE,
          ! BBL         a partir de la valeur TTT (avec l'expression analytique de la
          ! BBL         trajectoire dans la direction zonale)
          ! BBL

          CALL sub_dont_reachside(hi(n), du(n), gi(n), u(n), ttt(n), i1(n), i2(n))

        ENDIF

        IF (ty(n) > ttt(n)) THEN

          ! BBL __________________________________________________________________________
          ! BBL - BBL: On applique un traitement equivalent dans la direction meridienne
          ! BBL
          CALL sub_dont_reachside(hj(n), dv(n), gj(n), v(n), ttt(n), j1(n), j2(n))

        ENDIF


        IF (TRIM(mode) == 'quantitative') THEN
          !!NG: BBL
          IF ((zbuoy(n) > 1.e8_rprec).AND.(.NOT.key_2dquant)) THEN

            IF (tz(n) > ttt(n)) THEN
              ! BBL __________________________________________________________________________
              ! BBL -- BBL: On applique un traitement equivalent dans la direction verticale,
              ! BBL          si la particule n'est pas isobare
              ! BBL
              CALL sub_dont_reachside(hk(n), dw(n), gk(n), w(n), ttt(n), -k1(n), -k2(n))

            ENDIF
          ELSEIF(key_2dquant) THEN

            ! BBL __________________________________________________________________________
            ! BBL - BBL: Si la particule est isobare, on ne fait aucun calcul et on force la
            ! BBL         conservation de l'indice de maillage vertical
            ! BBL
            hk(n)=gk(n)

          ENDIF

        ELSE !! QUALITATIVE !!
          IF (zbuoy(n) > 1.e8_rprec) THEN

            IF (tz(n) > ttt(n)) THEN
              ! BBL __________________________________________________________________________
              ! BBL -- BBL: On applique un traitement equivalent dans la direction verticale,
              ! BBL          si la particule n'est pas isobare
              ! BBL
              CALL sub_dont_reachside(hk(n), dw(n), gk(n), w(n), ttt(n), -k1(n), -k2(n))

            ENDIF
          ELSE

            ! BBL __________________________________________________________________________
            ! BBL - BBL: Si la particule est isobare, on ne fait aucun calcul et on force la
            ! BBL         conservation de l'indice de maillage vertical
            ! BBL
            hk(n)=gk(n)

          ENDIF
        ENDIF

        !====================================================================
        ! case 2: we DO reach the side of the cell (for a given direction)
        ! no need for time update within the gridcell
        ! - we record the exiting transport
        ! - we update the positions index (hi, hj, hk)
        !====================================================================
        IF (key_2D) THEN
          logic_flag = ((tx(n) <= ttt(n)).AND.(tx(n) <= ty(n)))
        ELSE
          logic_flag = (tx(n) <= ttt(n))
        ENDIF
        IF (logic_flag) THEN

          ! BBL __________________________________________________________________________
          ! BBL -- BBL: On traite le cas ou le temps d'acces vers la section de sortie
          ! BBL         *zonale* (TX) est inferieur ou egal au temps TTT retenu pour
          ! BBL         l'integration de la particules
          ! BBL         NOTE: TX ne peut etre en fait STRICTEMENT inferieur a TTT (par
          ! BBL          construction), et le test capture en fait les cas TX=TTT
          ! BBL __________________________________________________________________________
          ! BBL -- BBL: Le nouvel indice de positionnement zonal est EXACTEMENT l'indice
          ! BBL          de la facette de sortie
          ! BBL
          hi(n)=REAL(i2(n),kind=rprec)

          IF (TRIM(mode) == 'quantitative') THEN

            ! BBL __________________________________________________________________________
            ! BBL BBL BBL: Dans le cas QUANTITATIVE, on memorise le passage de la particule
            ! BBL           dans les champs de transport lagrangien 2D appropries
            ! BBL

            CALL sub_quant_reachside_u( &
                 i2(n), j0(n), k0(n), 1, &
                 hi(n), hj(n), hk(n), 1._rprec, &
                 ttr(n))

          ENDIF

          !
          ! BBL __________________________________________________________________________
          ! BBL -- BBL: On effectue la mise a jour des indices des facettes *ENTREE* et
          ! BBL          *SORTIE* dans la direction zonale; cette mise a jour depend bien
          ! BBL           sur du sens du deplacement
          ! BBL
          IF (i2(n) > i1(n)) THEN
            i1(n) = i2(n)
            i2(n) = i2(n)+1
          ELSE
            i1(n) = i2(n)
            i2(n) = i2(n)-1
          ENDIF

        ENDIF

        IF (key_2D) THEN
          logic_flag = ((ty(n) <= ttt(n)).AND.(ty(n) <= tx(n)))
        ELSE
          logic_flag = (ty(n) <= ttt(n))
        ENDIF
        IF (logic_flag) THEN

          ! BBL __________________________________________________________________________
          ! BBL -- BBL: On applique un traitement equivalent dans la direction meridienne
          ! BBL
          hj(n) = REAL(j2(n),kind=rprec)

          IF (TRIM(mode) == 'quantitative') THEN

            CALL sub_quant_reachside_v(         &
                 i0(n), j2(n), k0(n), 1,        &
                 hi(n), hj(n), hk(n), 1._rprec, &
                 ttr(n))

          ENDIF

          IF (j2(n) > j1(n)) THEN
            j1(n) = j2(n)
            j2(n) = j2(n) + 1
          ELSE
            j1(n) = j2(n)
            j2(n) = j2(n) - 1
          ENDIF

        ENDIF


        IF (tz(n) <= ttt(n)) THEN

          ! BBL __________________________________________________________________________
          ! BBL -- BBL: On applique un traitement equivalent dans la direction verticale,
          ! BBL          sauf pour les particules isobares pour lesquelles aucun
          ! BBL          franchissement vertical de maille n'est necessaire
          ! BBL
          IF (TRIM(mode) == 'quantitative') THEN

            CALL sub_quant_reachside_w(i0(n), j0(n), k2(n), 1, ttr(n))

          ENDIF

          !!NG: BBL
          IF (.NOT.key_2dquant) THEN

            IF (zbuoy(n) > 1.e8_rprec) THEN

              hk(n) = -REAL(k2(n),kind=rprec)

              IF (k2(n) > k1(n)) THEN
                k1(n) = k2(n)
                k2(n) = k2(n) + 1
              ELSE
                k1(n) = k2(n)
                k2(n) = k2(n) - 1
              ENDIF

            ENDIF

          ENDIF

        ENDIF


        ! BBL __________________________________________________________________________
        ! BBL - BBL: Le cas des particules isobares necessite un traitement particulier
        ! BBL        En effet, pour un systeme de coordonnees suivant la topographie, si
        ! BBL         l'on se deplace dans une couche donnee sans modifier l'indice
        ! BBL         vertical (HK), on ne conserve pas forcement l'immersion de
        ! BBL         reference
        ! BBL        Ainsi, il convient de verifier que la maille TEMPERATURE dans
        ! BBL         laquelle la particule evolue couvre toujours un domaine de
        ! BBL         profondeur incluant l'immersion de reference a laquelle on
        ! BBL         souhaite faire voyager la particule isobare
        ! BBL
        IF (zbuoy(n) < 1.e8_rprec) THEN
          iflag(n)=0
          ! BBL __________________________________________________________________________
          ! BBL -- BBL: Dans le contexte le plus simple, l'extension verticale de la
          ! BBL          maille ambiante inclut bien l'immersion de reference, et il n'est
          ! BBL          pas necessaire de modifier l'indice de positionnement vertical
          ! BBL 

          !NG: potential bug : change <0 by <=rZero (15/05/2009)

          IF (                                                                      &
               (zbuoy(n) - fz( gi=hi(n), gj=hj(n), gk=-AINT(-hk(n)), il=1))       * &
               (zbuoy(n) - fz( gi=hi(n), gj=hj(n), gk=-1._rprec-AINT(-hk(n)), il=1))&
               <= rZero) THEN
            hk(n)=gk(n)
            iflag(n)=1
          ELSE
            ! BBL __________________________________________________________________________
            ! BBL -- BBL: Dans le cas contraire, il convient de scanner la verticale pour
            ! BBL          trouver la nouvelle maille TEMPERATURE dont l'extension verticale
            ! BBL          inclut l'immersion de reference de la particule isobare
            ! BBL

            DO kk0=1,kmt-1

              !NG: potential bug : change <0 by <=rZero (15/05/2009)

              IF (                                                                           &
                   (zbuoy(n) - fz(gi=hi(n),gj=hj(n),gk=-REAL(kk0, kind=rprec),il=1))       * &
                   (zbuoy(n) - fz(gi=hi(n),gj=hj(n),gk=-1._rprec-REAL(kk0,kind=rprec),il=1)) &
                   <= rZero) THEN
                hk(n)=-REAL(kk0, kind=rprec)-.5_rprec
                k1(n)=kk0
                k2(n)=kk0+1
                iflag(n)=1
              ENDIF
            END DO

          ENDIF

          IF (iflag(n) == 0) THEN
            ! BBL __________________________________________________________________________
            ! BBL -- BBL: La particule qui ne peut plus etre positionnee dans une maille
            ! BBL          TEMPERATURE appropriee (incluant son immersion de reference) voit
            ! BBL           son integration se finir, et son etat passe a IGO=3
            ! BBL
            WRITE(lun_error,7000)'out of zdomain', &
                 n,i0(n),j0(n),k0(n),l0(n),tfi(n),tfj(n),tfk(n),tfl(n)
            igo(n) = 3
          ENDIF
        ENDIF
        !
      ENDIF

    ENDDO

    !!----------------------------------------------------------------------------!!
    !!---5555555555555555555---- DO LOOP : ON TRAJECTORIES ---5555555555555555----!!
    !!----------------------------------------------------------------------------!!

    DO n=1,ntraj
      IF (igo(n) > 0) THEN
        ! BBL __________________________________________________________________________
        ! BBL - BBL: On ne fait les calculs suivants que pour les particules *actives*
        ! BBL         (etat IGO=1,2,3)
        !
        ! update of the T-grid index of the new embedding gridcell
        !
        ! BBL __________________________________________________________________________
        ! BBL - BBL: Il convient de mettre a jour les indices (I0, J0, K0) de la
        ! BBL         nouvelle maille TEMPERATURE contenant la particule, a partir de
        ! BBL         la connaissance des nouveaux indices des facettes d'entree et de
        ! BBL         sortie
        ! BBL
        i0(n)    = max0(i1(n), i2(n))
        j0(n)    = max0(j1(n), j2(n))
        k0(n)    = min0(k1(n), k2(n))

        !===========================!
        !- Periodicity in OPA-ORCA -!
        !===========================!
        IF (.NOT.(key_roms.OR.key_mars.OR.key_symphonie)) THEN

          !!NG: 8 march 2013 CALL sub_orca_north_pole(i1(n), j1(n), j0(n)+1)
          !!NG: 8 march 2013 CALL sub_orca_north_pole(i2(n), j2(n), j0(n)+1)
          !!NG: 8 march 2013 CALL sub_orca_north_pole(hi(n), hj(n), j0(n)+1)
          !!NG: 8 march 2013 CALL sub_orca_north_pole(i0(n), j0(n), j0(n)+1, jperfl(n))

          CALL sub_orca_north_pole(i1(n), j1(n), j0(n))
          CALL sub_orca_north_pole(i2(n), j2(n), j0(n))
          CALL sub_orca_north_pole(hi(n), hj(n), j0(n))
          CALL sub_orca_north_pole(i0(n), j0(n), j0(n), jperfl(n))

          !!NG: 26/06/2009: BUGGED CALL sub_orca_east_west_periodic(i1(n), min0(i1(n),i2(n))+1, i0(n))
          !!NG: 26/06/2009: BUGGED CALL sub_orca_east_west_periodic(i2(n), min0(i1(n),i2(n))+1, i0(n))
          !!NG: 26/06/2009: BUGGED CALL sub_orca_east_west_periodic(hi(n), min0(i1(n),i2(n))+1, i0(n))
          !!NG: 26/06/2009: BUGGED CALL sub_orca_east_west_periodic(i0(n), min0(i1(n),i2(n))+1, i0(n), iperfl(n))

          !!NG: 26/06/2009: Replace by:
          imin_tmp = min0(i1(n),i2(n))
          imax_tmp = max0(i1(n),i2(n))
          CALL sub_orca_east_west_periodic(i1(n), imin_tmp+1, imax_tmp)
          CALL sub_orca_east_west_periodic(i2(n), imin_tmp+1, imax_tmp)
          CALL sub_orca_east_west_periodic(hi(n), imin_tmp+1, imax_tmp)
          CALL sub_orca_east_west_periodic(i0(n), imin_tmp+1, imax_tmp, iperfl(n))

        ENDIF

        !====================================================================
        ! age update
        !====================================================================
        !
        ! BBL __________________________________________________________________________
        ! BBL - BBL: L'age des particules (en secondes) est mis a jour a partir de la
        ! BBL         valeur de TTT
        ! BBL
        agenew(n)=age(n)+ttt(n)*tfac(n)

        ! BBL __________________________________________________________________________
        ! BBL - BBL: On met a jour l'indice temporel des particules, selon le sens de
        ! BBL         l'integration (BACKWARD ou FORWARD)
        ! BBL
        IF (TRIM(forback) == 'forward') THEN
          hl(n) = gl(n) + ttt(n) * tfac(n) / tfic
        ELSE
          hl(n) = gl(n) - ttt(n) * tfac(n) / tfic
        ENDIF

        ! BBL __________________________________________________________________________
        ! BBL - BBL: Si une particule a ete stoppee dans l'attente d'une mise a jour du
        ! BBL         champ de transport charge en memoire, on force la valeur de
        ! BBL         l'indice temporel a l'instant exact de la transition
        ! BBL
        IF (igo(n) == 2)  THEN
          hl(n)     = NINT(hl(n)+.5_rprec) - .5_rprec
          agenew(n) = tswap(n)
        ENDIF

        !!NG: 20 may 2010
        !!NG: to correct some precision problem due to float operations with
        !!NG: high number. We force the new age to be exactly at the tfin value
        IF (igo(n) == 3)  THEN
          agenew(n) = tfin
        ENDIF
      ENDIF

    ENDDO


    !!================================================================================!!
    !!============================== QUALITATIVE MODE ================================!!
    !!================================================================================!!
    !
    !====================================================================
    ! we stop the integration if we reach a limit of the domain (it may
    !  happen in QUALITATIVE experiments) [QUALI]
    !====================================================================
    IF (TRIM(mode) == 'qualitative') THEN

      !!-------------------------------------!!
      !!------ DO LOOP: ON TRAJECTORIES -----!!
      !!-------------------------------------!!
      DO  n=1,ntraj

        IF (igo(n) >= 1) THEN
          ! BBL _________________________________________________________________
          ! BBL - BBL: On ne fait les calculs suivants que 
          ! BBL        pour les particules *actives* (etat IGO=1)
          !

          IF (k0(n) == 0) THEN
            ! BBL ______________________________________________________________
            ! BBL -- BBL: Traitement des sorties de domaine *par la surface*
            ! BBL
            WRITE(lun_error,7000)'out of zdomain', &
                 n,i0(n),j0(n),k0(n),l0(n),tfi(n),tfj(n),tfk(n),tfl(n)
            igo(n) = 3
          ENDIF

          !! NG 10 mai 2010: BUG IF ((igo(n) /= 3).AND.((j0(n) == dims_reg(2,1)).OR.(j0(n) == dims_reg(2,2)))) THEN
          IF (igo(n) /= 3) THEN
            IF ( (j0(n) == dims_reg(2,1)).OR.&
                 ( (j0(n) == dims_reg(2,2)).AND.(.NOT.key_jfold) ) ) THEN

              ! BBL _______________________________________________________________
              ! BBL -- BBL: Traitement des sorties de domaine 
              ! BBL         *par une frontiere meridienne*
              ! BBL
              WRITE(lun_error,7000)'out of ydomain', &
                   n,i0(n),j0(n),k0(n),l0(n),tfi(n),tfj(n),tfk(n),tfl(n)
              igo(n)=3
            END IF
          ENDIF

          !!NG          IF (.NOT.key_jfold) THEN
          !!NG            IF ((igo(n) /= 3).AND.(j0(n) == dims_reg(2,2))) THEN
          !!NG              ! BBL ______________________________________________________________
          !!NG              ! BBL -- BBL: Il y a un BUG ici, car le test suivant n'est jamais
          !!NG              ! BBL         verifie (le test precedent a deja fait le travail)
          !!NG              ! BBL         Il faudrait reecrire les tests de maniere plus propre 
          !!NG              ! BBL         dans le cas d'une condition de repliement meridien
          !!NG              ! BBL
          !!NG              WRITE(lun_error,7000)'out of Ydomain', &
          !!NG                   n,i0(n),j0(n),k0(n),l0(n),tfi(n),tfj(n),tfk(n),tfl(n)
          !!NG              igo(n)=3
          !!NG            ENDIF
          !!NG          ENDIF ! .NOT.key_jfold

          IF (.NOT.key_periodic) THEN
            IF ((igo(n) /= 3).AND.((i0(n) == dims_reg(1,1)).OR.(i0(n) == dims_reg(1,2)))) THEN
              ! BBL _________________________________________________________
              ! BBL -- BBL: Traitement des sorties de domaine 
              ! BBL         *par une frontiere zonale*
              ! BBL
              WRITE(lun_error,7000)'out of xdomain', &
                   n,i0(n),j0(n),k0(n),l0(n),tfi(n),tfj(n),tfk(n),tfl(n)
              igo(n)=3
            ENDIF
          ENDIF ! .NOT.key_periodic

        ENDIF
      ENDDO
      !
    ENDIF ! (TRIM(mode) == 'qualitative') => QUALITATIVE MODE

    !!================================================================================!!
    !!============================== QUANTITATIVE MODE ===============================!!
    !!================================================================================!!
    !
    !====================================================================
    ! quantitative analysis
    !
    ! in QUANTITATIVE experiments, the 2D projections of the transport
    !  field associated with 1 particle is summed depending on the final
    !  section [QUANT]
    !====================================================================
    ! BBL __________________________________________________________________________
    ! BBL BBL: En mode QUANTITATIVE, on teste l'atteinte des sections d'interception
    ! BBL
    IF (TRIM(mode) == 'quantitative') THEN

      !!-------------------------------------!!
      !!------ DO LOOP: ON TRAJECTORIES -----!!
      !!-------------------------------------!!
      DO n=1,ntraj

        IF (igo(n) > 0) THEN
          ! BBL __________________________________________________________________________
          ! BBL - BBL: On ne fait les calculs suivants que pour les particules *actives*
          ! BBL         (etat IGO=1,2,3)
          !------------------------------------------------------------------
          ! mtfin informs about the sections used in quantitative experiments
          !------------------------------------------------------------------
          nfin(n) = mtfin(i0(n),j0(n),k0(n))

          !!NG: REMARK - nothing here about criter2 like in trajec.f90 ?????

          !===============================================================
          ! a final criterion may define a "final section" (criter1)
          !===============================================================
          IF ((nfin(n) > 0).OR. &
               (criter1(tfi(n),tfj(n),tfk(n),tfl(n), &
               hi(n),hj(n),hk(n),hl(n),i0(n),j0(n),k0(n),l0(n)))) THEN

            !----------------------------------------------------------
            ! a final criterion uses "1 + nsect" as the section index
            ! we flag the existence of such a criterion
            !----------------------------------------------------------
            ! BBL __________________________________________________________________________
            ! BBL - BBL: Quand une particule atteint une section d'interception ou verifie
            ! BBL        critere d'interception (CRITER1), son integration s'arrete et son
            ! BBL        etat passe a IGO=3
            ! BBL
            igo(n)=3
            IF (nfin(n) <= 0) THEN
              nfin(n) = 1 + nsect
              icrit1 = 1
            ENDIF
            !-----------------------------------------------------------------
            ! we distinguish meandering and recirculating particles (criter0)
            !-----------------------------------------------------------------
            IF ((TRIM(bin) == 'nobin').AND.(nfin(n) == 1) .AND. &
                 criter0(hi(n),hj(n),hk(n),hl(n),i0(n),j0(n),k0(n),l0(n))) THEN

              nfin(n) = 0

            ENDIF
            !
            ! hereafter,
            ! nfin=0, if the initial section (with criter0) has been reached
            ! nfin=1, if the initial section (no criter0) or a final section has
            !  been reached
            ! nfin=1+nsect, if a final criterion (criter1) is satisfied
            !
            sectrans(nfin(n))=sectrans(nfin(n))+ttr(n)

            !???????????????????????????????????????????????????????????????????
            ! there IS NOW storage of transport projections for "criter1" cases
            !  (but never try to compute 2D streamfunctions with these fields !)
            !???????????????????????????????????????????????????????????????????

            !==================================================================
            ! statistics...
            !==================================================================

            IF (key_alltracers) THEN
              CALL sub_quant_statistics(               &
                   nfin(n), ttr(n), agenew(n), tcyc  , &
                   hi(n)  , hj(n) , hk(n) , 1._rprec , &
                   tfi(n) , tfj(n), tfk(n), 1._rprec , &
                   i0(n)  , j0(n) , k0(n) , 1        , &
                   tf(n)  , sf(n) , rf(n) , truage(n), &
                   ti(n), si(n), ri(n))
            ELSE
              CALL sub_quant_statistics(               &
                   nfin(n), ttr(n), agenew(n), tcyc  , &
                   hi(n)  , hj(n) , hk(n) , 1._rprec , &
                   tfi(n) , tfj(n), tfk(n), 1._rprec , &
                   i0(n)  , j0(n) , k0(n) , 1        , &
                   tf(n)  , sf(n) , rf(n) , truage(n)  )
            endif

            !! NG: 15_09_2008
            !! NG: Here we store the final depth which depend of time for ROMS
            !! NG: and symphonie.
            IF (key_roms.OR.key_mars.OR.key_symphonie) THEN
              prof=fz(gi=hi(n), gj=hj(n), gk=hk(n), il=1)
            ELSE
              prof=fz(gi=hi(n), gj=hj(n), gk=hk(n))
            ENDIF
            final_depth(n)=prof

          ENDIF

        ENDIF

      ENDDO

    ENDIF ! end of quantitative analysis

    !!----------------------------------------------------------------------------!!
    !!---6666666666666666666---- DO LOOP : ON TRAJECTORIES ---6666666666666666----!!
    !!----------------------------------------------------------------------------!!
    DO n=1,ntraj

      IF (igo(n) == 2) THEN
        ! BBL __________________________________________________________________________
        ! BBL - BBL: On applique ici un traitement specifique aux particules dans
        ! BBL         l'attente d'une mise a jour en memoire du champ de transport
        ! BBL __________________________________________________________________________
        ! BBL - BBL: On incremente (ou decremente) d'une unite l'indice temporel
        ! BBL
        IF (TRIM(forback) == 'forward') THEN
          l0(n) = l0(n) + 1
        ELSE
          l0(n) = l0(n) - 1
        ENDIF
        ! BBL __________________________________________________________________________
        ! BBL - BBL: Si necessaire, on boucle sur le cycle periodique
        ! BBL
        IF (l0(n) < 1)   l0(n) = l0(n) + lmt
        IF (l0(n) > lmt) l0(n) = l0(n) - lmt
        ! BBL __________________________________________________________________________
        ! BBL - BBL: On met a jour pour chaque particule l'instant a utiliser pour le
        ! BBL         prochain swap temporel en memoire (prochaine mise a jour du champ
        ! BBL         de transport)
        ! BBL
        tswap(n)=tswap(n)+tfic
      ENDIF

    ENDDO

    !!----------------------------------------------------------------------------!!
    !!---7777777777777777777---- DO LOOP : ON TRAJECTORIES ---7777777777777777----!!
    !!----------------------------------------------------------------------------!!
    DO  n=1,ntraj

      IF ((igo(n) > 0).AND.(igo(n) /= 3)) THEN
        ! BBL __________________________________________________________________________
        ! BBL - BBL: On applique ici un traitement specifique aux particules dont
        ! BBL         l'integration n'est pas encore terminee (mais deja commencee)
        ! BBL        On verifie que l'on n'a pas atteint un point TERRE

        CALL sub_reducmem_shift_or_not_ind(i0(n),j0(n),k0(n),is,js,ks)

        IF (tmask(is,js,ks,1) == 0) THEN
          WRITE(lun_error,7000)'   coast crash', &
               n,i0(n),j0(n),k0(n),l0(n),tfi(n),tfj(n),tfk(n),tfl(n)
          WRITE(lun_output,*)' '
          WRITE(lun_output,*)'coast crash occurs in some cases where it seems it should not '
          WRITE(lun_output,*)'    further investigation necessary'
          WRITE(lun_output,*)'coast crash for particle # ',n
          WRITE(lun_output,*)'i j k: ',i0(n),' ',j0(n),' ',k0(n)
          WRITE(lun_output,*)'is js ks: ',is,' ',js,' ',ks
          WRITE(lun_output,*)'gi gj gk:',gi(n),' ',gj(n),' ',gk(n)
          WRITE(lun_output,*)'hi hj hk:',hi(n),' ',hj(n),' ',hk(n)
          !! NG: Bug          WRITE(lun_output,*)'u: ',uu(i0(n)-1,j0(n),k0(n),1),uu(i0(n),j0(n),k0(n),1)
          !! NG: Bug          WRITE(lun_output,*)'v: ',vv(i0(n),j0(n)-1,k0(n),1),vv(i0(n),j0(n),k0(n),1)
          !! NG: Bug          WRITE(lun_output,*)'w: ',ww(i0(n),j0(n),k0(n)+1,1),ww(i0(n),j0(n),k0(n),1)
          WRITE(lun_output,*)'u: ',uu(is-1,js,ks,1),uu(is,js,ks,1)
          WRITE(lun_output,*)'v: ',vv(is,js-1,ks,1),vv(is,js,ks,1)
          WRITE(lun_output,*)'w: ',ww(is,js,ks+1,1),ww(is,js,ks,1)
          WRITE(lun_output,*)'tx ty tz ttt: ',tx(n),ty(n),tz(n),ttt(n)
          ! Do not actually stop- just write warning, because indexing seems sometimes wrong (maybe)
       ENDIF
        !
      ENDIF

    ENDDO !  n=1,ntraj

    !!================================================================================!!
    !!============================== QUALITATIVE MODE ===============================!!
    !!================================================================================!!
    !
    ! do we need an output of the position (qualitative experiment) ?
    ! (linear interpolation)
    ! - age: previous age
    ! - agenew: current age
    ! - tnext: output time
    !
    IF (TRIM(mode) == 'qualitative') THEN
      ! BBL __________________________________________________________________________
      ! BBL BBL: En mode QUALITATIVE, on s'interesse ici au calcul des positions
      ! BBL        GEOGRAPHIQUES des particules
      ! BBL       ARIANE a en memoire une position ancienne (GI, GJ, GK, GL) et une
      ! BBL        position nouvelle (HI, HJ, HK, HL) pour chaque particule active, et
      ! BBL         on sait qu'il faut calculer une position geographique toutes les
      ! BBL        TNEXT secondes
      ! BBL       Les valeurs de GL et HL n'etant pas les memes pour toutes les
      ! BBL        particules, on travaille ici de maniere incrementale, et les
      ! BBL        calculs se font par incrementation temporelle, *tant que cela
      ! BBL        s'avere necessaire*
      ! BBL       Par defaut, toutes les particules sont *bonnes pour le service*
      ! BBL        (IFLAG=1)
      ! BBL

      iflag(:) = 1 ! array notation

      !----------
      ! BBL __________________________________________________________________________
      ! BBL BBL: La variable IFL sera mise a la valeur 1 si une particule au moins a
      ! BBL       besoin de voir echantillonner sa trajectoire GEOGRAPHIQUE
      ! BBL
1111  ifl=0

      !! NG 29 april 2010
      flag_write_position = 0_iprec

      !!-------------------------------------!!
      !!------ DO LOOP: ON TRAJECTORIES -----!!
      !!-------------------------------------!!
      DO n=1,ntraj

        IF (igo(n) > 0) THEN
          ! BBL __________________________________________________________________________
          ! BBL - BBL: On ne fait les calculs suivants que pour les particules *actives*
          ! BBL         (etat IGO=1,2,3)
          !
          IF ((iflag(n) == 1).AND.(agenew(n) >= tnext(n))) THEN

            ! BBL __________________________________________________________________________
            ! BBL - BBL: Si l'age *nouveau* AGENEW d'une particule excede l'instant TNEXT
            ! BBL         correspondant a la demande d'une sortie GEOGRAPHIQUE, on calcule
            ! BBL         une position par interpolation entre les positions (GI, GJ, GK) et
            ! BBL         (HI, HJ, HK), avec un facteur de ponderation FT lie a
            ! BBL         l'emplacement temporel de TNEXT entre AGE et AGENEW
            ! BBL        NOTE: XT, YT et ZT ne sont pas des positions GEOGRAPHIQUES, mais
            ! BBL         des indices sur la grille du modele, tout comme GI, GJ et GK, et
            ! BBL         HI, HJ et HK
            ! BBL
            ft=(tnext(n)-age(n))/(agenew(n)-age(n))
            xt=gi(n)+(hi(n)-gi(n))*ft
            yt=gj(n)+(hj(n)-gj(n))*ft
            zt=gk(n)+(hk(n)-gk(n))*ft

            IF (TRIM(forback) == 'forward') THEN
              st=gl(n)+(tnext(n)-age(n))/tfic
            ELSE
              st=gl(n)-(tnext(n)-age(n))/tfic
            ENDIF

            IF (key_periodic) THEN
              ! BBL _____________________________________________________________________
              ! BBL -- BBL: La variable IPERFL est utilisee pour corriger le cas echeant
              ! BBL          l'interpolation lineaire precedente (donnant un resultat 
              ! BBL          faux si les deux positions de reference, GI et HI, sont 
              ! BBL          situes de part et d'autre de la periodicite
              ! BBL
              IF (iperfl(n) == 1) THEN
                xt=gi(n)+(hi(n)-REAL(imt-2,kind=rprec)-gi(n))*ft
              ENDIF
              IF (iperfl(n) == 2) THEN
                xt=gi(n)+(hi(n)+REAL(imt-2,kind=rprec)-gi(n))*ft
              ENDIF

            ENDIF


            IF (key_jfold) THEN

              ! extrapolation not yet implemented
              IF (jperfl(n) == 1) THEN
                xt=gi(n)
                yt=gj(n)
                zt=gk(n)
              ENDIF

            ENDIF


            !---------------------------!
            !- position output on file -!
            !---------------------------!
            ! BBL __________________________________________________________________________
            ! BBL -- BBL: On effectue ici la conversion entre indice de grille et
            ! BBL          coordonnees geographiques
            ! BBL
            x0=fx(xt,yt)
            y0=fy(xt,yt)

            IF (zbuoy(n) > 1.e8_rprec) THEN
              IF (key_roms.OR.key_mars.OR.key_symphonie) THEN
                z0 = fz(gi=xt, gj=yt, gk=zt, il = 1)
              ELSE
                z0 = fz(gi=xt, gj=yt, gk=zt)
              ENDIF
            ELSE

              ! BBL ________________________________________________________________________
              ! BBL -- BBL: Pour les particules isobares, aucun calcul n'est necessaire pour
              ! BBL          l'immersion...
              ! BBL
              z0=zbuoy(n)
            ENDIF

            IF (tnext(n) <= tfin) THEN

              !! NG 29 april 2010
              IF (small_virtual_memory) THEN

                flag_write_position = 1_iprec

                ! BBL _______________________________________________________________________
                ! BBL -- BBL: L'ecriture sur fichier se fait si on n'a pas atteint le temps
                ! BBL          limite TFIN

                !! NG 29 april 2010 ind_traj(n) = ind_traj(n) + 1

                !! NG 29 april 2010 array_qual_traj(n,ind_traj(n),1) = x0
                !! NG 29 april 2010 array_qual_traj(n,ind_traj(n),2) = y0
                !! NG 29 april 2010 array_qual_traj(n,ind_traj(n),3) = z0
                !! NG 29 april 2010 array_qual_traj(n,ind_traj(n),4) = tnext(n)/tcyc
                array_qual_traj_2D(n,1) = x0
                array_qual_traj_2D(n,2) = y0
                array_qual_traj_2D(n,3) = z0
                array_qual_traj_2D(n,4) = tnext(n)/tcyc

                !!NG              To print the vertical velocity
                !!NG              IF (n == 1) THEN
                !!NG              CALL sub_reducmem_shift_or_not_ind(i0(n),j0(n),k0(n),is,js,ks)
                !!NG              WRITE(0,*) z0, ( ww(is,js,ks,1)   - &
                !!NG                   ( gk(n) + REAL(k0(n), kind = rprec))   * & 
                !!NG                   (ww(is,js,ks+1,1) - ww(is,js,ks,1) )) / &
                !!NG                   (e1t(is,js,1,1)*e2t(is,js,1,1))
                !!NG              ENDIF

                IF (key_alltracers) THEN

                  tf(n) = zinter(tt,xt,yt,zt, 1._rprec)
                  sf(n) = zinter(ss,xt,yt,zt, 1._rprec)

                  IF (key_approximatesigma) THEN

                    rf(n) = zinter(rr,xt,yt,zt, 1._rprec)

                  ELSE 

                    rf(n)=sigma(zsigma,sf(n),tf(n))

                  ENDIF

                  !! NG 29 april 2010 array_qual_traj(n,ind_traj(n),5) = tf(n)
                  !! NG 29 april 2010 array_qual_traj(n,ind_traj(n),6) = sf(n)
                  !! NG 29 april 2010 array_qual_traj(n,ind_traj(n),7) = rf(n)
                  array_qual_traj_2D(n,5) = tf(n)
                  array_qual_traj_2D(n,6) = sf(n)
                  array_qual_traj_2D(n,7) = rf(n)

                  IF (key_iU_jV_kW) THEN
                    !! NG 29 april 2010 array_qual_traj(n,ind_traj(n),8)  = xt
                    !! NG 29 april 2010 array_qual_traj(n,ind_traj(n),9)  = yt
                    !! NG 29 april 2010 array_qual_traj(n,ind_traj(n),10) = zt
                    array_qual_traj_2D(n,8)  = xt
                    array_qual_traj_2D(n,9)  = yt
                    array_qual_traj_2D(n,10) = zt
                  END IF

                ELSE

                  IF (key_iU_jV_kW) THEN
                    !! NG 29 april 2010 array_qual_traj(n,ind_traj(n),5) = xt
                    !! NG 29 april 2010 array_qual_traj(n,ind_traj(n),6) = yt
                    !! NG 29 april 2010 array_qual_traj(n,ind_traj(n),7) = zt
                    array_qual_traj_2D(n,5) = xt
                    array_qual_traj_2D(n,6) = yt
                    array_qual_traj_2D(n,7) = zt
                  END IF

                ENDIF ! key_alltracers

              ELSE

                ! BBL _______________________________________________________________________
                ! BBL -- BBL: L'ecriture sur fichier se fait si on n'a pas atteint le temps
                ! BBL          limite TFIN

                ind_traj(n) = ind_traj(n) + iOne

                array_qual_traj(n,ind_traj(n),1) = x0
                array_qual_traj(n,ind_traj(n),2) = y0
                array_qual_traj(n,ind_traj(n),3) = z0
                array_qual_traj(n,ind_traj(n),4) = tnext(n)/tcyc

                !!NG              To print the vertical velocity
                !!NG              IF (n == 1) THEN
                !!NG              CALL sub_reducmem_shift_or_not_ind(i0(n),j0(n),k0(n),is,js,ks)
                !!NG              WRITE(0,*) z0, ( ww(is,js,ks,1)   - &
                !!NG                   ( gk(n) + REAL(k0(n), kind = rprec))   * & 
                !!NG                   (ww(is,js,ks+1,1) - ww(is,js,ks,1) )) / &
                !!NG                   (e1t(is,js,1,1)*e2t(is,js,1,1))
                !!NG              ENDIF

                IF (key_alltracers) THEN

                  tf(n) = zinter(tt,xt,yt,zt, 1._rprec)
                  sf(n) = zinter(ss,xt,yt,zt, 1._rprec)

                  IF (key_approximatesigma) THEN

                    rf(n) = zinter(rr,xt,yt,zt, 1._rprec)

                  ELSE 

                    rf(n)=sigma(zsigma,sf(n),tf(n))

                  ENDIF

                  array_qual_traj(n,ind_traj(n),5) = tf(n)
                  array_qual_traj(n,ind_traj(n),6) = sf(n)
                  array_qual_traj(n,ind_traj(n),7) = rf(n)

                  IF (key_iU_jV_kW) THEN

                    array_qual_traj(n,ind_traj(n),8)  = xt
                    array_qual_traj(n,ind_traj(n),9)  = yt
                    array_qual_traj(n,ind_traj(n),10) = zt

                  END IF

                ELSE

                  IF (key_iU_jV_kW) THEN

                    array_qual_traj(n,ind_traj(n),5) = xt
                    array_qual_traj(n,ind_traj(n),6) = yt
                    array_qual_traj(n,ind_traj(n),7) = zt

                  END IF

                ENDIF ! key_alltracers

              ENDIF
              !! NG 29 april 2010 

            ENDIF
            !
            ! BBL __________________________________________________________________________
            ! BBL -- BBL: On met a jour la date d'ecriture de la position suivante
            ! BBL

            !!NG: 9 march 2011: a Bug here because tlap + tlap +... + tlap is
            !!NG: 9 march 2011:  different that tlap * x
            !!NG: 9 march 2011: tnext(n) = tnext(n) + tlap

            inext(n) = inext(n) + iOne
            tnext(n) = tlap * REAL(inext(n), kind = rprec)

            ! BBL __________________________________________________________________________
            ! BBL -- BBL: Puisqu'une particule est *passee par ce morceau de code* il faut
            ! BBL          continuer les calculs et on utilise IFL=1
            ! BBL
            ifl = iOne
          ELSE
            ! BBL __________________________________________________________________________
            ! BBL -- BBL: Cette particule n'a plus de raison de voir un calcul de positions
            ! BBL          geographiques
            ! BBL
            iflag(n)= iZero
          ENDIF
          !
        ENDIF
      ENDDO

      !! NG 29 april 2010
      !! Write output trajectories: The initial positions

      IF (small_virtual_memory) THEN

        IF (flag_write_position  == iOne) THEN

          flag_write_position  = iZero

          IF (key_alltracers.AND.key_iU_jV_kW) THEN

            CALL sub_save_netcdf_trajectories_generic(             &
                 indice_traj                                     , &
                 array_qual_traj_2D(:,1), array_qual_traj_2D(:,2), &
                 array_qual_traj_2D(:,3), array_qual_traj_2D(:,4), &
                 array_qual_traj_2D(:,5), array_qual_traj_2D(:,6), &
                 array_qual_traj_2D(:,7), array_qual_traj_2D(:,8), &
                 array_qual_traj_2D(:,9), array_qual_traj_2D(:,10) )

          ELSEIF (key_alltracers.AND.(.NOT.key_iU_jV_kW)) THEN

            CALL sub_save_netcdf_trajectories_generic(             &
                 indice_traj                                     , &
                 array_qual_traj_2D(:,1), array_qual_traj_2D(:,2), &
                 array_qual_traj_2D(:,3), array_qual_traj_2D(:,4), &
                 traj_temp=array_qual_traj_2D(:,5), &
                 traj_salt=array_qual_traj_2D(:,6), &
                 traj_dens=array_qual_traj_2D(:,7))

          ELSEIF ((.NOT.key_alltracers).AND.key_iU_jV_kW) THEN
            CALL sub_save_netcdf_trajectories_generic(             &
                 indice_traj                                     , &
                 array_qual_traj_2D(:,1), array_qual_traj_2D(:,2), &
                 array_qual_traj_2D(:,3), array_qual_traj_2D(:,4), &
                 traj_iU=array_qual_traj_2D(:,5), &
                 traj_jV=array_qual_traj_2D(:,6), &
                 traj_kW=array_qual_traj_2D(:,7))  
          ELSE
            CALL sub_save_netcdf_trajectories_generic(             &
                 indice_traj                                     , &
                 array_qual_traj_2D(:,1), array_qual_traj_2D(:,2), &
                 array_qual_traj_2D(:,3), array_qual_traj_2D(:,4)  )
          ENDIF

          array_qual_traj_2D(:,:) = mask_value

          ind_traj = ind_traj + iOne

        END IF
      ENDIF
      !! NG 29 april 2010

      !
      ! BBL __________________________________________________________________________
      ! BBL BBL: Si IFL=1, une particule au moins necessite toujours un calcul de
      ! BBL       positions geographiques
      ! BBL
      IF (ifl == iOne) GOTO 1111
      !
    ENDIF

    !!----------------------------------------------------------------------------!!
    !!--88888888888888888888---- DO LOOP : ON TRAJECTORIES ---8888888888888888----!!
    !!----------------------------------------------------------------------------!!
    !
    ! swap for space and time index positions
    !
    ! BBL __________________________________________________________________________
    ! BBL BBL: Mise a jour des parametres de position spatiale et temporelle
    ! BBL
    DO n=1,ntraj

      IF (igo(n) > iZero) THEN
        ! BBL __________________________________________________________________________
        ! BBL - BBL: On ne fait les calculs suivants que pour les particules *actives*
        ! BBL         (etat IGO=1,2,3)
        ! BBL
        gi(n)  = hi(n)
        gj(n)  = hj(n)
        gk(n)  = hk(n)
        gl(n)  = hl(n)
        age(n) = agenew(n)
      ENDIF

    ENDDO

    !!===================================!!
    !!= TESTS TO VERIFY IF GAME IS OVER =!!
    !!===================================!!

    !!----------------------------------------------------------------------------!!
    !!--88888888888888888888---- DO LOOP : ON TRAJECTORIES ---8888888888888888----!!
    !!----------------------------------------------------------------------------!!
    ! BBL __________________________________________________________________________
    ! BBL BBL: IFL controle ici la necessite (IFL=1) ou non (IFL=0) de continuer les
    ! BBL       calculs avec le MEME champ de transport charge en memoire
    ! BBL
    ifl=0

    DO  n=1,ntraj
      ! BBL __________________________________________________________________________
      ! BBL BBL: On *masque* d'abord les particules dans l'attente de la mise a jour
      ! BBL       du champ de transport (IGO=-2)
      ! BBL
      IF (igo(n) == 2) igo(n) = -2
      ! BBL __________________________________________________________________________
      ! BBL BBL: On *neutralise* definitivement les particules ayant atteint le terme
      ! BBL       de leur trajectoire (IGO=-3)
      ! BBL
      IF (igo(n) == 3) igo(n) = -3
      ! BBL __________________________________________________________________________
      ! BBL BBL: Si une particule est toujours active (IGO=1), il faudra continuer
      ! BBL       (IFL=1)
      ! BBL
      IF (igo(n) == 1) ifl = 1
    ENDDO
    !
    ! BBL __________________________________________________________________________
    ! BBL BBL: On continue donc avec le meme champ de transport (IFL=1)
    ! BBL
    IF (ifl == 1) THEN
      GOTO 2222 !! GAME IS NOT OVER !!
    ENDIF

    !!----------------------------------------------------------------------------!!
    !!--99999999999999999999---- DO LOOP : ON TRAJECTORIES ---9999999999999999----!!
    !!----------------------------------------------------------------------------!!
    ! BBL __________________________________________________________________________
    ! BBL BBL: IFL controle maintenant la necessite (IFL=1) ou non (IFL=0) de
    ! BBL       continuer les calculs, en chargeant en memoire un nouveau champ de
    ! BBL       transport
    ! BBL
    ifl=0

    DO  n=1,ntraj
      ! BBL __________________________________________________________________________
      ! BBL BBL: Si une particule n'a pas atteint le terme de sa trajectoire
      ! BBL        (IGO NE -3), il faudra continuer (IFL=1)
      ! BBL
      IF (igo(n) /= -3) ifl=1
    ENDDO
    !
    ! BBL __________________________________________________________________________
    ! BBL BBL: Pas la peine de continuer...
    ! BBL
    IF (ifl == 0) THEN
      GOTO 3333 !! GAME IS OVER !!
    ENDIF
    !
    ! BBL __________________________________________________________________________
    ! BBL BBL: Il faut continuer, et mettre a jour l'indice temporel du champ de
    ! BBL       transport a lire, en prenant garde au sens d'integration (FORWARD ou
    ! BBL        BACKWARD)
    ! BBL      Tant que l'on reste dans le meme cycle temporel, on relance vers le
    ! BBL       label 9998, sinon on enchaine sur le label 9999 qui correspond a
    ! BBL       l'incrementation d'un nouveau cycle temporel
    ! BBL
    IF (TRIM(forback) == 'forward') THEN
      lref = lref + 1
      IF (lref <= lmt) THEN
        GOTO 9998
      ENDIF

    ELSE ! backward
      lref = lref - 1
      IF (lref >= 1) THEN
        GOTO 9998
      ENDIF

    ENDIF

  ENDDO ! DO LOOP: CYCLES (iter=1,maxcycles)

  !!==================================================================================!!
  !!============================ MAIN LOOP IS FINISHED ===============================!!
  !!==================================================================================!!

3333 CONTINUE

  IF (TRIM(mode) == 'quantitative') THEN
    ndone=0

    DO  n = 1, ntraj

      ! BBL ___________________________________________________________________
      ! BBL -- BBL: IGO=3 ou -3 signifie que l'integration de la particule est
      ! BBL          terminee
      ! BBL
      IF (ABS(igo(n)) >= 3) ndone = ndone + 1

    ENDDO

    WRITE(lun_standard,*)
    WRITE(lun_standard,*)'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
    WRITE(lun_standard,*)'HOW MANY PARTICLES REACH A SECTION: ', &
         INT(100._rprec*REAL(ndone,kind=rprec)/REAL(ntraj,kind=rprec)), '%'
    WRITE(lun_standard,*)'    - number of particules         : ', ntraj
    WRITE(lun_standard,*)'    - done                         : ', ndone
    WRITE(lun_standard,*)'    - rest                         : ', ntraj - ndone
    WRITE(lun_standard,*)'    - cycle                        : ', iter-1
    WRITE(lun_standard,*)'<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
  ENDIF

  !- Comments -!
  WRITE(lun_standard,*)''
  WRITE(lun_standard,*)'=================================='
  WRITE(lun_standard,*)'= MAIN LOOP ON CYCLE >>>>>> EXIT ='
  WRITE(lun_standard,*)'=================================='

  !!----------------------------------------------------------------------------!!
  !!--10-10-10-10-10-10-10---- DO LOOP : ON TRAJECTORIES ---10-10-10-10-10-10---!!
  !!----------------------------------------------------------------------------!!

  WRITE(lun_standard,*)''
  WRITE(lun_standard,*)'-----------------------'
  WRITE(lun_standard,*)'= Writing Output Data ='
  WRITE(lun_standard,*)'-----------------------'

  !!===============================================================================!!
  !!============================== QUALITATIVE MODE ===============================!!
  !!===============================================================================!!
  IF (TRIM(mode) == 'qualitative' ) THEN

    !----------------!
    !- Trajectories -!
    !----------------!
    !- Netcdf -!
    !----------!

    IF (small_virtual_memory) THEN

    ELSE

      IF (key_alltracers.AND.key_iU_jV_kW) THEN

        CALL sub_save_netcdf_trajectories_generic(                   &
             array_qual_traj(:,:,1), array_qual_traj(:,:,2), &
             array_qual_traj(:,:,3), array_qual_traj(:,:,4), &
             array_qual_traj(:,:,5), array_qual_traj(:,:,6), &
             array_qual_traj(:,:,7), array_qual_traj(:,:,8), &
             array_qual_traj(:,:,9), array_qual_traj(:,:,10) )

      ELSEIF (key_alltracers.AND.(.NOT.key_iU_jV_kW)) THEN

        CALL sub_save_netcdf_trajectories_generic(                   &
             array_qual_traj(:,:,1), array_qual_traj(:,:,2), &
             array_qual_traj(:,:,3), array_qual_traj(:,:,4), &
             traj_temp=array_qual_traj(:,:,5), &
             traj_salt=array_qual_traj(:,:,6), &
             traj_dens=array_qual_traj(:,:,7))

      ELSEIF ((.NOT.key_alltracers).AND.key_iU_jV_kW) THEN
        CALL sub_save_netcdf_trajectories_generic(                   &
             array_qual_traj(:,:,1), array_qual_traj(:,:,2), &
             array_qual_traj(:,:,3), array_qual_traj(:,:,4), &
             traj_iU=array_qual_traj(:,:,5), &
             traj_jV=array_qual_traj(:,:,6), &
             traj_kW=array_qual_traj(:,:,7))  
      ELSE
        CALL sub_save_netcdf_trajectories_generic(                   &
             array_qual_traj(:,:,1), array_qual_traj(:,:,2), &
             array_qual_traj(:,:,3), array_qual_traj(:,:,4)  )
      ENDIF

      !---------!
      !- ASCII -!
      !---------!
      IF (key_ascii_outputs) THEN

        DO n = 1, ntraj
          DO i = 1, ind_traj(n)
            WRITE(lun_traj,7055)n,array_qual_traj(n,i,1:nstat)
          ENDDO
        ENDDO

      ENDIF
    ENDIF

7055 FORMAT (i0,3(1x,f0.5),1x,f0.5,1x,3(1x,f0.5))

    !!================================================================================!!
    !!============================== QUANTITATIVE MODE ===============================!!
    !!================================================================================!!
    !
    ! BBL __________________________________________________________________________
    ! BBL BBL: On trouve ensuite la suite logique de l'archivage (mode QUANTITATIVE)
    ! BBL
  ELSE

    !-------------------!
    !- Final positions -!
    !-------------------! 
    !---------!
    !- ASCII -!
    !---------!
    IF (key_ascii_outputs) THEN
      DO  n=1,ntraj

        IF (nfin(n) >= 0) THEN
          WRITE(lun_output,7000)secname(nfin(n)), &
               n,i0(n),j0(n),k0(n),l0(n),tfi(n),tfj(n),tfk(n),tfl(n)
        ELSE
          WRITE(lun_output,7000)'meanders (-1)', &
               n,i0(n),j0(n),k0(n),l0(n),tfi(n),tfj(n),tfk(n),tfl(n)
        ENDIF

        IF (key_alltracers) THEN
          WRITE(lun_fin_pos,7059)hi(n),hj(n),-hk(n),hl(n), &
               ttr(n),truage(n),nfin(n),tf(n),sf(n),rf(n)
        ELSE
          WRITE(lun_fin_pos,7059)hi(n),hj(n),-hk(n),hl(n), &
               ttr(n),truage(n),nfin(n)
        ENDIF

      ENDDO

    ENDIF

7000 FORMAT(a14,', #',i0,' en:', 1x,i0,1x,i0,1x,i0,1x,i0,',  init.=',4(1x,f0.2),': ',f0.2)

7059 FORMAT(3(1x,f0.3),1x,f0.3,1x,f0.3,1x,f0.3,1x,i0,3(1x,f0.3))

    !----------------------------!
    !- statistics file in ASCII -!
    !----------------------------!
    WRITE(lun_stats,*)' '

    IF (key_unitm3) THEN
      WRITE(lun_stats,*)'total transport (in m3/s): ', trtot/(1+lmax-lmin), &
           ' ( x [1+lmax-lmin] =' ,trtot, ')'
      WRITE(lun_stats,*)'max_transport (in m3/s)  : ', max_transport
    ELSE
      WRITE(lun_stats,*)'total transport (in Sv)  : ', (trtot/(1+lmax-lmin))/1.e6_rprec , &
           ' ( x [1+lmax-lmin] =' ,trtot/1.e6_rprec, ')'
      WRITE(lun_stats,*)'max_transport (in Sv)    : ', max_transport/1.e6_rprec
    ENDIF
    WRITE(lun_stats,*)'# particles              : '  , ntraj


    IF ((TRIM(bin) == 'nobin').AND.(isn(-1) > 0)) THEN

      !------------------------------------------
      ! we mask uninitialized min. and max. ages
      !------------------------------------------
      DO i=1,nstat
        IF (ABS(s1min(i,-1)) > 1.e10_rprec) s1min(i,-1)=rZero
        IF (ABS(s1max(i,-1)) > 1.e10_rprec) s1max(i,-1)=rZero
        IF (ABS(s0min(i,-1)) > 1.e10_rprec) s0min(i,-1)=rZero
        IF (ABS(s0max(i,-1)) > 1.e10_rprec) s0max(i,-1)=rZero
      END DO
      !
      WRITE(lun_stats,*)' '
      WRITE(lun_stats,7057)'initial state ',' ','#',isn(-1)
      WRITE(lun_stats,7157)'stats. for: ',(chpstat(i),i=1,nstat)
      WRITE(lun_stats,7257)' min: ',(s0min(i,-1),i=1,nstat)
      WRITE(lun_stats,7257)' max: ',(s0max(i,-1),i=1,nstat)
      WRITE(lun_stats,7257)'mean: ',(s0x(i,-1)/sn(-1),i=1,nstat)
      WRITE(lun_stats,7257)'std. dev.: ', &
           ( SQRT(ABS((s0x2(i,-1)-s0x(i,-1)*s0x(i,-1)/sn(-1))/(sn(-1)-1._rprec))),i=1,nstat)

    ENDIF

    IF (TRIM(bin) /= 'nobin') THEN
      WRITE(lun_stats,*)' '
      WRITE(lun_stats,*)'no stats for initial state available since'
      WRITE(lun_stats,*)'   automatic positioning was by-passed'
    ENDIF

    !================================
    ! 2D flow projection and z stats
    !================================
    subtime = REAL(lmt, kind=rprec) / REAL(1+lmax-lmin, kind =rprec)

    !NG Performances (sub_quant_store_transp => * r_lmt)
    r_lmt =  1._rprec / REAL(lmt, kind=rprec)

    !!NG    WRITE(lun_xy_zonal)((subtime*uxy(i,j,1)*r_lmt,i=1,imt),j=1,jmt)
    !!NG    WRITE(lun_xy_merid)((subtime*vxy(i,j,1)*r_lmt,i=1,imt),j=1,jmt)
    !!NG    WRITE(lun_yz_merid)((subtime*vyz(j,k,1)*r_lmt,j=1,jmt),k=1,kmt)
    !!NG    WRITE(lun_yz_verti)((subtime*wyz(j,k,1)*r_lmt,j=1,jmt),k=1,kmt)
    !!NG    WRITE(lun_xz_zonal)((subtime*uxz(i,k,1)*r_lmt,i=1,imt),k=1,kmt)
    !!NG    WRITE(lun_xz_verti)((subtime*wxz(i,k,1)*r_lmt,i=1,imt),k=1,kmt)
    !!NG
    !!NG    IF (key_alltracers) THEN
    !!NG: seq-roms       WRITE(lun_yz_merid)((subtime*vyr(j,k)*r_lmt,j=1,jmt),k=1,nrclas)
    !!NG: seq-roms       WRITE(lun_yz_verti)((subtime*wyr(j,k)*r_lmt,j=1,jmt),k=1,nrclas)
    !!NG: seq-roms       WRITE(lun_xz_zonal)((subtime*uxr(i,k)*r_lmt,i=1,imt),k=1,nrclas)
    !!NG: seq-roms       WRITE(lun_xz_verti)((subtime*wxr(i,k)*r_lmt,i=1,imt),k=1,nrclas)
    !!NG: #endif
    !!NG    ENDIF
    !!NG
    !!NG    WRITE(lun_xy_stats)((uh(i,j,1)*r_lmt,i=1,imt),j=1,jmt)
    !!NG    WRITE(lun_xy_stats)((vh(i,j,1)*r_lmt,i=1,imt),j=1,jmt)
    !!NG    WRITE(lun_xy_stats)((zuh(i,j,1)*r_lmt,i=1,imt),j=1,jmt)
    !!NG    WRITE(lun_xy_stats)((zvh(i,j,1)*r_lmt,i=1,imt),j=1,jmt)
    !!NG    WRITE(lun_xy_stats)((z2uh(i,j,1)*r_lmt,i=1,imt),j=1,jmt)
    !!NG    WRITE(lun_xy_stats)((z2vh(i,j,1)*r_lmt,i=1,imt),j=1,jmt)
    !!NG
    !!NG    IF (key_alltracers) THEN
    !!NG      WRITE(lun_xy_stats)((tuh(i,j,1)*r_lmt,i=1,imt),j=1,jmt)
    !!NG      WRITE(lun_xy_stats)((tvh(i,j,1)*r_lmt,i=1,imt),j=1,jmt)
    !!NG      WRITE(lun_xy_stats)((t2uh(i,j,1)*r_lmt,i=1,imt),j=1,jmt)
    !!NG      WRITE(lun_xy_stats)((t2vh(i,j,1)*r_lmt,i=1,imt),j=1,jmt)
    !!NG      WRITE(lun_xy_stats)((suh(i,j,1)*r_lmt,i=1,imt),j=1,jmt)
    !!NG      WRITE(lun_xy_stats)((svh(i,j,1)*r_lmt,i=1,imt),j=1,jmt)
    !!NG      WRITE(lun_xy_stats)((s2uh(i,j,1)*r_lmt,i=1,imt),j=1,jmt)
    !!NG      WRITE(lun_xy_stats)((s2vh(i,j,1)*r_lmt,i=1,imt),j=1,jmt)
    !!NG      WRITE(lun_xy_stats)((ruh(i,j,1)*r_lmt,i=1,imt),j=1,jmt)
    !!NG      WRITE(lun_xy_stats)((rvh(i,j,1)*r_lmt,i=1,imt),j=1,jmt)
    !!NG      WRITE(lun_xy_stats)((r2uh(i,j,1)*r_lmt,i=1,imt),j=1,jmt)
    !!NG      WRITE(lun_xy_stats)((r2vh(i,j,1)*r_lmt,i=1,imt),j=1,jmt)
    !!NG
    !!NG    ENDIF

    uxy(:,:,:)  = uxy(:,:,:) * subtime * r_lmt
    vxy(:,:,:)  = vxy(:,:,:) * subtime * r_lmt
    vyz(:,:,:)  = vyz(:,:,:) * subtime * r_lmt
    wyz(:,:,:)  = wyz(:,:,:) * subtime * r_lmt
    uxz(:,:,:)  = uxz(:,:,:) * subtime * r_lmt
    wxz(:,:,:)  = wxz(:,:,:) * subtime * r_lmt
    uh(:,:,:)   = uh(:,:,:)            * r_lmt
    vh(:,:,:)   = vh(:,:,:)            * r_lmt
    zuh(:,:,:)  = zuh(:,:,:)           * r_lmt
    zvh(:,:,:)  = zvh(:,:,:)           * r_lmt
    z2uh(:,:,:) = z2uh(:,:,:)          * r_lmt
    z2vh(:,:,:) = z2vh(:,:,:)          * r_lmt

    IF (key_alltracers) THEN
      tuh(:,:,:)  = tuh(:,:,:)         * r_lmt
      tvh(:,:,:)  = tvh(:,:,:)         * r_lmt
      t2uh(:,:,:) = t2uh(:,:,:)        * r_lmt
      t2vh(:,:,:) = t2vh(:,:,:)        * r_lmt
      suh(:,:,:)  = suh(:,:,:)         * r_lmt
      svh(:,:,:)  = svh(:,:,:)         * r_lmt
      s2uh(:,:,:) = s2uh(:,:,:)        * r_lmt
      s2vh(:,:,:) = s2vh(:,:,:)        * r_lmt
      ruh(:,:,:)  = ruh(:,:,:)         * r_lmt
      rvh(:,:,:)  = rvh(:,:,:)         * r_lmt
      r2uh(:,:,:) = r2uh(:,:,:)        * r_lmt
      r2vh(:,:,:) = r2vh(:,:,:)        * r_lmt
    ENDIF

    !!---------------------------!
    !!= Quantitative statistics =!
    !!---------------------------!
    trtot1=rZero

    WRITE(lun_stats,*)' '
    nsect2 = nsect
    IF (icrit1 == 1) THEN
      nsect2 = 1 + nsect
    ENDIF

    DO n = 0, nsect2
      IF (key_unitm3) THEN
        WRITE(lun_stats,7357)secname(n),subtime*sectrans(n)*r_lmt,n
      ELSE
        WRITE(lun_stats,7357)secname(n),subtime*sectrans(n)/1.e6_rprec*r_lmt,n
      ENDIF
      trtot1=trtot1+sectrans(n)
    END DO

    IF (key_unitm3) THEN
      WRITE(lun_stats,7357)'total',subtime*trtot*r_lmt
      WRITE(lun_stats,7357)'except mnds', &
           subtime*(trtot-sectrans(0))*r_lmt
      WRITE(lun_stats,7357)'lost',subtime*(trtot-trtot1)*r_lmt
    ELSE
      WRITE(lun_stats,7357)'total',subtime*trtot/1.e6_rprec*r_lmt
      WRITE(lun_stats,7357)'except mnds', &
           subtime*(trtot-sectrans(0))/1.e6_rprec*r_lmt
      WRITE(lun_stats,7357)'lost',subtime*(trtot-trtot1)/1.e6_rprec*r_lmt
    ENDIF

7357 FORMAT(a15,1x,f0.4,1x,i0)

    DO n=0,nsect2
      IF (isn(n) > 0) THEN
        WRITE(lun_stats,*)' '
        WRITE(lun_stats,7057)'final state ',secname(n),' #',isn(n)
        WRITE(lun_stats,7157)'stats. ini: ',(chpstat(i),i=1,nstat)
        WRITE(lun_stats,7257)' min: ',(s0min(i,n),i=1,nstat)
        WRITE(lun_stats,7257)' max: ',(s0max(i,n),i=1,nstat)
        WRITE(lun_stats,7257)'mean: ',(s0x(i,n)/sn(n),i=1,nstat)
        WRITE(lun_stats,7257)'std. dev.: ',( SQRT( &
             ABS((s0x2(i,n)-s0x(i,n)*s0x(i,n)/sn(n))/(sn(n)-1._rprec))),i=1,nstat)
        WRITE(lun_stats,7157)'stats. fin: ',(chpstat(i),i=1,nstat)
        WRITE(lun_stats,7257)' min: ',(s1min(i,n),i=1,nstat)
        WRITE(lun_stats,7257)' max: ',(s1max(i,n),i=1,nstat)
        WRITE(lun_stats,7257)'mean: ',(s1x(i,n)/sn(n),i=1,nstat)
        WRITE(lun_stats,7257)'std. dev.: ',( SQRT( &
             ABS((s1x2(i,n)-s1x(i,n)*s1x(i,n)/sn(n))/(sn(n)-1._rprec))),i=1,nstat)
      ENDIF
    END DO

    !----------------------------------------------------------------------!
    !- Define statistic variables and their attributes in quantative mode -!
    !----------------------------------------------------------------------!
    CALL sub_save_netcdf_data_init_stats()

    !------------------------------------!
    !- Save statistics in Netcdf Format -!
    !------------------------------------!
    Call sub_save_netcdf_stats()

    !----------------------------!
    !- Close Netcdf output file -!
    !----------------------------!
    WRITE(lun_standard,*)''
    WRITE(lun_standard,*)'--------------------------------'
    WRITE(lun_standard,*)'= Close Statistics NetCDF File ='
    WRITE(lun_standard,*)'--------------------------------'
    CALL sub_save_netcdf_data_close(ncid_data_stats)

  ENDIF

7057 FORMAT(a13,a15,a2,i7)
7157 FORMAT(a13,7a10)
7257 FORMAT(a13,7f10.3)


  !=====================================================!
  !- OLD binary storage of initial and final positions -!
  !=====================================================!
  IF (key_ascii_outputs) THEN
    DO n=1,ntraj
      ! BBL __________________________________________________________________________
      ! BBL BBL: Il est temps de sortir sur fichier les positions initiales et
      ! BBL        finales de toutes les particules
      ! BBL
      WRITE(lun_final)hi(n),hj(n),-hk(n),hl(n),ttr(n),agenew(n)
      WRITE(lun_init)tfi(n),tfj(n),tfk(n),tfl(n),ttr(n), tage(n)
    ENDDO
  ENDIF

  !! NG: 15_09_2008
  IF (TRIM(mode) == 'quantitative') THEN 
    !-----------------!
    !- Initial Depth -!
    !-----------------!
    !- Netcdf -!
    !----------!

    CALL sub_save_netcdf_init_depth(init_depth(:))

  ENDIF

  !-------------------!
  !- Final positions -!
  !-------------------! 
  !- Netcdf -!
  !----------!
  !! NG: 15_09_2008
  IF (key_alltracers) THEN
    CALL sub_save_netcdf_final_pos( &
         hi(:),hj(:),-hk(:),hl(:) , &
         agenew(:), nfin(:)       , &
         final_depth(:)           , &
         tf(:),sf(:),rf(:)          )
  ELSE
    CALL sub_save_netcdf_final_pos(                           &
         hi(1:ntraj),hj(1:ntraj),-hk(1:ntraj),hl(1:ntraj)   , &
         agenew(1:ntraj), nfin(1:ntraj),final_depth(1:ntraj)  )
  ENDIF


  !----------------------------!
  !- Close Netcdf output file -!
  !----------------------------!
  WRITE(lun_standard,*)''
  WRITE(lun_standard,*)'-------------------------------'
  WRITE(lun_standard,*)'= Close Positions NetCDF File ='
  WRITE(lun_standard,*)'-------------------------------'
  CALL sub_save_netcdf_data_close(ncid_data_pos)

  !=====================!
  !- DEALLOCATE memory -!
  !=====================!

  WRITE(lun_standard,*)''
  WRITE(lun_standard,*)'---------------------'
  WRITE(lun_standard,*)'= Deallocate Memory ='
  WRITE(lun_standard,*)'---------------------'

  IF (ALLOCATED(igo)) THEN
    CALL sub_memory(-size(igo),'i','igo','trajec_seq')
    DEALLOCATE(igo)
  END IF
  IF (ALLOCATED(icolim)) THEN
    CALL sub_memory(-size(icolim),'i','icolim','trajec_seq')
    DEALLOCATE(icolim)
  END IF
  IF (ALLOCATED(i0)) THEN
    CALL sub_memory(-size(i0),'i','i0','trajec_seq')
    DEALLOCATE(i0)
  END IF
  IF (ALLOCATED(j0)) THEN
    CALL sub_memory(-size(j0),'i','j0','trajec_seq')
    DEALLOCATE(j0)
  END IF
  IF (ALLOCATED(k0)) THEN
    CALL sub_memory(-size(k0),'i','k0','trajec_seq')
    DEALLOCATE(k0)
  END IF
  IF (ALLOCATED(l0)) THEN
    CALL sub_memory(-size(l0),'i','l0','trajec_seq')
    DEALLOCATE(l0)
  END IF
  IF (ALLOCATED(i1)) THEN
    CALL sub_memory(-size(i1),'i','i1','trajec_seq')
    DEALLOCATE(i1)
  END IF
  IF (ALLOCATED(i2)) THEN
    CALL sub_memory(-size(i2),'i','i2','trajec_seq')
    DEALLOCATE(i2)
  END IF
  IF (ALLOCATED(j1)) THEN
    CALL sub_memory(-size(j1),'i','j1','trajec_seq')
    DEALLOCATE(j1)
  END IF
  IF (ALLOCATED(j2)) THEN
    CALL sub_memory(-size(j2),'i','j2','trajec_seq')
    DEALLOCATE(j2)
  END IF
  IF (ALLOCATED(k1)) THEN
    CALL sub_memory(-size(k1),'i','k1','trajec_seq')
    DEALLOCATE(k1)
  END IF
  IF (ALLOCATED(k2)) THEN
    CALL sub_memory(-size(k2),'i','k2','trajec_seq')
    DEALLOCATE(k2)
  END IF
  IF (ALLOCATED(iflag)) THEN
    CALL sub_memory(-size(iflag),'i','iflag','trajec_seq')
    DEALLOCATE(iflag)
  END IF
  IF (ALLOCATED(iperfl)) THEN
    CALL sub_memory(-size(iperfl),'i','iperfl','trajec_seq')
    DEALLOCATE(iperfl)
  END IF
  IF (ALLOCATED(jperfl)) THEN
    CALL sub_memory(-size(jperfl),'i','jperfl','trajec_seq')
    DEALLOCATE(jperfl)
  END IF
  IF (ALLOCATED(nfin)) THEN
    CALL sub_memory(-size(nfin),'i','nfin','trajec_seq')
    DEALLOCATE(nfin)
  END IF

  IF (ALLOCATED(age)) THEN
    CALL sub_memory(-size(age),'r','age','trajec_seq')
    DEALLOCATE(age)
  END IF
  IF (ALLOCATED(agenew)) THEN
    CALL sub_memory(-size(agenew),'r','agenew','trajec_seq')
    DEALLOCATE(agenew)
  END IF
  IF (ALLOCATED(truage)) THEN
    CALL sub_memory(-size(truage),'r','truage','trajec_seq')
    DEALLOCATE(truage)

  END IF
  IF (ALLOCATED(gi)) THEN
    CALL sub_memory(-size(gi),'r','gi','trajec_seq')
    DEALLOCATE(gi)
  END IF
  IF (ALLOCATED(gj)) THEN
    CALL sub_memory(-size(gj),'r','gj','trajec_seq')
    DEALLOCATE(gj)
  END IF
  IF (ALLOCATED(gk)) THEN
    CALL sub_memory(-size(gk),'r','gk','trajec_seq')
    DEALLOCATE(gk)
  END IF
  IF (ALLOCATED(gl)) THEN
    CALL sub_memory(-size(gl),'r','gl','trajec_seq')
    DEALLOCATE(gl)
  END IF
  IF (ALLOCATED(hi)) THEN
    CALL sub_memory(-size(hi),'r','hi','trajec_seq')
    DEALLOCATE(hi)
  END IF
  IF (ALLOCATED(hj)) THEN
    CALL sub_memory(-size(hj),'r','hj','trajec_seq')
    DEALLOCATE(hj)
  END IF
  IF (ALLOCATED(hk)) THEN
    CALL sub_memory(-size(hk),'r','hk','trajec_seq')
    DEALLOCATE(hk)
  END IF
  IF (ALLOCATED(hl)) THEN
    CALL sub_memory(-size(hl),'r','hl','trajec_seq')
    DEALLOCATE(hl)
  END IF
  IF (ALLOCATED(inext)) THEN
    CALL sub_memory(-size(inext),'i','inext','trajec_seq')
    DEALLOCATE(inext)
  END IF
  IF (ALLOCATED(tnext)) THEN
    CALL sub_memory(-size(tnext),'r','tnext','trajec_seq')
    DEALLOCATE(tnext)
  END IF
  IF (ALLOCATED(tswap)) THEN
    CALL sub_memory(-size(tswap),'r','tswap','trajec_seq')
    DEALLOCATE(tswap)
  END IF
  IF (ALLOCATED(zbuoy)) THEN
    CALL sub_memory(-size(zbuoy),'r','zbuoy','trajec_seq')
    DEALLOCATE(zbuoy)
  END IF
  IF (ALLOCATED(u)) THEN
    CALL sub_memory(-size(u),'r','u','trajec_seq')
    DEALLOCATE(u)
  END IF
  IF (ALLOCATED(v)) THEN
    CALL sub_memory(-size(v),'r','v','trajec_seq')
    DEALLOCATE(v)
  END IF
  IF (ALLOCATED(w)) THEN
    CALL sub_memory(-size(w),'r','w','trajec_seq')
    DEALLOCATE(w)
  END IF
  IF (ALLOCATED(du)) THEN
    CALL sub_memory(-size(du),'r','du','trajec_seq')
    DEALLOCATE(du)
  END IF
  IF (ALLOCATED(dv)) THEN
    CALL sub_memory(-size(dv),'r','dv','trajec_seq')
    DEALLOCATE(dv)
  END IF
  IF (ALLOCATED(dw)) THEN
    CALL sub_memory(-size(dw),'r','dw','trajec_seq')
    DEALLOCATE(dw)
  END IF
  IF (ALLOCATED(tx)) THEN
    CALL sub_memory(-size(tx),'r','tx','trajec_seq')
    DEALLOCATE(tx)
  END IF
  IF (ALLOCATED(ty)) THEN
    CALL sub_memory(-size(ty),'r','ty','trajec_seq')
    DEALLOCATE(ty)
  END IF
  IF (ALLOCATED(tz)) THEN
    CALL sub_memory(-size(tz),'r','tz','trajec_seq')
    DEALLOCATE(tz)
  END IF
  IF (ALLOCATED(tf)) THEN
    CALL sub_memory(-size(tf),'r','tf','trajec_seq')
    DEALLOCATE(tf)
  END IF
  IF (ALLOCATED(sf)) THEN
    CALL sub_memory(-size(sf),'r','sf','trajec_seq')
    DEALLOCATE(sf)
  END IF
  IF (ALLOCATED(rf)) THEN
    CALL sub_memory(-size(rf),'r','rf','trajec_seq')
    DEALLOCATE(rf)
  END IF
  IF (ALLOCATED(tfac)) THEN
    CALL sub_memory(-size(tfac),'r','tfac','trajec_seq')
    DEALLOCATE(tfac)
  END IF
  IF (ALLOCATED(ttt)) THEN
    CALL sub_memory(-size(ttt),'r','ttt','trajec_seq')
    DEALLOCATE(ttt)
  END IF

  !! NG: 15_09_2008
  IF (ALLOCATED(init_depth)) THEN
    CALL sub_memory(-size(init_depth),'r','init_depth','trajec_seq')
    DEALLOCATE(init_depth)
  END IF
  IF (ALLOCATED(final_depth)) THEN
    CALL sub_memory(-size(final_depth),'r','final_depth','trajec_seq')
    DEALLOCATE(final_depth)
  END IF

  IF (ALLOCATED(ti)) THEN
    CALL sub_memory(-size(ti),'r','ti','trajec_seq')
    DEALLOCATE(ti)
  END IF
  IF (ALLOCATED(si)) THEN
    CALL sub_memory(-size(si),'r','si','trajec_seq')
    DEALLOCATE(si)
  END IF
  IF (ALLOCATED(ri)) THEN
    CALL sub_memory(-size(ri),'r','ri','trajec_seq')
    DEALLOCATE(ri)
  END IF

  IF (small_virtual_memory) THEN
    IF (ALLOCATED(array_qual_traj_2D)) THEN
      CALL sub_memory(-size(array_qual_traj_2D),'r','array_qual_traj_2D','trajec_seq')
      DEALLOCATE(array_qual_traj_2D)
    END IF
  ELSE
    IF (ALLOCATED(ind_traj)) THEN
      CALL sub_memory(-size(ind_traj),'i','ind_traj','trajec_seq')
      DEALLOCATE(ind_traj)
    END IF
    IF (ALLOCATED(array_qual_traj)) THEN
      CALL sub_memory(-size(array_qual_traj),'r','array_qual_traj','trajec_seq')
      DEALLOCATE(array_qual_traj)
    END IF
  END IF

  CALL sub_netcdf_dealloc_mem()
  CALL sub_input_data_dealloc_mem()
  CALL sub_coord_dealloc()
  CALL sub_scalef_dealloc()
  IF (key_roms)  CALL sub_h_roms_dealloc()
  CALL sub_transp_dealloc()
  IF (key_alltracers) THEN
    CALL sub_tracer_dealloc()
  ENDIF
  if (key_vvl) then
     call sub_ssh_dealloc()
  endif
  CALL sub_posin_dealloc()
  CALL sub_tmask_dealloc()
  CALL sub_seq_dealloc()

  IF (TRIM(mode) == 'quantitative') THEN ! [QUANT]
    CALL sub_txt_dealloc()
    CALL sub_flags_dealloc()
    CALL sub_stats_dealloc()
    CALL sub_stati_dealloc()
    CALL sub_quant_dealloc()
  ENDIF

END SUBROUTINE trajec_seq
