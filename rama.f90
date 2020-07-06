!-----------------------------------------------------!
!          Racking And Monotonic Analysis             !
!                 of CFS wall panlels                 !
!_____________________________________________________!

PROGRAM RAMA
    USE Hanuman
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(3)

 INTEGER::fixed_freedoms,i,iel,k,loaded_nodes,ndim,ndofb,ndofp,nelsb,nelsp,neq,nodb=2,nodp, &
   nodofp=2,nodofb=3,nn,nnp,nnb,nprops=1,bnp_types,pnp_types,nr,nlen,count=0,j,eq_num=1,load_count=1,strt

 REAL(iwp)::axial,penalty=1.0e20_iwp,zero=0.0_iwp
!-------------------Dynamic Arrays---------------!
 INTEGER,ALLOCATABLE::betype(:),pe_typr(:),e_st_vb(:),e_st_vp(:),g_e_st_mb(:,:),g_e_st_mp(:,:),g_nd_numb(:,:),g_nd_nump(:,:),eqnb(:,:),eqnp(:,:), &
   node(:),e_nd_numb(:),e_nd_nump(:),sense(:),bndry(:,:),res_eq(:),load(:,:),load_eq(:)

 REAL(iwp),ALLOCATABLE::actionb(:),actionp(:),e_coordb(:,:),e_coordp(:,:),e_disp(:),e_disb(:),g_coordb(:,:),g_coordp(:,:),keleb(:,:),kelep(:,:),  &
   kstruct(:,:),loads(:),pprop(:,:),x_coords(:),y_coords(:),points(:,:),bprop(:,:),&
   value(:),kv(:,:),ksub(:,:),load_val(:),full_load(:),load_sub(:),shp_fun(:),jac(:,:), &
   loc_sp_der(:,:),g_sp_der(:,:),str_disp_m(:,:)

 CHARACTER(LEN=15)::filename
!_________________________________________________!
!                                                 !
!-------------PROGRAM STARTS HERE-----------------!
!_________________________________________________!

 PRINT*,"ENTER THE INPUT FILENAME"
 READ*,filename
 OPEN(UNIT=10,FILE=filename)
 OPEN(UNIT=11,FILE='RESULT.dat',STATUS='NEW')

 !------------------------------------TOTAL MODEL DATA-----------------------------!

 READ(10,*)nelsb,nnb,ndim,bnp_types               !beam
 nodofb=ndim
 ndofb=nodb*nodofb
 print*,ndofb
! READ(10,*)element,nodp,dir,nxe,nye,nip,pnp_types  !plate
! CALL mesh_size(element,nodp,nelsp,nnp,nxe,nye)
! ndofp=nodp*nodofp

 !------------------------------------BEAM----------------------------------------------!

 ALLOCATE(eqnb(nodofb,nnb),bndry(nodofb,nnb),load(nodofb,nnb),keleb(ndofb,ndofb),e_coordb(nodb,ndim),g_coordb(ndim,nnb),    &
   e_disb(ndofb),actionb(ndofb),g_nd_numb(nodb,nelsb),e_nd_numb(nodb),e_st_vb(ndofb),g_e_st_mb(ndofb,nelsb),&
   betype(nelsb),bprop(nprops,bnp_types))

 !------------------------------------------------------------------------------------------
 READ(10,*)bprop
 betype=1
 IF(bnp_types>1)READ(10,*)betype
 READ(10,*)g_coordb
 READ(10,*)g_nd_numb
 print*,g_nd_numb 
 eqnb=1
 strt=1
 CALL num_eq(eqnb,strt) 
 neq=MAXVAL(eqnb)
print*,neq
 !---------------------------------------------------PLATE-----------------------------------------

 !allocate(neqp(nodofp,nnp))
 !eqnp=1
 !strt=eqnb
 !call num_eq(eqnp,strt)
 !neq=maxval(eqnp)

!--------------------------------------------------------------
!------------------------------GLOBAL--------------------------

 ALLOCATE(kstruct(neq,neq),kv(neq,neq),full_load(neq))


!________________________________________________________________________!
!                                                                        !
!                          BEAM COMPUTATIONS                             !
!________________________________________________________________________!

!-----------------------Loop over elements to form global steering matrix--------------!
 elements_1: DO iel=1,nelsb
   e_nd_numb=g_nd_numb(:,iel)
   CALL end_to_est(e_nd_numb,eqnb,e_st_vb)
   g_e_st_mb(:,iel)=e_st_vb   
  END DO elements_1
print*,"reached here"
print*,g_e_st_mb
print*,e_st_vb
!-----------------------global stiffness matrix assembly------------------
  kstruct=0
  elements_2: DO iel=1,nelsb
    e_nd_numb=g_nd_numb(:,iel)
    e_coordb=TRANSPOSE(g_coordb(:,e_nd_numb))
    CALL truss(keleb,bprop(1,betype(iel)),e_coordb)
    kv=0
    print*,"reached2"
    e_st_vb=g_e_st_mb(:,iel)
    CALL assemble(kv,keleb,e_st_vb)
    kstruct=kv+kstruct
  END DO elements_2
  write(unit=11,fmt=245)kstruct
  245 format(12f12.3)
 ! write(11,*)" "





END PROGRAM RAMA

