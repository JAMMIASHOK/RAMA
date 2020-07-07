!-----------------------------------------------------!
!          Racking And Monotonic Analysis             !
!                 of CFS wall panlels                 !
!_____________________________________________________!

PROGRAM RAMA
    USE Hanuman
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(3)

 INTEGER::fixed_freedoms,i,iel,k,loaded_nodes,ndim,ndofb,ndofp,nelsb,nelsp,neq,nodb=2,nodp,nxe,nye,nip,nels, &
   nodofp=2,nodofb=3,nn,nnp,nnb,nprops=1,npprops=2,bnp_types,pnp_types,nr,nlen,count=0,j,eq_num=1,load_count=1,strt,nst=3

 REAL(iwp)::axial,penalty=1.0e20_iwp,zero=0.0_iwp,det
!-------------------Dynamic Arrays---------------!
 INTEGER,ALLOCATABLE::betype(:),petype(:),e_st_vb(:),e_st_vp(:),g_e_st_mb(:,:),g_e_st_mp(:,:),g_nd_numb(:,:),&
  g_nd_nump(:,:),eqnb(:,:),eqnp(:,:), &
   node(:),e_nd_numb(:),e_nd_nump(:),sense(:),bndry(:,:),res_eq(:),load(:,:),load_eq(:)

 REAL(iwp),ALLOCATABLE::actionb(:),actionp(:),e_coordb(:,:),e_coordp(:,:),e_disp(:),e_disb(:),g_coordb(:,:),&
  g_coordp(:,:),keleb(:,:),kelep(:,:),  &
   kstruct(:,:),loads(:),pprop(:,:),x_coords(:),y_coords(:),points(:,:),bprop(:,:),&
   value(:),kv(:,:),ksub(:,:),load_val(:),full_load(:),load_sub(:),shp_fun(:),jac(:,:), &
   loc_sp_der(:,:),sp_der(:,:),bee(:,:),constm(:,:),weights(:)

 CHARACTER(LEN=15)::filename,dir,element
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
 
 READ(10,*)element,nodp,dir,nxe,nye,nip,pnp_types  !plate
 CALL mesh_size(element,nodp,nelsp,nnp,nxe,nye)
 ndofp=nodp*nodofp

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
 
 eqnb=1
 strt=0
 CALL num_eq(eqnb,strt) 
 neq=MAXVAL(eqnb)
 
 !---------------------------------------------------PLATE-----------------------------------------

 allocate(eqnp(nodofp,nnp))
 eqnp=1
 strt=neq+1
 call num_eq(eqnp,strt)
 neq=maxval(eqnp)

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
!-----------------------global stiffness matrix assembly------------------
  kstruct=0
  elements_2: DO iel=1,nelsb
    e_nd_numb=g_nd_numb(:,iel)
    e_coordb=TRANSPOSE(g_coordb(:,e_nd_numb))
    CALL beam(keleb,bprop,betype,iel,e_coordb)
    kv=0
    e_st_vb=g_e_st_mb(:,iel)
    CALL assemble(kv,keleb,e_st_vb)
    kstruct=kv+kstruct
  END DO elements_2
  
 ! write(11,*)" "
!_______------------------_______________--------- Till here fine-------_____________-----------------______________!

  !___________________________________________________________________________!
!                                                                           !
!                           PLATE COMPUTATIONS                              !
!___________________________________________________________________________!

allocate(kelep(ndofp,ndofp),e_coordp(nodp,ndim),g_coordp(ndim,nnp),    &
e_disp(ndofp),actionp(ndofp),g_nd_nump(nodp,nelsp),e_nd_nump(nodp),e_st_vp(ndofp),g_e_st_mp(ndofp,nelsp),&
petype(nelsp),pprop(npprops,pnp_types),x_coords(nxe+1),y_coords(nye+1),points(nip,ndim),shp_fun(nodp),jac(ndim,ndim),&
loc_sp_der(ndim,nodp),sp_der(ndim,nodp),bee(nst,ndofp),constm(nst,nst),weights(nip))
READ(10,*)pprop
print*,pprop
 petype=1
 IF(pnp_types>1)read(10,*)petype
 READ(10,*)x_coords,y_coords

!-----------------------loop the elements to find global steering matrix-----
 
 elements_3: DO iel=1,nelsp
   CALL geom_rect(element,iel,x_coords,y_coords,e_coordp,e_nd_nump,dir)
   CALL end_to_est(e_nd_nump,eqnp,e_st_vp)
   g_nd_nump(:,iel)=e_nd_nump
   g_coordp(:,e_nd_nump)=TRANSPOSE(e_coordp)
   g_e_st_mp(:,iel)=e_st_vp
 END DO elements_3
 
 !-----------------------element stiffness integration and assembly--------
 CALL sample(points,weights)
 
 print*,"reached here"
 elements_4: DO iel=1,nelsp
  
   CALL cmat(constm,pprop(1,petype(iel)),pprop(2,petype(iel)))
   
   e_nd_nump=g_nd_nump(:,iel)
   e_st_vp=g_e_st_mp(:,iel)
   e_coordp=TRANSPOSE(g_coordp(:,e_nd_nump))
   
   kv=0
   int_pts_1: DO i=1,nip
     CALL shape_fun(shp_fun,points,i)
     CALL shape_der(loc_sp_der,points,i)
     jac=MATMUL(loc_sp_der,e_coordp)
     det=determinant(jac)
     CALL invert(jac)
     sp_der=MATMUL(jac,loc_sp_der)
     CALL beemat(bee,sp_der)
     kelep=kelep+MATMUL(MATMUL(TRANSPOSE(bee),constm),bee)*det*weights(i)
   END DO int_pts_1
   CALL assemble(kv,kelep,e_st_vp)
   
   kstruct=kv+kstruct
  END DO elements_4
  write(unit=11,fmt=245)kstruct
  245 format(12f12.3)
!-----------------------------------------------Compiles well--------------------------




END PROGRAM RAMA

