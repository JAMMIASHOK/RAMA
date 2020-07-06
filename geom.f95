MODULE geom
CONTAINS

SUBROUTINE emb_2d_bc(nx1,nx2,ny1,ny2,nf)
!
! This subroutine generates the nf array for a 2-d slope geometry.
!
 IMPLICIT NONE
 INTEGER,INTENT(IN)::nx1,nx2,ny1,ny2
 INTEGER,INTENT(OUT)::nf(:,:)
 INTEGER::nm,ic,i,j,nye
 nye=ny1+ny2
 nm=0
 ic=0
 DO i=1,2*nye
   nm=nm+1
   nf(1,nm)=0
   ic=ic+1
   nf(2,nm)=ic
 END DO
 nm=nm+1
 nf(1,nm)=0
 nf(2,nm)=0
 DO j=1,nx1-1
   DO i=1,nye
     nm=nm+1
     ic=ic+1
     nf(1,nm)=ic
     ic=ic+1
     nf(2,nm)=ic
   END DO
   nm=nm+1
   nf(1,nm)=0
   nf(2,nm)=0
   DO i=1,2*nye
     nm=nm+1
     ic=ic+1
     nf(1,nm)=ic
     ic=ic+1
     nf(2,nm)=ic
   END DO
   nm=nm+1
   nf(1,nm)=0
   nf(2,nm)=0
 END DO
   DO i=1,nye
     nm=nm+1
     ic=ic+1
     nf(1,nm)=ic
     ic=ic+1
     nf(2,nm)=ic
   END DO
   nm=nm+1
   nf(1,nm)=0
   nf(2,nm)=0
   DO i=1,2*ny1
     nm=nm+1
     ic=ic+1
     nf(1,nm)=ic
     ic=ic+1
     nf(2,nm)=ic
   END DO 
   IF(nx2==0)THEN
   DO i=1,2*ny2
     nm=nm+1
     nf(1,nm)=0
     ic=ic+1
     nf(2,nm)=ic
   END DO
   nm=nm+1
   nf(1,nm)=0
   nf(2,nm)=0
 ELSE
   DO i=1,2*ny2
     nm=nm+1
     ic=ic+1
     nf(1,nm)=ic
     ic=ic+1
     nf(2,nm)=ic
   END DO
   nm=nm+1
   nf(1,nm)=0
   nf(2,nm)=0
   DO j=1,nx2
     DO i=1,ny2
       nm=nm+1
       ic=ic+1
       nf(1,nm)=ic
       ic=ic+1
       nf(2,nm)=ic
     END DO
     nm=nm+1
     nf(1,nm)=0
     nf(2,nm)=0
     DO i=1,2*ny2
       nm=nm+1
       IF(j==nx2)THEN
         nf(1,nm)=0
       ELSE
         ic=ic+1
         nf(1,nm)=ic
       END IF
       ic=ic+1
       nf(2,nm)=ic
     END DO
     nm=nm+1
     nf(1,nm)=0
     nf(2,nm)=0
   END DO
 END IF
RETURN
END SUBROUTINE emb_2d_bc






















SUBROUTINE emb_2d_geom(iel,nx1,nx2,ny1,ny2,w1,s1,w2,h1,h2,coord,num)
!
! This subroutine forms the nodal coordinates and numbering for a 2-d
! slope of 8-node quadrilaterals. Nodes numbering in the y-direction,
! elements numbered in the x-direction.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::w1,s1,w2,h1,h2
 INTEGER,INTENT(IN)::iel,nx1,nx2,ny1,ny2
 REAL(iwp),INTENT(OUT)::coord(:,:)
 INTEGER,INTENT(OUT)::num(:)
 REAL(iwp)::facx,facy,facb,facs,frh,zero=0.0_iwp,pt5=0.5_iwp,one=1.0_iwp
 INTEGER::nxe,nye,nc,nt,ip,iq
 nxe=nx1+nx2
 nye=ny1+ny2
 nt=nx1*ny1
 nc=(3*nye+2)*nx1+2*ny1
 facx=s1/ny1
 facy=h1/ny1
 facb=(w1+s1)/nx1
 facs=zero
 IF(ny2/=0)facs=h2/ny2
 frh=zero
 IF(nx2/=0)frh=w2/nx2
 IF(iel<=nt)THEN
   iq=(iel-1)/nx1+1
   ip=iel-(iq-1)*nx1
 ELSE
   iq=(iel-nt-1)/nxe+ny1+1
   ip=iel-nt-(iq-ny1-1)*nxe
 END IF
 IF(ip<=nx1)THEN
   num(1)=(ip-1)*(3*nye+2)+2*iq+1
   num(2)=num(1)-1
   num(3)=num(1)-2
   num(4)=(ip-1)*(3*nye+2)+2*nye+iq+1
   num(5)=ip*(3*nye+2)+2*iq-1
   num(6)=num(5)+1
   num(7)=num(5)+2
   num(8)=num(4)+1
   IF(iq<=ny1)THEN
     coord(1,1)=(ip-one)*(w1+iq*facx)/nx1
     coord(3,1)=(ip-one)*(w1+(iq-1)*facx)/nx1
     coord(5,1)=ip*(w1+(iq-1)*facx)/nx1
     coord(7,1)=ip*(w1+iq*facx)/nx1
     coord(1,2)=-iq*facy
     coord(3,2)=-(iq-1)*facy
     coord(5,2)=-(iq-1)*facy
     coord(7,2)=-iq*facy
   ELSE
     coord(1,1)=(ip-one)*facb
     coord(3,1)=(ip-one)*facb
     coord(5,1)=ip*facb
     coord(7,1)=ip*facb
     coord(1,2)=-h1-(iq-ny1)*facs
     coord(3,2)=-h1-(iq-ny1-1)*facs
     coord(5,2)=-h1-(iq-ny1-1)*facs
     coord(7,2)=-h1-(iq-ny1)*facs
   END IF
 ELSE
   num(1)=nc+(ip-nx1-1)*(3*ny2+2)+2*(iq-ny1)+1
   num(2)=num(1)-1
   num(3)=num(1)-2
   num(4)=nc+(ip-nx1-1)*(3*ny2+2)+2*ny2+iq-ny1+1
   num(5)=nc+(ip-nx1)*(3*ny2+2)+2*(iq-ny1)-1
   num(6)=num(5)+1
   num(7)=num(5)+2
   num(8)=num(4)+1
   coord(1,1)=w1+s1+(ip-nx1-1)*frh
   coord(3,1)=coord(1,1)
   coord(5,1)=w1+s1+(ip-nx1)*frh
   coord(7,1)=coord(5,1)
   coord(1,2)=-h1-(iq-ny1)*facs
   coord(3,2)=-h1-(iq-ny1-1)*facs
   coord(5,2)=-h1-(iq-ny1-1)*facs
   coord(7,2)=-h1-(iq-ny1)*facs
 END IF
 coord(2:6:2,:)=pt5*(coord(1:5:2,:)+coord(3:7:2,:))
 coord(8,:)=pt5*(coord(7,:)+coord(1,:))
RETURN
END SUBROUTINE emb_2d_geom
SUBROUTINE emb_3d_bc(ifix,nx1,nx2,ny1,ny2,nze,nf)
!
! This subroutine generates the nf array for a 3-d slope geometry.
! Side boundary conditions:
! smooth-smooth  ifix=1  (for checking against plane-strain)
! rough -smooth  ifix=2  (for centerline symmetry)
! rough -rough   ifix=3
!
 IMPLICIT NONE
 INTEGER,INTENT(IN)::nx1,nx2,ny1,ny2,nze,ifix
 INTEGER,INTENT(OUT)::nf(:,:)
 INTEGER::nm,ic,i,j,k,nye
 nye=ny1+ny2
 nm=0
 ic=0
!
!                       boundary condition for the front face
!
 DO i=1,2*nye
   nm=nm+1
   nf(:,nm)=0
   IF(ifix==1)THEN
     ic=ic+1
     nf(2,nm)=ic
   END IF
 END DO
 nm=nm+1
 nf(:,nm)=0
 DO j=1,nx1
   DO i=1,nye
     nm=nm+1
     nf(:,nm)=0
     IF(ifix==1)THEN
       ic=ic+1
       nf(1,nm)=ic
       ic=ic+1
       nf(2,nm)=ic
     END IF
   END DO
   nm=nm+1
   nf(:,nm)=0
   DO i=1,2*nye
     nm=nm+1
     nf(:,nm)=0
     IF(ifix==1)THEN
       ic=ic+1
       nf(1,nm)=ic
       ic=ic+1
       nf(2,nm)=ic
     END IF
   END DO
   nm=nm+1
   nf(:,nm)=0
 END DO
 DO j=1,nx2
   DO i=1,ny2
     nm=nm+1
     nf(:,nm)=0
     IF(ifix==1)THEN
       ic=ic+1
       nf(1,nm)=ic
       ic=ic+1
       nf(2,nm)=ic
     END IF
   END DO
   nm=nm+1
   nf(:,nm)=0
   DO i=1,2*ny2
     nm=nm+1
     nf(:,nm)=0
     IF(ifix==1)THEN
       IF(j<nx2)THEN
         ic=ic+1
         nf(1,nm)=ic
       END IF 
       ic=ic+1
       nf(2,nm)=ic
     END IF
   END DO
   nm=nm+1
   nf(:,nm)=0
 END DO    
!
!                       boundary conditions for the middle layer
!
 DO i=1,nye
   nm=nm+1
   nf(1,nm)=0
   ic=ic+1
   nf(3,nm)=ic
   ic=ic+1
   nf(2,nm)=ic
 END DO
 nm=nm+1
 nf(:,nm)=0
 DO j=1,nx1
   DO i=1,nye
     nm=nm+1
     ic=ic+1
     nf(1,nm)=ic
     ic=ic+1
     nf(3,nm)=ic
     ic=ic+1
     nf(2,nm)=ic
   END DO
   nm=nm+1
   nf(:,nm)=0
 END DO
 DO j=1,nx2
   DO i=1,ny2
     nm=nm+1
     IF(j==nx2)THEN
       nf(1,nm)=0
     ELSE
       ic=ic+1
       nf(1,nm)=ic
     END IF
     ic=ic+1
     nf(3,nm)=ic
     ic=ic+1
     nf(2,nm)=ic
   END DO
   nm=nm+1
   nf(:,nm)=0
 END DO
!all other layers
DO k=1,nze-1
  DO i=1,2*nye
    nm=nm+1
    nf(1,nm)=0
    ic=ic+1
    nf(3,nm)=ic
    ic=ic+1
    nf(2,nm)=ic
  END DO
  nm=nm+1
  nf(:,nm)=0
  DO j=1,nx1
    DO i=1,nye
      nm=nm+1
      ic=ic+1
      nf(1,nm)=ic
      ic=ic+1
      nf(3,nm)=ic
      ic=ic+1
      nf(2,nm)=ic
    END DO
    nm=nm+1
    nf(:,nm)=0
    DO i=1,2*nye
      nm=nm+1
      ic=ic+1
      nf(1,nm)=ic
      ic=ic+1
      nf(3,nm)=ic
      ic=ic+1
      nf(2,nm)=ic
    END DO
    nm=nm+1
    nf(:,nm)=0
  END DO
  DO j=1,nx2
    DO i=1,ny2
      nm=nm+1
      ic=ic+1
      nf(1,nm)=ic
      ic=ic+1
      nf(3,nm)=ic
      ic=ic+1
      nf(2,nm)=ic
    END DO
    nm=nm+1
    nf(:,nm)=0
    DO i=1,2*ny2
      nm=nm+1
      IF(j==nx2)THEN
        nf(1,nm)=0
      ELSE
        ic=ic+1
        nf(1,nm)=ic
      END IF
      ic=ic+1
      nf(3,nm)=ic
      ic=ic+1
      nf(2,nm)=ic
    END DO
    nm=nm+1
    nf(:,nm)=0
  END DO
! middle layers
  DO i=1,nye
    nm=nm+1
    nf(1,nm)=0
    ic=ic+1
    nf(3,nm)=ic
    ic=ic+1
    nf(2,nm)=ic
  END DO
  nm=nm+1
  nf(:,nm)=0
  DO j=1,nx1
    DO i=1,nye
      nm=nm+1
      ic=ic+1
      nf(1,nm)=ic
      ic=ic+1
      nf(3,nm)=ic
      ic=ic+1
      nf(2,nm)=ic
    END DO
    nm=nm+1
    nf(:,nm)=0
  END DO
  DO j=1,nx2
    DO i=1,ny2
      nm=nm+1
      IF(j==nx2)THEN
        nf(1,nm)=0
      ELSE
        ic=ic+1
        nf(1,nm)=ic
      END IF
      ic=ic+1
      nf(3,nm)=ic
      ic=ic+1
      nf(2,nm)=ic
    END DO
    nm=nm+1
    nf(:,nm)=0
  END DO
END DO
!
!                       boundary conditions for back face
!
 DO i=1,2*nye
   nm=nm+1
   nf(:,nm)=0
   IF(ifix==1.OR.ifix==2)THEN
     ic=ic+1
     nf(2,nm)=ic
   END IF
 END DO
 nm=nm+1
 nf(:,nm)=0
 DO j=1,nx1
   DO i=1,nye
     nm=nm+1
     nf(:,nm)=0
     IF(ifix==1.OR.ifix==2)THEN
       ic=ic+1
       nf(1,nm)=ic
       ic=ic+1
       nf(2,nm)=ic
     END IF
   END DO
   nm=nm+1
   nf(:,nm)=0
   DO i=1,2*nye
     nm=nm+1
     nf(:,nm)=0
     IF(ifix==1.OR.ifix==2)THEN
       ic=ic+1
       nf(1,nm)=ic
       ic=ic+1
       nf(2,nm)=ic
     END IF
   END DO
   nm=nm+1
   nf(:,nm)=0
 END DO
 DO j=1,nx2
   DO i=1,ny2
     nm=nm+1
     nf(:,nm)=0
     IF(ifix==1.OR.ifix==2)THEN
       ic=ic+1
       nf(1,nm)=ic
       ic=ic+1
       nf(2,nm)=ic
     END IF
   END DO
   nm=nm+1
   nf(:,nm)=0
   DO i=1,2*ny2
     nm=nm+1
     nf(:,nm)=0
     IF(ifix==1.OR.ifix==2)THEN
       IF(j<nx2)THEN
         ic=ic+1
         nf(1,nm)=ic
       END IF 
       ic=ic+1
       nf(2,nm)=ic
     END IF
   END DO
   nm=nm+1
   nf(:,nm)=0
 END DO    
RETURN
END SUBROUTINE emb_3d_bc
SUBROUTINE emb_3d_geom(iel,nx1,nx2,ny1,ny2,nze,w1,s1,w2,h1,h2,d1,coord,num)
!
! This subroutine forms the nodal coordinates and numbering for a 3-d
! slope of 20-node hexahedra. Nodes and elements numbered in xz planes
! going in the y-direction.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::w1,s1,w2,h1,h2, d1
 INTEGER,INTENT(IN)::iel,nx1,nx2,ny1,ny2,nze
 REAL(iwp),INTENT(OUT)::coord(:,:)
 INTEGER,intent(OUT)::num(:)
 INTEGER::nxe,nye,nc,nt,nn_mid,nn_main,ip,iq,is,nd,ns
 REAL(iwp)::facx,facy,facz,facb,facs,frh,zero=0.0_iwp,pt5=0.5_iwp
!
!                       basic variable initialization
!
 nxe=nx1+nx2
 nye=ny1+ny2
!number of elements in a slice
 ns=nx1*nye+nx2*ny2
!number of elements in sloped region
 nt=nx1*ny1
!number of nodes in main sloped region
 nc=(3*nye+2)*nx1+2*ny1
!number of nodes in mid sloped region
 nd=(nye+1)*nx1+ny1
!number of nodes in a middle layer
 nn_mid=(nye+1)*(nx1+1)+(ny2+1)*nx2
!number of nodes in a main layer
 nn_main=(3*nye+2)*nx1+2*nye+1+(3*ny2+2)*nx2
!
 facx=s1/ny1
 facy=h1/ny1
 facz=d1/nze
 facb=(w1+s1)/nx1
 facs=zero
 IF(ny2/=0)facs=h2/ny2
 frh=zero
 IF(nx2/=0)frh=w2/nx2
!
!                       determine element index values (ip, iq, is)
!
 is=(iel-1)/ns+1
 IF(iel-(is-1)*ns<=nt)THEN !if it is in sloped region
   iq=(iel-(is-1)*ns-1)/nx1+1
   ip=iel-(is-1)*ns-(iq-1)*nx1
 ELSE                      !if it is in unsloped region
   iq=(iel-(is-1)*ns-nt-1)/nxe+ny1+1
   ip=iel-(is-1)*ns-nt-(iq-ny1-1)*nxe
 END IF
!
! Calculate node indices (num()) and coordinates (coord(,))
 IF(ip<=nx1)THEN !if iel is in sloped region
   num(1)=(ip-1)*(3*nye+2)+2*iq+1+(is-1)*(nn_mid+nn_main)
   num(2)=num(1)-1
   num(3)=num(1)-2
   num(4)=(ip-1)*(3*nye+2)+2*nye+iq+1+(is-1)*(nn_mid+nn_main)
   num(5)=ip*(3*nye+2)+2*iq-1+(is-1)*(nn_mid+nn_main)
   num(6)=num(5)+1
   num(7)=num(5)+2
   num(8)=num(4)+1
   num(9)=(nn_main+(ip-1)*(nye+1)+iq+1)*is+(is-1)*                        &
     (nn_mid-(ip-1)*(nye+1)-iq-1)
   num(10)=num(9)-1
   num(11)=1+(nn_main+ip*(nye+1)+iq-1)*is+(is-1)*(nn_mid-ip*(nye+1)-iq+1)
   num(12)=num(11)+1
   num(13)=(ip-1)*(3*nye+2)+2*iq+1+is*(nn_mid+nn_main)
   num(14)=num(13)-1
   num(15)=num(13)-2
   num(16)=(ip-1)*(3*nye+2)+2*nye+iq+1+is*(nn_mid+nn_main)
   num(17)=ip*(3*nye+2)+2*iq-1+is*(nn_mid+nn_main)
   num(18)=num(17)+1
   num(19)=num(17)+2
   num(20)=num(16)+1
   IF(iq<=ny1)THEN ! if iel is in trapezoidal region
     coord(1,1)=(ip-1)*(w1+iq*facx)/nx1
     coord(3,1)=(ip-1)*(w1+(iq-1)*facx)/nx1
     coord(5,1)=ip*(w1+(iq-1)*facx)/nx1
     coord(7,1)=ip*(w1+iq*facx)/nx1
     coord(1,2)=-iq*facy
     coord(3,2)=-(iq-1)*facy
     coord(5,2)=-(iq-1)*facy
     coord(7,2)=-iq*facy
   ELSE            ! if iel is in rectangular region
     coord(1,1)=(ip-1)*facb
     coord(3,1)=(ip-1)*facb
     coord(5,1)=ip*facb
     coord(7,1)=ip*facb
     coord(1,2)=-h1-(iq-ny1)*facs
     coord(3,2)=-h1-(iq-ny1-1)*facs
     coord(5,2)=-h1-(iq-ny1-1)*facs
     coord(7,2)=-h1-(iq-ny1)*facs
   END IF
 ELSE              ! if iel is in unsloped region
   num(1)=(nn_mid+nn_main)*(is-1)+nc+(3*ny2+2)*(ip-nx1-1)+(iq-ny1)*2+1
   num(2)=num(1)-1
   num(3)=num(1)-2
   num(4)=(nn_mid+nn_main)*(is-1)+nc+(3*ny2+2)*(ip-nx1-1)+2*ny2+1+iq-ny1
   num(5)=(nn_mid+nn_main)*(is-1)+nc+(3*ny2+2)*(ip-nx1)+(iq-ny1)*2-1
   num(6)=num(5)+1
   num(7)=num(5)+2
   num(8)=num(4)+1
   num(9)=(nn_mid+nn_main)*(is-1)+nn_main+nd+(ny2+1)*(ip-nx1-1)+iq-ny1+1
   num(10)=num(9)-1
   num(11)=(nn_mid+nn_main)*(is-1)+nn_main+nd+(ny2+1)*(ip-nx1)+iq-ny1
   num(12)=num(11)+1
   num(13)=(nn_mid+nn_main)*is+nc+(3*ny2+2)*(ip-nx1-1)+(iq-ny1)*2+1
   num(14)=num(13)-1
   num(15)=num(13)-2
   num(16)=(nn_mid+nn_main)*is+nc+(3*ny2+2)*(ip-nx1-1)+2*ny2+1+iq-ny1
   num(17)=(nn_mid+nn_main)*is+nc+(3*ny2+2)*(ip-nx1)+(iq-ny1)*2-1
   num(18)=num(17)+1
   num(19)=num(17)+2
   num(20)=num(16)+1
!
   coord(1,1)=w1+s1+(ip-nx1-1)*frh
   coord(3,1)=coord(1,1)
   coord(5,1)=w1+s1+(ip-nx1)*frh
   coord(7,1)=coord(5,1)
   coord(1,2)=-h1-(iq-ny1)*facs
   coord(3,2)=-h1-(iq-ny1-1)*facs
   coord(5,2)=-h1-(iq-ny1-1)*facs
   coord(7,2)=-h1-(iq-ny1)*facs
 END IF
 coord(9,1)=coord(1,1)
 coord(10,1)=coord(3,1)
 coord(11,1)=coord(5,1)
 coord(12,1)=coord(7,1)
 coord(13,1)=coord(1,1)
 coord(15,1)=coord(3,1)
 coord(17,1)=coord(5,1)
 coord(19,1)=coord(7,1)
 coord(9,2) =coord(1,2)
 coord(10,2)=coord(3,2)
 coord(11,2)=coord(5,2)
 coord(12,2)=coord(7,2)
 coord(13,2)=coord(1,2)
 coord(15,2)=coord(3,2)
 coord(17,2)=coord(5,2)
 coord(19,2)=coord(7,2)
 coord(1:8:1,3)=-(is-1)*facz
 coord(9:12:1,3)=-(is*facz-facz*pt5)
 coord(13:20:1,3)=-is*facz    
 coord(2:6:2,:)=pt5*(coord(1:5:2,:)+coord(3:7:2,:))
 coord(14:18:2,:)=pt5*(coord(13:17:2,:)+coord(15:19:2,:))
 coord(8,:)=pt5*(coord(7,:)+coord(1,:))
 coord(20,:)=pt5*(coord(19,:)+coord(13,:))
RETURN
END SUBROUTINE emb_3d_geom






















SUBROUTINE fmcoem(g_num,g_coord,fwidth,fdepth,width,depth,                &
  lnxe,lifts,fnxe,fnye,itype)
!
! This subroutine returns g_coord for the mesh.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 INTEGER,INTENT(IN)::g_num(:,:),lnxe,lifts,fnxe,fnye,itype
 REAL(iwp),INTENT(IN)::fwidth(:),fdepth(:),width(:),depth(:)
 REAL(iwp),INTENT(OUT)::g_coord(:,:)
 INTEGER::ig,i,j,ii,lnye
 REAL(iwp)::pt5=0.5_iwp
 lnye=1
 ig=0
 DO i=1,fnye
   DO j=1,fnxe
     ig=ig+1
     g_coord(1,g_num(1,ig))=fwidth(j)
     g_coord(1,g_num(2,ig))=fwidth(j)
     g_coord(1,g_num(3,ig))=fwidth(j)
     g_coord(1,g_num(5,ig))=fwidth(j+1)
     g_coord(1,g_num(6,ig))=fwidth(j+1)
     g_coord(1,g_num(7,ig))=fwidth(j+1)
     g_coord(1,g_num(4,ig))=                                              &
       pt5*(g_coord(1,g_num(3,ig))+g_coord(1,g_num(5,ig)))
     g_coord(1,g_num(8,ig))=                                              &
       pt5*(g_coord(1,g_num(1,ig))+g_coord(1,g_num(7,ig)))
     g_coord(2,g_num(1,ig))=fdepth(i)
     g_coord(2,g_num(8,ig))=fdepth(i)
     g_coord(2,g_num(7,ig))=fdepth(i)
     g_coord(2,g_num(3,ig))=fdepth(i+1)
     g_coord(2,g_num(4,ig))=fdepth(i+1)
     g_coord(2,g_num(5,ig))=fdepth(i+1)
     g_coord(2,g_num(2,ig))=                                              &
       pt5*(g_coord(2,g_num(1,ig))+g_coord(2,g_num(3,ig)))
     g_coord(2,g_num(6,ig))=                                              &
       pt5*(g_coord(2,g_num(7,ig))+g_coord(2,g_num(5,ig)))
   END DO
 END DO
 IF(itype==1)THEN
   DO ii=1,lifts-1
     DO i=1,lnye
       DO j=1,lnxe-(ii-1)*1
         ig=ig+1
         IF(j==1)THEN
           g_coord(1,g_num(1,ig))=width((ii-1)+j)
           g_coord(1,g_num(5,ig))=width((ii-1)+j+1)
           g_coord(1,g_num(6,ig))=g_coord(1,g_num(5,ig))
           g_coord(1,g_num(7,ig))=g_coord(1,g_num(5,ig))
           g_coord(1,g_num(3,ig))=                                        &
             pt5*(g_coord(1,g_num(1,ig))+g_coord(1,g_num(5,ig)))
           g_coord(1,g_num(2,ig))=                                        &
             pt5*(g_coord(1,g_num(3,ig))+g_coord(1,g_num(1,ig)))
           g_coord(1,g_num(4,ig))=                                        &
             pt5*(g_coord(1,g_num(3,ig))+g_coord(1,g_num(5,ig)))
           g_coord(1,g_num(8,ig))=                                        &
             pt5*(g_coord(1,g_num(1,ig))+g_coord(1,g_num(7,ig)))
           g_coord(2,g_num(1,ig))=depth((ii-1)+i)
           g_coord(2,g_num(8,ig))=depth((ii-1)+i)
           g_coord(2,g_num(7,ig))=depth((ii-1)+i)
           g_coord(2,g_num(5,ig))=depth((ii-1)+i+1)
           g_coord(2,g_num(3,ig))=                                        &
             pt5*(g_coord(2,g_num(1,ig))+g_coord(2,g_num(5,ig)))
           g_coord(2,g_num(2,ig))=                                        &
             pt5*(g_coord(2,g_num(1,ig))+g_coord(2,g_num(3,ig)))
           g_coord(2,g_num(4,ig))=                                        &
             pt5*(g_coord(2,g_num(3,ig))+g_coord(2,g_num(5,ig)))
           g_coord(2,g_num(6,ig))=                                        &
             pt5*(g_coord(2,g_num(5,ig))+g_coord(2,g_num(7,ig)))
         ELSE
           g_coord(1,g_num(1,ig))=width((ii-1)+j)
           g_coord(1,g_num(2,ig))=width((ii-1)+j)
           g_coord(1,g_num(3,ig))=width((ii-1)+j)
           g_coord(1,g_num(5,ig))=width((ii-1)+j+1)
           g_coord(1,g_num(6,ig))=width((ii-1)+j+1)
           g_coord(1,g_num(7,ig))=width((ii-1)+j+1)
           g_coord(1,g_num(4,ig))=                                        &
             pt5*(g_coord(1,g_num(3,ig))+g_coord(1,g_num(5,ig)))
           g_coord(1,g_num(8,ig))=                                        &
             pt5*(g_coord(1,g_num(1,ig))+g_coord(1,g_num(7,ig)))
           g_coord(2,g_num(1,ig))=depth((ii-1)+i)
           g_coord(2,g_num(8,ig))=depth((ii-1)+i)
           g_coord(2,g_num(7,ig))=depth((ii-1)+i)
           g_coord(2,g_num(3,ig))=depth((ii-1)+i+1)
           g_coord(2,g_num(4,ig))=depth((ii-1)+i+1)
           g_coord(2,g_num(5,ig))=depth((ii-1)+i+1)
           g_coord(2,g_num(2,ig))=                                        &
             pt5*(g_coord(2,g_num(1,ig))+g_coord(2,g_num(3,ig)))
           g_coord(2,g_num(6,ig))=                                        &
             pt5*(g_coord(2,g_num(5,ig))+g_coord(2,g_num(7,ig)))
         END IF
       END DO
     END DO
   END DO
 ELSE
   DO ii=1,lifts-1
     DO i=1,lnye
       DO j=1,lnxe-(ii-1)*1
         ig=ig+1
         if(j==1)then
           g_coord(1,g_num(1,ig))=width((ii-1)+j)
           g_coord(1,g_num(5,ig))=width((ii-1)+j+1)
           g_coord(1,g_num(6,ig))=g_coord(1,g_num(5,ig))
           g_coord(1,g_num(7,ig))=g_coord(1,g_num(5,ig))
           g_coord(1,g_num(3,ig))=g_coord(1,g_num(5,ig))
           g_coord(1,g_num(4,ig))=g_coord(1,g_num(5,ig))
           g_coord(1,g_num(8,ig))=                                        &
             pt5*(g_coord(1,g_num(1,ig))+g_coord(1,g_num(7,ig)))
           g_coord(1,g_num(2,ig))=                                        &
             pt5*(g_coord(1,g_num(1,ig))+g_coord(1,g_num(5,ig)))
           g_coord(2,g_num(1,ig))=depth((ii-1)+i)
           g_coord(2,g_num(8,ig))=depth((ii-1)+i)
           g_coord(2,g_num(7,ig))=depth((ii-1)+i)
           g_coord(2,g_num(5,ig))=depth((ii-1)+i+1)
           g_coord(2,g_num(2,ig))=                                        &
             pt5*(g_coord(2,g_num(1,ig))+g_coord(2,g_num(5,ig)))
           g_coord(2,g_num(6,ig))=                                        &
             pt5*(g_coord(2,g_num(5,ig))+g_coord(2,g_num(7,ig)))
           g_coord(2,g_num(3,ig))=g_coord(2,g_num(5,ig))
           g_coord(2,g_num(4,ig))=g_coord(2,g_num(5,ig))
         ELSE
           g_coord(1,g_num(1,ig))=width((ii-1)+j)
           g_coord(1,g_num(2,ig))=width((ii-1)+j)
           g_coord(1,g_num(3,ig))=width((ii-1)+j)
           g_coord(1,g_num(5,ig))=width((ii-1)+j+1)
           g_coord(1,g_num(6,ig))=width((ii-1)+j+1)
           g_coord(1,g_num(7,ig))=width((ii-1)+j+1)
           g_coord(1,g_num(4,ig))=                                        &
             pt5*(g_coord(1,g_num(3,ig))+g_coord(1,g_num(5,ig)))
           g_coord(1,g_num(8,ig))=                                        &
             pt5*(g_coord(1,g_num(1,ig))+g_coord(1,g_num(7,ig)))
           g_coord(2,g_num(1,ig))=depth((ii-1)+i)
           g_coord(2,g_num(8,ig))=depth((ii-1)+i)
           g_coord(2,g_num(7,ig))=depth((ii-1)+i)
           g_coord(2,g_num(3,ig))=depth((ii-1)+i+1)
           g_coord(2,g_num(4,ig))=depth((ii-1)+i+1)
           g_coord(2,g_num(5,ig))=depth((ii-1)+i+1)
           g_coord(2,g_num(2,ig))=                                        &
             pt5*(g_coord(2,g_num(1,ig))+g_coord(2,g_num(3,ig)))
           g_coord(2,g_num(6,ig))=                                        &
             pt5*(g_coord(2,g_num(5,ig))+g_coord(2,g_num(7,ig)))
         END IF
       END DO
     END DO
   END DO
 END IF
RETURN
END SUBROUTINE fmcoem                       
SUBROUTINE fmglem(fnxe,fnye,lnxe,g_num,lifts)
!
! This subroutine returns g_num for the mesh.
!
 IMPLICIT NONE
 INTEGER,INTENT(IN)::fnxe,fnye,lnxe,lifts
 INTEGER,INTENT(OUT)::g_num(:,:)
 INTEGER::ig,i,j,ii,ilast       
 ig=0
 DO i=1,fnye
   DO j=1,fnxe
     ig=ig+1
     g_num(1,ig)=(j*2-1)+(i-1)*(fnxe+1)+(i-1)*(2*fnxe+1)
     g_num(8,ig)=g_num(1,ig)+1
     g_num(7,ig)=g_num(1,ig)+2
     g_num(2,ig)=i*(2*fnxe+1)+j+(i-1)*(fnxe+1)
     g_num(6,ig)=g_num(2,ig)+1
     g_num(3,ig)=i*(2*fnxe+1)+i*(fnxe+1)+(j*2-1)
     g_num(4,ig)=g_num(3,ig)+1
     g_num(5,ig)=g_num(3,ig)+2
   END DO
 END DO
 ilast=g_num(5,ig)
 DO ii=1,lifts-1
   DO j=1,lnxe-(ii-1)*1
     ig=ig+1
     g_num(1,ig)=ilast-(lnxe-ii+1)*2+(j-1)*2
     g_num(8,ig)=g_num(1,ig)+1
     g_num(7,ig)=g_num(1,ig)+2
     g_num(2,ig)=ilast+j
     g_num(6,ig)=g_num(2,ig)+1
     g_num(3,ig)=ilast+(lnxe-ii+1)+1+(j*2-1)
     g_num(4,ig)=g_num(3,ig)+1
     g_num(5,ig)=g_num(3,ig)+2
   END DO
   ilast=g_num(5,ig)
 END DO
RETURN
END SUBROUTINE fmglem
SUBROUTINE geom_freesurf(iel,nxe,fixed_seep,fixed_down,down,              &
  width,angs,surf,coord,num)
!
! This subroutine forms the coordinates and steering vector
! for 4-node quads numbering in the x-direction
! (Laplace's equation, variable mesh, 1-freedom per node).
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::width(:),angs(:),surf(:),down 
 REAL(iwp),INTENT(OUT)::coord(:,:)
 INTEGER,INTENT(IN)::iel,nxe,fixed_seep,fixed_down
 INTEGER,INTENT(OUT)::num(:)
 REAL(iwp)::angr(SIZE(angs)),pi,b1,b2,tan1,tan2,fac1,fac2,zero=0.0_iwp,   &
   pt5=0.5_iwp,one=1.0_iwp,small=0.0001_iwp,d180=180.0_iwp
 INTEGER::ip,iq
 pi=ACOS(-one)
 angr=angs*pi/d180   
 iq=(iel-1)/nxe+1
 ip=iel-(iq-1)*nxe
 num(1)=iq*(nxe+1)+ip
 num(2)=(iq-1)*(nxe+1)+ip
 num(3)=num(2)+1
 num(4)=num(1)+1
 IF(iq<=fixed_seep+1)THEN
   b1=(surf(ip)-down)/real(fixed_seep+1)
   b2=(surf(ip+1)-down)/real(fixed_seep+1)
   coord(1,2)=down+(fixed_seep+1-iq)*b1 
   coord(2,2)=down+(fixed_seep+2-iq)*b1
   coord(3,2)=down+(fixed_seep+2-iq)*b2 
   coord(4,2)=down+(fixed_seep+1-iq)*b2
 ELSE
   b1=real(fixed_down+fixed_seep-iq)/real(fixed_down-1) 
   b2=real(fixed_down+fixed_seep-iq+1)/real(fixed_down-1) 
   coord(1,2)=down*b1
   coord(2,2)=down*b2
   coord(3,2)=coord(2,2)
   coord(4,2)=coord(1,2)
 END IF
 IF(ABS(angr(ip)-pi*pt5)<small)THEN
   fac1=zero
 ELSE
   tan1=TAN(angr(ip))      
   fac1=one/tan1
 END IF
 IF(ABS(angr(ip+1)-pi*pt5)<small)THEN
   fac2=zero
 ELSE
   tan2=TAN(angr(ip+1))    
   fac2=one/tan2 
 END IF
 coord(1,1)=width(ip)+coord(1,2)*fac1 
 coord(2,1)=width(ip)+coord(2,2)*fac1
 coord(3,1)=width(ip+1)+coord(3,2)*fac2 
 coord(4,1)=width(ip+1)+coord(4,2)*fac2
RETURN
END SUBROUTINE geom_freesurf                                               

SUBROUTINE geom_rect(element,iel,x_coords,y_coords,e_nd_coord,num,dir)
!
! This subroutine forms the coordinates and connectivity for a
! rectangular mesh of rt. angled triangular elements (3, 6, 10 or 15-node)
! or quadrilateral elements (4, 8 or 9-node) counting in the
! x- or y-dir. 
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::x_coords(:),y_coords(:)
 REAL(iwp),INTENT(OUT)::coord(:,:)
 CHARACTER(LEN=15),INTENT(IN)::element
 CHARACTER(LEN=1),INTENT(IN)::dir
 INTEGER,INTENT(IN)::iel
 INTEGER,INTENT(OUT)::num(:)
 INTEGER::ip,iq,jel,fac1,nod,nxe,nye
 REAL(iwp)::pt5=0.5_iwp,two=2.0_iwp,d3=3.0_iwp 
 nxe=UBOUND(x_coords,1)-1
 nod=UBOUND(num,1)
 IF(element=='triangle')THEN
   nye=(UBOUND(y_coords,1)-1)*2
   IF(dir=='x'.OR.dir=='r')THEN
     jel=2*nxe*((iel-1)/(2*nxe))
     ip=(iel-jel+1)/2
     iq=2*((iel-1)/(2*nxe)+1)-1+((iel/2)*2)/iel
   ELSE  
     jel=(iel-1)/nye
     ip=jel+1
     iq=iel-nye*jel
   END IF
   SELECT CASE(nod)
   CASE(3)
     IF(MOD(iq,2)/=0)THEN
       IF(dir=='x'.OR.dir=='r')THEN
         num(1)=(nxe+1)*(iq-1)/2+ip
         num(2)=num(1)+1              
         num(3)=(nxe+1)*(iq+1)/2+ip
       ELSE
         num(1)=(ip-1)*(nye+2)/2+(iq+1)/2
         num(2)=num(1)+(nye+2)/2
         num(3)=num(1)+1
       END IF
!
       coord(1,1)=x_coords(ip)
       coord(1,2)=y_coords((iq+1)/2)
       coord(2,1)=x_coords(ip+1)   
       coord(2,2)=y_coords((iq+1)/2)
       coord(3,1)=x_coords(ip)   
       coord(3,2)=y_coords((iq+3)/2)
     ELSE
       IF(dir=='x'.OR.dir=='r')THEN
         num(1)=(nxe+1)*iq/2+ip+1     
         num(2)=num(1)-1               
         num(3)=(nxe+1)*(iq-2)/2+ip+1
       ELSE
         num(1)=ip*(nye+2)/2+(iq+2)/2
         num(2)=(ip-1)*(nye+2)/2+(iq+1)/2+1
         num(3)=num(1)-1
       END IF
!
       coord(1,1)=x_coords(ip+1)
       coord(1,2)=y_coords((iq+2)/2)
       coord(2,1)=x_coords(ip)   
       coord(2,2)=y_coords((iq+2)/2)
       coord(3,1)=x_coords(ip+1) 
       coord(3,2)=y_coords(iq/2)
     END IF
   CASE(6)
     IF(MOD(iq,2)/=0)THEN
       IF(dir=='x'.OR.dir=='r')THEN
         num(1)=(iq-1)*(2*nxe+1)+2*ip-1
         num(2)=num(1)+1 
         num(3)=num(1)+2 
         num(4)=(iq-1)*(2*nxe+1)+2*nxe+2*ip+1
         num(5)=(iq+1)*(2*nxe+1)+2*ip-1
         num(6)=num(4)-1 
       ELSE
         num(1)=2*(nye+1)*(ip-1)+iq
         num(2)=2*(nye+1)*(ip-1)+nye+1+iq
         num(3)=2*(nye+1)*ip+iq
         num(4)=num(2)+1
         num(5)=num(1)+2 
         num(6)=num(1)+1
       END IF
!
       coord(1,1)=x_coords(ip)
       coord(1,2)=y_coords((iq+1)/2)
       coord(3,1)=x_coords(ip+1)   
       coord(3,2)=y_coords((iq+1)/2)
       coord(5,1)=x_coords(ip)   
       coord(5,2)=y_coords((iq+3)/2)
     ELSE
       IF(dir=='x'.OR.dir=='r')THEN
         num(1)=iq*(2*nxe+1)+2*ip+1
         num(2)=num(1)-1 
         num(3)=num(1)-2 
         num(4)=(iq-2)*(2*nxe+1)+2*nxe+2*ip+1
         num(5)=(iq-2)*(2*nxe+1)+2*ip+1
         num(6)=num(4)+1 
       ELSE 
         num(1)=2*(nye+1)*ip+iq+1 
         num(2)=2*(nye+1)*(ip-1)+nye+iq+2
         num(3)=2*(nye+1)*(ip-1)+iq+1
         num(4)=num(2)-1 
         num(5)=num(1)-2
         num(6)=num(1)-1
       END IF
!
       coord(1,1)=x_coords(ip+1)
       coord(1,2)=y_coords((iq+2)/2)
       coord(3,1)=x_coords(ip)   
       coord(3,2)=y_coords((iq+2)/2)
       coord(5,1)=x_coords(ip+1) 
       coord(5,2)=y_coords(iq/2)
     END IF
     coord(2,:)=pt5*(coord(1,:)+coord(3,:))
     coord(4,:)=pt5*(coord(3,:)+coord(5,:))
     coord(6,:)=pt5*(coord(5,:)+coord(1,:))
   CASE(10)
     IF(MOD(iq,2)/=0)THEN
       IF(dir=='x'.OR.dir=='r')THEN
         num(1)=(iq-1)/2*(3*nxe+1)*3+3*ip-2
         num(2)=num(1)+1
         num(3)=num(1)+2
         num(4)=num(1)+3
         num(5)=(iq-1)/2*(3*nxe+1)*3+3*nxe+1+3*ip
         num(6)=(iq-1)/2*(3*nxe+1)*3+6*nxe+2+3*ip-1
         num(7)=(iq-1)/2*(3*nxe+1)*3+9*nxe+3+3*ip-2
         num(8)=num(6)-1
         num(9)=num(5)-2
         num(10)=num(9)+1
       ELSE
         num(1)=(9*(nye-2)/2+12)*(ip-1)+3*(iq-1)/2+1
         num(2)=(9*(nye-2)/2+12)*(ip-1)+3*(nye-2)/2+4+3*(iq-1)/2+1
         num(3)=(9*(nye-2)/2+12)*(ip-1)+3*(nye-2)+8+3*(iq-1)/2+1
         num(4)=(9*(nye-2)/2+12)*(ip-1)+9*(nye-2)/2+12+3*(iq-1)/2+1
         num(5)=num(3)+1 
         num(6)=num(2)+2
         num(7)=num(1)+3
         num(8)=num(1)+2
         num(9)=num(1)+1
         num(10)=num(2)+1
       END IF
!
       coord(1,1)=x_coords(ip)
       coord(2,1)=x_coords(ip)+(x_coords(ip+1)-x_coords(ip))/d3
       coord(3,1)=x_coords(ip)+two*(x_coords(ip+1)-x_coords(ip))/d3
       coord(4,1)=x_coords(ip+1)
       coord(4,2)=y_coords((iq+1)/2)
       coord(5,2)=y_coords((iq+1)/2)+                                     &
         (y_coords((iq+3)/2)-y_coords((iq+1)/2))/d3
       coord(6,2)=y_coords((iq+1)/2)+                                     &
         two*(y_coords((iq+3)/2)-y_coords((iq+1)/2))/d3
       coord(7,2)=y_coords((iq+3)/2)
     ELSE
       IF(dir=='x'.OR.dir=='r')THEN
         num(1)=(iq-2)/2*(3*nxe+1)*3+9*nxe+3+3*ip+1
         num(2)=num(1)-1
         num(3)=num(1)-2
         num(4)=num(1)-3
         num(5)=(iq-2)/2*(3*nxe+1)*3+6*nxe+2+3*ip-1
         num(6)=(iq-2)/2*(3*nxe+1)*3+3*nxe+1+3*ip
         num(7)=(iq-2)/2*(3*nxe+1)*3+3*ip+1
         num(8)=num(6)+1
         num(9)=num(5)+2
         num(10)=num(9)-1
       ELSE
         num(1)=(9*(nye-2)/2+12)*(ip-1)+9*(nye-2)/2+12+3*iq/2+1
         num(2)=(9*(nye-2)/2+12)*(ip-1)+3*(nye-2)+8+3*iq/2+1
         num(3)=(9*(nye-2)/2+12)*(ip-1)+3*(nye-2)/2+4+3*iq/2+1
         num(4)=(9*(nye-2)/2+12)*(ip-1)+3*iq/2+1
         num(5)=num(3)-1
         num(6)=num(2)-2
         num(7)=num(1)-3
         num(8)=num(1)-2
         num(9)=num(1)-1
         num(10)=num(2)-1
       END IF
!
       coord(1,1)=x_coords(ip+1)
       coord(2,1)=x_coords(ip+1)-(x_coords(ip+1)-x_coords(ip))/d3
       coord(3,1)=x_coords(ip+1)-two*(x_coords(ip+1)-x_coords(ip))/d3
       coord(4,1)=x_coords(ip)
       coord(4,2)=y_coords((iq+2)/2)
       coord(5,2)=y_coords((iq+2)/2)-(y_coords((iq+2)/2)-y_coords(iq/2))/d3
       coord(6,2)=y_coords((iq+2)/2)-                                     &
         two*(y_coords((iq+2)/2)-y_coords(iq/2))/d3
       coord(7,2) =y_coords(iq/2)
     END IF
     coord(5,1)=coord(3,1)
     coord(6,1)=coord(2,1)
     coord(7,1)=coord(1,1)
     coord(8,1)=coord(1,1)
     coord(9,1)=coord(1,1)
     coord(10,1)=coord(2,1)
     coord(1,2)=coord(4,2)
     coord(2,2)=coord(4,2)
     coord(3,2)=coord(4,2)
     coord(8,2)=coord(6,2)
     coord(9,2)=coord(5,2)
     coord(10,2)=coord(5,2)
   CASE(15)
     IF(MOD(iq,2)/=0)THEN
       IF(dir=='x'.OR.dir=='r')THEN
       fac1=4*(4*nxe+1)*(iq-1)/2
         num(1)=fac1+4*ip-3
         num(2)=num(1)+1
         num(3)=num(1)+2
         num(4)=num(1)+3
         num(5)=num(1)+4
         num(6)=fac1+ 4*nxe+1+4*ip
         num(7)=fac1+ 8*nxe+1+4*ip
         num(8)=fac1+12*nxe+1+4*ip
         num(9)=fac1+16*nxe+1+4*ip
         num(10)=num(8)-1
         num(11)=num(7)-2
         num(12)=num(6)-3
         num(13)=num(12)+1
         num(14)=num(12)+2
         num(15)=num(11)+1
       ELSE
         fac1=4*(2*nye+1)*(ip-1)+2*iq-1 
         num(1)=fac1
         num(2)=fac1+2*nye+1
         num(3)=fac1+4*nye+2 
         num(4)=fac1+6*nye+3 
         num(5)=fac1+8*nye+4
         num(6)=fac1+6*nye+4 
         num(7)=fac1+4*nye+4 
         num(8)=fac1+2*nye+4
         num(9)=fac1+4 
         num(10)=fac1+3 
         num(11)=fac1+2 
         num(12)=fac1+1
         num(13)=fac1+2*nye+2 
         num(14)=fac1+4*nye+3
         num(15)=fac1+2*nye+3  
       END IF
!
       coord(1,1)=x_coords(ip)
       coord(1,2)=y_coords((iq+1)/2)
       coord(5,1)=x_coords(ip+1)   
       coord(5,2)=y_coords((iq+1)/2)
       coord(9,1)=x_coords(ip)   
       coord(9,2)=y_coords((iq+3)/2)
     ELSE
       IF(dir=='x'.OR.dir=='r')THEN
         fac1=4*(4*nxe+1)*(iq-2)/2
         num(1)=fac1+16*nxe+5+4*ip
         num(2)=num(1)-1
         num(3)=num(1)-2
         num(4)=num(1)-3
         num(5)=num(1)-4
         num(6)=fac1+12*nxe+1+4*ip
         num(7)=fac1+8*nxe+1+4*ip
         num(8)=fac1+4*nxe+1+4*ip
         num(9)=fac1+4*ip+1
         num(10)=num(8)+1
         num(11)=num(7)+2
         num(12)=num(6)+3
         num(13)=num(12)-1
         num(14)=num(12)-2
         num(15)=num(11)-1
       ELSE
         fac1=4*(2*nye+1)*(ip-1)+2*iq+8*nye+5 
         num(1)=fac1 
         num(2)=fac1-2*nye-1
         num(3)=fac1-4*nye-2 
         num(4)=fac1-6*nye-3 
         num(5)=fac1-8*nye-4
         num(6)=fac1-6*nye-4  
         num(7)=fac1-4*nye-4 
         num(8)=fac1-2*nye-4
         num(9)=fac1-4
         num(10)=fac1-3 
         num(11)=fac1-2 
         num(12)=fac1-1
         num(13)=fac1-2*nye-2  
         num(14)=fac1-4*nye-3
         num(15)=fac1-2*nye-3 
       END IF
!
       coord(1,1)=x_coords(ip+1)
       coord(1,2)=y_coords((iq+2)/2)
       coord(5,1)=x_coords(ip)   
       coord(5,2)=y_coords((iq+2)/2)
       coord(9,1)=x_coords(ip+1) 
       coord(9,2)=y_coords(iq/2)
     END IF
     coord(3,:)=pt5*(coord(1,:)+coord(5,:))
     coord(7,:)=pt5*(coord(5,:)+coord(9,:))
     coord(11,:)=pt5*(coord(9,:)+coord(1,:))
     coord(2,:)=pt5*(coord(1,:)+coord(3,:))
     coord(4,:)=pt5*(coord(3,:)+coord(5,:))
     coord(6,:)=pt5*(coord(5,:)+coord(7,:))
     coord(8,:)=pt5*(coord(7,:)+coord(9,:))
     coord(10,:)=pt5*(coord(9,:)+coord(11,:))
     coord(12,:)=pt5*(coord(11,:)+coord(1,:))
     coord(15,:)=pt5*(coord(7,:)+coord(11,:))
     coord(14,:)=pt5*(coord(3,:)+coord(7,:))
     coord(13,:)=pt5*(coord(2,:)+coord(15,:))
   CASE DEFAULT
     WRITE(11,'(a)')"Wrong number of nodes for triangular element"
     STOP
   END SELECT
 ELSE
   nye=UBOUND(y_coords,1)-1
   IF(dir=='x'.OR.dir=='r')THEN
     iq=(iel-1)/nxe+1
     ip=iel-(iq-1)*nxe
   ELSE
     ip=(iel-1)/nye+1
     iq=iel-(ip-1)*nye
   END IF
   SELECT CASE(nod)
   CASE(4)
     IF(dir=='x'.OR.dir=='r')THEN
       num(1)=iq*(nxe+1)+ip		        		
       num(2)=(iq-1)*(nxe+1)+ip				
       num(3)=num(2)+1					
       num(4)=num(1)+1					
     ELSE
       num(1)=(ip-1)*(nye+1)+iq+1
       num(2)=num(1)-1
       num(3)=ip*(nye+1)+iq
       num(4)=num(3)+1
     END IF
!
     coord(1:2,1)=x_coords(ip)
     coord(3:4,1)=x_coords(ip+1)
     coord(1,2)=y_coords(iq+1)
     coord(2:3,2)=y_coords(iq)
     coord(4,2)=coord(1,2)
   CASE(5)
     IF(dir=='x'.OR.dir=='r')THEN
       num(1)=iq*(2*nxe+1)+ip
       num(2)=(iq-1)*(2*nxe+1)+ip
       num(3)=num(2)+1
       num(4)=num(1)+1
       num(5)=iq*(2*nxe+1)+ip-nxe
     ELSE
       num(1)=(ip-1)*(2*nye+1)+iq+1
       num(2)=num(1)-1
       num(3)=ip*(2*nye+1)+iq
       num(4)=num(3)+1
       num(5)=ip*(2*nye+1)+iq-nye
     END IF
!
     coord(1:2,1)=x_coords(ip)
     coord(3:4,1)=x_coords(ip+1)
     coord(1,2)=y_coords(iq+1)
     coord(2:3,2)=y_coords(iq)
     coord(4,2)=coord(1,2)
     coord(5,:)=0.25_iwp*(coord(1,:)+coord(2,:)+coord(3,:)+coord(4,:))
   CASE(8)
     IF(dir=='x'.OR.dir=='r')THEN
       num(1)=iq*(3*nxe+2)+2*ip-1                 
       num(2)=iq*(3*nxe+2)+ip-nxe-1		  
       num(3)=(iq-1)*(3*nxe+2)+2*ip-1		   
       num(4)=num(3)+1
       num(5)=num(4)+1
       num(6)=num(2)+1
       num(7)=num(1)+2
       num(8)=num(1)+1
     ELSE
       num(1)=(ip-1)*(3*nye+2)+2*iq+1
       num(2)=num(1)-1
       num(3)=num(1)-2
       num(4)=(ip-1)*(3*nye+2)+2*nye+iq+1
       num(5)=ip*(3*nye+2)+2*iq-1
       num(6)=num(5)+1
       num(7)=num(5)+2
       num(8)=num(4)+1
     END IF
!
     coord(1:3,1)=x_coords(ip)
     coord(5:7,1)=x_coords(ip+1)
     coord(4,1)=pt5*(coord(3,1)+coord(5,1))
     coord(8,1)=pt5*(coord(7,1)+coord(1,1))
     coord(1,2)=y_coords(iq+1)
     coord(7:8,2)=y_coords(iq+1)
     coord(3:5,2)=y_coords(iq)
     coord(2,2)=pt5*(coord(1,2)+coord(3,2))
     coord(6,2)=pt5*(coord(5,2)+coord(7,2))
   CASE(9)
     IF(dir=='x'.OR.dir=='r')THEN
       num(1)=iq*(4*nxe+2)+2*ip-1
       num(2)=iq*(4*nxe+2)+2*ip-nxe-4
       num(3)= (iq-1)*(4*nxe+2)+2*ip-1
       num(4)=num(3)+1
       num(5)=num(4)+1
       num(6)=num(2)+2
       num(7)=num(1)+2
       num(8)=num(1)+1
       num(9)=num(2)+1
     ELSE
       num(1)=(ip-1)*2*(2*nye+1)+2*iq+1
       num(2)=num(1)-1
       num(3)=num(1)-2
       num(4)=(ip-1)*2*(2*nye+1)+2*nye+2*iq
       num(5)=ip*2*(2*nye+1)+2*iq-1
       num(6)=num(5)+1
       num(7)=num(5)+2
       num(8)=num(4)+2
       num(9)=num(4)+1
     END IF
!
     coord(1:3,1)=x_coords(ip)
     coord(5:7,1)=x_coords(ip+1)
     coord(4,1)=pt5*(coord(3,1)+coord(5,1))
     coord(8,1)=pt5*(coord(7,1)+coord(1,1))
     coord(1,2)=y_coords(iq+1)
     coord(7:8,2)=y_coords(iq+1)
     coord(3:5,2)=y_coords(iq)
     coord(2,2)=pt5*(coord(1,2)+coord(3,2))
     coord(6,2)=pt5*(coord(5,2)+coord(7,2))
     coord(9,:)=pt5*(coord(4,:)+coord(8,:))
   CASE DEFAULT
     WRITE(11,'(a)')"Wrong number of nodes for quadrilateral element"
     STOP
   END SELECT
 END IF
RETURN
END SUBROUTINE geom_rect
SUBROUTINE hexahedron_xz(iel,x_coords,y_coords,z_coords,coord,num)
!
! This subroutine generates nodal coordinates and numbering for
! 8, 14 or 20-node "bricks" counting x-z planes in the y-direction.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 INTEGER,INTENT(IN)::iel
 REAL(iwp),INTENT(IN)::x_coords(:),y_coords(:),z_coords(:)
 REAL(iwp),INTENT(OUT)::coord(:,:)
 INTEGER,INTENT(OUT)::num(:)
 REAL(iwp)::pt5=0.5_iwp
 INTEGER::fac1,fac2,ip,iq,is,iplane,nod,nxe,nze
 nxe=UBOUND(x_coords,1)-1
 nze=UBOUND(z_coords,1)-1
 nod=UBOUND(num,1)
 iq=(iel-1)/(nxe*nze)+1
 iplane=iel-(iq-1)*nxe*nze
 is=(iplane-1)/nxe+1 
 ip=iplane-(is-1)*nxe
 SELECT CASE(nod)
 CASE(8)
   fac1=(nxe+1)*(nze+1)*(iq-1)
   num(1)=fac1+is*(nxe+1)+ip
   num(2)=fac1+(is-1)*(nxe+1)+ip
   num(3)=num(2)+1
   num(4)=num(1)+1
   num(5)=(nxe+1)*(nze+1)*iq+is*(nxe+1)+ip
   num(6)=(nxe+1)*(nze+1)*iq+(is-1)*(nxe+1)+ip
   num(7)=num(6)+1
   num(8)=num(5)+1
!
   coord(1:2,1)=x_coords(ip)
   coord(5:6,1)=x_coords(ip)
   coord(3:4,1)=x_coords(ip+1)
   coord(7:8,1)=x_coords(ip+1)
!
   coord(1:4,2)=y_coords(iq)
   coord(5:8,2)=y_coords(iq+1)
!
   coord(2:3,3)=z_coords(is)
   coord(6:7,3)=z_coords(is)
   coord(1:4:3,3)=z_coords(is+1)
   coord(5:8:3,3)=z_coords(is+1)
!
 CASE(14)
   fac1=(2*nxe+1)*(2*nze+1)*(iq-1)
   fac2=(2*nxe+1)*(2*nze+1)*iq
   num(1)=fac1+is*(2*nxe+1)+ip
   num(2)=num(1)-(2*nxe+1)
   num(3)=num(2)+1
   num(4)=num(1)+1
   num(5)=num(2)+nxe+1
   num(6)=fac1+(nxe+1)*(nze+1)+nxe*nze+(is-1)*(2*nxe+1)+nxe+ip
   num(7)=num(6)-nxe
   num(8)=num(6)+1
   num(9)=num(8)+nxe
   num(10)=fac2+is*(2*nxe+1)+ip
   num(11)=num(10)-(2*nxe+1)
   num(12)=num(11)+1
   num(13)=num(10)+1
   num(14)=num(11)+nxe+1
!
   coord(1:2,1)=x_coords(ip)
   coord(6,1)=x_coords(ip)
   coord(10:11,1)=x_coords(ip)
   coord(5:9:2,1)=pt5*(x_coords(ip)+x_coords(ip+1))
   coord(14,1)=pt5*(x_coords(ip)+x_coords(ip+1))
   coord(3:4,1)=x_coords(ip+1)
   coord(8,1)=x_coords(ip+1)
   coord(12:13,1)=x_coords(ip+1)
!
   coord(1:5,2)=y_coords(iq)
   coord(6:9,2)=pt5*(y_coords(iq)+y_coords(iq+1))
   coord(10:14,2)=y_coords(iq+1)
!
   coord(2:3,3)=z_coords(is)
   coord(7,3)=z_coords(is)
   coord(11:12,3)=z_coords(is)
   coord(5:6,3)=pt5*(z_coords(is)+z_coords(is+1))
   coord(8:14:6,3)=pt5*(z_coords(is)+z_coords(is+1))
   coord(1:4:3,3)=z_coords(is+1)
   coord(9,3)=z_coords(is+1)
   coord(10:13:3,3)=z_coords(is+1)
!
 CASE(20)
   fac1=((2*nxe+1)*(nze+1)+(2*nze+1)*(nxe+1))*(iq-1)
   fac2=((2*nxe+1)*(nze+1)+(2*nze+1)*(nxe+1))*iq
   num(1)=fac1+(3*nxe+2)*is+2*ip-1
   num(2)=fac1+(3*nxe+2)*is-nxe+ip-1
   num(3)=num(1)-3*nxe-2
   num(4)=num(3)+1
   num(5)=num(4)+1
   num(6)=num(2)+1
   num(7)=num(1)+2
   num(8)=num(1)+1
   num(9)=fac2-(nxe+1)*(nze+1)+(nxe+1)*is+ip
   num(10)=num(9)-nxe-1
   num(11)=num(10)+1
   num(12)=num(9)+1
   num(13)=fac2+(3*nxe+2)*is+2*ip-1
   num(14)=fac2+(3*nxe+2)*is-nxe+ip-1
   num(15)=num(13)-3*nxe-2
   num(16)=num(15)+1
   num(17)=num(16)+1
   num(18)=num(14)+1
   num(19)=num(13)+2
   num(20)=num(13)+1 
!
   coord(1:3,1)=x_coords(ip)
   coord(9:10,1)=x_coords(ip)
   coord(13:15,1)=x_coords(ip)
   coord(5:7,1)=x_coords(ip+1)
   coord(11:12,1)=x_coords(ip+1)
   coord(17:19,1)=x_coords(ip+1)
   coord(4:8:4,1)=pt5*(x_coords(ip)+x_coords(ip+1))
   coord(16:20:4,1)=pt5*(x_coords(ip)+x_coords(ip+1))
!
   coord(1:8,2)=y_coords(iq)
   coord(13:20,2)=y_coords(iq+1)
   coord(9:12,2)=pt5*(y_coords(iq)+y_coords(iq+1))
!
   coord(1,3)=z_coords(is+1)
   coord(7:9,3)=z_coords(is+1)
   coord(12:13,3)=z_coords(is+1)
   coord(19:20,3)=z_coords(is+1)
   coord(3:5,3)=z_coords(is)
   coord(10:11,3)=z_coords(is)
   coord(15:17,3)=z_coords(is)
   coord(2:6:4,3)=pt5*(z_coords(is)+z_coords(is+1))
   coord(14:18:4,3)=pt5*(z_coords(is)+z_coords(is+1))
!
 CASE DEFAULT
   WRITE(11,'(a)')"Wrong number of nodes for hexahedral element"
   STOP
 END SELECT
RETURN
END SUBROUTINE hexahedron_xz
SUBROUTINE mesh_size(element,nod,nels,nn,nxe,nye,nze)
!
!  This subroutine returns the number of elements (nels) and the number
!  of nodes (nn) in a 2-d geometry-created mesh.
!
 IMPLICIT NONE
 CHARACTER(LEN=15),INTENT(IN)::element
 INTEGER,INTENT(IN)::nod,nxe,nye
 INTEGER,INTENT(IN),OPTIONAL::nze
 INTEGER,INTENT(OUT)::nels,nn
 IF(element=="triangle")THEN
   nels=nxe*nye*2
   IF(nod==3)nn=(nxe+1)*(nye+1)
   IF(nod==6)nn=(2*nxe+1)*(2*nye+1)
   IF(nod==10)nn=(3*nxe+1)*(3*nye+1)
   IF(nod==15)nn=(4*nxe+1)*(4*nye+1)
 ELSE IF(element=="quadrilateral")THEN
   nels=nxe*nye
   IF(nod==4)nn=(nxe+1)*(nye+1)
   IF(nod==5)nn=(nxe+1)*(nye+1)+nxe*nye
   IF(nod==8)nn=(2*nxe+1)*(nye+1)+(nxe+1)*nye
   IF(nod==9)nn=(2*nxe+1)*(2*nye+1)
 ELSE IF(element=="hexahedron")THEN
   nels=nxe*nye*nze
   IF(nod==8)nn=(nxe+1)*(nye+1)*(nze+1)
   IF(nod==14)nn=4*nxe*nye*nze+2*(nxe*nye+nye*nze+nze*nxe)+nxe+nye+nze+1
   IF(nod==20)nn=((2*nxe+1)*(nze+1)+(nxe+1)*nze)*(nye+1)+                 &
     (nxe+1)*(nze+1)*nye
 END IF
RETURN
END SUBROUTINE mesh_size

END MODULE geom