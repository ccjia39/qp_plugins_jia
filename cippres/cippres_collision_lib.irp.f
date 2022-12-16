use bitmasks ! you need to include the bitmasks_module.f90 features
use general
!  If an entity is declared with a BEGIN_PROVIDER ... END_PROVIDER block, 
! then it is an IRP entity and it will behave as a global variable in the whole program. 
 BEGIN_PROVIDER [integer, n_sta_coll_max]
  implicit none
     if(n_csf_max<5000) then
        n_sta_coll_max = n_csf_max
     else
        n_sta_coll_max = 5000
     endif
 END_PROVIDER 

 BEGIN_PROVIDER [integer, ib_coll]
  implicit none
   ib_coll = 1
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, b_coll]
  implicit none
   b_coll = 0d0
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, v_coll]
  implicit none
     call ezfio_get_cippres_v_coll(v_coll)
 END_PROVIDER 

 BEGIN_PROVIDER [integer, i_state_coll]
  implicit none
     call ezfio_get_cippres_i_state_coll(i_state_coll)
 END_PROVIDER 

 BEGIN_PROVIDER [integer, stamin_bound]
  implicit none
     call ezfio_get_cippres_stamin_bound(stamin_bound)
 END_PROVIDER 

 BEGIN_PROVIDER [integer, stamax_bound]
  implicit none
     call ezfio_get_cippres_stamax_bound(stamax_bound)
 END_PROVIDER 

 BEGIN_PROVIDER [integer, stamin_si]
  implicit none
     call ezfio_get_cippres_stamin_si(stamin_si)
 END_PROVIDER 

 BEGIN_PROVIDER [integer, stamax_si]
  implicit none
     call ezfio_get_cippres_stamax_si(stamax_si)
 END_PROVIDER 

 BEGIN_PROVIDER [integer, stamin_di]
  implicit none
     call ezfio_get_cippres_stamin_di(stamin_di)
 END_PROVIDER 

 BEGIN_PROVIDER [integer, stamax_di]
  implicit none
     call ezfio_get_cippres_stamax_di(stamax_di)
 END_PROVIDER 

 BEGIN_PROVIDER [integer, n_time]
  implicit none
     call ezfio_get_cippres_n_time(n_time)
 END_PROVIDER 

 BEGIN_PROVIDER [integer, n_bimp]
  implicit none
     call ezfio_get_cippres_n_bimp(n_bimp)
 END_PROVIDER 

 BEGIN_PROVIDER [integer, n_pcenter]
  implicit none
     call ezfio_get_cippres_n_pcenter(n_pcenter)
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, charge_pcenter, (n_pcenter)]
  implicit none
     call ezfio_get_cippres_charge_pcenter(charge_pcenter)
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, bgrid, (n_bimp)]
  implicit none
     call ezfio_get_cippres_bgrid(bgrid)
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, zgrid, (n_time)]
  implicit none
     call ezfio_get_cippres_zgrid(zgrid)
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, tgrid, (n_time)]
  implicit none
     call ezfio_get_cippres_tgrid(tgrid)
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, coll_w1e_mo, (mo_num,mo_num,n_time,n_bimp)]
  implicit none

  double precision:: zg, bg
  integer :: i, j, k, l, ib, iz

  open(unit=10,file='ints/onee_int_tt.txt')
  do ib = 1, n_bimp
     read(10,*)bg
!     print*,bg
    do iz = 1, n_time
     read(10,*)zg
       do i = 1, mo_num
         do j = 1, mo_num
             read(10,*)k,l,coll_w1e_mo(j,i,iz,ib)  !!mcoup(izgrid,j,i),movl(izgrid,j,i)
!             print*,k,l,coll_w1e_mo(j,i,iz,ib)  !!mcoup(izgrid,j,i),movl(izgrid,j,i)
         enddo
       enddo
    enddo
  enddo
  print*,'read ok'
  close(10)
 END_PROVIDER 


 BEGIN_PROVIDER [integer, ntdet]
&BEGIN_PROVIDER [integer, ntsta]
&BEGIN_PROVIDER [integer, npdet]
&BEGIN_PROVIDER [integer, npsta]

  integer :: i, j
  ntdet = 0
  npdet = 0
 
  open(unit=10,file='ints/tcistates_det.txt')
   read(10,*)i,ntdet,ntsta
  close(10)
  open(unit=10,file='ints/pcistates_det.txt')
   read(10,*)i,npdet,npsta
  close(10)

 END_PROVIDER 

 BEGIN_PROVIDER [integer, tdeta, (elec_alpha_num,ntdet)]
&BEGIN_PROVIDER [integer, tdetb, (elec_beta_num,ntdet)]
&BEGIN_PROVIDER [double precision, ctdet, (ntdet)]
&BEGIN_PROVIDER [double precision, tci_e, (ntsta)]
&BEGIN_PROVIDER [double precision, tci_sta, (ntdet,ntsta)]
&BEGIN_PROVIDER [integer, pdeta, (elec_alpha_num,npdet)]
&BEGIN_PROVIDER [integer, pdetb, (elec_beta_num,npdet)]
&BEGIN_PROVIDER [double precision, cpdet, (npdet)]
&BEGIN_PROVIDER [double precision, pci_e, (npsta)]
&BEGIN_PROVIDER [double precision, pci_sta, (npdet,npsta)]

 integer :: i, j ,k, l

 tdeta(:,:) = 0
 tdetb(:,:) = 0
 ctdet(:) = 0d0
 pdeta(:,:) = 0
 pdetb(:,:) = 0
 cpdet(:) = 0d0

 open(unit=10,file='ints/tcistates_det.txt')
    read(10,*)i,j,k
  do i = 1, ntdet
    read(10,*)ctdet(i), (tdeta(j,i),j=1,elec_alpha_num), (tdetb(j,i),j=1,elec_beta_num)
  enddo 
  do i = 1, ntsta
    read(10,*)tci_e(i)
    read(10,*)(tci_sta(j,i),j=1,ntdet)
  enddo 
 close(10)

 open(unit=10,file='ints/pcistates_det.txt')
    read(10,*)i,j,k
  do i = 1, npdet
    read(10,*)cpdet(i), (pdeta(j,i),j=1,elec_alpha_num), (pdetb(j,i),j=1,elec_beta_num)
  enddo 
  do i = 1, npsta
    read(10,*)pci_e(i)
    read(10,*)(pci_sta(j,i),j=1,npdet)
  enddo 
 close(10)

 END_PROVIDER 

 BEGIN_PROVIDER [double precision, coll_couplings, (n_sta_coll_max,n_sta_coll_max,n_time)]
 use general
 use SlaterDeterminant
 implicit none
 integer :: i, j, k, l, imo
 integer :: ib, ic, it
 double precision, allocatable :: eigval1(:),eigvec1(:,:),eigval2(:),eigvec2(:,:),coll_csf_mat(:,:),coll_mat(:,:), mattmp(:,:)
 double precision, dimension(mo_num,mo_num) :: w1e

 integer :: ne, nea, neb, n_mo
 double complex, dimension(mo_num,mo_num) :: ovmo
 double complex, dimension(mo_num,mo_num,mo_num,mo_num) :: r12mo

 double complex :: ov, h1e, r12

 integer :: nsta, ncsf
 double precision :: hij

 double precision :: t1, t2

 integer :: ni, nf

 logical :: exists

 PROVIDE ezfio_filename !HF_bitmask mo_coef

!! TO BE DONE:  read all integrals and store them in single matrices 1 -> (ntdet+npdet)
!!              read eigevec from cistates_det.txt which should be written manually to incorporate both target and proj. MOs
!!              initiate psi with target CI coeff.
!!              test lowdin rules code (ok for tt-1e ints.)

   print*,'Computing coll_couplings', b_coll 
   call cpu_time(t1)
   coll_couplings(:,:,:) = 0d0

   allocate(coll_csf_mat(ntdet,ntdet))
 
   allocate(eigval1(ntsta),eigval2(ntsta))
   eigval1(:) = tci_e(:)
   eigval2(:) = tci_e(:)

   allocate(eigvec1(ntdet,ntsta))
   allocate(eigvec2(ntdet,ntsta))
   eigvec1(:,:) = tci_sta(:,:)
   eigvec2(:,:) = tci_sta(:,:)

   nsta = ntsta
   if(nsta>n_sta_coll_max) then
    print*, "nsta > n_sta_coll_max, I stop"
    stop
   endif

   ncsf = ntdet
   allocate(coll_mat(nsta,nsta))
   coll_mat(:,:) = 0d0
   allocate(mattmp(ncsf,nsta))

  do it = 1, n_time
   w1e(:,:) = coll_w1e_mo(:,:,it,ib_coll)
   coll_csf_mat(:,:) = 0d0

   nea = elec_alpha_num 
   neb = elec_beta_num 
   ne = nea + neb
   n_mo = mo_num
   ovmo(:,:) = 0d0
   do imo= 1, n_mo
     ovmo(imo,imo) = 1d0
   enddo

   h1e = 0d0
   r12mo(:,:,:,:) = 0d0 

!$OMP PARALLEL DO PRIVATE(i,j,k,l,hij)
! SCHEDULE(DYNAMIC) 
    do i = 1, ntdet
      do j = 1, ntdet
!        print*,tdeta(:,i),tdeta(:,j),tdetb(:,i),tdetb(:,j)
        call lowdin(ne,nea,neb,n_mo,ovmo,dcmplx(w1e,0d0),r12mo,ctdet(i),ctdet(j),tdeta(:,i),tdeta(:,j),tdetb(:,i),tdetb(:,j),ov,h1e,r12)
!        print*, j,i,h1e,ov
!        stop
        coll_csf_mat(j,i) = real(h1e)
      enddo
    enddo
!$OMP END PARALLEL DO

    coll_mat(:,:) = 0d0
    CALL DGEMM('N','N',ncsf,nsta,ncsf,1.d0,coll_csf_mat(1:ncsf,1:ncsf),ncsf,eigvec1(1:ncsf,1:nsta),ncsf,0.d0,mattmp,ncsf)
    CALL DGEMM('N','N',nsta,nsta,ncsf,1.d0,transpose(eigvec2(1:ncsf,1:nsta)),nsta,mattmp,ncsf,0.d0,coll_mat(1:nsta,1:nsta),nsta)

    coll_couplings(1:nsta,1:nsta,it) = coll_mat(1:nsta,1:nsta)
   enddo

   call cpu_time(t2)
   print*,t2-t1
   print*,' '

 deallocate(coll_csf_mat,eigval1,eigval2,eigvec1,eigvec2,coll_mat,mattmp) 

 END_PROVIDER

