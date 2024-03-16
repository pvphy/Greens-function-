module  global
    integer(kind=16)::i,i1,j,j1,k,p,get_fact,n,r,ncr_i,temp,value1,k1,dr,dj,p1,p2,p3,flag,a,v_l1,v_r1,n_r,n_l,dim_rpob,q,s,x,y,y_1
    integer(kind=16)::val2,val1,u_mult,nos,nop
    integer,allocatable::b(:),c(:),c_1(:),b_comp(:,:)
    integer(kind=16)::counter,l,l3,l4,count1,count2,count3
    double precision::u,eta,spec_f,pi,t_s,t_f,t1_s,t1_f,t,v
    double precision,allocatable,dimension(:)::v1,v2
    complex*16::trace,trace1
    complex*16,allocatable::d_mat(:,:),tau01(:,:),tau10(:,:),work(:),matrix(:,:),I_mat(:,:)
    integer(kind=16),allocatable,dimension(:)::IPIV  
    !complex*16,allocatable::a11(:,:),a12(:,:),a21(:,:),a22(:,:)
    complex*16,target,allocatable::d_eff(:,:)
    !complex*16,allocatable,dimension(:,:):: a1_ia3_a2,ia3,b0,b1,b2,b3,a1_ia3,b2_a1_ia3,a2_b0

endmodule global

program main
    use global
    implicit none
    integer(kind=16)::l2,nwp,m5,m4
    double precision::w,dw,wmin,wmax,wmax1,inpt,t21,t22,t11,t12,t31
    complex*16,allocatable::a1(:,:),a2(:,:),a0(:,:),a3(:,:),new_D(:,:)
    complex*16,pointer::a22(:,:)
    complex*16,allocatable,dimension(:,:):: a1_ia3_a2,ia3,b0,b1,b2,b3,a1_ia3,b2_a1_ia3,a2_b0
    complex*16,allocatable::a_1(:,:),a_2(:,:),a_0(:,:),a_3(:,:)
    complex*16,allocatable::b_1(:,:),b_2(:,:),b_0(:,:),b_3(:,:),ib0(:,:)

    open(8,file='input.dat',status='unknown')
    do i=1,9
        read(8,*)inpt
        if (i.eq.1) nos=int(inpt)
        if (i.eq.2) nop=int(inpt)
        if (i.eq.3) t=dble(inpt)
        if (i.eq.4) u=dble(inpt)
        if (i.eq.5) eta=dble(inpt)
        if (i.eq.6) wmin=dble(inpt)
        if (i.eq.7) wmax1=dble(inpt)
        if (i.eq.8) nwp=int(inpt)
        if (i.eq.9) v=dble(inpt)
        
    end do
    close(8)
 
    
    print*,'no of site______=',nos
    print*,'no of particle__=',nop
    print*,'t_______________=',t
    print*,'u_______________=',u
    print*,'eta_____________=',eta
    print*,'wmin____________=',wmin
    print*,'wmax____________=',wmax1+u*(nop-1)
    print*,'no of w points__=',nwp
    print*,'V_______________=',v
    

    allocate(b(nos),c(nos),v1(nos))
	
    pi=4.0*atan(1.0)

    do i=1,nos
        if(i.le.(nos/2)) v1(i)=dble((nos/2)-i)
        if(i.gt.(nos/2)) v1(i)=dble(i-((nos/2)+1))
    
    enddo

    call diag_pob
    m5=p2/2
    print*,m5
    allocate(a1_ia3_a2(m5,m5),ia3(m5,m5),ib0(m5,m5),b0(m5,m5),b1(m5,m5),b2(m5,m5),b3(m5,m5)&
    &,a1_ia3(m5,m5),b2_a1_ia3(m5,m5),a2_b0(m5,m5))
    allocate(a0(m5,m5),a1(m5,m5),a2(m5,m5),a3(m5,m5))
    allocate(a_0(m5/2,m5/2),a_1(m5/2,m5/2),a_2(m5/2,m5/2),a_3(m5/2,m5/2))
    allocate(b_0(m5/2,m5/2),b_1(m5/2,m5/2),b_2(m5/2,m5/2),b_3(m5/2,m5/2))

    allocate(d_mat(p2,p2),I_mat(p2,p2)) 
    allocate(d_eff(p2,p2))  
    call matgen(p2)
    call I_matrix
    wmax=wmax1+u*(nop-1)
    w=wmin
    dw=abs(wmax-wmin)/dble(nwp)
    do l2=1,nwp
        call cpu_time(t1_s)
        w=w+dw
        
        d_eff=((complex(w,eta)*(I_mat))-d_mat)



        a0=0.0d0;a1=0.0d0;a2=0.0d0;a3=0.0d0
        b0=0.0d0;b1=0.0d0;b2=0.0d0;b3=0.0d0
        a1_ia3=0.0d0
        a2_b0=0.0d0
        a1_ia3_a2=0.0d0
        b2_a1_ia3=0.0d0
        do i=1,m5
            do j=1,m5
                a0(i,j)=d_eff(i,j)
            enddo
            do j=m5+1,p2
                a1(i,j-m5)=d_eff(i,j)
            enddo
        enddo
        do i=m5+1,p2
            do j=1,m5
                a2(i-m5,j)=d_eff(i,j)
            enddo
            do j=m5+1,p2
                a3(i-m5,j-m5)=d_eff(i,j)
            enddo
        enddo
        m4=m5/2


        call cpu_time(t11)
        call euler_if(a3,m5,m4,a_0,a_1,a_2,a_3)
   
        !call inverse(m5,a3)
        call cpu_time(t12)
        write(20,*) t12-t11
        flush(20)

        do i=1,m4
            do j=1,m4
                ia3(i,j)=a_0(i,j)
            enddo
            do j=m4+1,m5
                ia3(i,j)=a_1(i,j-m4)
            enddo
        enddo
        do i=m4+1,m5
            do j=1,m4
                ia3(i,j)=a_2(i-m4,j)
            enddo
            do j=m4+1,m5
                ia3(i,j)=a_3(i-m4,j-m4)
            enddo
        enddo
        !********
        
        a1_ia3=matmul(a1,ia3)
        a1_ia3_a2=matmul(a1_ia3,a2)
        b0=a0-a1_ia3_a2

        call cpu_time(t21)

        call euler_if(b0,m5,m4,b_0,b_1,b_2,b_3)
        call cpu_time(t22)
        write(21,*) t22-t21
        flush(21)
        
        do i=1,m4
            do j=1,m4
                ib0(i,j)=b_0(i,j)
            enddo
            do j=m4+1,m5
                ib0(i,j)=b_1(i,j-m4)
            enddo
        enddo
        do i=m4+1,m5
            do j=1,m4
                ib0(i,j)=b_2(i-m4,j)
            enddo
            do j=m4+1,m5
                ib0(i,j)=b_3(i-m4,j-m4)
            enddo
        enddo
        !call inverse(m5,b0)


        !****
        a2_b0=matmul(a2,ib0)
        b2=-matmul(ia3,a2_b0)
        !****

        b1=-matmul(ib0,a1_ia3)
        
        !******
        b2_a1_ia3=matmul(b2,a1_ia3)
        b3=ia3-b2_a1_ia3

        call cpu_time(t31)
        write(22,*) t31-t1_s
        flush(22)
       

        trace=complex(0.0d0,0.0d0)
        do l3=1,m5  
            trace=trace+ib0(l3,l3)+b3(l3,l3)      
        enddo
       
        write(900,*)w,-imag(trace)/pi
        flush(900)
        d_eff=0.0d0
       call cpu_time(t1_f)  
       write(23,*)t1_f-t1_s,'sec'
       flush(23)
       
   enddo   
      
endprogram main

subroutine diag_pob
    use global
    implicit none 
    integer(kind=16)::m11
    call ncr(nos,nop,p2)

    m11=p2

    allocate(b_comp(m11,nos))

    do value1=0,m11-1
        r=nop
        call binary_gen
        do k=1,nos
            b_comp(value1+1,k)=c(k)
            
        enddo

    enddo

    ! do i=1,m11
	! 	write(104,*)(b_comp(i,j),j=1,nos)
	! enddo
    ! stop
endsubroutine diag_pob

subroutine matgen(dim1)
    use global
    implicit none
    integer(kind=16)::i5,dim1,k2,b_1(nos),i4
    double precision:: i7
        
    d_mat=complex(0.0d0,0.0d0)
    do k2=1,nos-1
    do i=1,dim1 
        if((b_comp(i,k2).eq.0).and.(b_comp(i,k2+1).eq.1))then
            do l=1,nos
                b_1(l)=b_comp(i,l)
            enddo
            !print*,b_1
            b_1(k2)=1
            b_1(k2+1)=0
            !print*,b_1
        do j=1,dim1
            if((b_comp(j,k2).eq.1).and.(b_comp(j,k2+1).eq.0))then
            i5=0
            do k=1,nos
                if(b_1(k).eq.b_comp(j,k))then
                    i5=i5+1
                endif
            enddo
            if (i5==(nos))  d_mat(i,j)=t   
            
            endif
        enddo
        endif
    enddo
    enddo

    do k2=1,nos-1
        do i=1,dim1 
            if((b_comp(i,k2).eq.1).and.(b_comp(i,k2+1).eq.0))then
                do l=1,nos
                    b_1(l)=b_comp(i,l)
                enddo
                !print*,b_1
                b_1(k2)=0
                b_1(k2+1)=1
                
            do j=1,dim1
                if((b_comp(j,k2).eq.0).and.(b_comp(j,k2+1).eq.1))then
                i5=0
                do k=1,nos
                    if(b_1(k).eq.b_comp(j,k))then
                        i5=i5+1
                    endif
                enddo
                if (i5==(nos))  d_mat(i,j)=t   
               
                endif
            enddo
            endif
        enddo
    enddo

    do i=1,dim1
        i4=0
        do k=1,nos-1
                    
            if((b_comp(i,k).eq.1).and.((b_comp(i,k+1).eq.1)))then
                i4=i4+1
            endif
        enddo
        i7=0.0d0
        do k=1,nos   
            i7=i7+v1(k)*dble(b_comp(i,k))
        enddo

        d_mat(i,i)=i4*u+i7*v
    
    enddo




    do i=1,dim1
        do j=1,dim1
            if(d_mat(i,j).ne.d_mat(j,i)) print*,i,j,'not hermitian'
            !print*,i,j,real(d_mat(i,j))
        enddo
    enddo

endsubroutine matgen

subroutine pob(pob2)
    use global
    implicit none
    integer(kind=16)::jc_pj,pj,pob1,pob2,j_1,i2,j2
    pob2=0
    do j_1=0,nos-1
        j2=nos-j_1
       
        pj=0
        ! print*,'j',j
        do i2=j2,nos
            pj=pj+c(i2)
        enddo

        call ncr(j_1,pj,jc_pj)
    
        pob1=c(j2)*jc_Pj
        pob2=pob2+pob1
        ! print*,j_1,pj,c(j),jc_pj,pob1,pob2

    enddo
    !  print*,pob2
endsubroutine pob


subroutine binary_gen
    use global
    implicit none
	! print*,p2
	!do value1=0,(ncr_i-1)	
	!b=0  
    temp=value1
    k=r+1
    b=0
    c=0
    !a=0
    do k1=1,r
        dr=-1
        k=k+dr
        j=nos
        flag=0
        do while (flag.eq.0)
            dj=-1
            j=j+dj
            call ncr(j,k,p1)  
            p=p1
            if((temp.gt.p).or.(temp.eq.p))then
                temp=temp-p
                b(nos-j)=1
                flag=1       
            else
                flag=0
            endif 
       
        enddo 
    enddo 
    c=b

endsubroutine binary_gen

subroutine ncr(n2,r2,ncr2)  
    implicit none 
    integer(kind=16)::get_fact,ncr2,n2,r2
    if(n2.lt.r2)then
        ncr2=0
    else
        ncr2=(get_fact(n2)/((get_fact(n2-r2)*get_fact(r2))))
    endif
endsubroutine ncr


integer(kind=16) function get_fact(inte)
    implicit none 
    integer(kind=16)::fact1,inte,i1
    !if (inte.lt.0) then
    !get_fact=0
    !else
    fact1=1 
    do i1=1,inte
        fact1=fact1*i1
    end do
    get_fact=fact1
    !endif
 
end function get_fact

subroutine inverse(m,mat0)
    use global
    implicit none
    integer::info,error
    integer(kind=16)::m  
    complex*16::mat0(m,m) 
    allocate(WORK(m),IPIV(m),stat=error)
    if (error.ne.0)then
        write(84,*)"error:not enough memory"
        flush(84)
        stop
    end if
    !  print*,shape(d_eff)
    call zgetrf(m,m,mat0,m,IPIV,info)
    if(info.ne.0) then
        write(84,*)"failed"
        flush(84)
        stop
    else
    ! print*,"failed"
    end if
    ! print*,info
    call zgetri(m,mat0,m,IPIV,WORK,m,info)
    if(info .ne. 0) then
        write(84,*)"failed"
        flush(84)
        stop
    end if
    deallocate(IPIV,WORK,stat=error)
    if (error.ne.0)then
        write(84,*)"error:fail to release"
        flush(84)
        stop
    end if

endsubroutine inverse

subroutine I_matrix
    use global 
    implicit none
    I_mat=complex(0.0d0,0.0d0)
    do i=1,p2
       I_mat(i,i)=complex(1.0d0,0.0d0)
    enddo
endsubroutine I_matrix





subroutine euler_if(d_pass,p_1,m5,b0,b1,b2,b3)
    use global
    implicit none
    integer(kind=16)::m5,p_1
    double precision:: t22,t21,t11,t12
    complex*16::a00(m5,m5),a11(m5,m5),a22(m5,m5),a33(m5,m5)
    complex*16::a1_ia3_a2(m5,m5),ia3(m5,m5),b0(m5,m5),b1(m5,m5),b2(m5,m5),b3(m5,m5)&
    &,a1_ia3(m5,m5),b2_a1_ia3(m5,m5),a2_b0(m5,m5),d_pass(p_1,p_1)
    !complex*16::a0(m5,m5),a1(m5,m5),a2(m5,m5),a3(m5,m5)



    do i=1,m5
        do j=1,m5
            a00(i,j)=d_pass(i,j)
        enddo
        do j=m5+1,p_1
            a11(i,j-m5)=d_pass(i,j)
        enddo
    enddo
    do i=m5+1,p_1
        do j=1,m5
            a22(i-m5,j)=d_pass(i,j)
        enddo
        do j=m5+1,p_1
            a33(i-m5,j-m5)=d_pass(i,j)
            
        enddo
    enddo

    b0=0.0d0;b1=0.0d0;b2=0.0d0;b3=0.0d0
    a1_ia3=0.0d0
    a2_b0=0.0d0
    a1_ia3_a2=0.0d0
    b2_a1_ia3=0.0d0
    call cpu_time(t11)
    call inverse(m5,a33)
    call cpu_time(t12)
    write(20,*) t12-t11
    flush(20)
    !********
    ia3=a33
    a1_ia3=matmul(a11,ia3)
    a1_ia3_a2=matmul(a1_ia3,a22)
    b0=a00-a1_ia3_a2
    call cpu_time(t21)
    call inverse(m5,b0)
    call cpu_time(t22)
    write(21,*) t22-t21
    flush(21)

    !****
    a2_b0=matmul(a22,b0)
    b2=-matmul(ia3,a2_b0)
    !****

    b1=-matmul(b0,a1_ia3)
    
    !******
    b2_a1_ia3=matmul(b2,a1_ia3)
    b3=ia3-b2_a1_ia3
    
    
endsubroutine euler_if

