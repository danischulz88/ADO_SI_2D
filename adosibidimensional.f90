program ado_bidimensional
implicit none

integer:: n, kk, jj, nado, nit, l,ll, ndx, ndy,  i, j, k, prop
double precision, dimension (:,:), allocatable::  s
double precision::  hx, hy, compx, compy, tol
double precision, dimension(:,:,:), allocatable :: psi_medioxy, psi_medioy,psi_mediox
double precision, dimension (:,:), allocatable ::sigma_t, sigma_s0, sigma_s1, q, R
double precision, dimension(:), allocatable:: as, bs
double precision, dimension(882):: dadoscamila
! variaveis para medicao de tempo de execucao
double precision :: t
integer, dimension(8) :: t0,t1
!variaveis de automatizacao
integer, dimension(6)::direcoes, discret
    !prop eh o numero de vezes que o numero de pontos fornecidos pelo ado forcado  cabe na malha
    !Varredura 2D

    ! Considera o particionamento de uma regiao de comprimento LxM em J celulas ao longo do eixo X,
    !com J+1 paredes;e K celulas no eixo Y com K+1 paredes;fluxos angulares sao calculados em cada
    !parede (logo, psi tem N direcoes em J+1 posicoes em X e K+1 posicoes em Y);fluxos angulares
    !medios sao calculados no meio de cada celula (logo,psi_medio tem N direcoes em J posicoes em
    !X e K posicoes em Y); o tamanho de cada celula pode ser constante ou nao.
    open(unit=10,file='fluxbidi.txt')
    open(unit=20,file='regioessi.txt')
    open(unit=30,file='regioessiado.txt')
    open(unit=40,file='dadosproblema3.dat',status='old')
    100 format(' ',1I20, 1I20, 1E20.6)
    !inicializacoes do problema
    do i=1, 882
      read(40, *)dadoscamila(i)
    end do
    ! print *,'Digite a ordem da quadratura (ex. S_4 digite 4)'
    !read *, n
    direcoes(1)=2
    direcoes(2)=4
    direcoes(3)=6
    direcoes(4)=8
    direcoes(5)=12
    direcoes(6)=16
!$$$$$$     discret(1)=2
!$$$$$$     discret(2)=4
!$$$$$$     discret(3)=8
!$$$$$$     discret(4)=10
!$$$$$$     discret(5)=50
!$$$$$$     discret(6)=100
! $$$$$$     discret(1)=3
! $$$$$$     discret(2)=6
! $$$$$$     discret(3)=15
! $$$$$$     discret(4)=30
! $$$$$$     discret(5)=60
! $$$$$$     discret(6)=120
    ! discret(1)=50
    ! discret(2)=100
    ! discret(3)=150
    ! discret(4)=200
    ! discret(5)=250
    ! discret(6)=300
    discret(1)=21
    discret(2)=42
    discret(3)=63
    discret(4)=84
    discret(5)=105
    discret(6)=126
    ndx=2!numero de regioes no eixo x
    ndy=2!numero de regioes no eixo y

	nado=4!numero de pontos em x e/ou y para as direcoes

    allocate(R(ndx, ndy), as(ndx), bs(ndy))
    do l=1, 6
      do ll=1, 6

            n=direcoes(l)!numero de direcoes para a varredura
            ! definicoes das constantes do problema
            compx=30.d0;compy=30.d0; jj=discret(ll);kk=jj; hx=compx/jj; hy=compy/kk;tol=1.0d-6
             !compx=1.d0;compy=1.d0; jj=10;kk=jj; hx=compx/jj; hy=compy/kk;tol=1.0d-6
            print*, 'n ; jxk', n, jj,'x',kk
            prop=jj/21

            allocate(psi_medioxy(n*(n+2)/2,jj, kk))
            allocate(psi_mediox(n*(n+2)/2, jj, kk+1))
            allocate(psi_medioy(n*(n+2)/2, jj+1, kk))
            allocate(sigma_t(jj,kk))
            allocate(sigma_s0(jj,kk))
            allocate(sigma_s1(jj,kk))
            allocate(s(jj,kk),q(jj,kk) )
            !pontos de divisao das regioes para obter fluxos medios
            as(1)=10.d0
            bs(1)=10.d0
            as(2)=30.d0
            bs(2)=30.d0
!$$$$$$             as(3)=50.d0
!$$$$$$             bs(3)=50.d0
            !parametros gerais do programa
			sigma_t=2.d0
			sigma_s0=0.1d0
            sigma_s1=0.0d0
            !parametros na regiao da fonte
            sigma_t(1:nint(as(1)/hx), 1:nint(bs(1)/hy))=1.0d0
            sigma_s0(1:nint(as(1)/hx), 1:nint(bs(1)/hy))=0.5d0
            sigma_s1(1:nint(as(1)/hx), 1:nint(bs(1)/hy))=0.0d0

			!inicializacao da fonte
            q=0.
            q(1:nint(as(1)/hx), 1:nint(bs(1)/hy)) = 1.d0 !magnitude da fonte em [0,as]x[0,bs]
            !inicializacao dos parametros de entrada
            psi_mediox = 0.d0
            psi_medioy = 0.d0
            psi_medioxy= 0.d0
            s=0.d0


            ! inicializa psi na fronteira
            psi_mediox(n*(n+2)/4+1:n*(n+2)/2,1:jj, kk+1) = 0.d0
            psi_medioy(n*(n+2)/4+1:n*(n+2)/2,jj+1,1: kk) = 0.d0


            ! chama SI
            call date_and_time(values=t0)
            call si_bi(n, jj, kk, hx, hy, tol, sigma_t, sigma_s0,sigma_s1, q, psi_mediox, psi_medioy, psi_medioxy, s, nit, as, bs,&
                         ndx, ndy,R)
            call date_and_time(values=t1)
            t = diferenca_tempo(t0,t1)
            print *,'numero de iteracoes: nit=', nit
            print *,'SI: Tempo de execucao: t=',t,' s'
            do i=1, ndx
              do j=1, ndy
                print *,'SI:i, j,  fluxos nas regioes R(i,j)=',i, j, R(i,j)
              enddo
            end do
             write(20, *)n, jj, nit
             do i=1, ndx
              do j=i , ndy
                write(20, 100)i, j, R(i,j)
              enddo
           end do


!$$$$$$             ! inicializa psi na fronteira
!$$$$$$             psi_mediox(n*(n+2)/4+1:n*(n+2)/2,1:jj, kk+1) = 0.d0
!$$$$$$             psi_medioy(n*(n+2)/4+1:n*(n+2)/2,jj+1,1: kk) = 0.d0
!$$$$$$             s=0.d0
!$$$$$$
!$$$$$$             !Chama SI+ado
!$$$$$$             call date_and_time(values=t0)
!$$$$$$             call ado_bi(nado, jj, kk, sigma_t(1,1), sigma_s0(1,1), compx,compy, as, bs, hx, hy,q(1,1), s, R1, R2, R3, R4)
!$$$$$$             print *,' fluxos nas regioes R1=', R1, 'R2=', R2, 'R3=', R3, 'R4=', R4
!$$$$$$             !n= numero de divisoes por eixo
!$$$$$$             !M= numero de direcoes totais
!$$$$$$             !j=M/2
!$$$$$$             !jj= numero de celulas na discretizacao espacial x
!$$$$$$             !kk= numero de celulas na discretizacao espacial y
!$$$$$$             call si_bi(n, jj, kk, hx, hy, tol, sigma_t, sigma_s0,sigma_s1, q, psi_mediox, psi_medioy, psi_medioxy, s, nit, as, bs,&
!$$$$$$                         R1, R2, R3, R4)
!$$$$$$             call date_and_time(values=t1)
!$$$$$$             t = diferenca_tempo(t0,t1)
!$$$$$$             print *,'numero de iteracoes: nit=', nit
!$$$$$$             print *,'SI+ADO: Tempo de execucao: t=',t,' s'
!$$$$$$             print *,' fluxos nas regioes R1=', R1, 'R2=', R2, 'R3=', R3, 'R4=', R4

!$$$$$$             ! inicializa psi na fronteira
!$$$$$$             psi_mediox(n*(n+2)/4+1:n*(n+2)/2,1:jj, kk+1) = 0.d0
!$$$$$$             psi_medioy(n*(n+2)/4+1:n*(n+2)/2,jj+1,1: kk) = 0.d0
!$$$$$$             s=0.d0
!$$$$$$
!$$$$$$             ! chama ado
!$$$$$$             call date_and_time(values=t0)
!$$$$$$             call ado_bi(n, jj, kk, sigma_t(1,1), sigma_s0(1,1), compx,compy, as, bs, hx, hy,q(1,1), s, R1, R2, R3, R4)
!$$$$$$             call date_and_time(values=t1)
!$$$$$$             t = diferenca_tempo(t0,t1)
!$$$$$$             !print *,'ADO: Tempo de execucao: t=',t,' s'
!$$$$$$             !print *,' fluxos nas regioes R1=', R1, 'R2=', R2, 'R3=', R3, 'R4=', R4
!$$$$$$
!$$$$$$             ! inicializa psi na fronteira

            !Chama SI+ado(forcado
            !inicializacao dos parametros
            psi_mediox = 0.d0
            psi_medioy = 0.d0
            psi_medioxy= 0.d0
            s=0.d0
            call date_and_time(values=t0)
            do i=1, 21
              do j=1, 21
                k=2*j+42*(i-1)-1
                S((i-1)*prop+1:prop*i, (j-1)*prop+1:prop*j)=(dadoscamila(k)+dadoscamila(k+1))/2
              end do
            end do

            call si_bi(n, jj, kk, hx, hy, tol, sigma_t, sigma_s0,sigma_s1, q, psi_mediox, psi_medioy, psi_medioxy, s, nit, as, bs,&
                        ndx, ndy,R)
            call date_and_time(values=t1)
            t = diferenca_tempo(t0,t1)
            print *,'SI+ADO(forcado):numero de iteracoes: nit=', nit
           ! print *,'SI+ADO(forcado): Tempo de execucao: t=',t,' s'
            do i=1, ndx
              do j=1, ndy
                print *,'SI+ADO(forcado):i, j,  fluxos nas regioes R(i,j)=',i, j, R(i,j)
              enddo
            end do
             write(30, *)n, jj, nit
             do i=1, ndx
              do j=1, ndy
                write(30, 100)i, j, R(i,j)
              enddo
           end do
           deallocate (psi_medioxy, psi_mediox,psi_medioy,sigma_t, sigma_s0,sigma_s1, s,q)
        enddo
   end do

    contains

        !=====================================================================================================
        subroutine ado_bi(n, jj, kk, sigma_t, sigma_s, a, b, as, bs, hx, hy,q, s, R)
        !=====================================================================================================
        implicit none
        integer, intent(in)::n, jj, kk
        integer:: M, j
        double precision, intent(in):: sigma_s, sigma_t, a, b, hx, hy, q
        double precision, dimension(:), intent(in):: as, bs
        double precision:: aas, bbs
        double precision, dimension(:,:), intent(out)::s, R
        double precision, dimension(:), allocatable:: peso,Omegamu, Omegaeta, ni, gama, VetSol
        double precision, dimension (:,:), allocatable::  U, V, phiy, U2, V2, phix

        	M=(n*(n+2))/2!numero total de direcoes
    		j=M/2 !numero de direcoes em cada sentido
            open(unit=10,file='adobi.txt')
            allocate(peso(M), Omegamu(j), Omegaeta(j), U(j,j), V(j,j), phiy(M, j), U2(j,j), V2(j,j), phix(M,j), &
                        ni(j), gama(j), VetSol(8*M))
            aas=as(1)
            bbs=bs(1)
            call nosepesosado(n,M, Omegamu, Omegaeta, peso)

            call eigens(j, Omegamu, sigma_t, sigma_s, peso, U, V,ni)

            call phis(M, j, U, V, phiy)

            call eigens(j, Omegaeta, sigma_t, sigma_s, peso, U2, V2, gama)

            call phis(M, j, U2, V2, phix)

            Call matsist(M,j, Omegamu, a, phiy, aas, ni, Omegaeta, b, phix, bbs, gama, sigma_s, sigma_t,q, peso, VetSol)

            call solutions(M, j,kk,jj, phiy, ni, a, phix, gama, b, aas, bbs, peso,VetSol, hx, hy, s)
            call fluxo_regioes(jj,kk, hx, hy,ndx, ndy,as, bs, s, R)

        end subroutine ado_bi

        !=================================================================================================
        subroutine si_bi(n, j, k, hx, hy, tol, sigma_t, sigma_s0,sigma_s1, q, psi_mediox, psi_medioy, psi_medioxy, s, nit,&
        					as, bs,ndx, ndy,R)
        !=================================================================================================
        implicit none
        integer, intent(in)::n, j, k, ndx, ndy
        double precision::erro
        double precision, intent(in)::tol, hx, hy
        double precision, dimension(:), intent(in)::as, bs
        double precision, dimension(:,:), intent(in)::sigma_t, sigma_s0,sigma_s1, q
        double precision, dimension(:,:,:), intent(inout)::psi_mediox, psi_medioy, psi_medioxy
        double precision, dimension(:,:), intent(inout):: s
        integer, intent(out):: nit
        double precision, dimension(:,:), intent(out):: R
        double precision, dimension(:), allocatable::mu, peso
        double precision, dimension(:,:), allocatable:: Omega, J1, J2
!        integer:: m, y, x
    	!print*, 'n/2, (n*(n+2))/2,(n*(n+2))/8, j, k',n/2, (n*(n+2))/2,(n*(n+2))/8, j, k
        allocate (mu(n/2),peso((n*(n+2))/2),Omega((n*(n+2))/8,2), J1(j,k), J2(j,k))

    	J1=0.d0
        J2=0.d0

        call nosepesossi( n, peso, omega)
            do nit=1,5000
               	call  varredura(n, j, k, hx, hy, Omega, sigma_t, sigma_s0,sigma_s1, s,J1, J2, q, psi_mediox, psi_medioy,&
                	 			psi_medioxy)

!$$$$$$                 do m=1, n*(n+2)/2
!$$$$$$                     do y=1, k
!$$$$$$                         do x=1, j
!$$$$$$                             write(10, *)m, x, y, psi_medioxy(m, x, y)
!$$$$$$                         enddo
!$$$$$$                     end do
!$$$$$$                 end do


               	call calcula_erro(n, j, k, psi_medioxy, erro,peso, s)
                !$$$$$$       print*, 'nit, fluxo', nit, s(1,1)
                !$$$$$$     pause
                if (erro<tol) then
                    exit
                end if

            	call corrente(n, j, k, peso, Omega, psi_medioxy, J1, J2)
            end do

            call fluxo_regioes(j,k, hx, hy,ndx, ndy,as, bs, s,R)

        end subroutine si_bi

        !=================================================================================================
        subroutine nosepesossi( n, peso, omega)
        !=================================================================================================
        implicit none
        integer, intent(in)::n
        double precision, dimension(:), intent(out) :: peso
        double precision, dimension(:), allocatable :: w, mu
        double precision, dimension(:,:), intent(out):: Omega
        integer :: i, l, x
       		allocate (mu(n/2))
            do
                select case (n)
                    case (2)
                        allocate (w(1))
                        mu(1)=0.5773502; w(1)=1.d0; peso=w(1)
                        exit
                    case (4)
                        allocate (w(1))
                        mu(1)=0.3500212; mu(2)=0.8688903; w(1)=0.3333333; peso=w(1)
                        exit
                    case (6)
                        allocate (w(2))
                        mu(1)=0.2666355; mu(2)=0.6815076; mu(3)=0.9261808; w(1)=0.1761263; w(2)=0.1572071
                        peso(1)=w(1); peso(2)=w(2);peso(3)=w(2);peso(4)=w(1);peso(5)=w(2);peso(6)=w(1);
                        do i=1,6
                        peso(n*(n+2)/8+i)=peso(i)
                        peso(n*(n+2)/4+i)=peso(i)
                        peso(3*n*(n+2)/8+i)=peso(i)
                        end do
                        exit
                    case (8)
                        allocate (w(3))
                        mu(1)=0.2182179; mu(2)=0.5773503; mu(3)=0.7867958; mu(4)=0.9511897; w(1)=0.1209877;
                        w(2)=0.0907407; w(3)=0.0925926
                        peso(1)=w(1);peso(2)=w(2);peso(3)=w(2);peso(4)=w(2);peso(5)=w(3);peso(6)=w(2);
                        peso(7)=w(1);peso(8)=w(2);peso(9)=w(2); peso(10)=w(1);
                        do i=1,10
                        peso(n*(n+2)/8+i)=peso(i)
                        peso(n*(n+2)/4+i)=peso(i)
                        peso(3*n*(n+2)/8+i)=peso(i)
                        end do
                        exit
                    case (12)
                        allocate (w(5))
                        mu(1)=0.1672126; mu(2)=0.4595476; mu(3)=0.6280191; mu(4)=0.7600210; mu(5)=0.8722706;
                        mu(6)=0.9716377; w(1)=0.0707626; w(2)=0.0558811; w(3)=0.0373377; w(4)=0.0502819;
                        w(5)=0.0258513
                        peso(1)=w(1);peso(16)=w(1);peso(21)=w(1); peso(2)=w(2);peso(11)=w(2);peso(15)=w(2);
                        peso(17)=w(2);peso(20)=w(2);peso(3)=w(2); peso(4)=w(3);peso(6)=w(3);peso(7)=w(3);
                        peso(10)=w(3);peso(18)=w(3);peso(19)=w(3); peso(5)=w(4);peso(12)=w(4);peso(14)=w(4);
                        peso(8)=w(5);peso(9)=w(5);peso(13)=w(5);
                        do i=1,21
                        peso(n*(n+2)/8+i)=peso(i)
                        peso(n*(n+2)/4+i)=peso(i)
                        peso(3*n*(n+2)/8+i)=peso(i)
                        end do
                        exit
                    case (16)
                        allocate (w(8))
                        mu(1)=0.1389568; mu(2)=0.3922893; mu(3)=0.5370966; mu(4)=0.6504264; mu(5)=0.7467506;
                        mu(6)=0.8319966;mu(7)=0.9092855;mu(8)=0.9805009; w(1)=0.0489872; w(2)=0.0413296;
                        w(3)=0.0212326; w(4)=0.0256207; w(5)=0.0360486;w(6)=0.0144589;w(7)=0.0344958; w(8)=0.0085179
                        peso(1)=w(1);peso(2)=w(2);peso(3)=w(2);peso(4)=w(3);peso(5)=w(5);peso(6)=w(3);
                        peso(7)=w(4);peso(8)=w(6);peso(9)=w(6);peso(10)=w(4);peso(11)=w(4);peso(12)=w(7);
                        peso(13)=w(8);peso(14)=w(7);peso(15)=w(4);peso(16)=w(3);peso(17)=w(6);peso(18)=w(8);
                        peso(19)=w(8);peso(20)=w(6);peso(21)=w(3);peso(22)=w(2);peso(23)=w(5);peso(24)=w(6);
                        peso(25)=w(7);peso(26)=w(6);peso(27)=w(5);peso(28)=w(2);peso(29)=w(1);peso(30)=w(2);
                        peso(31)=w(3);peso(32)=w(4);peso(33)=w(4);peso(34)=w(3);peso(35)=w(2);peso(36)=w(1);
                        do i=1,36
                        peso(n*(n+2)/8+i)=peso(i)
                        peso(n*(n+2)/4+i)=peso(i)
                        peso(3*n*(n+2)/8+i)=peso(i)
                        end do
                        exit
                    end select
                end do

                !preenche o vetor de posicoes mu e eta no primeiro quadrante
                l=0
                do i=n/2,1, -1
                    do x=1, n/2-i+1
                        l=l+1
                        Omega(l,1)=mu(i)
                        Omega(l,2)=mu(x)
                        !print*, 'omega', Omega(l,1), Omega(l,2)
                    end do
            end do

        end subroutine nosepesossi

        !=================================================================================================
        subroutine varredura(n, j, k,hx, hy, Omega, sigma_t, sigma_s0, sigma_s1, s,J1, J2, q, psi_mediox, psi_medioy, psi_medioxy)
        !=================================================================================================
        implicit none
        integer, intent(in)::n, j, k
        double precision, intent(in):: hx, hy
        double precision, dimension(:,:), intent(in):: Omega, sigma_t,sigma_s0,sigma_s1, s, J1, J2, q
        double precision, dimension(:,:,:), intent(inout) :: psi_medioy,psi_mediox
        double precision, dimension(:,:,:), intent(inout) :: psi_medioxy
        integer :: y, x, m

            !varredura do terceiro quadrante
            do y=k,1,-1
                do x=j,1,-1
                    do m=1,n*(n+2)/8
                        !print*, 'psis contorno',x,y,
                        !$$$$$$     psi_medioxy(n*(n+2)/4+m,x,y)=((2*Omega(m,1)/hx)*psi_medioy(n*(n+2)/4+m,x+1,y)+&
                        !$$$$$$              (2*Omega(m,2)/hy)*psi_mediox(n*(n+2)/4+m, x,y+1)+(sigma_s0(x,y)/4)*s(x,y)+q(x,y))/&
                        !$$$$$$                     (sigma_t(x,y)+(2*Omega(m,1)/hx)+(2*Omega(m,2)/hy))
                        psi_medioxy(n*(n+2)/4+m,x,y)=((2*Omega(m,1)/hx)*psi_medioy(n*(n+2)/4+m,x+1,y)+&
                                                     (2*Omega(m,2)/hy)*psi_mediox(n*(n+2)/4+m, x,y+1)+sigma_s0(x,y)*s(x,y)+&
                                                     0.75*sigma_s1(x,y)*(-Omega(m,1)*J1(x,y)-Omega(m,2)*J2(x,y))+q(x,y))/&
                                                     (sigma_t(x,y)+(2*Omega(m,1)/hx)+(2*Omega(m,2)/hy))
                        psi_medioy(n*(n+2)/4+m,x,y)=2*psi_medioxy(n*(n+2)/4+m,x,y)-psi_medioy(n*(n+2)/4+m,x+1,y)
                        psi_mediox(n*(n+2)/4+m,x,y)=2*psi_medioxy(n*(n+2)/4+m,x,y)-psi_mediox(n*(n+2)/4+m,x,y+1)
                        !$$$$$$             print*, 'psi_medio 3q', n*(n+2)/4+m,x,y,psi_medioxy(n*(n+2)/4+m, x, y)
                        psi_medioy(3*n*(n+2)/8+m,1, y) = psi_medioy(n*(n+2)/4+m,1,y)
                        psi_mediox(n*(n+2)/8+m,x,1) = psi_mediox(n*(n+2)/4+m,x,1)
                    end do
                end do
            end do

            !varredura do quarto quadrante
            do y=k,1,-1
                do x=1,j
                    do m=1,n*(n+2)/8
                        psi_medioxy(3*n*(n+2)/8+m,x,y)=((2*Omega(m,1)/hx)*psi_medioy(3*n*(n+2)/8+m,x,y)+&
                                                       (2*Omega(m,2)/hy)*psi_mediox(3*n*(n+2)/8+m, x,y+1)+&
                                                       sigma_s0(x,y)*s(x,y)+0.75*sigma_s1(x,y)*(Omega(m,1)*J1(x,y)&
                                                       -Omega(m,2)*J2(x,y))+q(x,y))/&
                                                       (sigma_t(x,y)+(2*Omega(m,1)/hx)+(2*Omega(m,2)/hy))
                        psi_medioy(3*n*(n+2)/8+m,x+1,y)=2*psi_medioxy(3*n*(n+2)/8+m,x,y)-psi_medioy(3*n*(n+2)/8+m,x,y)
                        psi_mediox(3*n*(n+2)/8+m,x,y)=2*psi_medioxy(3*n*(n+2)/8+m,x,y)-psi_mediox(3*n*(n+2)/8+m,x,y+1)
                        !$$$$$$             print*, 'psi_medio 4q', 3*n*(n+2)/8+m,x,y,psi_medioxy(3*n*(n+2)/8+m, x, y)
                        psi_mediox(m,x,1) = psi_mediox(3*n*(n+2)/8+m,x,1)
                    end do
                end do
            end do

            !varredura do segundo quadrante
            do y=1,k
                do x=j,1,-1
                    do m=1,n*(n+2)/8
                    psi_medioxy(n*(n+2)/8+m,x,y)=((2*Omega(m,1)/hx)*psi_medioy(n*(n+2)/8+m, x+1,y)+&
                                                 (2*Omega(m,2)/hy)*psi_mediox(n*(n+2)/8+m, x,y)+sigma_s0(x,y)*s(x,y)+&
                                                 0.75*sigma_s1(x,y)*(-Omega(m,1)*J1(x,y)+Omega(m,2)*J2(x,y))+q(x,y))/&
                                                 (sigma_t(x,y)+(2*Omega(m,1)/hx)+(2*Omega(m,2)/hy))
                    psi_medioy(n*(n+2)/8+m,x,y)=2*psi_medioxy(n*(n+2)/8+m,x,y)-psi_medioy(n*(n+2)/8+m,x+1,y)
                    psi_mediox(n*(n+2)/8+m,x,y+1)=2*psi_medioxy(n*(n+2)/8+m,x,y)-psi_mediox(n*(n+2)/8+m,x,y)
                    !$$$$$$            print*, 'psi_medio 2q', n*(n+2)/8+m,x,y,psi_medioxy(n*(n+2)/8+m, x, y)
                    psi_medioy(m,1, y) = psi_medioy(n*(n+2)/8+m,1,y)
                    end do
                end do
            end do

            !varredura no primeiro quadrante
            do y=1,k
                do x=1,j
                    do m=1,n*(n+2)/8
                    psi_medioxy(m,x,y)=((2*Omega(m,1)/hx)*psi_medioy(m,x,y)+(2*Omega(m,2)/hy)*psi_mediox(m, x,y)+&
                                        sigma_s0(x,y)*s(x,y)+0.75*sigma_s1(x,y)*(Omega(m,1)*J1(x,y)+Omega(m,2)*J2(x,y))&
                                        +q(x,y))/(sigma_t(x,y)+(2*Omega(m,1)/hx)+(2*Omega(m,2)/hy))
                    psi_medioy(m,x+1,y)=2*psi_medioxy(m,x,y)-psi_medioy(m,x,y)
                    psi_mediox(m,x,y+1)=2*psi_medioxy(m,x,y)-psi_mediox(m,x,y)
                    !$$$$$$             print*, 'psi_medioxy 1q', m,x,y,psi_medioxy(m, x, y)
                    end do
                end do
            end do
        end subroutine varredura
        !==============================================================================================================
        subroutine calcula_erro(n, j, k, psi_medioxy, erro,peso, s)
        !==============================================================================================================
        implicit none
        integer, intent(in):: k, j, n
        double precision, intent(out) :: erro
        double precision, dimension(:), intent(in):: peso
        double precision, dimension(:,:, :), intent(in) ::psi_medioxy
        double precision, dimension(:,:), intent(out)::s
        integer:: x, y, m
        double precision:: phi_ant
        erro = -huge(erro)
               do y=1,k
                   do x=1,j
                       phi_ant=s(x,y)
                       s(x,y)=0.d0
                       do m=1,n*(n+2)/2
                           s(x,y)=s(x,y)+psi_medioxy(m,x,y)*peso(m)
                       end do
                       s(x,y)=s(x,y)/4.d0
                       !print*,'x, y, s(x,y)',x, y, s(x,y)
                       !erro=max(erro,abs(s(x,y)-phi_ant))
                       !print*, 'x, y, max, s', x, y, max(erro,abs(s(x,y)-phi_ant)), abs(s(x,y))
                       erro=max(erro,abs(s(x,y)-phi_ant)/abs(s(x,y)))
                   end do
               enddo
               !pause
        end subroutine calcula_erro


        !==============================================================================================================
        subroutine corrente(n, j, k, peso, Omega, psi_medioxy, J1, J2)
        !==============================================================================================================
        implicit none
        integer, intent(in):: n, j, k
        double precision, dimension(:), intent(in)::peso
        double precision, dimension(:, :), intent(in)::Omega
        double precision, dimension(:,:, :), intent(in) ::psi_medioxy
        double precision, dimension(:,:), intent(out)::J1, J2
        integer:: y, x, m
           do y=1,k
               do x=1,j
                   J1(x,y)=0.d0
                   J2(x,y)=0.d0
                   do m=1,n*(n+2)/8
                       J1(x,y)=J1(x,y)+Omega(m, 1)*psi_medioxy(m,x,y)*peso(m)-Omega(m, 1)*psi_medioxy(n*(n+2)/8+m,x,y)*&
                       			peso(n*(n+2)/8+m)-Omega(m, 1)*psi_medioxy(n*(n+2)/4+m,x,y)*peso(n*(n+2)/4+m)+Omega(m, 1)*&
                                psi_medioxy(3*n*(n+2)/8+m,x,y)*peso(3*n*(n+2)/8+m)
                       J2(x,y)=J2(x,y)+Omega(m, 2)*psi_medioxy(m,x,y)*peso(m)+Omega(m, 2)*psi_medioxy(n*(n+2)/8+m,x,y)*&
                       			peso(n*(n+2)/8+m)-Omega(m, 2)*psi_medioxy(n*(n+2)/4+m,x,y)*peso(n*(n+2)/4+m)-Omega(m, 2)*&
                                psi_medioxy(3*n*(n+2)/8+m,x,y)*peso(3*n*(n+2)/8+m)
                   end do

               end do
           enddo
        end subroutine corrente

        !==============================================================================================================
        subroutine fluxo_regioes(j,k, hx, hy, ndx, ndy, as, bs, s, R)
        !==============================================================================================================
        implicit none
        integer, intent(in):: j, k, ndx, ndy
        double precision, intent(in)::hx, hy
        double precision, dimension(:), intent(in):: as, bs
        double precision, dimension(:,:), intent(in):: s
        double precision, dimension(:,:), intent(out):: R
        integer, dimension(:), allocatable:: x,y
        integer::i, l, m, n
        double precision:: A

        allocate(x(0:ndx), y(0:ndy))
                R=0.d0
                x(0)=0
                y(0)=0
                do i=1, ndx
                  x(i)=nint(as(i)/hx)
               	end do
                do i=1, ndy
                  y(i)=nint(bs(i)/hy)
               	end do

				!i, l correm nas regioes
                !m, n sao indices para o somatorio em cada celula da malha em cada regiao

                do i=1, ndx
                  do l=1, ndy
                    A=0
                    do m=y(l-1)+1,y(l)
                      do n=x(i-1)+1,x(i)
                        A=A+hx*hy
                        R(i, l)=R(i, l)+s(n, m)*hx*hy
                      end do
                    end do
                    R(i, l)=R(i, l)/A
                  end do
                end do

    	end subroutine fluxo_regioes
        !====================================================
        subroutine nosepesosado(n,M,Omegamu, Omegaeta, peso)
        !====================================================
        implicit none
        integer, intent(in):: n,M
        double precision, dimension(:), intent(out)::peso
        double precision, dimension(:), intent(out):: Omegamu, Omegaeta
        double precision, dimension (:), allocatable:: w, mu
        integer:: l,i,x

            allocate(mu(n/2))

            do
              select case (n)
                case (2)
                    allocate (w(1))
                    mu(1)=0.5773502d0; w(1)=1.d0; peso=w(1)
                    exit
                case (4)
                    allocate (w(1))
                    mu(1)=0.3500212d0; mu(2)=0.8688903d0;w(1)=0.333333d0; peso=w(1);! w(1)=1.d0/3; peso=w(1)
                    exit
                case (6)
                    allocate (w(2))
                    mu(1)=0.2666355d0; mu(2)=0.6815076d0; mu(3)=0.9261808d0; w(1)=0.1761263d0; w(2)=0.1572071d0
                    peso(1)=w(1); peso(2)=w(2);peso(3)=w(2);peso(4)=w(1);peso(5)=w(2);peso(6)=w(1);
                    do i=1,6
                        peso(M/4+i)=peso(i)
                        peso(M/2+i)=peso(i)
                        peso(3*M/4+i)=peso(i)
                    end do
                    exit
                case (8)
                    allocate (w(3))
                    mu(1)=0.2182179; mu(2)=0.5773503; mu(3)=0.7867958; mu(4)=0.9511897; w(1)=0.1209877;
                    w(2)=0.0907407; w(3)=0.0925926
                    peso(1)=w(1);peso(2)=w(2);peso(3)=w(2);peso(4)=w(2);peso(5)=w(3);peso(6)=w(2);
                    peso(7)=w(1);peso(8)=w(2);peso(9)=w(2); peso(10)=w(1);
                    do i=1,10
                      peso(M/4+i)=peso(i)
                      peso(M/2+i)=peso(i)
                      peso(3*M/4+i)=peso(i)
                    end do
                    exit
                case (12)
                    allocate (w(5))
                    mu(1)=0.1672126; mu(2)=0.4595476; mu(3)=0.6280191; mu(4)=0.7600210; mu(5)=0.8722706;
                    mu(6)=0.9716377; w(1)=0.0707626; w(2)=0.0558811; w(3)=0.0373377; w(4)=0.0502819;
                    w(5)=0.0258513
                    peso(1)=w(1);peso(16)=w(1);peso(21)=w(1); peso(2)=w(2);peso(11)=w(2);peso(15)=w(2);
                    peso(17)=w(2);peso(20)=w(2);peso(3)=w(2); peso(4)=w(3);peso(6)=w(3);peso(7)=w(3);
                    peso(10)=w(3);peso(18)=w(3);peso(19)=w(3); peso(5)=w(4);peso(12)=w(4);peso(14)=w(4);
                    peso(8)=w(5);peso(9)=w(5);peso(13)=w(5);
                    do i=1,21
                      peso(M/4+i)=peso(i)
                      peso(M/2+i)=peso(i)
                      peso(3*M/4+i)=peso(i)
                    end do
                    exit
                case (16)
                    allocate (w(8))
                    mu(1)=0.1389568; mu(2)=0.3922893; mu(3)=0.5370966; mu(4)=0.6504264; mu(5)=0.7467506;
                    mu(6)=0.8319966;mu(7)=0.9092855;mu(8)=0.9805009; w(1)=0.0489872; w(2)=0.0413296;
                    w(3)=0.0212326; w(4)=0.0256207; w(5)=0.0360486;w(6)=0.0144589;w(7)=0.0344958; w(8)=0.0085179
                    peso(1)=w(1);peso(2)=w(2);peso(3)=w(2);peso(4)=w(3);peso(5)=w(5);peso(6)=w(3);
                    peso(7)=w(4);peso(8)=w(6);peso(9)=w(6);peso(10)=w(4);peso(11)=w(4);peso(12)=w(7);
                    peso(13)=w(8);peso(14)=w(7);peso(15)=w(4);peso(16)=w(3);peso(17)=w(6);peso(18)=w(8);
                    peso(19)=w(8);peso(20)=w(6);peso(21)=w(3);peso(22)=w(2);peso(23)=w(5);peso(24)=w(6);
                    peso(25)=w(7);peso(26)=w(6);peso(27)=w(5);peso(28)=w(2);peso(29)=w(1);peso(30)=w(2);
                    peso(31)=w(3);peso(32)=w(4);peso(33)=w(4);peso(34)=w(3);peso(35)=w(2);peso(36)=w(1);
                    do i=1,36
                      peso(M/4+i)=peso(i)
                      peso(M/2+i)=peso(i)
                      peso(3*M/4+i)=peso(i)
                    end do
                    exit
                end select
            end do


            !preenche o vetor de posicoes mu e eta no primeiro quadrante
            l=0
            do i=n/2,1, -1
              do x=1, n/2-i+1
                l=l+1
                Omegamu(l)=mu(i)
                Omegaeta(l)=mu(x)
                !print*, 'l', l
                !print*, 'omega', Omega(l,1), Omega(l,2)
              end do
            end do
            do i=1, M/4
             omegamu(i+M/4)=omegamu(i)
             omegaeta(i+M/4)=omegaeta(i)!pensando no outro ordenamento
            end do
            deallocate(mu, w)
        end subroutine nosepesosado

        !========================================================
        subroutine eigens(j, Omegamu, sigma_t, sigma_s, peso, U, V, ni)
        !========================================================
        implicit none
        integer, intent(in)::j
        double precision, intent(in)::sigma_t, sigma_s
        double precision, dimension(:), intent(in)::peso, Omegamu
        double precision, dimension(:), intent(out):: ni
        double precision, dimension(:,:), intent(out)::U, V
        double precision, dimension(:), allocatable:: D, wi,work,scale, rconde, rcondv,wr
        double precision, dimension (:,:), allocatable:: A, vr
        double precision, dimension(1)::aux
        integer, dimension(1):: iwork
        logical, dimension(:), allocatable :: bwork
        integer::i,k,ilo, ihi, lwork,info
        double precision::abnrm


            allocate(D(j), A(j,j), wi(j),vr(j,j), wr(j))

            D=0.d0
            A=0.d0
            do i=1, j
                D(i)=(sigma_t/Omegamu(i))**2
            end do

    !$$$$$$         do i=1, j
    !$$$$$$             print*, 'D', D(i)
    !$$$$$$         end do

            !constroi matriz -A (por conveniencia)
            do i=1,j
                do k=1, j
                    A(i,k)=-(sigma_t*sigma_s*peso(k))/(2.d0*Omegamu(i)**2)
                end do
            end do

            !altera matriz A pra D-A
            do i=1,j
                A(i,i)=D(i)+A(i,i)
            end do

    !$$$$$$         do i=1,j
    !$$$$$$             do k=1, j
    !$$$$$$                 print*, 'A, i, j',i,k, A(i,k)
    !$$$$$$             enddo
    !$$$$$$         end do

            lwork = 10*j
            allocate(work(lwork),bwork(n),scale(j),rconde(j),rcondv(j))
            ilo = 1
            ihi = j


            call dgeevx('N','N','V','N',j,A,j,wr,wi,aux,1,vr,j,ilo,ihi,scale,abnrm,rconde,rcondv,work,lwork,iwork,info)
    !$$$$$$             if (info==0) then
    !$$$$$$                 !print *,'DGEES terminou com sucesso...'
    !$$$$$$                 print *,'autovalores:'
    !$$$$$$                 do i=1,j
    !$$$$$$                     print*,wr(i),wi(i)
    !$$$$$$                 end do
    !$$$$$$                 pause
    !$$$$$$                 print *,'autovetores:'
    !$$$$$$                     do i=1,j
    !$$$$$$                         print *,vr(i,:)
    !$$$$$$                     end do
    !$$$$$$             end if
    !$$$$$$             pause
            deallocate(work,bwork,scale,rconde,rcondv)

            do i=1,j
                U(i,:)=vr(i,:)
            end do

    !$$$$$$         do i=1,j
    !$$$$$$           do k=1, j
    !$$$$$$             print*,'i,k, U',i,k, U(i,k)
    !$$$$$$           end do
    !$$$$$$         end do
    !$$$$$$           pause

             do i=1,j
                 if(wr(i).LE.0.d0) then
                     wr(i)=dabs(wr(i))
                 end if
                 ni(i)=1.d0/dsqrt(wr(i))
             end do

    !$$$$$$          print *,'ni=',ni
    !$$$$$$          pause

            do i=1, j
              do k=1,j
                V(i,k)=Omegamu(i)*U(i,k)/(ni(k)*sigma_t)
              end do
            end do
    !$$$$$$         do i=1, j
    !$$$$$$           do k=1,j
    !$$$$$$              print*, 'V, i, k', i,k, V(i,k)
    !$$$$$$            end do
    !$$$$$$         end do
            deallocate(D, A, wi,vr, wr)
         end subroutine eigens

        !=============================================================
        subroutine phis(M, j, U, V, phiy)
        !=============================================================
        implicit none
        integer, intent (in)::M,j
        double precision, dimension(:,:), intent(in)::U, V
        double precision, dimension(:,:), intent(out)::phiy
        integer:: i,k

            do i=1, j
                do k=1, j
                    phiy(i,k)=0.5d0*(U(i,k)+V(i,k))
                    phiy(j+i, k)=0.5d0*(U(i,k)-V(i,k))
                end do
            end do

    !$$$$$$         do i=1, M
    !$$$$$$             do k=1, j
    !$$$$$$                 print*, 'i,j, phi', i, k, phiy(i,k)
    !$$$$$$             end do
    !$$$$$$         end do
        end subroutine phis

        !=============================================================
        subroutine matsist(M,j, Omegamu, a, phiy, as, ni, Omegaeta, b, phix, bs, gama, sigma_s, sigma_t,&
                            q,peso, VetSol)
        !=============================================================
        implicit none
        integer, intent(in)::M,j
        double precision, intent(in)::a, as, b, bs, sigma_s, sigma_t,q
        double precision, dimension(:), intent(in):: Omegamu,ni, Omegaeta, gama, peso
        double precision, dimension(:,:), intent(in):: phiy, phix
        integer::i,k,z,info
        integer, dimension(:), allocatable :: ipiv
        double precision, dimension(:), intent(out):: VetSol
        double precision, dimension(:,:), allocatable:: MAT
            z=M/4
            allocate (MAT(8*M, 8*M))
            MAT=0.d0
            VetSol=0.d0

            !matrizes m/4*m/2
            do i=1, z
                do k=1, j
                    !A
                    MAT(i,k)=-Omegamu(i)/a*phiy(i, k)
                    MAT(z+i, k)=Omegamu(i)/a*phiy(j+i, k)
                    MAT(j+i, k)=-Omegamu(i)/a*phiy(z+i, k)
                    MAT(3*z+i, k)=Omegamu(i)/a*phiy(3*z+i,k)
                    MAT(M+i, k)=Mat(i,k)
                    MAT(M+z+i, k)= MAT(z+i, k)
                    MAT(M+j+i, k)=MAT(j+i, k)
                    MAT(M+3*z+i, k)=MAT(3*z+i, k)
                    !B
                    MAT(i,j+k)=-Omegamu(i)/a*phiy(j+i, k)*exp(-as/ni(k))
                    MAT(z+i, j+k)=Omegamu(i)/a*phiy(i, k)*exp(-as/ni(k))
                    MAT(j+i, j+k)=-Omegamu(i)/a*phiy(3*z+i, k)*exp(-as/ni(k))
                    MAT(3*z+i, j+k)=Omegamu(i)/a*phiy(z+i,k)*exp(-as/ni(k))
                    MAT(M+i, j+k)=Mat(i,j+k)
                    MAT(M+z+i, j+k)= MAT(z+i, j+k)
                    MAT(M+j+i, j+k)=MAT(j+i, j+k)
                    MAT(M+3*z+i, j+k)=MAT(3*z+i, j+k)
                    !C
                    MAT(i,m+k)=Omegamu(i)/a*phiy(i, k)*exp(-(a-as)/ni(k))
                    MAT(j+i, m+k)=Omegamu(i)/a*phiy(z+i, k)*exp(-(a-as)/ni(k))
                    MAT(M+i, m+k)=Mat(i,m+k)
                    MAT(M+j+i, m+k)=MAT(j+i, m+k)
                    !D
                    MAT(i,m+j+k)=Omegamu(i)/a*phiy(j+i, k)
                    MAT(j+i, m+j+k)=Omegamu(i)/a*phiy(j+z+i, k)
                    MAT(M+i, m+j+k)=Mat(i,m+j+k)
                    MAT(M+j+i, m+j+k)=MAT(j+i, m+j+k)
                    !E
                    MAT(2*m+i,2*m+k)=-Omegaeta(i)/b*phix(i, k)
                    MAT(2*m+z+i, 2*m+k)=Omegaeta(i)/b*phix(j+i, k)
                    MAT(2*m+j+i, 2*m+k)=-Omegaeta(i)/b*phix(z+i, k)
                    MAT(2*m+3*z+i, 2*m+k)=Omegaeta(i)/b*phix(3*z+i,k)
                    MAT(3*m+i, 2*m+k)=Mat(2*m+i,2*m+k)
                    MAT(3*m+z+i, 2*m+k)= MAT(2*m+z+i, 2*m+k)
                    MAT(3*m+j+i, 2*m+k)=MAT(2*m+j+i, 2*m+k)
                    MAT(3*m+3*z+i, 2*m+k)=MAT(2*m+3*z+i, 2*m+k)
                    !F
                    MAT(2*m+i,2*m+j+k)=-Omegaeta(i)/b*phix(j+i, k)*exp(-bs/gama(k))
                    MAT(2*m+z+i, 2*m+j+k)=Omegaeta(i)/b*phix(i, k)*exp(-bs/gama(k))!tinha um menos na frente?tirei ps camila
                    MAT(2*m+j+i, 2*m+j+k)=-Omegaeta(i)/b*phix(3*z+i, k)*exp(-bs/gama(k))
                    MAT(2*m+3*z+i, 2*m+j+k)=Omegaeta(i)/b*phix(z+i,k)*exp(-bs/gama(k))
                    MAT(3*M+i, 2*m+j+k)=Mat(2*m+i,2*m+j+k)
                    MAT(3*M+z+i, 2*m+j+k)= MAT(2*m+z+i, 2*m+j+k)
                    MAT(3*M+j+i, 2*m+j+k)=MAT(2*m+j+i, 2*m+j+k)
                    MAT(3*M+3*z+i, 2*m+j+k)=MAT(2*m+3*z+i, 2*m+j+k)
                    !G
                    MAT(2*m+i,3*m+k)=Omegaeta(i)/b*phix(i, k)*exp(-(b-bs)/gama(k))
                    MAT(2*m+j+i, 3*m+k)=Omegaeta(i)/b*phix(z+i, k)*exp(-(b-bs)/gama(k))
                    MAT(3*M+i, 3*m+k)=MAT(2*m+i,3*m+k)
                    MAT(3*M+j+i, 3*m+k)=MAT(2*m+j+i, 3*m+k)
                    !H
                    MAT(2*m+i,3*m+j+k)=Omegaeta(i)/b*phix(j+i, k)
                    MAT(2*m+j+i, 3*m+j+k)=Omegaeta(i)/b*phix(j+z+i, k)
                    MAT(3*M+i, 3*m+j+k)=Mat(2*m+i,3*m+j+k)
                    MAT(3*M+j+i, 3*m+j+k)=MAT(2*m+j+i, 3*m+j+k)
                end do
            enddo
            !matrizes m/2*m/2
            do i=1, j
                do k=1, j
                    !A
                    MAT(5*m+i,k)=phiy(i, k)-phiy(j+i, k)
                    Mat(6*M+i,k)=phiy(i, k)*exp(-as/ni(k))
                    Mat(6*M+j+i, k)=phiy(j+i, k)*exp(-as/ni(k))
                    !B
                    MAT(5*m+i,j+k)=-Mat(5*m+i,k)*exp(-as/ni(k))
                    Mat(6*M+i,j+k)=phiy(j+i, k)
                    Mat(6*M+j+i, j+k)=phiy(i, k)
                    !C
                    MAT(4*M+i, m+k)=phiy(j+i, k)*exp(-(a-as)/ni(k))
                    Mat(6*M+i,m+k)=-phiy(i, k)
                    Mat(6*M+j+i, m+k)=-phiy(j+i, k)
                    !D
                    MAT(4*M+i, m+j+k)=phiy(i, k)
                    Mat(6*M+i,m+j+k)=-phiy(j+i, k)*exp(-(a-as)/ni(k))
                    Mat(6*M+j+i, m+j+k)=-phiy(i, k)*exp(-(a-as)/ni(k))
                    !E
                    MAT(5*m+j+i,2*m+k)=phix(i, k)-phix(j+i, k)
                    Mat(7*M+i,2*m+k)=phix(i, k)*exp(-bs/gama(k))
                    Mat(7*M+j+i, 2*m+k)=phix(j+i, k)*exp(-bs/gama(k))
                    !F
                    MAT(5*m+j+i, 2*m+j+k)=-MAT(5*m+j+i,2*m+k)*exp(-bs/gama(k))
                    Mat(7*M+i, 2*m+j+k)=phix(j+i, k)
                    Mat(7*M+j+i,  2*m+j+k)=phix(i, k)
                    !G
                    MAT(4*M+j+i, 3*m+k)=phix(j+i, k)*exp(-(b-bs)/gama(k))
                    Mat(7*M+i,3*m+k)=-phix(i, k)
                    Mat(7*M+j+i, 3*m+k)=-phix(j+i, k)
                    !H
                    MAT(4*M+j+i, 3*m+j+k)=phix(i, k)
                    Mat(7*M+i,3*m+j+k)=-phix(j+i, k)*exp(-(b-bs)/gama(k))
                    Mat(7*M+j+i, 3*m+j+k)=-phix(i, k)*exp(-(b-bs)/gama(k))
                end do
            enddo

            !matrizes m*m/4
            do i=1, z
                !O
                Mat(i, 4*m+i)=-Omegamu(i)/a
                Mat(z+i, 4*m+j+i)=Omegamu(i)/a
                Mat(j+i, 4*m+z+i)=-Omegamu(i)/a
                Mat(3*z+i, 4*m+3*z+i)=Omegamu(i)/a
                Mat(m+i, 4*m+i)=Mat(i, 4*m+i)
                Mat(m+z+i, 4*m+j+i)=Mat(z+i, 4*m+j+i)
                Mat(m+j+i, 4*m+z+i)=Mat(j+i, 4*m+z+i)
                Mat(m+3*z+i, 4*m+3*z+i)=Mat(3*z+i, 4*m+3*z+i)
                !P
                Mat(i, 5*m+i)=Omegamu(i)/a
                Mat(j+i, 5*m+z+i)=Omegamu(i)/a
                Mat(m+i, 5*m+i)=Mat(i, 5*m+i)
                Mat(m+j+i, 5*m+z+i)=Mat(j+i, 5*m+z+i)
                !R
                Mat(2*m+i, 6*m+i)=-Omegaeta(i)/b
                Mat(2*m+z+i, 6*m+j+i)=Omegaeta(i)/b
                Mat(2*m+j+i, 6*m+z+i)=-Omegaeta(i)/b
                Mat(2*m+3*z+i, 6*m+3*z+i)=Omegaeta(i)/b
                Mat(3*m+i, 6*m+i)=Mat(2*m+i, 6*m+i)
                Mat(3*m+z+i, 6*m+j+i)=Mat(2*m+z+i, 6*m+j+i)
                Mat(3*m+j+i, 6*m+z+i)=Mat(2*m+j+i, 6*m+z+i)
                Mat(3*m+3*z+i, 6*m+3*z+i)=Mat(2*m+3*z+i, 6*m+3*z+i)
                !S
                Mat(2*m+i, 7*m+i)=Omegaeta(i)/b
                Mat(2*m+j+i, 7*m+z+i)=Omegaeta(i)/b
                Mat(3*m+i, 7*m+i)=Mat(2*m+i, 7*m+i)
                Mat(3*m+j+i, 7*m+z+i)=Mat(2*m+j+i, 7*m+z+i)
            end do
            !matriz mxm
            do i=1, m
                !O
                Mat(2*m+1:3*m, 4*m+i)=-(sigma_s/4)*peso(i)
                Mat(2*m+i, 4*m+i)=Mat(2*m+i, 4*m+i)+sigma_t
                mat(6*m+i, 4*m+i)=1.d0
                !P
                Mat(3*m+1:4*m, 5*m+i)=-(sigma_s/4)*peso(i)
                Mat(3*m+i, 5*m+i)=Mat(3*m+i, 5*m+i)+sigma_t
                mat(6*m+i, 5*m+i)=-1.d0
                !R
                Mat(1:m, 6*m+i)=-(sigma_s/4)*peso(i)
                Mat(i, 6*m+i)=Mat(i, 6*m+i)+sigma_t
                mat(7*m+i, 6*m+i)=1.d0
                !S
                Mat(m+1:2*m, 7*m+i)=-(sigma_s/4)*peso(i)
                Mat(m+i, 7*m+i)=Mat(m+i, 7*m+i)+sigma_t
                mat(7*m+i, 7*m+i)=-1.d0
            end do

            !matrizes m/2xm
            do i=1,j
                !O
                mat(5*m+i, 4*m+i)=1.d0
                mat(5*m+i, 4*m+j+i)=-1.d0
                !P
                mat(4*m+i, 5*m+j+i)=1.d0
                !R
                mat(5*m+j+i, 6*m+i)=1.d0
                mat(5*m+j+i, 6*m+j+i)=-1.d0
                !S
                mat(4*m+j+i, 7*m+j+i)=1.d0
            enddo
    !$$$$$$         do i=1, 8*m
    !$$$$$$             do k=1, 8*m
    !$$$$$$                 write(10, *)i, k, mat(i,k)
    !$$$$$$             end do
    !$$$$$$         end do


            !Vetor Independente
            VetSol(1:M)=q*as/a
            VetSol(2*M+1:3*M)=q*bs/b


            !dgetrf(M,N,A,lda, ipiv, info) calcula a fatoracao LU de uma matriz A , mxn usando pivotamento parcial
            !sem troca de linhas onde:
            !M(input), integer, numero de linhas da matriz
            !N(input), integer, numero de colunas da matriz A
            !A(input/output), real/complex matriz que deve ser fatorada, com saida na forma A=P*L*U, com L com diagonal de um
            !lda(input) integer, dimens�o maxima do vetor A. Lda>=max(1,M)
            !ipiv(output)integer, vetor de dimensao(min(M,N)).
            !info(output) integer, =0 successful exit, ou se ~=0 deu algum erro
                allocate (ipiv(8*M))
                call dgetrf (8*M,8*M,Mat,8*M,ipiv,info)
    !$$$$$$             if (info==0) then
    !$$$$$$                     print *,'DGETRF terminou com sucesso...'
    !$$$$$$             end if

            !dgetrs(trans,N,nrhs,A,lda, ipiv, B, LDB, info) resolve um sistema linear AX=B, A'X=B ou A*X=b
            !usando a fatoracao LU onde:
            !trans(input, caracter, 'N',AX=B (no transpose), 'T', A'X=B(transpose), 'C', A*X=B (conjugate transpose)
            !N(input), integer, ordem da matriz A
            !nrhs(input), integer, numero de colunas de B
            !A(input), fatores L e U da fatoracao A=P*L*U, com L com diagonal de um
            !lda(input) integer, dimensao maxima do vetor A. Lda>=max(1,M)
            !ipiv(input)integer, vetor de dimensao(min(M,N)).
            !B(input/output) real/complex array dimension(M, nrhs), entra B e sai X, solucao do sistema
            !ldb(input) integer, numero de linhas de B
            !info(output) integer, =0 successful exit, ou se <0 deu algum erro
              call dgetrs('N', 8*M , 1, MAT , 8*M, ipiv, VetSol, 8*M, info)
               deallocate (ipiv)
    !$$$$$$            do i=1, 8*m
    !$$$$$$              print*, 'i, VetSol(i)',i, VetSol(i)
    !$$$$$$            end do
    !$$$$$$             pause
            deallocate(MAT)

        end subroutine matsist

        !==============================================================================================================================
        subroutine solutions(M, j,kk,jj, phiy, ni, a, phix, gama, b, as, bs, peso,VetSol, hx, hy, s)
        !==========================================================================================================================
        implicit none
        integer, intent(in):: M, j, kk, jj
        double precision, intent(in):: a, b, as, bs, hx, hy
        double precision, dimension(:), intent(in)::  ni, gama, peso, VetSol
        double precision, dimension (:,:), intent(in)::phiy, phix
        double precision, dimension (:,:), intent(out):: s
        double precision :: x,y, psiy1, psiy2, psix1, psix2, fluxescx, fluxescy, r1, r2, r3, r4
        double precision, dimension(:), allocatable:: psiy, psix,psixsource, psiysource
        double precision, dimension(:,:), allocatable:: phixx, phiyy
        integer:: i, k,l, p

            allocate(psiy(M), psix(M), psiysource(M), psixsource(M), phixx(jj,kk), phiyy(jj,kk))
            x=0.d0
            y=0.d0
            s=0.d0
            do  l=1, kk
                do p=1, jj
                    x=hx*(2*p-1)/2.d0
                    y=hy*(2*l-1)/2.d0

                    psiy=0.d0
                    psix=0.d0

                    !solu�cao na fonte
                    if (x<=as) then
                        do i=1, j
                            psix1=0.d0
                            psix2=0.d0
                            do k=1,j
                               psix1=VetSol(2*M+k)*phix(i, k)*exp(-y/gama(k))+VetSol(2*M+j+k)*phix(j+i, k)*exp(-(bs-y)/gama(k))&
                               +psix1
                               psix2=VetSol(2*M+k)*phix(j+i, k)*exp(-y/gama(k))+VetSol(2*m+j+k)*phix(i, k)*exp(-(bs-y)/gama(k))&
                               +psix2
                            end do
                            psix(i)=psix1+VetSol(6*M+i)
                            psix(j+i)=psix2+VetSol(6*M+j+i)
                        end do
                    end if
                    if (y<=bs) then
                        do i=1, j
                            psiy1=0.d0
                            psiy2=0.d0
                            do k=1,j
                                psiy1=VetSol(k)*phiy(i, k)*exp(-x/ni(k))+VetSol(j+k)*phiy(j+i, k)*exp(-(as-x)/ni(k))+psiy1
                                psiy2=VetSol(k)*phiy(j+i, k)*exp(-x/ni(k))+VetSol(j+k)*phiy(i, k)*exp(-(as-x)/ni(k))+psiy2
                            end do
                            psiy(i)=psiy1+VetSol(4*m+i)
                            psiy(j+i)=psiy2+VetSol(4*m+j+i)
                        end do
                    end if
                    if (x>=as) then
                    !solucao sem fonte
                        do i=1, j
                            psix1=0.d0
                            psix2=0.d0
                            do k=1,j
                                psix1=VetSol(3*M+k)*phix(i, k)*exp(-(y-bs)/gama(k))+&
                                        &VetSol(3*M+j+k)*phix(j+i, k)*exp(-(b-y)/gama(k))+psix1
                                psix2=VetSol(3*M+k)*phix(j+i, k)*exp(-(y-bs)/gama(k))+&
                                        &VetSol(3*M+j+k)*phix(i, k)*exp(-(b-y)/gama(k))+psix2
                            end do
                            psix(i)=psix1+VetSol(7*M+i)
                            psix(j+i)=psix2+VetSol(7*M+j+i)
                        end do
                    end if
                    if (y>=bs) then
                    !solucao sem fonte
                        do i=1, j
                            psiy1=0.d0
                            psiy2=0.d0
                            do k=1,j
                                psiy1=VetSol(M+k)*phiy(i, k)*exp(-(x-as)/ni(k))+VetSol(M+j+k)*phiy(j+i, k)*exp(-(a-x)/ni(k))+psiy1
                                psiy2=VetSol(M+k)*phiy(j+i, k)*exp(-(x-as)/ni(k))+VetSol(M+j+k)*phiy(i, k)*exp(-(a-x)/ni(k))+psiy2
                            end do
                            psiy(i)=psiy1+VetSol(5*m+i)
                            psiy(j+i)=psiy2+VetSol(5*m+j+i)
                        end do
                    end if
            !$$$$$$     do i=1, m
            !$$$$$$       print*, 'i, psiysource, psiy', i, psiysource(i), psiy(i)
            !$$$$$$     end do
            !$$$$$$     pause
            !$$$$$$     do i=1, m
            !$$$$$$       print*, 'i, psixsource, psix', i, psixsource(i), psix(i)
            !$$$$$$     end do
                    fluxescx=0.d0
                    fluxescy=0.d0
                    !fluxos escalares
                    do i=1, m
                        fluxescx=.25*psix(i)*peso(i)+fluxescx
                        fluxescy=.25*psiy(i)*peso(i)+fluxescy
                    end do
!$$$$$$                     phixx(p,l)=fluxescx
!$$$$$$                     phiyy(p,l)=fluxescy
                    s(p,l)=(fluxescx+fluxescy)/2.d0
    !                print*, 'fluxx, fluxy, x, y, s', fluxescx,  fluxescy,x, y, s(p,l)
                end do
            end do
!$$$$$$             call fluxo_regioes(jj,kk, hx, hy,as, bs, phixx, R1, R2, R3, R4)
!$$$$$$             print*, 'phix:r1, r2, r3, r4', r1, r2, r3, r4
!$$$$$$             call fluxo_regioes(jj,kk, hx, hy,as, bs, phiyy, R1, R2, R3, R4)
!$$$$$$             print*, 'phiy:r1, r2, r3, r4', r1, r2, r3, r4
        end subroutine solutions

         !============================================================================================
        double precision function diferenca_tempo(t0,t1) result (valor)
        !============================================================================================
        implicit none
        integer, dimension(8), intent(in) :: t0,t1
        ! t0 e t1 sao valores retornados pela chamada do intr�nseco DATE_AND_TIME()
        ! valor e a diferenca de tempo, em segundos
        ! somente pode ser usada para medicao de tempo de trechos de programa que executem
        ! em menos de um mes
        double precision :: t0_ms,t1_ms
            t0_ms = t0(3)*24*3600.0D3+t0(5)*3600.0D3+t0(6)*60.0D3+t0(7)*1.0D3+t0(8)
            t1_ms = t1(3)*24*3600.0D3+t1(5)*3600.0D3+t1(6)*60.0D3+t1(7)*1.0D3+t1(8)
            valor = (t1_ms-t0_ms)*1.0E-3
        end function diferenca_tempo

end program ado_bidimensional
