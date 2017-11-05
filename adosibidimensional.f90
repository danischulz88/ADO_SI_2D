PROGRAM ado_bidimensional
  IMPLICIT NONE

  INTEGER:: n, kk, jj, nado, nit, l,ll, ndx, ndy,  i, j, k, prop, ncell
  DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE::  s
  DOUBLE PRECISION::  hx, hy, compx, compy, tol
  DOUBLE PRECISION:: in_sigma_t, in_sigma_s0, in_sigma_s1, source_sigma_t, source_sigma_s0, source_sigma_s1, insource
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: psi_medioxy, psi_medioy,psi_mediox
  DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE ::sigma_t, sigma_s0, sigma_s1, q, R, C, dig
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: as, bs
  DOUBLE PRECISION, DIMENSION(5000):: dadoscamila
  CHARACTER(len=80) :: nome_arquivo,buffer
  LOGICAL :: existe
  ! variaveis para medicao de tempo de execucao
  DOUBLE PRECISION :: t
  INTEGER, DIMENSION(8) :: t0,t1
  !variaveis de automatizacao
  INTEGER, DIMENSION(6)::direcoes, discret
  !prop eh o numero de vezes que o numero de pontos fornecidos pelo ado forcado  cabe na malha
  !Varredura 2D

  ! Considera o particionamento de uma regiao de comprimento LxM em J celulas ao longo do eixo X,
  !com J+1 paredes;e K celulas no eixo Y com K+1 paredes;fluxos angulares sao calculados em cada
  !parede (logo, psi tem N direcoes em J+1 posicoes em X e K+1 posicoes em Y);fluxos angulares
  !medios sao calculados no meio de cada celula (logo,psi_medio tem N direcoes em J posicoes em
  !X e K posicoes em Y); o tamanho de cada celula pode ser constante ou nao.
100 FORMAT(' ',1I20, 1I20, 1E20.6)

  OPEN(unit=10,file='fluxbidi.txt')
  OPEN(unit=20,file='regioessi.txt')
  OPEN(unit=30,file='regioessiado.txt')
  OPEN(unit=50,file='digsesi.txt')
  OPEN(unit=60,file='digseado.txt')

  !numero de células no eixo x para a descretização forçada
  OPEN(unit=40,file='dadosproblema50.dat',status='old')
  ncell=0
  DO
     READ(40, *, END=99)
     ncell=ncell+1
  END DO
99 PRINT*, 'numero de linhas no arquvo', ncell
  ncell=SQRT(ncell/2.0)
  PRINT*, 'ncell', ncell
  CLOSE(40)

  OPEN(unit=40,file='dadosproblema50.dat',status='old')
  !inicializacoes do problema
  DO i=1, 2*ncell*ncell
     READ(40, *)dadoscamila(i)
  END DO

  ! print *,'Digite a ordem da quadratura (ex. S_4 digite 4)'
  !read *, n
  direcoes(1)=2
  direcoes(2)=4
  direcoes(3)=6
  direcoes(4)=8
  direcoes(5)=12
  direcoes(6)=16
  PRINT *,'Entre o nome do arquivo de dados:'
  READ *,nome_arquivo
  INQUIRE(file=nome_arquivo,exist=existe)
  IF (.NOT.existe) THEN
     PRINT *,'Arquivo inexistente, terminando...'
     STOP
  END IF


  OPEN(unit=90,file=nome_arquivo,status='old', action='read')
  !numero de celulas
  DO i=1,6
     READ(90, *),discret(i)
  END DO
  PRINT*, 'discret', discret
  !numero de regioes no eixo x
  READ(90, '(/A80)')buffer ; READ(buffer, *, err=900)ndx
  PRINT*, 'ndx', ndx
  !numero de regioes no eixo y
  READ(90, '(/A80)')buffer ;PRINT*, 'ok', buffer; READ(buffer, *, err=900)ndy
  PRINT*, 'ndy', ndy
  ALLOCATE(R(ndx, ndy), C(ndx, ndy),dig(ndx, ndy), as(ndx), bs(ndy))
  nado=4!numero de pontos em x e/ou y para as direcoes

  ! definicoes das constantes do problema
  READ(90, '(/A80)')buffer ; READ(buffer, *, err=900)compx
  READ(90, '(/A80)')buffer ; READ(buffer, *, err=900)compy
  READ(90, '(/A80)')buffer ; READ(buffer, *, err=900)tol
  PRINT*, 'compx, compy, tol', compx, compy, tol

  DO i=1, ndx
     READ(90, *), as(i)
  END DO
  PRINT*, 'as', as

  DO i=1, ndy
     READ(90, '(/A80)')buffer ; READ(buffer, *, err=900)bs(i)
  END DO
  PRINT*, 'bs',bs

  !parametros gerais do programa
  READ(90, '(/A80)')buffer ; READ(buffer, *, err=900)in_sigma_t,in_sigma_s0,in_sigma_s1
  PRINT*, 'sigmas', in_sigma_t,in_sigma_s0,in_sigma_s1
  !parametros na regiao da fonte
  READ(90, '(/A80)')buffer ; READ(buffer, *, err=900)source_sigma_t,source_sigma_S0,source_sigma_S1
  PRINT*, 'sigmas fonte', source_sigma_t,source_sigma_S0,source_sigma_S1
  !fonte
  READ(90, '(/A80)')buffer ; READ(buffer, *, err=900)insource
  PRINT*, 'fonte', insource
  ! Fecha o arquivo
  CLOSE(unit=90,status='keep')
  GOTO 200
900 PRINT *,'Erro na leitura do arquivo, terminando...'
  CLOSE(unit=90,status='keep')
  STOP
200 CONTINUE

  ALLOCATE(s(ncell, ncell))
  s=0
  DO i=1, ncell
     DO j=1, ncell
        k=2*j+2*ncell*(i-1)-1
        S(i,j)=(dadoscamila(k)+dadoscamila(k+1))/2
     END DO
  END DO

  CALL fluxo_regioes(ncell,ncell, 1.d0, 1.d0, ndx, ndy, as, bs, s, C)

  DO i=1, ndx
     DO j=1, ndy
        PRINT *,'SI+ADO(forcado Camila):i, j,  fluxos nas regioes C(i,j)=',i, j, C(i,j)
     ENDDO
  END DO


  DEALLOCATE(s)

  DO ll=1, 6
     DO l=1, 6

        n=direcoes(l)!numero de direcoes para a varredura

        jj=discret(ll);kk=jj; hx=compx/jj; hy=compy/kk

        !compx=1.d0;compy=1.d0; jj=10;kk=jj; hx=compx/jj; hy=compy/kk;tol=1.0d-6
        PRINT*, 'n ; jxk', n, jj,'x',kk
        prop=jj/50!jj pelo número de pontos que contem a discretização de dados externos no eixo x

        ALLOCATE(psi_medioxy(n*(n+2)/2,jj, kk))
        ALLOCATE(psi_mediox(n*(n+2)/2, jj, kk+1))
        ALLOCATE(psi_medioy(n*(n+2)/2, jj+1, kk))
        ALLOCATE(sigma_t(jj,kk))
        ALLOCATE(sigma_s0(jj,kk))
        ALLOCATE(sigma_s1(jj,kk))
        ALLOCATE(s(jj,kk),q(jj,kk) )
        !pontos de divisao das regioes para obter fluxos medios


        !parametros gerais do programa
        sigma_t(1:jj, 1:kk)=in_sigma_t
        sigma_s0(1:jj, 1:kk)=in_sigma_s0
        sigma_s1(1:jj, 1:kk)=in_sigma_s1
        !parametros na regiao da fonte
        sigma_t(1:NINT(as(1)/hx), 1:NINT(bs(1)/hy))=source_sigma_t
        sigma_s0(1:NINT(as(1)/hx), 1:NINT(bs(1)/hy))=source_sigma_S0
        sigma_s1(1:NINT(as(1)/hx), 1:NINT(bs(1)/hy))=source_sigma_S1


        !inicializacao da fonte
        q=0.

        q(1:NINT(as(1)/hx), 1:NINT(bs(1)/hy))=insource!magnitude da fonte em [0,as]x[0,bs]
        ! s=0
        ! DO i=1, ncell
        !    DO j=1, ncell
        !       k=2*j+2*ncell*(i-1)-1
        !       S((i-1)*prop+1:prop*i, (j-1)*prop+1:prop*j)=(dadoscamila(k)+dadoscamila(k+1))/2
        !    END DO
        ! END DO
        !
        !
        ! CALL fluxo_regioes(jj,kk, hx, hy, ndx, ndy, as, bs, s, C)
        !
        ! DO i=1, ndx
        !    DO j=1, ndy
        !       PRINT *,'SI+ADO(forcado Camila):i, j,  fluxos nas regioes C(i,j)=',i, j, C(i,j)
        !    ENDDO
        ! END DO

        !inicializacao dos parametros de entrada
        psi_mediox = 0.d0
        psi_medioy = 0.d0
        psi_medioxy= 0.d0
        s=0.d0


        ! inicializa psi na fronteira
        psi_mediox(n*(n+2)/4+1:n*(n+2)/2,1:jj, kk+1) = 0.d0
        psi_medioy(n*(n+2)/4+1:n*(n+2)/2,jj+1,1: kk) = 0.d0

        PRINT*, 'chama SI'
        ! chama SI
        CALL date_and_TIME(values=t0)
        CALL si_bi(n, jj, kk, hx, hy, tol, sigma_t, sigma_s0,sigma_s1, q, psi_mediox, psi_medioy, psi_medioxy, s, nit, as, bs,&
             ndx, ndy,R)
        CALL date_and_TIME(values=t1)
        t = diferenca_tempo(t0,t1)
        PRINT *,'numero de iteracoes: nit=', nit
        PRINT *,'SI: Tempo de execucao: t=',t,' s'
        DO i=1, ndx
           DO j=1, ndy
              PRINT *,'SI:i, j,  fluxos nas regioes R(i,j)=',i, j, R(i,j)
           ENDDO
        END DO
        WRITE(20, *)n, jj, nit
        DO i=1, ndx
           DO j=i , ndy
              WRITE(20, 100)i, j, R(i,j)
           ENDDO
        END DO
        CALL digse(ndx, ndy, C, R, dig)
        DO i=1, ndx
           DO j=1, ndy
              PRINT *,'SI:i, j,  digtos exatos nas regioes R(i,j)=',i, j, dig(i,j)
           ENDDO
        END DO
        WRITE(50, *)n, jj, nit
        DO i=1, ndx
           DO j=i , ndy
              WRITE(50, 100)i, j, dig(i,j)
           ENDDO
        END DO

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
        PRINT*, 'chama si+ado'
        CALL date_and_TIME(values=t0)
        DO i=1, ncell
           DO j=1, ncell
              k=2*j+2*ncell*(i-1)-1
              S((i-1)*prop+1:prop*i, (j-1)*prop+1:prop*j)=(dadoscamila(k)+dadoscamila(k+1))/2
           END DO
        END DO
        ! S(1:nint(as(1)/hx), 1:nint(bs(1)/hy))=1.8362d0
        ! S(1:nint(as(1)/hx), nint(bs(1)/hy)+1:kk)=0.010658d0
        ! S(nint(as(1)/hx)+1:jj, 1:nint(bs(1)/hy))=0.010658d0
        ! S(nint(as(1)/hx)+1:jj, nint(bs(1)/hy)+1:kk)=0.00011739d0

        CALL si_bi(n, jj, kk, hx, hy, tol, sigma_t, sigma_s0,sigma_s1, q, psi_mediox, psi_medioy, psi_medioxy, s, nit, as, bs,&
             ndx, ndy,R)
        CALL date_and_TIME(values=t1)
        t = diferenca_tempo(t0,t1)
        PRINT *,'SI+ADO(forcado):numero de iteracoes: nit=', nit
        ! print *,'SI+ADO(forcado): Tempo de execucao: t=',t,' s'
        DO i=1, ndx
           DO j=1, ndy
              PRINT *,'SI+ADO(forcado):i, j,  fluxos nas regioes R(i,j)=',i, j, R(i,j)
           ENDDO
        END DO
        WRITE(30, *)n, jj, nit
        DO i=1, ndx
           DO j=1, ndy
              WRITE(30, 100)i, j, R(i,j)
           ENDDO
        END DO
        CALL digse(ndx, ndy, C, R, dig)
        DO i=1, ndx
           DO j=1, ndy
              PRINT *,'SI+ado:i, j,  digtos exatos nas regioes R(i,j)=',i, j, dig(i,j)
           ENDDO
        END DO
        WRITE(60, *)n, jj, nit
        DO i=1, ndx
           DO j=i , ndy
              WRITE(60, 100)i, j, dig(i,j)
           ENDDO
        END DO
        DEALLOCATE (psi_medioxy, psi_mediox,psi_medioy,sigma_t, sigma_s0,sigma_s1, s,q)
     ENDDO
  END DO

CONTAINS

  !=====================================================================================================
  SUBROUTINE ado_bi(n, jj, kk, sigma_t, sigma_s, a, b, as, bs, hx, hy,q, s, R)
    !=====================================================================================================
    IMPLICIT NONE
    INTEGER, INTENT(in)::n, jj, kk
    INTEGER:: M, j
    DOUBLE PRECISION, INTENT(in):: sigma_s, sigma_t, a, b, hx, hy, q
    DOUBLE PRECISION, DIMENSION(:), INTENT(in):: as, bs
    DOUBLE PRECISION:: aas, bbs
    DOUBLE PRECISION, DIMENSION(:,:), INTENT(out)::s, R
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: peso,Omegamu, Omegaeta, ni, gama, VetSol
    DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE::  U, V, phiy, U2, V2, phix

    M=(n*(n+2))/2!numero total de direcoes
    j=M/2 !numero de direcoes em cada sentido
    OPEN(unit=10,file='adobi.txt')
    ALLOCATE(peso(M), Omegamu(j), Omegaeta(j), U(j,j), V(j,j), phiy(M, j), U2(j,j), V2(j,j), phix(M,j), &
         ni(j), gama(j), VetSol(8*M))
    aas=as(1)
    bbs=bs(1)
    CALL nosepesosado(n,M, Omegamu, Omegaeta, peso)

    CALL eigens(j, Omegamu, sigma_t, sigma_s, peso, U, V,ni)

    CALL phis(M, j, U, V, phiy)

    CALL eigens(j, Omegaeta, sigma_t, sigma_s, peso, U2, V2, gama)

    CALL phis(M, j, U2, V2, phix)

    CALL matsist(M,j, Omegamu, a, phiy, aas, ni, Omegaeta, b, phix, bbs, gama, sigma_s, sigma_t,q, peso, VetSol)

    CALL solutions(M, j,kk,jj, phiy, ni, a, phix, gama, b, aas, bbs, peso,VetSol, hx, hy, s)
    CALL fluxo_regioes(jj,kk, hx, hy,ndx, ndy,as, bs, s, R)

  END SUBROUTINE ado_bi

  !=================================================================================================
  SUBROUTINE si_bi(n, j, k, hx, hy, tol, sigma_t, sigma_s0,sigma_s1, q, psi_mediox, psi_medioy, psi_medioxy, s, nit,&
       as, bs,ndx, ndy,R)
    !=================================================================================================
    IMPLICIT NONE
    INTEGER, INTENT(in)::n, j, k, ndx, ndy
    DOUBLE PRECISION::erro
    DOUBLE PRECISION, INTENT(in)::tol, hx, hy
    DOUBLE PRECISION, DIMENSION(:), INTENT(in)::as, bs
    DOUBLE PRECISION, DIMENSION(:,:), INTENT(in)::sigma_t, sigma_s0,sigma_s1, q
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(inout)::psi_mediox, psi_medioy, psi_medioxy
    DOUBLE PRECISION, DIMENSION(:,:), INTENT(inout):: s
    INTEGER, INTENT(out):: nit
    DOUBLE PRECISION, DIMENSION(:,:), INTENT(out):: R
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE::mu, peso
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: Omega, J1, J2
    !        integer:: m, y, x
    !print*, 'n/2, (n*(n+2))/2,(n*(n+2))/8, j, k',n/2, (n*(n+2))/2,(n*(n+2))/8, j, k
    ALLOCATE (mu(n/2),peso((n*(n+2))/2),Omega((n*(n+2))/8,2), J1(j,k), J2(j,k))

    J1=0.d0
    J2=0.d0

    CALL nosepesossi( n, peso, omega)
    DO nit=1,5000
       CALL  varredura(n, j, k, hx, hy, Omega, sigma_t, sigma_s0,sigma_s1, s,J1, J2, q, psi_mediox, psi_medioy,&
            psi_medioxy)

       !$$$$$$                 do m=1, n*(n+2)/2
       !$$$$$$                     do y=1, k
       !$$$$$$                         do x=1, j
       !$$$$$$                             write(10, *)m, x, y, psi_medioxy(m, x, y)
       !$$$$$$                         enddo
       !$$$$$$                     end do
       !$$$$$$                 end do


       CALL calcula_erro(n, j, k, psi_medioxy, erro,peso, s)
       !$$$$$$       print*, 'nit, fluxo', nit, s(1,1)
       !$$$$$$     pause
       IF (erro<tol) THEN
          EXIT
       END IF

       CALL corrente(n, j, k, peso, Omega, psi_medioxy, J1, J2)
    END DO

    CALL fluxo_regioes(j,k, hx, hy,ndx, ndy,as, bs, s,R)

  END SUBROUTINE si_bi

  !=================================================================================================
  SUBROUTINE nosepesossi( n, peso, omega)
    !=================================================================================================
    IMPLICIT NONE
    INTEGER, INTENT(in)::n
    DOUBLE PRECISION, DIMENSION(:), INTENT(out) :: peso
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: w, mu
    DOUBLE PRECISION, DIMENSION(:,:), INTENT(out):: Omega
    INTEGER :: i, l, x
    ALLOCATE (mu(n/2))
    DO
       SELECT CASE (n)
       CASE (2)
          ALLOCATE (w(1))
          mu(1)=0.5773502; w(1)=1.d0; peso=w(1)
          EXIT
       CASE (4)
          ALLOCATE (w(1))
          mu(1)=0.3500212; mu(2)=0.8688903; w(1)=0.3333333; peso=w(1)
          EXIT
       CASE (6)
          ALLOCATE (w(2))
          mu(1)=0.2666355; mu(2)=0.6815076; mu(3)=0.9261808; w(1)=0.1761263; w(2)=0.1572071
          peso(1)=w(1); peso(2)=w(2);peso(3)=w(2);peso(4)=w(1);peso(5)=w(2);peso(6)=w(1);
          DO i=1,6
             peso(n*(n+2)/8+i)=peso(i)
             peso(n*(n+2)/4+i)=peso(i)
             peso(3*n*(n+2)/8+i)=peso(i)
          END DO
          EXIT
       CASE (8)
          ALLOCATE (w(3))
          mu(1)=0.2182179; mu(2)=0.5773503; mu(3)=0.7867958; mu(4)=0.9511897; w(1)=0.1209877;
          w(2)=0.0907407; w(3)=0.0925926
          peso(1)=w(1);peso(2)=w(2);peso(3)=w(2);peso(4)=w(2);peso(5)=w(3);peso(6)=w(2);
          peso(7)=w(1);peso(8)=w(2);peso(9)=w(2); peso(10)=w(1);
          DO i=1,10
             peso(n*(n+2)/8+i)=peso(i)
             peso(n*(n+2)/4+i)=peso(i)
             peso(3*n*(n+2)/8+i)=peso(i)
          END DO
          EXIT
       CASE (12)
          ALLOCATE (w(5))
          mu(1)=0.1672126; mu(2)=0.4595476; mu(3)=0.6280191; mu(4)=0.7600210; mu(5)=0.8722706;
          mu(6)=0.9716377; w(1)=0.0707626; w(2)=0.0558811; w(3)=0.0373377; w(4)=0.0502819;
          w(5)=0.0258513
          peso(1)=w(1);peso(16)=w(1);peso(21)=w(1); peso(2)=w(2);peso(11)=w(2);peso(15)=w(2);
          peso(17)=w(2);peso(20)=w(2);peso(3)=w(2); peso(4)=w(3);peso(6)=w(3);peso(7)=w(3);
          peso(10)=w(3);peso(18)=w(3);peso(19)=w(3); peso(5)=w(4);peso(12)=w(4);peso(14)=w(4);
          peso(8)=w(5);peso(9)=w(5);peso(13)=w(5);
          DO i=1,21
             peso(n*(n+2)/8+i)=peso(i)
             peso(n*(n+2)/4+i)=peso(i)
             peso(3*n*(n+2)/8+i)=peso(i)
          END DO
          EXIT
       CASE (16)
          ALLOCATE (w(8))
          mu(1)=0.1389568; mu(2)=0.3922893; mu(3)=0.5370966; mu(4)=0.6504264; mu(5)=0.7467506;
          mu(6)=0.8319966;mu(7)=0.9092855;mu(8)=0.9805009; w(1)=0.0489872; w(2)=0.0413296;
          w(3)=0.0212326; w(4)=0.0256207; w(5)=0.0360486;w(6)=0.0144589;w(7)=0.0344958; w(8)=0.0085179
          peso(1)=w(1);peso(2)=w(2);peso(3)=w(2);peso(4)=w(3);peso(5)=w(5);peso(6)=w(3);
          peso(7)=w(4);peso(8)=w(6);peso(9)=w(6);peso(10)=w(4);peso(11)=w(4);peso(12)=w(7);
          peso(13)=w(8);peso(14)=w(7);peso(15)=w(4);peso(16)=w(3);peso(17)=w(6);peso(18)=w(8);
          peso(19)=w(8);peso(20)=w(6);peso(21)=w(3);peso(22)=w(2);peso(23)=w(5);peso(24)=w(6);
          peso(25)=w(7);peso(26)=w(6);peso(27)=w(5);peso(28)=w(2);peso(29)=w(1);peso(30)=w(2);
          peso(31)=w(3);peso(32)=w(4);peso(33)=w(4);peso(34)=w(3);peso(35)=w(2);peso(36)=w(1);
          DO i=1,36
             peso(n*(n+2)/8+i)=peso(i)
             peso(n*(n+2)/4+i)=peso(i)
             peso(3*n*(n+2)/8+i)=peso(i)
          END DO
          EXIT
       END SELECT
    END DO

    !preenche o vetor de posicoes mu e eta no primeiro quadrante
    l=0
    DO i=n/2,1, -1
       DO x=1, n/2-i+1
          l=l+1
          Omega(l,1)=mu(i)
          Omega(l,2)=mu(x)
          !print*, 'omega', Omega(l,1), Omega(l,2)
       END DO
    END DO

  END SUBROUTINE nosepesossi

  !=================================================================================================
  SUBROUTINE varredura(n, j, k,hx, hy, Omega, sigma_t, sigma_s0, sigma_s1, s,J1, J2, q, psi_mediox, psi_medioy, psi_medioxy)
    !=================================================================================================
    IMPLICIT NONE
    INTEGER, INTENT(in)::n, j, k
    DOUBLE PRECISION, INTENT(in):: hx, hy
    DOUBLE PRECISION, DIMENSION(:,:), INTENT(in):: Omega, sigma_t,sigma_s0,sigma_s1, s, J1, J2, q
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(inout) :: psi_medioy,psi_mediox
    DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(inout) :: psi_medioxy
    INTEGER :: y, x, m

    !varredura do terceiro quadrante
    DO y=k,1,-1
       DO x=j,1,-1
          DO m=1,n*(n+2)/8
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
          END DO
       END DO
    END DO

    !varredura do quarto quadrante
    DO y=k,1,-1
       DO x=1,j
          DO m=1,n*(n+2)/8
             psi_medioxy(3*n*(n+2)/8+m,x,y)=((2*Omega(m,1)/hx)*psi_medioy(3*n*(n+2)/8+m,x,y)+&
                  (2*Omega(m,2)/hy)*psi_mediox(3*n*(n+2)/8+m, x,y+1)+&
                  sigma_s0(x,y)*s(x,y)+0.75*sigma_s1(x,y)*(Omega(m,1)*J1(x,y)&
                  -Omega(m,2)*J2(x,y))+q(x,y))/&
                  (sigma_t(x,y)+(2*Omega(m,1)/hx)+(2*Omega(m,2)/hy))
             psi_medioy(3*n*(n+2)/8+m,x+1,y)=2*psi_medioxy(3*n*(n+2)/8+m,x,y)-psi_medioy(3*n*(n+2)/8+m,x,y)
             psi_mediox(3*n*(n+2)/8+m,x,y)=2*psi_medioxy(3*n*(n+2)/8+m,x,y)-psi_mediox(3*n*(n+2)/8+m,x,y+1)
             !$$$$$$             print*, 'psi_medio 4q', 3*n*(n+2)/8+m,x,y,psi_medioxy(3*n*(n+2)/8+m, x, y)
             psi_mediox(m,x,1) = psi_mediox(3*n*(n+2)/8+m,x,1)
          END DO
       END DO
    END DO

    !varredura do segundo quadrante
    DO y=1,k
       DO x=j,1,-1
          DO m=1,n*(n+2)/8
             psi_medioxy(n*(n+2)/8+m,x,y)=((2*Omega(m,1)/hx)*psi_medioy(n*(n+2)/8+m, x+1,y)+&
                  (2*Omega(m,2)/hy)*psi_mediox(n*(n+2)/8+m, x,y)+sigma_s0(x,y)*s(x,y)+&
                  0.75*sigma_s1(x,y)*(-Omega(m,1)*J1(x,y)+Omega(m,2)*J2(x,y))+q(x,y))/&
                  (sigma_t(x,y)+(2*Omega(m,1)/hx)+(2*Omega(m,2)/hy))
             psi_medioy(n*(n+2)/8+m,x,y)=2*psi_medioxy(n*(n+2)/8+m,x,y)-psi_medioy(n*(n+2)/8+m,x+1,y)
             psi_mediox(n*(n+2)/8+m,x,y+1)=2*psi_medioxy(n*(n+2)/8+m,x,y)-psi_mediox(n*(n+2)/8+m,x,y)
             !$$$$$$            print*, 'psi_medio 2q', n*(n+2)/8+m,x,y,psi_medioxy(n*(n+2)/8+m, x, y)
             psi_medioy(m,1, y) = psi_medioy(n*(n+2)/8+m,1,y)
          END DO
       END DO
    END DO

    !varredura no primeiro quadrante
    DO y=1,k
       DO x=1,j
          DO m=1,n*(n+2)/8
             psi_medioxy(m,x,y)=((2*Omega(m,1)/hx)*psi_medioy(m,x,y)+(2*Omega(m,2)/hy)*psi_mediox(m, x,y)+&
                  sigma_s0(x,y)*s(x,y)+0.75*sigma_s1(x,y)*(Omega(m,1)*J1(x,y)+Omega(m,2)*J2(x,y))&
                  +q(x,y))/(sigma_t(x,y)+(2*Omega(m,1)/hx)+(2*Omega(m,2)/hy))
             psi_medioy(m,x+1,y)=2*psi_medioxy(m,x,y)-psi_medioy(m,x,y)
             psi_mediox(m,x,y+1)=2*psi_medioxy(m,x,y)-psi_mediox(m,x,y)
             !$$$$$$             print*, 'psi_medioxy 1q', m,x,y,psi_medioxy(m, x, y)
          END DO
       END DO
    END DO
  END SUBROUTINE varredura
  !==============================================================================================================
  SUBROUTINE calcula_erro(n, j, k, psi_medioxy, erro,peso, s)
    !==============================================================================================================
    IMPLICIT NONE
    INTEGER, INTENT(in):: k, j, n
    DOUBLE PRECISION, INTENT(out) :: erro
    DOUBLE PRECISION, DIMENSION(:), INTENT(in):: peso
    DOUBLE PRECISION, DIMENSION(:,:, :), INTENT(in) ::psi_medioxy
    DOUBLE PRECISION, DIMENSION(:,:), INTENT(out)::s
    INTEGER:: x, y, m
    DOUBLE PRECISION:: phi_ant
    erro = -HUGE(erro)
    DO y=1,k
       DO x=1,j
          phi_ant=s(x,y)
          s(x,y)=0.d0
          DO m=1,n*(n+2)/2
             s(x,y)=s(x,y)+psi_medioxy(m,x,y)*peso(m)
          END DO
          s(x,y)=s(x,y)/4.d0
          !print*,'x, y, s(x,y)',x, y, s(x,y)
          !erro=max(erro,abs(s(x,y)-phi_ant))
          !print*, 'x, y, max, s', x, y, max(erro,abs(s(x,y)-phi_ant)), abs(s(x,y))
          erro=MAX(erro,ABS(s(x,y)-phi_ant)/ABS(s(x,y)))
       END DO
    ENDDO
    !pause
  END SUBROUTINE calcula_erro


  !==============================================================================================================
  SUBROUTINE corrente(n, j, k, peso, Omega, psi_medioxy, J1, J2)
    !==============================================================================================================
    IMPLICIT NONE
    INTEGER, INTENT(in):: n, j, k
    DOUBLE PRECISION, DIMENSION(:), INTENT(in)::peso
    DOUBLE PRECISION, DIMENSION(:, :), INTENT(in)::Omega
    DOUBLE PRECISION, DIMENSION(:,:, :), INTENT(in) ::psi_medioxy
    DOUBLE PRECISION, DIMENSION(:,:), INTENT(out)::J1, J2
    INTEGER:: y, x, m
    DO y=1,k
       DO x=1,j
          J1(x,y)=0.d0
          J2(x,y)=0.d0
          DO m=1,n*(n+2)/8
             J1(x,y)=J1(x,y)+Omega(m, 1)*psi_medioxy(m,x,y)*peso(m)-Omega(m, 1)*psi_medioxy(n*(n+2)/8+m,x,y)*&
                  peso(n*(n+2)/8+m)-Omega(m, 1)*psi_medioxy(n*(n+2)/4+m,x,y)*peso(n*(n+2)/4+m)+Omega(m, 1)*&
                  psi_medioxy(3*n*(n+2)/8+m,x,y)*peso(3*n*(n+2)/8+m)
             J2(x,y)=J2(x,y)+Omega(m, 2)*psi_medioxy(m,x,y)*peso(m)+Omega(m, 2)*psi_medioxy(n*(n+2)/8+m,x,y)*&
                  peso(n*(n+2)/8+m)-Omega(m, 2)*psi_medioxy(n*(n+2)/4+m,x,y)*peso(n*(n+2)/4+m)-Omega(m, 2)*&
                  psi_medioxy(3*n*(n+2)/8+m,x,y)*peso(3*n*(n+2)/8+m)
          END DO

       END DO
    ENDDO
  END SUBROUTINE corrente

  !==============================================================================================================
  SUBROUTINE fluxo_regioes(j,k, hx, hy, ndx, ndy, as, bs, s, R)
    !==============================================================================================================
    IMPLICIT NONE
    INTEGER, INTENT(in):: j, k, ndx, ndy
    DOUBLE PRECISION, INTENT(in)::hx, hy
    DOUBLE PRECISION, DIMENSION(:), INTENT(in):: as, bs
    DOUBLE PRECISION, DIMENSION(:,:), INTENT(in):: s
    DOUBLE PRECISION, DIMENSION(:,:), INTENT(out):: R
    INTEGER, DIMENSION(:), ALLOCATABLE:: x,y
    INTEGER::i, l, m, n
    DOUBLE PRECISION:: A

    ALLOCATE(x(0:ndx), y(0:ndy))
    R=0.d0
    x(0)=0
    y(0)=0
    DO i=1, ndx
       x(i)=NINT(as(i)/hx)
    END DO
    DO i=1, ndy
       y(i)=NINT(bs(i)/hy)
    END DO

    !i, l correm nas regioes
    !m, n sao indices para o somatorio em cada celula da malha em cada regiao

    DO i=1, ndx
       DO l=1, ndy
          A=0
          DO m=y(l-1)+1,y(l)
             DO n=x(i-1)+1,x(i)
                A=A+hx*hy
                R(i, l)=R(i, l)+s(n, m)*hx*hy
             END DO
          END DO
          R(i, l)=R(i, l)/A
       END DO
    END DO

  END SUBROUTINE fluxo_regioes

  !==========================================================
  SUBROUTINE digse(ndx, ndy, C, R, dig)
    !============================================================
    IMPLICIT NONE
    INTEGER, INTENT(in)::ndx, ndy
    DOUBLE PRECISION, DIMENSION(:,:), INTENT(in):: C, R
    DOUBLE PRECISION, DIMENSION(:,:), INTENT(out)::dig
    INTEGER::i, j
    dig=0

    DO i=1, ndx
       DO j=1, ndy
          dig(i,j)=-LOG(ABS(C(i,j)-R(i,j))/ABS(C(i,j)))
       END DO
    END DO

  END SUBROUTINE digse
  !====================================================
  SUBROUTINE nosepesosado(n,M,Omegamu, Omegaeta, peso)
    !====================================================
    IMPLICIT NONE
    INTEGER, INTENT(in):: n,M
    DOUBLE PRECISION, DIMENSION(:), INTENT(out)::peso
    DOUBLE PRECISION, DIMENSION(:), INTENT(out):: Omegamu, Omegaeta
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE:: w, mu
    INTEGER:: l,i,x

    ALLOCATE(mu(n/2))

    DO
       SELECT CASE (n)
       CASE (2)
          ALLOCATE (w(1))
          mu(1)=0.5773502d0; w(1)=1.d0; peso=w(1)
          EXIT
       CASE (4)
          ALLOCATE (w(1))
          mu(1)=0.3500212d0; mu(2)=0.8688903d0;w(1)=0.333333d0; peso=w(1);! w(1)=1.d0/3; peso=w(1)
          EXIT
       CASE (6)
          ALLOCATE (w(2))
          mu(1)=0.2666355d0; mu(2)=0.6815076d0; mu(3)=0.9261808d0; w(1)=0.1761263d0; w(2)=0.1572071d0
          peso(1)=w(1); peso(2)=w(2);peso(3)=w(2);peso(4)=w(1);peso(5)=w(2);peso(6)=w(1);
          DO i=1,6
             peso(M/4+i)=peso(i)
             peso(M/2+i)=peso(i)
             peso(3*M/4+i)=peso(i)
          END DO
          EXIT
       CASE (8)
          ALLOCATE (w(3))
          mu(1)=0.2182179; mu(2)=0.5773503; mu(3)=0.7867958; mu(4)=0.9511897; w(1)=0.1209877;
          w(2)=0.0907407; w(3)=0.0925926
          peso(1)=w(1);peso(2)=w(2);peso(3)=w(2);peso(4)=w(2);peso(5)=w(3);peso(6)=w(2);
          peso(7)=w(1);peso(8)=w(2);peso(9)=w(2); peso(10)=w(1);
          DO i=1,10
             peso(M/4+i)=peso(i)
             peso(M/2+i)=peso(i)
             peso(3*M/4+i)=peso(i)
          END DO
          EXIT
       CASE (12)
          ALLOCATE (w(5))
          mu(1)=0.1672126; mu(2)=0.4595476; mu(3)=0.6280191; mu(4)=0.7600210; mu(5)=0.8722706;
          mu(6)=0.9716377; w(1)=0.0707626; w(2)=0.0558811; w(3)=0.0373377; w(4)=0.0502819;
          w(5)=0.0258513
          peso(1)=w(1);peso(16)=w(1);peso(21)=w(1); peso(2)=w(2);peso(11)=w(2);peso(15)=w(2);
          peso(17)=w(2);peso(20)=w(2);peso(3)=w(2); peso(4)=w(3);peso(6)=w(3);peso(7)=w(3);
          peso(10)=w(3);peso(18)=w(3);peso(19)=w(3); peso(5)=w(4);peso(12)=w(4);peso(14)=w(4);
          peso(8)=w(5);peso(9)=w(5);peso(13)=w(5);
          DO i=1,21
             peso(M/4+i)=peso(i)
             peso(M/2+i)=peso(i)
             peso(3*M/4+i)=peso(i)
          END DO
          EXIT
       CASE (16)
          ALLOCATE (w(8))
          mu(1)=0.1389568; mu(2)=0.3922893; mu(3)=0.5370966; mu(4)=0.6504264; mu(5)=0.7467506;
          mu(6)=0.8319966;mu(7)=0.9092855;mu(8)=0.9805009; w(1)=0.0489872; w(2)=0.0413296;
          w(3)=0.0212326; w(4)=0.0256207; w(5)=0.0360486;w(6)=0.0144589;w(7)=0.0344958; w(8)=0.0085179
          peso(1)=w(1);peso(2)=w(2);peso(3)=w(2);peso(4)=w(3);peso(5)=w(5);peso(6)=w(3);
          peso(7)=w(4);peso(8)=w(6);peso(9)=w(6);peso(10)=w(4);peso(11)=w(4);peso(12)=w(7);
          peso(13)=w(8);peso(14)=w(7);peso(15)=w(4);peso(16)=w(3);peso(17)=w(6);peso(18)=w(8);
          peso(19)=w(8);peso(20)=w(6);peso(21)=w(3);peso(22)=w(2);peso(23)=w(5);peso(24)=w(6);
          peso(25)=w(7);peso(26)=w(6);peso(27)=w(5);peso(28)=w(2);peso(29)=w(1);peso(30)=w(2);
          peso(31)=w(3);peso(32)=w(4);peso(33)=w(4);peso(34)=w(3);peso(35)=w(2);peso(36)=w(1);
          DO i=1,36
             peso(M/4+i)=peso(i)
             peso(M/2+i)=peso(i)
             peso(3*M/4+i)=peso(i)
          END DO
          EXIT
       END SELECT
    END DO


    !preenche o vetor de posicoes mu e eta no primeiro quadrante
    l=0
    DO i=n/2,1, -1
       DO x=1, n/2-i+1
          l=l+1
          Omegamu(l)=mu(i)
          Omegaeta(l)=mu(x)
          !print*, 'l', l
          !print*, 'omega', Omega(l,1), Omega(l,2)
       END DO
    END DO
    DO i=1, M/4
       omegamu(i+M/4)=omegamu(i)
       omegaeta(i+M/4)=omegaeta(i)!pensando no outro ordenamento
    END DO
    DEALLOCATE(mu, w)
  END SUBROUTINE nosepesosado

  !========================================================
  SUBROUTINE eigens(j, Omegamu, sigma_t, sigma_s, peso, U, V, ni)
    !========================================================
    IMPLICIT NONE
    INTEGER, INTENT(in)::j
    DOUBLE PRECISION, INTENT(in)::sigma_t, sigma_s
    DOUBLE PRECISION, DIMENSION(:), INTENT(in)::peso, Omegamu
    DOUBLE PRECISION, DIMENSION(:), INTENT(out):: ni
    DOUBLE PRECISION, DIMENSION(:,:), INTENT(out)::U, V
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: D, wi,work,scale, rconde, rcondv,wr
    DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE:: A, vr
    DOUBLE PRECISION, DIMENSION(1)::aux
    INTEGER, DIMENSION(1):: iwork
    LOGICAL, DIMENSION(:), ALLOCATABLE :: bwork
    INTEGER::i,k,ilo, ihi, lwork,info
    DOUBLE PRECISION::abnrm


    ALLOCATE(D(j), A(j,j), wi(j),vr(j,j), wr(j))

    D=0.d0
    A=0.d0
    DO i=1, j
       D(i)=(sigma_t/Omegamu(i))**2
    END DO

    !$$$$$$         do i=1, j
    !$$$$$$             print*, 'D', D(i)
    !$$$$$$         end do

    !constroi matriz -A (por conveniencia)
    DO i=1,j
       DO k=1, j
          A(i,k)=-(sigma_t*sigma_s*peso(k))/(2.d0*Omegamu(i)**2)
       END DO
    END DO

    !altera matriz A pra D-A
    DO i=1,j
       A(i,i)=D(i)+A(i,i)
    END DO

    !$$$$$$         do i=1,j
    !$$$$$$             do k=1, j
    !$$$$$$                 print*, 'A, i, j',i,k, A(i,k)
    !$$$$$$             enddo
    !$$$$$$         end do

    lwork = 10*j
    ALLOCATE(work(lwork),bwork(n),SCALE(j),rconde(j),rcondv(j))
    ilo = 1
    ihi = j


    CALL dgeevx('N','N','V','N',j,A,j,wr,wi,aux,1,vr,j,ilo,ihi,scale,abnrm,rconde,rcondv,work,lwork,iwork,info)
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
    DEALLOCATE(work,bwork,scale,rconde,rcondv)

    DO i=1,j
       U(i,:)=vr(i,:)
    END DO

    !$$$$$$         do i=1,j
    !$$$$$$           do k=1, j
    !$$$$$$             print*,'i,k, U',i,k, U(i,k)
    !$$$$$$           end do
    !$$$$$$         end do
    !$$$$$$           pause

    DO i=1,j
       IF(wr(i).LE.0.d0) THEN
          wr(i)=dabs(wr(i))
       END IF
       ni(i)=1.d0/dsqrt(wr(i))
    END DO

    !$$$$$$          print *,'ni=',ni
    !$$$$$$          pause

    DO i=1, j
       DO k=1,j
          V(i,k)=Omegamu(i)*U(i,k)/(ni(k)*sigma_t)
       END DO
    END DO
    !$$$$$$         do i=1, j
    !$$$$$$           do k=1,j
    !$$$$$$              print*, 'V, i, k', i,k, V(i,k)
    !$$$$$$            end do
    !$$$$$$         end do
    DEALLOCATE(D, A, wi,vr, wr)
  END SUBROUTINE eigens

  !=============================================================
  SUBROUTINE phis(M, j, U, V, phiy)
    !=============================================================
    IMPLICIT NONE
    INTEGER, INTENT (in)::M,j
    DOUBLE PRECISION, DIMENSION(:,:), INTENT(in)::U, V
    DOUBLE PRECISION, DIMENSION(:,:), INTENT(out)::phiy
    INTEGER:: i,k

    DO i=1, j
       DO k=1, j
          phiy(i,k)=0.5d0*(U(i,k)+V(i,k))
          phiy(j+i, k)=0.5d0*(U(i,k)-V(i,k))
       END DO
    END DO

    !$$$$$$         do i=1, M
    !$$$$$$             do k=1, j
    !$$$$$$                 print*, 'i,j, phi', i, k, phiy(i,k)
    !$$$$$$             end do
    !$$$$$$         end do
  END SUBROUTINE phis

  !=============================================================
  SUBROUTINE matsist(M,j, Omegamu, a, phiy, as, ni, Omegaeta, b, phix, bs, gama, sigma_s, sigma_t,&
       q,peso, VetSol)
    !=============================================================
    IMPLICIT NONE
    INTEGER, INTENT(in)::M,j
    DOUBLE PRECISION, INTENT(in)::a, as, b, bs, sigma_s, sigma_t,q
    DOUBLE PRECISION, DIMENSION(:), INTENT(in):: Omegamu,ni, Omegaeta, gama, peso
    DOUBLE PRECISION, DIMENSION(:,:), INTENT(in):: phiy, phix
    INTEGER::i,k,z,info
    INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv
    DOUBLE PRECISION, DIMENSION(:), INTENT(out):: VetSol
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: MAT
    z=M/4
    ALLOCATE (MAT(8*M, 8*M))
    MAT=0.d0
    VetSol=0.d0

    !matrizes m/4*m/2
    DO i=1, z
       DO k=1, j
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
          MAT(i,j+k)=-Omegamu(i)/a*phiy(j+i, k)*EXP(-as/ni(k))
          MAT(z+i, j+k)=Omegamu(i)/a*phiy(i, k)*EXP(-as/ni(k))
          MAT(j+i, j+k)=-Omegamu(i)/a*phiy(3*z+i, k)*EXP(-as/ni(k))
          MAT(3*z+i, j+k)=Omegamu(i)/a*phiy(z+i,k)*EXP(-as/ni(k))
          MAT(M+i, j+k)=Mat(i,j+k)
          MAT(M+z+i, j+k)= MAT(z+i, j+k)
          MAT(M+j+i, j+k)=MAT(j+i, j+k)
          MAT(M+3*z+i, j+k)=MAT(3*z+i, j+k)
          !C
          MAT(i,m+k)=Omegamu(i)/a*phiy(i, k)*EXP(-(a-as)/ni(k))
          MAT(j+i, m+k)=Omegamu(i)/a*phiy(z+i, k)*EXP(-(a-as)/ni(k))
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
          MAT(2*m+i,2*m+j+k)=-Omegaeta(i)/b*phix(j+i, k)*EXP(-bs/gama(k))
          MAT(2*m+z+i, 2*m+j+k)=Omegaeta(i)/b*phix(i, k)*EXP(-bs/gama(k))!tinha um menos na frente?tirei ps camila
          MAT(2*m+j+i, 2*m+j+k)=-Omegaeta(i)/b*phix(3*z+i, k)*EXP(-bs/gama(k))
          MAT(2*m+3*z+i, 2*m+j+k)=Omegaeta(i)/b*phix(z+i,k)*EXP(-bs/gama(k))
          MAT(3*M+i, 2*m+j+k)=Mat(2*m+i,2*m+j+k)
          MAT(3*M+z+i, 2*m+j+k)= MAT(2*m+z+i, 2*m+j+k)
          MAT(3*M+j+i, 2*m+j+k)=MAT(2*m+j+i, 2*m+j+k)
          MAT(3*M+3*z+i, 2*m+j+k)=MAT(2*m+3*z+i, 2*m+j+k)
          !G
          MAT(2*m+i,3*m+k)=Omegaeta(i)/b*phix(i, k)*EXP(-(b-bs)/gama(k))
          MAT(2*m+j+i, 3*m+k)=Omegaeta(i)/b*phix(z+i, k)*EXP(-(b-bs)/gama(k))
          MAT(3*M+i, 3*m+k)=MAT(2*m+i,3*m+k)
          MAT(3*M+j+i, 3*m+k)=MAT(2*m+j+i, 3*m+k)
          !H
          MAT(2*m+i,3*m+j+k)=Omegaeta(i)/b*phix(j+i, k)
          MAT(2*m+j+i, 3*m+j+k)=Omegaeta(i)/b*phix(j+z+i, k)
          MAT(3*M+i, 3*m+j+k)=Mat(2*m+i,3*m+j+k)
          MAT(3*M+j+i, 3*m+j+k)=MAT(2*m+j+i, 3*m+j+k)
       END DO
    ENDDO
    !matrizes m/2*m/2
    DO i=1, j
       DO k=1, j
          !A
          MAT(5*m+i,k)=phiy(i, k)-phiy(j+i, k)
          Mat(6*M+i,k)=phiy(i, k)*EXP(-as/ni(k))
          Mat(6*M+j+i, k)=phiy(j+i, k)*EXP(-as/ni(k))
          !B
          MAT(5*m+i,j+k)=-Mat(5*m+i,k)*EXP(-as/ni(k))
          Mat(6*M+i,j+k)=phiy(j+i, k)
          Mat(6*M+j+i, j+k)=phiy(i, k)
          !C
          MAT(4*M+i, m+k)=phiy(j+i, k)*EXP(-(a-as)/ni(k))
          Mat(6*M+i,m+k)=-phiy(i, k)
          Mat(6*M+j+i, m+k)=-phiy(j+i, k)
          !D
          MAT(4*M+i, m+j+k)=phiy(i, k)
          Mat(6*M+i,m+j+k)=-phiy(j+i, k)*EXP(-(a-as)/ni(k))
          Mat(6*M+j+i, m+j+k)=-phiy(i, k)*EXP(-(a-as)/ni(k))
          !E
          MAT(5*m+j+i,2*m+k)=phix(i, k)-phix(j+i, k)
          Mat(7*M+i,2*m+k)=phix(i, k)*EXP(-bs/gama(k))
          Mat(7*M+j+i, 2*m+k)=phix(j+i, k)*EXP(-bs/gama(k))
          !F
          MAT(5*m+j+i, 2*m+j+k)=-MAT(5*m+j+i,2*m+k)*EXP(-bs/gama(k))
          Mat(7*M+i, 2*m+j+k)=phix(j+i, k)
          Mat(7*M+j+i,  2*m+j+k)=phix(i, k)
          !G
          MAT(4*M+j+i, 3*m+k)=phix(j+i, k)*EXP(-(b-bs)/gama(k))
          Mat(7*M+i,3*m+k)=-phix(i, k)
          Mat(7*M+j+i, 3*m+k)=-phix(j+i, k)
          !H
          MAT(4*M+j+i, 3*m+j+k)=phix(i, k)
          Mat(7*M+i,3*m+j+k)=-phix(j+i, k)*EXP(-(b-bs)/gama(k))
          Mat(7*M+j+i, 3*m+j+k)=-phix(i, k)*EXP(-(b-bs)/gama(k))
       END DO
    ENDDO

    !matrizes m*m/4
    DO i=1, z
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
    END DO
    !matriz mxm
    DO i=1, m
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
    END DO

    !matrizes m/2xm
    DO i=1,j
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
    ENDDO
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
    ALLOCATE (ipiv(8*M))
    CALL dgetrf (8*M,8*M,Mat,8*M,ipiv,info)
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
    CALL dgetrs('N', 8*M , 1, MAT , 8*M, ipiv, VetSol, 8*M, info)
    DEALLOCATE (ipiv)
    !$$$$$$            do i=1, 8*m
    !$$$$$$              print*, 'i, VetSol(i)',i, VetSol(i)
    !$$$$$$            end do
    !$$$$$$             pause
    DEALLOCATE(MAT)

  END SUBROUTINE matsist

  !==============================================================================================================================
  SUBROUTINE solutions(M, j,kk,jj, phiy, ni, a, phix, gama, b, as, bs, peso,VetSol, hx, hy, s)
    !==========================================================================================================================
    IMPLICIT NONE
    INTEGER, INTENT(in):: M, j, kk, jj
    DOUBLE PRECISION, INTENT(in):: a, b, as, bs, hx, hy
    DOUBLE PRECISION, DIMENSION(:), INTENT(in)::  ni, gama, peso, VetSol
    DOUBLE PRECISION, DIMENSION (:,:), INTENT(in)::phiy, phix
    DOUBLE PRECISION, DIMENSION (:,:), INTENT(out):: s
    DOUBLE PRECISION :: x,y, psiy1, psiy2, psix1, psix2, fluxescx, fluxescy, r1, r2, r3, r4
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: psiy, psix,psixsource, psiysource
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: phixx, phiyy
    INTEGER:: i, k,l, p

    ALLOCATE(psiy(M), psix(M), psiysource(M), psixsource(M), phixx(jj,kk), phiyy(jj,kk))
    x=0.d0
    y=0.d0
    s=0.d0
    DO  l=1, kk
       DO p=1, jj
          x=hx*(2*p-1)/2.d0
          y=hy*(2*l-1)/2.d0

          psiy=0.d0
          psix=0.d0

          !solu�cao na fonte
          IF (x<=as) THEN
             DO i=1, j
                psix1=0.d0
                psix2=0.d0
                DO k=1,j
                   psix1=VetSol(2*M+k)*phix(i, k)*EXP(-y/gama(k))+VetSol(2*M+j+k)*phix(j+i, k)*EXP(-(bs-y)/gama(k))&
                        +psix1
                   psix2=VetSol(2*M+k)*phix(j+i, k)*EXP(-y/gama(k))+VetSol(2*m+j+k)*phix(i, k)*EXP(-(bs-y)/gama(k))&
                        +psix2
                END DO
                psix(i)=psix1+VetSol(6*M+i)
                psix(j+i)=psix2+VetSol(6*M+j+i)
             END DO
          END IF
          IF (y<=bs) THEN
             DO i=1, j
                psiy1=0.d0
                psiy2=0.d0
                DO k=1,j
                   psiy1=VetSol(k)*phiy(i, k)*EXP(-x/ni(k))+VetSol(j+k)*phiy(j+i, k)*EXP(-(as-x)/ni(k))+psiy1
                   psiy2=VetSol(k)*phiy(j+i, k)*EXP(-x/ni(k))+VetSol(j+k)*phiy(i, k)*EXP(-(as-x)/ni(k))+psiy2
                END DO
                psiy(i)=psiy1+VetSol(4*m+i)
                psiy(j+i)=psiy2+VetSol(4*m+j+i)
             END DO
          END IF
          IF (x>=as) THEN
             !solucao sem fonte
             DO i=1, j
                psix1=0.d0
                psix2=0.d0
                DO k=1,j
                   psix1=VetSol(3*M+k)*phix(i, k)*EXP(-(y-bs)/gama(k))+&
                        &VetSol(3*M+j+k)*phix(j+i, k)*EXP(-(b-y)/gama(k))+psix1
                   psix2=VetSol(3*M+k)*phix(j+i, k)*EXP(-(y-bs)/gama(k))+&
                        &VetSol(3*M+j+k)*phix(i, k)*EXP(-(b-y)/gama(k))+psix2
                END DO
                psix(i)=psix1+VetSol(7*M+i)
                psix(j+i)=psix2+VetSol(7*M+j+i)
             END DO
          END IF
          IF (y>=bs) THEN
             !solucao sem fonte
             DO i=1, j
                psiy1=0.d0
                psiy2=0.d0
                DO k=1,j
                   psiy1=VetSol(M+k)*phiy(i, k)*EXP(-(x-as)/ni(k))+VetSol(M+j+k)*phiy(j+i, k)*EXP(-(a-x)/ni(k))+psiy1
                   psiy2=VetSol(M+k)*phiy(j+i, k)*EXP(-(x-as)/ni(k))+VetSol(M+j+k)*phiy(i, k)*EXP(-(a-x)/ni(k))+psiy2
                END DO
                psiy(i)=psiy1+VetSol(5*m+i)
                psiy(j+i)=psiy2+VetSol(5*m+j+i)
             END DO
          END IF
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
          DO i=1, m
             fluxescx=.25*psix(i)*peso(i)+fluxescx
             fluxescy=.25*psiy(i)*peso(i)+fluxescy
          END DO
          !$$$$$$                     phixx(p,l)=fluxescx
          !$$$$$$                     phiyy(p,l)=fluxescy
          s(p,l)=(fluxescx+fluxescy)/2.d0
          !                print*, 'fluxx, fluxy, x, y, s', fluxescx,  fluxescy,x, y, s(p,l)
       END DO
    END DO
    !$$$$$$             call fluxo_regioes(jj,kk, hx, hy,as, bs, phixx, R1, R2, R3, R4)
    !$$$$$$             print*, 'phix:r1, r2, r3, r4', r1, r2, r3, r4
    !$$$$$$             call fluxo_regioes(jj,kk, hx, hy,as, bs, phiyy, R1, R2, R3, R4)
    !$$$$$$             print*, 'phiy:r1, r2, r3, r4', r1, r2, r3, r4
  END SUBROUTINE solutions

  !============================================================================================
  DOUBLE PRECISION FUNCTION diferenca_tempo(t0,t1) RESULT (valor)
    !============================================================================================
    IMPLICIT NONE
    INTEGER, DIMENSION(8), INTENT(in) :: t0,t1
    ! t0 e t1 sao valores retornados pela chamada do intr�nseco DATE_AND_TIME()
    ! valor e a diferenca de tempo, em segundos
    ! somente pode ser usada para medicao de tempo de trechos de programa que executem
    ! em menos de um mes
    DOUBLE PRECISION :: t0_ms,t1_ms
    t0_ms = t0(3)*24*3600.0D3+t0(5)*3600.0D3+t0(6)*60.0D3+t0(7)*1.0D3+t0(8)
    t1_ms = t1(3)*24*3600.0D3+t1(5)*3600.0D3+t1(6)*60.0D3+t1(7)*1.0D3+t1(8)
    valor = (t1_ms-t0_ms)*1.0E-3
  END FUNCTION diferenca_tempo

END PROGRAM ado_bidimensional
