!     Last change:  MC   30 Jul 2008   12:47 pm
module ShapeFile
USE ByteStreamer
implicit none
private

TYPE PolyLine
  REAL(KIND=8), DIMENSION(4)         :: Box         !Bounding Box Xmin, Ymin, Xmax, Ymax
  INTEGER     , DIMENSION(:), POINTER:: PartIndex => NULL()   !Index of first points in Part. Array indexes are with respect 0
  TYPE(point) , DIMENSION(:), POINTER:: Points => NULL()      !points for all Parts
end type

TYPE Point
  REAL(KIND=8):: x
  REAL(KIND=8):: y
end type



INTEGER, DIMENSION(:,:), ALLOCATABLE, private:: Indice

INTEGER:: FileLastByte     ! ultimo Byte del file
INTEGER:: NRecIndex        ! numero di oggetti nello SHP
INTEGER:: ShapeType        ! tipo di oggetti  1Point 3Polyline
INTEGER:: Versione         ! 1000

INTEGER:: SHPunit= 0
INTEGER:: SHXunit= 0
INTEGER:: NextRecordFileAddr=101      ! puntatore al prox record (espresso in bytes, primo byte del file= 1)

! usato solo in scrittura
INTEGER:: numRecord                   ! numero di record scritti (usato solo in scrittura)
REAL(KIND=8):: Xmin1, Xmax1, Ymin1, Ymax1

interface ReadShapeOgge
  module procedure ReadPoly
  module procedure ReadPoint
end interface

PUBLIC:: PolyLine, Point

PUBLIC:: Head
PUBLIC:: EstremiShape
PUBLIC:: IllustraSHP
PUBLIC:: NumOggetti
PUBLIC:: TipoShape
PUBLIC:: ReadShapeOgge
PUBLIC:: Termina

PUBLIC:: IniziaScrittura
PUBLIC:: ScriviPolyLine
PUBLIC:: ScriviPoint
PUBLIC:: ChiudiScrittura

contains

subroutine Termina()
CLOSE(SHPunit)
if (ALLOCATED(Indice) ) DEALLOCATE(Indice)
end subroutine


subroutine IllustraSHP()
REAL(KIND=8):: Xmin, Xmax, Ymin, Ymax
CHARACTER(LEN= 16):: tipo
call TipoShape(tipo)
WRITE(*,*) "Versione ", Versione
WRITE(*,*) "Numero di oggetti ", NRecIndex
WRITE(*,*) "File length ",  FileLastByte, " bytes (8bit words)"
WRITE(*,"(1X,A,I2,A)") "Shape type "//TRIM(tipo)//" (", ShapeType,")"
call EstremiShape(Xmin, Xmax, Ymin, Ymax)
WRITE(*,*) "Range coordinata X ", Xmin, Xmax
WRITE(*,*) "Range coordinata Y ", Ymin, Ymax
end subroutine

function NumOggetti()
INTEGER:: NumOggetti
NumOggetti = NRecIndex
end function


subroutine Head(nome, esito)
CHARACTER(LEN=*), INTENT(IN):: nome
INTEGER, INTENT(OUT), OPTIONAL:: esito !0 Ok
INTEGER(KIND=4):: i, k1, k2, ibyte
NextRecordFileAddr= 101

SHPunit = UnitaLibera()
if (.not. Stream_Open(SHPunit, FileName= TRIM(nome)//".SHP", scrittura= .false. , stato= 1)) then
   if (PRESENT(esito) ) esito= 1
   return
endif

SHXunit = UnitaLibera()
if (.not. Stream_Open(SHXunit, FileName= TRIM(nome)//".SHX", scrittura= .false. , stato= 1)) then
   if (PRESENT(esito) ) esito= 2
   return
endif

i = Stream_Read_Integer(SHPunit, nbytes=4, posizione=1)
if ( 9994 /= BigToLittle(i) ) then
   if (PRESENT(esito) ) esito= 3
   return
endif

i = Stream_Read_Integer(SHPunit, nbytes=4, posizione=25)   ! File length in 16bit words
FileLastByte=  BigToLittle(i)* 2                           ! File length in  8bit words (bytes)

Versione = Stream_Read_Integer(SHPunit, nbytes=4, posizione=29)   ! 1000


ShapeType = Stream_Read_Integer(SHPunit, nbytes=4, posizione= 33)   ! in uso 1 Point 3 Polyline

i = Stream_Read_Integer(SHXunit, nbytes=4, posizione=25)            ! File length in 16bit words
NRecIndex=  BigToLittle(i)* 2      ! File length in  8bit words (bytes)
NRecIndex= (NRecIndex- 100)/8      ! numero di Oggetti

if (ALLOCATED(Indice) ) DEALLOCATE(Indice)
ALLOCATE(Indice(NRecIndex, 2) )
do i= 1, NRecIndex
  iByte= 101+ (i-1)*8
  k1= Stream_Read_Integer(SHXunit, nbytes=4, posizione= iByte)
  k2= Stream_Read_Integer(SHXunit, nbytes=4)
  k1= BigToLittle(k1)* 2
  k2= BigToLittle(k2)
  Indice(i,:)= (/k1, k2/)
end do
call Stream_Close(SHXunit)
end subroutine Head


subroutine EstremiShape(Xmin, Xmax, Ymin, Ymax)
REAL(KIND=8), INTENT(OUT):: Xmin, Xmax, Ymin, Ymax
Xmin = Stream_Read_Real(ShpUnit, nbytes= 8, posizione= 37)
Ymin = Stream_Read_Real(ShpUnit, nbytes= 8, posizione= 45)
Xmax = Stream_Read_Real(ShpUnit, nbytes= 8, posizione= 53)
Ymax = Stream_Read_Real(ShpUnit, nbytes= 8, posizione= 61)
end subroutine

subroutine PosizionaSuOggetto(i)
INTEGER, INTENT(IN)::i
if (i<1 .or. i>NRecIndex) then
    NextRecordFileAddr = FileLastByte
  else
    NextRecordFileAddr = Indice(i, 1)+ 1
endif
end subroutine

subroutine LeggiRecordHeader( RecordNumber, ContentLength, ContentAddr, EOF)
INTEGER, INTENT(OUT), optional:: RecordNumber  ! record number
INTEGER, INTENT(OUT), optional:: ContentLength ! content length (in 16-bit words)
INTEGER, INTENT(OUT), optional:: ContentAddr   ! indirizzo del contenuto
LOGICAL, INTENT(OUT), optional:: EOF
INTEGER(KIND=4):: RecNum  ! record number
INTEGER(KIND=4):: ConLen ! content length (in 16-bit words)
INTEGER(KIND=4):: ConAddr   ! indirizzo del contenuto

if (NextRecordFileAddr < FileLastByte) then
   RecNum =  Stream_Read_Integer(SHPunit, nbytes=4, posizione= NextRecordFileAddr)
   RecNum = BigToLittle(RecNum)

   ConLen =  Stream_Read_Integer(SHPunit, nbytes=4)
   ConLen = BigToLittle(ConLen)

   ConAddr= Stream_Pos(ShpUnit)
!   INQUIRE(ShpUnit, NEXTREC= ConAddr)
   NextRecordFileAddr = NextRecordFileAddr+ (4+ ConLen )* 2
   if (PRESENT(RecordNumber )) RecordNumber   = RecNum
   if (PRESENT(ContentLength)) ContentLength  = ConLen
   if (PRESENT(ContentAddr  )) ContentAddr    = ConAddr
   if (PRESENT(EOF)) EOF= .false.
  else
   if (PRESENT(EOF)) EOF= .true.
endif
end subroutine


subroutine ReadPoly(P, posi, esito)
TYPE(PolyLine), INTENT(OUT):: P
INTEGER, INTENT(IN), OPTIONAL:: posi
integer, INTENT(OUT), OPTIONAL:: esito   ! 0 Ok      1 non trovato/EOF     2 shape non corrisponde

INTEGER(KIND=4):: NumParts, NumPoints, i
REAL(KIND=8):: x
!REAL(KIND=8), DIMENSION(2):: XY
INTEGER:: j, ipos, ilen
LOGICAL:: fine

if (PRESENT(posi) ) then
  call PosizionaSuOggetto(posi) ! posizionati sull'oggetto posi
 else
  ! oppure leggi l'oggetto sequenziale
endif
call LeggiRecordHeader(EOF= fine)

if (fine) then
  if (PRESENT(esito) ) esito= 1
else
  i = Stream_Read_Integer(ShpUnit, nbytes= 4)
  if (i==3) then
    do j=1, 4
      p%Box(j) = Stream_Read_Real(ShpUnit, nbytes= 8)
    enddo

    NumParts = Stream_Read_Integer(ShpUnit, nbytes= 4)
    NumPoints= Stream_Read_Integer(ShpUnit, nbytes= 4)

    if (ASSOCIATED(P%Points   ) ) DEALLOCATE(  P%Points    )
    if (ASSOCIATED(P%PartIndex) ) DEALLOCATE(  P%PartIndex )
    ALLOCATE(P%Points(NumPoints), P%PartIndex(NumParts) )

    do j=1, NumParts
      P%PartIndex(j) = Stream_Read_Integer(ShpUnit, nbytes= 4)
    enddo

    do j=1, NumPoints
      P%points (j)%x = Stream_Read_Real(ShpUnit, nbytes= 8)
      P%points (j)%y = Stream_Read_Real(ShpUnit, nbytes= 8)
    enddo

    if (PRESENT(esito) ) esito= 0
  else
    if (PRESENT(esito) ) esito= 2
  endif
endif
end subroutine


subroutine ReadPoint(p, posi, esito)
TYPE(point), INTENT(OUT):: p
INTEGER, INTENT(IN), OPTIONAL:: posi
integer, INTENT(OUT), OPTIONAL:: esito   ! 0 Ok      1 non trovato/EOF     2 shape non corrisponde
INTEGER(KIND=4):: i
INTEGER:: ipos, ilen
LOGICAL:: fine

if (PRESENT(posi) ) then
  call PosizionaSuOggetto(posi) ! posizionati sull'oggetto posi
 else
  ! oppure leggi l'oggetto sequenziale
endif
call LeggiRecordHeader(EOF= fine)

if (fine) then
  if (PRESENT(esito) ) esito= 1
else
  i= Stream_Read_Integer(ShpUnit, nbytes= 4)
  if (i==1) then
    p%x = Stream_Read_Real(ShpUnit, nbytes= 8)
    p%y = Stream_Read_Real(ShpUnit, nbytes= 8) 
    if (PRESENT(esito) ) esito= 0
  else
    if (PRESENT(esito) ) esito= 2
  endif
endif
end subroutine



function BigToLittle(i)
INTEGER(KIND=4), INTENT(IN):: i
INTEGER(KIND=4):: BigToLittle
INTEGER(KIND=4):: j
call MVBITS(i, 0, 8, j, 24)
call MVBITS(i, 8, 8, j, 16)
call MVBITS(i,16, 8, j,  8)
call MVBITS(i,24, 8, j,  0)
BigToLittle= j
end function



subroutine IniziaScrittura(nome, ITipoShape)
CHARACTER(LEN=*), INTENT(IN):: nome
INTEGER, INTENT(IN)::   ITipoShape  !1 3
INTEGER(KIND=4):: i, k1, k2
REAL(KIND=8):: x
NextRecordFileAddr= 101
numRecord = 0
Xmin1 = HUGE(xmin1)
Ymin1 = HUGE(xmin1)
Xmax1 =-HUGE(xmin1)
Ymax1 =-HUGE(xmin1)

SHPunit = UnitaLibera()
if (.not. Stream_Open(SHPunit, TRIM(nome)//".SHP", scrittura=.TRUE.) ) stop

SHXunit = UnitaLibera()
if (.not. Stream_Open(SHXunit, TRIM(nome)//".SHX", scrittura=.true.) ) stop

call Intesta(SHPunit)
call Intesta(SHXunit)
contains
 subroutine intesta(u)
 INTEGER, INTENT(IN):: u
      call Stream_Write(u,  BigToLittle( 9994), nbytes= 4 , posizione=  1)
      call Stream_Write(u,    0, nbytes= 4 , posizione=  5)   ! unused
      call Stream_Write(u,    0, nbytes= 4 , posizione=  9)   ! unused
      call Stream_Write(u,    0, nbytes= 4 , posizione= 13)   ! unused
      call Stream_Write(u,    0, nbytes= 4 , posizione= 17)   ! unused
      call Stream_Write(u,    0, nbytes= 4 , posizione= 21)   ! unused
      call Stream_Write(u,    0, nbytes= 4 , posizione= 25)   ! file lenght
      call Stream_Write(u, 1000, nbytes= 4 , posizione= 29)   ! version
      call Stream_Write(u, ITipoShape, nbytes= 4 , posizione= 33)  ! Tipo shape
      call Stream_Write(u, 0.0, nbytes= 8 , posizione= 37)  ! Xmin
      call Stream_Write(u, 0.0, nbytes= 8 , posizione= 45)  ! Ymin
      call Stream_Write(u, 0.0, nbytes= 8 , posizione= 53)  ! Xmax
      call Stream_Write(u, 0.0, nbytes= 8 , posizione= 61)  ! Ymax
      call Stream_Write(u, 0.0, nbytes= 8 , posizione= 69)  ! Zmin
      call Stream_Write(u, 0.0, nbytes= 8 , posizione= 77)  ! Zmax
      call Stream_Write(u, 0.0, nbytes= 8 , posizione= 85)  ! Mmin
      call Stream_Write(u, 0.0, nbytes= 8 , posizione= 93)  ! Mmax 
 end subroutine
end subroutine IniziaScrittura



subroutine ChiudiScrittura()
INTEGER:: inext, SHPlen, SHXlen
inext= Stream_Pos(SHPunit)
SHPlen= BigToLittle(INT( inext/2 , 4) )
SHXlen= BigToLittle(INT( 50+ numRecord* 4 , 4) )

call Stream_Write(SHPunit, SHPlen, nbytes= 4 , posizione= 25)  ! file lenght
call Stream_Write(SHPunit, Xmin1 , nbytes= 8 , posizione= 37)  ! Xmin
call Stream_Write(SHPunit, Ymin1 , nbytes= 8 , posizione= 45)  ! Ymin
call Stream_Write(SHPunit, Xmax1 , nbytes= 8 , posizione= 53)  ! Xmax
call Stream_Write(SHPunit, Ymax1 , nbytes= 8 , posizione= 61)  ! Ymax
call Stream_Close(SHPunit)

call Stream_Write(SHXunit, SHXlen, nbytes= 4 , posizione= 25)  ! file lenght
call Stream_Write(SHXunit, Xmin1 , nbytes= 8 , posizione= 37)  ! Xmin
call Stream_Write(SHXunit, Ymin1 , nbytes= 8 , posizione= 45)  ! Ymin
call Stream_Write(SHXunit, Xmax1 , nbytes= 8 , posizione= 53)  ! Xmax
call Stream_Write(SHXunit, Ymax1 , nbytes= 8 , posizione= 61)  ! Ymax
call Stream_Close(SHXunit)
end subroutine ChiudiScrittura


subroutine ScriviPolyline(P)
TYPE(PolyLine), INTENT(INOUT):: P
INTEGER:: ilen, inext
INTEGER:: numParts, nPunti
INTEGER:: i

if (.NOT. ASSOCIATED(P%points) ) return
if (.NOT. ASSOCIATED(P%PartIndex) ) then
  ALLOCATE(P%PartIndex(1) )
  P%PartIndex= 0
endif
numParts = UBOUND(P%PartIndex, 1)
nPunti   = UBOUND(P%points,    1)

P%box(1)= MINVAL(p%points%X)
P%box(2)= MINVAL(p%points%Y)
P%box(3)= MAXVAL(p%points%X)
P%box(4)= MAXVAL(p%points%Y)

Xmin1 =  MIN( Xmin1, P%box(1) )
Ymin1 =  MIN( Ymin1, P%box(2) )
Xmax1 =  MAX( Xmax1, P%box(3) )
Ymax1 =  MAX( Ymax1, P%box(4) )
numRecord= numRecord+1
inext= Stream_Pos(SHPunit)

ilen= 22 + 2* numParts+ 8* nPunti    ! lunghezza record in 16bit words

! ----------- record Header -------------------------------
call Stream_Write(SHPunit, BigToLittle(numRecord) , nbytes= 4 )
call Stream_Write(SHPunit, BigToLittle(ilen)      , nbytes= 4 )
! ---------------------------------------------------------

call Stream_Write(SHPunit,         3, nbytes= 4 )  ! ShapeType = 3 PolyLine                4
call Stream_Write(SHPunit, P%box(1) , nbytes= 8 )  !                                      32
call Stream_Write(SHPunit, P%box(2) , nbytes= 8 )  !
call Stream_Write(SHPunit, P%box(3) , nbytes= 8 )  !
call Stream_Write(SHPunit, P%box(4) , nbytes= 8 )  !
call Stream_Write(SHPunit, numParts , nbytes= 4 )  ! numParts                              4
call Stream_Write(SHPunit, nPunti   , nbytes= 4 )  ! numPoints                             4

do i= 1, numParts
  call Stream_Write(SHPunit, P%PartIndex(i), nbytes= 4 )  ! index to first points in parts    4* numParts
enddo

do i= 1, nPunti
  call Stream_Write(SHPunit, P%Points(i)%x , nbytes= 8 ) !                                 16* nPunti
  call Stream_Write(SHPunit, P%Points(i)%y , nbytes= 8 )
enddo

call Stream_Write(SHXunit, bigtolittle(inext/ 2) , nbytes= 4 )
call Stream_Write(SHXunit, bigtolittle(ilen)     , nbytes= 4 )
end subroutine


subroutine ScriviPoint(P)
TYPE(Point), INTENT(IN):: P
INTEGER:: ilen, inext

Xmin1 =  MIN( Xmin1, P%X )
Ymin1 =  MIN( Ymin1, P%Y )
Xmax1 =  MAX( Xmax1, P%X )
Ymax1 =  MAX( Ymax1, P%Y )
numRecord= numRecord+1
inext= Stream_Pos(SHPunit)

ilen= 10                             ! lunghezza record in 16bit words
! ----------- record Header -------------------------------
call Stream_Write(SHPunit, BigToLittle(numRecord) , nbytes= 4 )
call Stream_Write(SHPunit, BigToLittle(ilen)      , nbytes= 4 )
! ---------------------------------------------------------

call Stream_Write(SHPunit,         1, nbytes= 4 )  ! ShapeType = 1 Point                   4 bytes

call Stream_Write(SHPunit, P%x , nbytes= 8 ) ! X                                     8 bytes
call Stream_Write(SHPunit, P%y , nbytes= 8 ) ! Y                                     8 bytes

call Stream_Write(SHXunit, bigtolittle(inext/ 2) , nbytes= 4 )
call Stream_Write(SHXunit, bigtolittle(ilen)     , nbytes= 4 )
end subroutine

function UnitaLibera()
INTEGER:: UnitaLibera
INTEGER:: i
LOGICAL:: aperto
do i= 500, 1000
  INQUIRE(UNIT=i, OPENED= aperto)
  if (.not. aperto) exit
enddo
unitaLibera= i
end function


subroutine TipoShape(stri)
CHARACTER(LEN=*):: stri
select case (ShapeType)
    case ( 0)
      stri = "NullShape"
    case ( 1)
      stri = "Point"
    case ( 3)
      stri = "PolyLine"
    case ( 5)
      stri = "Polygon"
    case ( 8)
      stri = "MultiPoint"
    case (11)
      stri = "PointZ"
    case (13)
      stri = "PolyLineZ"
    case (15)
      stri = "PolygonZ"
    case (18)
      stri = "MultiPointZ"
    case (21)
      stri = "PointM"
    case (23)
      stri = "PolyLineM"
    case (25)
      stri = "PolygonM"
    case (28)
      stri = "MultiPointM"
    case (31)
      stri = "MultiPatch"
    case default
      stri = "Unknown"
end select
end subroutine
end module
