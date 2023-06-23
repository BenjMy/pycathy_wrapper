!     Last change:  MC    3 Sep 2008   11:34 am
module DBF
! gestione di files DBF sia in lettura che in scrittura
!
! specidfcihe dei files DBF  (from http://www.clicketyclick.dk/databases/xbase/format/index.html)
!
!00h /   0| Version number      *1|  ^
!         |-----------------------|  |
!01h /   1| Date of last update   |  |
!02h /   2|      YYMMDD        *21|  |
!03h /   3|                    *14|  |
!         |-----------------------|  |
!04h /   4| Number of records     | Record
!05h /   5| in data file          | header
!06h /   6| ( 32 bits )        *14|  |
!07h /   7|                       |  |
!         |-----------------------|  |
!08h /   8| Length of header   *14|  |
!09h /   9| structure ( 16 bits ) |  |
!         |-----------------------|  |
!0Ah /  10| Length of each record |  |
!0Bh /  11| ( 16 bits )     *2 *14|  |
!         |-----------------------|  |
!0Ch /  12| ( Reserved )        *3|  |
!0Dh /  13|                       |  |
!         |-----------------------|  |
!0Eh /  14| Incomplete transac.*12|  |
!         |-----------------------|  |
!0Fh /  15| Encryption flag    *13|  |
!         |-----------------------|  |
!10h /  16| Free record thread    |  |
!11h /  17| (reserved for LAN     |  |
!12h /  18|  only )               |  |
!13h /  19|                       |  |
!         |-----------------------|  |
!14h /  20| ( Reserved for        |  |            _        |=======================| ______
!         |   multi-user dBASE )  |  |           / 00h /  0| Field name in ASCII   |  ^
!         : ( dBASE III+ - )      :  |          /          : (terminated by 00h)   :  |
!         :                       :  |         |           |                       |  |
!1Bh /  27|                       |  |         |   0Ah / 10|                       |  |
!         |-----------------------|  |         |           |-----------------------| For
!1Ch /  28| MDX flag (dBASE IV)*14|  |         |   0Bh / 11| Field type (ASCII) *20| each
!         |-----------------------|  |         |           |-----------------------| field
!1Dh /  29| Language driver     *5|  |        /    0Ch / 12| Field data address    |  |
!         |-----------------------|  |       /             |                     *6|  |
!1Eh /  30| ( Reserved )          |  |      /              | (in memory !!!)       |  |
!1Fh /  31|                     *3|  |     /       0Fh / 15| (dBASE III+)          |  |
!         |=======================|__|____/                |-----------------------|  | <-
!20h /  32|                       |  |  ^          10h / 16| Field length       *22|  |   |
!         |- - - - - - - - - - - -|  |  |                  |-----------------------|  |   | *7
!         |                    *19|  |  |          11h / 17| Decimal count      *23|  |   |
!         |- - - - - - - - - - - -|  |  Field              |-----------------------|  | <-
!         |                       |  | Descriptor  12h / 18| ( Reserved for        |  |
!         :. . . . . . . . . . . .:  |  |array     13h / 19|   multi-user dBASE)*18|  |
!         :                       :  |  |                  |-----------------------|  |
!      n  |                       |__|__v_         14h / 20| Work area ID       *16|  |
!         |-----------------------|  |    \                |-----------------------|  |
!      n+1| Terminator (0Dh)      |  |     \       15h / 21| ( Reserved for        |  |
!         |=======================|  |      \      16h / 22|   multi-user dBASE )  |  |
!      m  | Database Container    |  |       \             |-----------------------|  |
!         :                    *15:  |        \    17h / 23| Flag for SET FIELDS   |  |
!         :                       :  |         |           |-----------------------|  |
!    / m+263                      |  |         |   18h / 24| ( Reserved )          |  |
!         |=======================|__v_ ___    |           :                       :  |
!         :                       :    ^       |           :                       :  |
!         :                       :    |       |           :                       :  |
!         :                       :    |       |   1Eh / 30|                       |  |
!         | Record structure      |    |       |           |-----------------------|  |
!         |                       |    |        \  1Fh / 31| Index field flag    *8|  |
!         |                       |    |         \_        |=======================| _v_____
!         |                       | Records
!         |-----------------------|    |
!         |                       |    |          _        |=======================| _______
!         |                       |    |         / 00h /  0| Record deleted flag *9|  ^
!         |                       |    |        /          |-----------------------|  |
!         |                       |    |       /           | Data               *10|  One
!         |                       |    |      /            : (ASCII)            *17: record
!         |                       |____|_____/             |                       |  |
!         :                       :    |                   |                       | _v_____
!         :                       :____|_____              |=======================|
!         :                       :    |     \             | Field deleted flag  *9|
!         |                       |    |      \            |-----------------------|
!         |                       |    |       \           |                       |
!         |                       |    |        \          |                       |
!         |                       |    |         \_        |-----------------------|
!         |                       |    |
!         |=======================|    |
!         |__End_of_File__________| ___v____  End of file ( 1Ah )  *11
!
!
USE ByteStreamer

implicit none
!private

PUBLIC:: DBFfile

PUBLIC:: dbase_lib_info

PUBLIC:: InitDbf, LoadRecord, CloseDbf, dbtest
PUBLIC:: DisplayRecord
PUBLIC:: DisplayRecordLinea, DisplayFieldNamesLinea

PUBLIC:: NumRecords, NumFields
PUBLIC:: FieldName, FieldType, FieldLen, FieldDec
PUBLIC:: DisplayStruct

PUBLIC:: FindRecord
PUBLIC:: GetField, iField, FieldPresent

PUBLIC:: CreaIndice, FindIndice


PUBLIC:: CreaDB ! crea un nuovo file di database vuoto
PUBLIC:: Add_Record
PUBLIC:: Scrivi_Record ! aggiorna il record attuale SCRIVENDOLO
PUBLIC:: Set_Field     ! non scrive, poi bisogna chiamare Scrivi_Record per rendere effettivo il cambiamento

interface Set_Field
    module procedure Set_Field_stri
    module procedure Set_Field_real
    module procedure Set_Field_integer
    module procedure Set_Field_logical
    module procedure Set_Field_data
end interface


interface GetField
    module procedure GetField1
    module procedure GetField2
end interface

CHARACTER(LEN=80):: dbase_lib_info= "Libreria di accesso DBF: versione 2.00 25/07/2008 - maurizio.cingi@eniaspa.it"


TYPE FieldDescr
!  sequence                             !                                     lunghezza (bytes)
  CHARACTER(LEN=11)::  FieldName       !                                     11
  CHARACTER(LEN= 1)::  FieldType       !                                      1
  INTEGER::            FieldAddress    !                                      4
  INTEGER::            FieldLength     !                                      1
  INTEGER::            Decimal         !                                      1
  INTEGER::            Reserv          !                                      2
  INTEGER::            WorkAreaID      !                                      1
  INTEGER::            SetField        !                                      1
  CHARACTER(LEN= 2)::  MultiUser       !                                      2
  character(LEN= 7)::  Reserv2         !                                      7
  INTEGER(KIND=1)  ::  MdxIndex        !                                      1
end type

TYPE Header
!  sequence              !                                     lunghezza (bytes)
  INTEGER:: byte00      !                                                   1
  INTEGER:: y, m, d     ! Last update                                       1
  INTEGER:: Nrecord=0   ! number of records in file                         4
  INTEGER:: HeaderSize  ! Header size                                       2
  INTEGER:: RecordSize  ! Record size (compreso il byte flag deleted *)     2
  INTEGER:: Reserved    !                                                   2
  INTEGER:: Transaction ! 01 begin transaction 00 end transaction           1
  INTEGER:: Encripted   ! 01 encripted  00 normal visible                   1
  character(len=12):: Multiuser               !                            12
  INTEGER:: IndexExists ! 01 production index exists 00 index upon demand   1
  INTEGER:: Language    ! language driver ID                                1
  INTEGER:: Reserved2   !                                                   2
end type

TYPE DBFfile
!  private
  CHARACTER(LEN=1024):: NomeFile= ""
  INTEGER            :: uni      = 0
  TYPE(Header)       :: Head
  INTEGER            :: NFields  = 0
  LOGICAL            :: Scrittura= .false.
  TYPE(FieldDescr), DIMENSION(:), POINTER:: F => NULL()
  INTEGER::                                 iRecord= 0          ! record corrente
  CHARACTER(LEN=1), DIMENSION(:), POINTER:: Record => NULL()    ! contenuto del record corrente
  LOGICAL::                                 del    = .false.    ! flag deleted del record corrente
  INTEGER, DIMENSION(:), POINTER         :: indice => NULL()
end type

CHARACTER(LEN=1), PARAMETER:: tab_char = CHAR(9)


contains

function GetField1(db,iField)
TYPE(DBFfile), INTENT(in):: db
INTEGER, INTENT(IN):: iField
CHARACTER(LEN=256):: GetField1
INTEGER:: inizio, j
GetField1=""
if (iField<1 .OR. iField>db%NFields) return

inizio= db%f(iField)%FieldAddress- 1
do j=1, db%f(iField)%FieldLength
  GetField1(j:j) = db%Record(inizio+j)
enddo
end function


function GetField2(db,nomeField)
TYPE(DBFfile), INTENT(in):: db
CHARACTER(LEN=*), INTENT(IN):: nomeField
CHARACTER(LEN=256):: GetField2
INTEGER:: ifi
ifi= iField(db, nomeField)
GetField2= GetField1(db, ifi)
end function

function FieldPresent(db,nomeField)
TYPE(DBFfile), INTENT(in):: db
CHARACTER(LEN=*), INTENT(IN):: nomeField
LOGICAL:: FieldPresent

FieldPresent= (iField(db, nomeField) /= 0)
end function


function NumRecords(db)
TYPE(DBFfile), INTENT(in):: db
INTEGER:: NumRecords
NumRecords= db%head%Nrecord
end function


function NumFields(db)
TYPE(DBFfile), INTENT(in):: db
INTEGER:: NumFields
NumFields= db%NFields
end function




subroutine InitDbf(nome, db, uni, ConsentiScrittura, esito)
CHARACTER(LEN=*), INTENT(IN):: nome
TYPE(DBFfile), INTENT(INOUT):: db
INTEGER, INTENT(IN) , OPTIONAL:: uni
logical, INTENT(IN) , OPTIONAL:: ConsentiScrittura
logical, INTENT(OUT), OPTIONAL:: esito
INTEGER(KIND=1):: i

CHARACTER(LEN=1):: aaa
INTEGER:: ibyte, nFields, Addr, ii

!call CloseDbf(db)

if (PRESENT(ConsentiScrittura) ) then
  db%scrittura= ConsentiScrittura
 else
  db%scrittura= .false.
end if

if (PRESENT(esito) ) esito= .false.

db%uni= 11
if (PRESENT(uni) ) then
   if (uni > 1) db%uni= uni
endif

db%NomeFile= nome
db%iRecord= 0
if (.not. Stream_Open(uni= db%uni, fileName=nome , scrittura= db%scrittura, stato= 1)) stop "apertura"

db%Head%byte00      = Stream_Read_Integer(db%uni,  1, posizione=1, unsigned= .TRUE.)
db%Head%y           = Stream_Read_Integer(db%uni,  1, unsigned= .TRUE.)
db%Head%m           = Stream_Read_Integer(db%uni,  1, unsigned= .TRUE.)
db%Head%d           = Stream_Read_Integer(db%uni,  1, unsigned= .TRUE.)
db%Head%Nrecord     = Stream_Read_Integer(db%uni,  4, unsigned= .TRUE.)
db%Head%HeaderSize  = Stream_Read_Integer(db%uni,  2, unsigned= .TRUE.)
db%Head%RecordSize  = Stream_Read_Integer(db%uni,  2, unsigned= .TRUE.)
db%Head%Reserved    = Stream_Read_Integer(db%uni,  2, unsigned= .TRUE.)
db%Head%Transaction = Stream_Read_Integer(db%uni,  1, unsigned= .TRUE.)
db%Head%Encripted   = Stream_Read_Integer(db%uni,  1, unsigned= .TRUE.)
db%Head%Multiuser   = Stream_Read_Char(db%uni, 12)
db%Head%IndexExists = Stream_Read_Integer(db%uni,  1, unsigned= .TRUE.)
db%Head%Language    = Stream_Read_Integer(db%uni,  1, unsigned= .TRUE.)
db%Head%Reserved2   = Stream_Read_Integer(db%uni,  2, unsigned= .TRUE.)

nFields=  (db%Head%HeaderSize - 1) / 32 - 1
db%NFields= nFields

ALLOCATE(db%F(nFields))
ALLOCATE(db%Record(db%Head%RecordSize-1))
db%Record=""
ibyte= 33
Addr= 1
do i=1, nFields
  db%f(i)%FieldName   = Stream_Read_Char(db%uni, 11, posizione= ibyte)
  db%f(i)%FieldType   = Stream_Read_Char(db%uni,  1)

  db%f(i)%FieldAddress= Stream_Read_Integer(db%uni,  4, unsigned= .TRUE.)
  db%f(i)%FieldLength = Stream_Read_Integer(db%uni,  1, unsigned= .TRUE.)
  db%f(i)%Decimal     = Stream_Read_Integer(db%uni,  1, unsigned= .TRUE.)
  db%f(i)%Reserv      = Stream_Read_Integer(db%uni,  2, unsigned= .TRUE.)
  db%f(i)%WorkAreaID  = Stream_Read_Integer(db%uni,  1, unsigned= .TRUE.)
  db%f(i)%MultiUser   = Stream_Read_Char(db%uni, 2  )
  db%f(i)%SetField    = Stream_Read_Integer(db%uni,  1, unsigned= .TRUE.)
  db%f(i)%Reserv2     = Stream_Read_Char(db%uni, 7)
  db%f(i)%MdxIndex    = Stream_Read_Integer(db%uni,  1, unsigned= .TRUE.)
  db%f(i)%FieldAddress= Addr
  Addr = Addr + db%f(i)%FieldLength
  ibyte= ibyte+ 32
  ii= INDEX(db%f(i)%FieldName,CHAR(0))
  if (ii>0) db%f(i)%FieldName(ii:)= ""
end do
ii= Stream_Pos(db%uni)
aaa= Stream_Read_Char(db%uni,1)
if ( aaa == CHAR(13) ) then ! controllo terminator
  if (PRESENT(esito) ) esito= .true.
 else
  WRITE(*,*) "terminatore", ICHAR(aaa), " in posizione ", ii
  READ(*,*)
!  call CloseDbf(db)
endif

end subroutine




subroutine CloseDbf(db)
TYPE(DBFfile), INTENT(INOUT):: db
if (db%uni/= 0) then
  if (db%scrittura)  call Stream_Write (db%uni, data_dbf() , nbytes=3, posizione= 2)
  call Stream_Close(db%uni)
  db%uni= 0
endif
db%Head%Nrecord = 0
db%nFields      = 0
db%iRecord      = 0
if (ASSOCIATED(db%F)     ) DEALLOCATE(db%F)
if (ASSOCIATED(db%Record)) DEALLOCATE(db%Record)
end subroutine




subroutine LoadRecord(db, irecord,esito)
TYPE(DBFfile), INTENT(INOUT):: db
INTEGER, OPTIONAL, INTENT(IN):: irecord
logical, OPTIONAL, INTENT(OUT)::esito
INTEGER:: irec
CHARACTER(LEN=1):: delFlag
INTEGER:: iad, Lrecord
if (PRESENT(esito)) esito= .false.

if (.not. ASSOCIATED(db%Record) ) return
if (PRESENT(irecord) ) then
 irec = irecord
else
 irec = db%iRecord + 1
endif

if (irec > 0 .AND. irec <= db%head%nRecord)  then
  iad= db%Head%HeaderSize+ (irec-1)* (db%Head%RecordSize)+ 1
  delFlag   = Stream_Read_Char(db%uni, 1, posizione= iad) ! Flag cancellazione
  Lrecord   = db%Head%RecordSize- 1                       ! RecordSize comprende Flag cancellazione
  db%Record = TRANSFER( Stream_Read_Char(db%uni, Lrecord), db%Record)
  db%del= (delFlag== "*")
  db%iRecord= irec
  if (PRESENT(esito)) esito= .true.
 else
  db%Record= ""
  db%del= .false.
  db%iRecord= 0
endif

end subroutine

subroutine DisplayRecord(db,stri)
TYPE(DBFfile), INTENT(IN):: db
CHARACTER(LEN=*), INTENT(OUT):: stri
CHARACTER(LEN=2048):: aa
INTEGER:: i
stri=""
if (db%iRecord<1 .OR. db%iRecord>db%head%nRecord) then
  stri= "-- nessun record caricato --"
 else
  WRITE(aa,"(A,I5.5)") "record #", db%iRecord
  if (db%del ) stri = TRIM(stri)// " **D E L E T E D**"
  call AttaccaLinea(stri, aa)
  do i=1, db%nFields
    WRITE(aa,"(A)") FieldName(db,i)//"-"//TRIM(GetField(db, i))//"-"
    call AttaccaLinea(stri, aa)
  end do
endif
end subroutine


subroutine DisplayRecordLinea(db, lin)
TYPE(DBFfile), INTENT(IN):: db
INTEGER:: i
CHARACTER(LEN=*), INTENT(OUT):: lin
lin=""
do i=1, db%nFields
  lin= TRIM(lin)//TRIM(adjustl(GetField(db, i)) ) // tab_char
end do
end subroutine

subroutine DisplayFieldNamesLinea(db, lin)
TYPE(DBFfile), INTENT(IN):: db
INTEGER:: i
CHARACTER(LEN=*), INTENT(OUT):: lin
lin=""
do i=1, db%nFields
  lin= TRIM(lin)//TRIM(adjustl(FieldName(db, i)) )// tab_char
end do
end subroutine


function FieldName(db, iField)
TYPE(DBFfile), INTENT(IN):: db
INTEGER, INTENT(IN):: iField
CHARACTER(LEN=11):: FieldName
if (iField>0 .AND. iField <= db%NFields) then
  FieldName= db%f(iField)%FieldName
 else
  FieldName= ""
endif
end function

function FieldLen(db, iField)
TYPE(DBFfile), INTENT(IN):: db
INTEGER, INTENT(IN):: iField
INTEGER:: FieldLen
FieldLen= db%f(iField)%FieldLength
end function

function FieldDec(db, iField)
TYPE(DBFfile), INTENT(IN):: db
INTEGER, INTENT(IN):: iField
INTEGER:: FieldDec
FieldDec= db%f(iField)%Decimal
end function




function FieldType(db, iField)
TYPE(DBFfile), INTENT(IN):: db
INTEGER, INTENT(IN):: iField
CHARACTER(LEN=20):: FieldType
if (iField>0 .AND. iField <= db%NFields) then
  select case (db%f(iField)%FieldType)
      case ("C")
         FieldType="ASCII"
      case ("D")
         FieldType="DATE"
      case ("F")
         FieldType="NUMERIC FLOAT"
      case ("N")
         FieldType="NUMERIC FIX"
      case ("L")
         FieldType="LOGICAL"
      case ("M")
         FieldType="MEMO"
      case ("V")
         FieldType="VARIABLE BIN"
      case ("P")
         FieldType="PICTURE"
      case ("B")
         FieldType="BINARY DATA"
      case ("G")
         FieldType="OLE object"
      case ("2")
         FieldType="Short INT"
      case ("4")
         FieldType="Long INT"
      case ("8")
         FieldType="double IEEE"
      case default
         FieldType="unknow"
  end select
 else
  FieldType= ""
endif
end function

function iField(db,nome)
INTEGER:: iField
TYPE(DBFfile),    INTENT(IN):: db
CHARACTER(LEN=*), INTENT(IN):: nome
INTEGER:: i
do i=1, db%nFields
  if (db%f(i)%FieldName==nome) then
     iField= i
     return
  end if
end do
iField = 0
end function

subroutine DisplayStruct(database,stri)
TYPE(DBFfile), INTENT(IN):: database
CHARACTER(LEN=*), INTENT(OUT):: stri
CHARACTER(LEN=1024):: LL
INTEGER:: i, itemp
if (.NOT. ASSOCIATED(database%F) ) then
  stri="- Database non valido -"
 else
  stri=""
  WRITE(LL,"(A)") "---- begin Display struttura DataBase --------------------------------- "
  call AttaccaLinea(stri, LL)
  WRITE(LL,"(A,1X,Z2.2)") "Type of file ", database%Head%byte00
  call AttaccaLinea(stri, LL)
  WRITE(LL,"(A,1X,2(I2.2,'/'),I4.4)") "Last updated ", database%Head%d, database%Head%m, database%Head%y+1900
  call AttaccaLinea(stri, LL)
  WRITE(LL,"(A, I10, A)") "il file contiene ",   NumRecords(database)," records "
  call AttaccaLinea(stri, LL)
  call AttaccaLinea(stri, "---------------------------------------------")
  WRITE(LL,"(A, I5, A,I8,A)") "il record contiene ", NumFields(database)," campi,", database%Head%RecordSize, " bytes"
  call AttaccaLinea(stri, LL)
  call AttaccaLinea(stri, "---------------------------------------------")
  itemp= 0
  do i=1, NumFields(database)
    WRITE(LL,"(1X,I2.2,2(1X,A),2I5,I10)") i, FieldName(database, i), FieldType(database, i), FieldLen(database, i),&
                                     & FieldDec(database, i), database%f(i)%FieldAddress
    call AttaccaLinea(stri, LL)
    itemp= itemp+ FieldLen(database, i)
  end do
  call AttaccaLinea(stri, "----end Display struttura DataBase --------------------------------- ")
endif
end subroutine

subroutine AttaccaLinea(stri, linea)
CHARACTER(LEN=*), INTENT(INOUT):: stri
CHARACTER(LEN=*), INTENT(IN)   :: linea
if (stri/="") stri= TRIM(stri)//CHAR(13)//CHAR(10)
stri= TRIM(stri)//linea
end subroutine


subroutine DBTest(nome)
TYPE(DBFfile):: database
INTEGER:: i
CHARACTER(LEN=1024):: linea
CHARACTER(LEN=20480):: stri
CHARACTER(LEN=*), INTENT(IN):: nome

call InitDbf(TRIM(nome)//".dbf", database,11)
OPEN(12,FILE=TRIM(nome)//".txt.",FORM="formatted")
call DisplayStruct(database, stri)
WRITE(12,"(A)") TRIM(stri)
call DisplayFieldNamesLinea(database, linea)
WRITE(12,"(A)") TRIM(linea)
do i=1, NumRecords(database)
  call LoadRecord(database,i)
  call DisplayRecordLinea(database, linea)
  WRITE(12,"(A)") TRIM(linea)
end do
call CloseDbf(database)
end subroutine

function FindRecord(db, nomeField, valore, next)
INTEGER:: FindRecord
TYPE(DBFfile):: db
CHARACTER(LEN=*),INTENT(IN):: nomeField, valore
LOGICAL, INTENT(IN), OPTIONAL:: next
INTEGER:: i, j

CHARACTER(LEN=1):: tipo
CHARACTER(LEN= 20):: for
CHARACTER(LEN=256):: val1
REAL(KIND=8):: xval
INTEGER:: inizio

inizio= 1
if (PRESENT(next)) then
  if (next) then
    inizio= db%iRecord + 1
  endif
endif

FindRecord= 0
j= iField(db, nomeField)
if (j==0) return
tipo= db%f(j)%FieldType

select case (tipo)
  case ("N")
    READ(valore, "(F20.0)") xval
    if (db%f(j)%Decimal > 0) then
      WRITE(for,"('(F',I2.2,'.',I2.2,')')") db%f(j)%FieldLength, db%f(j)%Decimal
      WRITE(val1,for) xval
     else
      WRITE(for,"('(I',I2.2,')')") db%f(j)%FieldLength
      WRITE(val1,for) INT(xval)
    endif
  case ("F")
    READ(valore, "(F20.0)") xval
    WRITE(for,"('(E',I2.2,'.',I2.2,')')") db%f(j)%FieldLength, db%f(j)%Decimal
    WRITE(val1,for) xval
  case default
    val1= valore
end select

do i= inizio, NumRecords(db)
   call LoadRecord(db, i)
   if ((GetField(db, j) == val1)) then
     FindRecord= i
     return
   endif
end do
end function

subroutine CreaDB(nome, nfields, NomiCampi, TipoCampi, lungCampi, decimCampi, unita, esito)
CHARACTER(LEN=*), INTENT(IN):: nome
INTEGER, INTENT(IN):: nfields
CHARACTER(LEN=*), DIMENSION(nfields), INTENT(IN):: NomiCampi
CHARACTER(LEN=1), DIMENSION(nfields), INTENT(IN):: TipoCampi
integer         , DIMENSION(nfields), INTENT(IN):: lungCampi
integer         , DIMENSION(nfields), INTENT(IN):: decimCampi
CHARACTER(LEN=11):: NomeField
INTEGER, INTENT(IN), OPTIONAL:: unita
LOGICAL, INTENT(OUT), OPTIONAL:: esito
!TipoCampi -> C ascii D date yyyymmdd F numeric N numeric (int?) L logical "YyNnTtFf ?"

INTEGER:: uni, Hsize, i

uni= 11
if (PRESENT(esito) ) esito= .false.
if (PRESENT(unita) ) then
   if (unita > 1) uni= unita
endif

if (.not. Stream_Open(uni= uni, fileName=nome , scrittura= .true. )) return


! dbf version number     plain dbf                    bytes  0 :  0
call Stream_Write (uni, 3 , nbytes=1, posizione= 1)

! date of last update                                 bytes  1 :  3      y-1900, m, d
call Stream_Write (uni, data_dbf(), nbytes=3)

! number of records in file       0 records nel file  bytes  4 :  7
call Stream_Write (uni, 0, nbytes=4)

! Header size in bytes                         bytes  8 :  9
HSize= 32* (nfields+1 ) + 1
call Stream_Write (uni, Hsize, nbytes=2, unsigned=.true.)

! Record size   (compreso deleted flag)        bytes 10 : 11
call Stream_Write (uni, SUM(lungCampi,1)+ 1, nbytes=2, unsigned=.true.)

! Reserved                                     bytes 12 : 13
call Stream_Write (uni, 0, nbytes=2)

! Incomplete transaction                       bytes 14 : 14
call Stream_Write (uni, 0, nbytes=1)

! Encrypted        normal visible              bytes 15 : 15
call Stream_Write (uni, 0, nbytes=1)

! Free record thread (reserved for LAN only)   bytes 16 : 19
call Stream_Write (uni, 0, nbytes=4)

! Reserved for multiuser dBASE)                bytes 20 : 27
call Stream_Write (uni, 0, nbytes=4)
call Stream_Write (uni, 0, nbytes=4)

! MDX flag (index)   index upon demand         bytes 28 : 28
call Stream_Write (uni, 0, nbytes=1)

! Language driver    windows ANSI              bytes 29 : 29
call Stream_Write (uni, 3, nbytes=1)

! Reserved                                     bytes 30 : 31
call Stream_Write (uni, 0, nbytes=2)



! Field descriptors -------------------------------------
do i=1, nfields
  nomefield= TRIM(NomiCampi(i) )//CHAR(0)
  ! field name
  call Stream_Write (uni, nomefield, nbytes=11)

  !field type
  call Stream_Write (uni, TipoCampi(i) , nbytes=1)

  ! field address in memory
  call Stream_Write (uni, 0, nbytes=4)

  ! field length
  call Stream_Write (uni, lungCampi(i), nbytes=1, unsigned=.true.)

  ! decimal count
  call Stream_Write (uni, decimCampi(i), nbytes=1, unsigned=.true.)

  ! reserved
  call Stream_Write (uni, 0, nbytes=2)

  ! work area ID
  call Stream_Write (uni, 0, nbytes=1)

  ! reserved for multi-user dB
  call Stream_Write (uni, 0, nbytes=2)

  ! reserved for SET FIELDS
  call Stream_Write (uni, 0, nbytes=1)

  ! reserved
  call Stream_Write (uni, 0, nbytes=4)
  call Stream_Write (uni, 0, nbytes=2)
  call Stream_Write (uni, 0, nbytes=1)

  ! index field flag
  call Stream_Write (uni, 0, nbytes=1)

end do
!--------------------------------------------------------

! Header Record Terminator
  call Stream_Write (uni, 13, nbytes=1)

! End_Of_File
  call Stream_Write (uni, 26, nbytes=1)
call Stream_Close(uni)
if (PRESENT(esito) ) esito= .true.

end subroutine


subroutine Add_Record(db)
TYPE(DBFfile), INTENT(inout):: db
INTEGER:: ibyte
if (.not. db%scrittura) return
ibyte= db%Head%HeaderSize+ (db%Head%Nrecord * db%Head%RecordSize)+ 1

db%Record=""
call Stream_Write (db%uni, CHAR(32)   , nbytes= 1, posizione= ibyte)
call WriteDBrecord(db)
call Stream_Write (db%uni, 26, nbytes=1)  ! End_Of_File

db%Head%Nrecord= db%Head%Nrecord+ 1
call Stream_Write (db%uni, db%Head%Nrecord , nbytes= 4, posizione= 5)
db%iRecord= db%Head%Nrecord
end subroutine

subroutine WriteDBrecord(db)
TYPE(DBFfile), INTENT(in):: db
INTEGER:: i
do i=1, db%Head%RecordSize- 1    ! RecordSize comprende Flag cancellazione
  call Stream_Write (db%uni, db%Record(i), nbytes= 1)
end do
end subroutine WriteDBrecord


subroutine Scrivi_Record(db)
TYPE(DBFfile), INTENT(inout):: db
INTEGER:: ibyte
CHARACTER(LEN=1):: flag
if (.not. db%scrittura) return
ibyte= db%Head%HeaderSize+ (db%iRecord- 1) * db%Head%RecordSize+ 1
if ( db%del ) then
  flag="*"
 else
  flag=" "
endif
call Stream_Write (db%uni, flag   , nbytes= 1, posizione= ibyte)
call WriteDBrecord(db)
end subroutine

subroutine Set_Field_stri(db, nomeField, stri)
TYPE(DBFfile), INTENT(inout):: db
CHARACTER(LEN=*), INTENT(IN):: nomeField
CHARACTER(LEN=*), INTENT(IN):: stri
CHARACTER(LEN=256):: ss
INTEGER:: i, inizio, fine

i= iField(db, nomeField)
if (i== 0) return
ss= stri

inizio= db%f(i)%FieldAddress
fine  = db%f(i)%FieldAddress+ db%f(i)%FieldLength- 1
db%Record(inizio: fine)= TRANSFER(ss,db%record, db%f(i)%FieldLength)
end subroutine

subroutine Set_Field_real(db, nomeField, x)
TYPE(DBFfile), INTENT(inout):: db
CHARACTER(LEN=*), INTENT(IN):: nomeField
real, INTENT(IN):: x
CHARACTER(LEN= 20):: stri
CHARACTER(LEN= 10):: formato
INTEGER:: i, inizio, fine
i= iField(db, nomeField)
if (i== 0) return
inizio= db%f(i)%FieldAddress
fine  = db%f(i)%FieldAddress+ db%f(i)%FieldLength- 1
if (db%f(i)%FieldType=="F") then
  WRITE(formato, "('(F',I3.3,'.',I3.3,')')")  db%f(i)%FieldLength, db%f(i)%Decimal
  WRITE(stri, formato) x
  if (stri(1:1)=="*") then
    WRITE(formato, "('(E',I3.3,'.',I3.3,')')")  db%f(i)%FieldLength, db%f(i)%FieldLength- 5
    WRITE(stri, formato) x
  endif
  db%Record(inizio: fine)= TRANSFER(stri,db%record, db%f(i)%FieldLength)
endif
if (db%f(i)%FieldType=="N") then
  WRITE(formato, "('(I',I3.3,')')")  db%f(i)%FieldLength
  WRITE(stri, formato) NINT(x)
  db%Record(inizio: fine)= TRANSFER(stri,db%record, db%f(i)%FieldLength)
endif
end subroutine

subroutine Set_Field_reald(db, nomeField, x)
TYPE(DBFfile), INTENT(inout):: db
CHARACTER(LEN=*), INTENT(IN):: nomeField
real(kind=8), INTENT(IN):: x
CHARACTER(LEN= 20):: stri
CHARACTER(LEN= 10):: formato
INTEGER:: i, inizio, fine
i= iField(db, nomeField)
if (i== 0) return
inizio= db%f(i)%FieldAddress
fine  = db%f(i)%FieldAddress+ db%f(i)%FieldLength- 1
if (db%f(i)%FieldType=="F") then
  WRITE(formato, "('(F',I3.3,'.',I3.3,')')")  db%f(i)%FieldLength, db%f(i)%Decimal
  WRITE(stri, formato) x
  if (stri(1:1)=="*") then
    WRITE(formato, "('(E',I3.3,'.',I3.3,')')")  db%f(i)%FieldLength, db%f(i)%FieldLength- 5
    WRITE(stri, formato) x
  endif
  db%Record(inizio: fine)= TRANSFER(stri,db%record, db%f(i)%FieldLength)
endif
if (db%f(i)%FieldType=="N") then
  WRITE(formato, "('(I',I3.3,')')")  db%f(i)%FieldLength
  WRITE(stri, formato) NINT(x)
  db%Record(inizio: fine)= TRANSFER(stri,db%record, db%f(i)%FieldLength)
endif
end subroutine

subroutine Set_Field_integer(db, nomeField, x)
TYPE(DBFfile), INTENT(inout):: db
CHARACTER(LEN=*), INTENT(IN):: nomeField
integer, INTENT(IN):: x
CHARACTER(LEN= 20):: stri
CHARACTER(LEN= 10):: formato
INTEGER:: i, inizio, fine
i= iField(db, nomeField)
if (i== 0) return
inizio= db%f(i)%FieldAddress
fine  = db%f(i)%FieldAddress+ db%f(i)%FieldLength- 1
if (db%f(i)%FieldType=="F") then
  WRITE(formato, "('(F',I3.3,'.',I3.3,')')")  db%f(i)%FieldLength, db%f(i)%Decimal
  WRITE(stri, formato) REAL(x)
  db%Record(inizio: fine)= TRANSFER(stri,db%record, db%f(i)%FieldLength)
endif
if (db%f(i)%FieldType=="N") then
  WRITE(formato, "('(I',I3.3,')')")  db%f(i)%FieldLength
  WRITE(stri, formato) x
  db%Record(inizio: fine)= TRANSFER(stri,db%record, db%f(i)%FieldLength)
endif
end subroutine

subroutine Set_Field_logical(db, nomeField, x)
TYPE(DBFfile), INTENT(inout):: db
CHARACTER(LEN=*), INTENT(IN):: nomeField
logical, INTENT(IN):: x
CHARACTER(LEN= 20):: stri

INTEGER:: i, inizio, fine
i= iField(db, nomeField)
if (i== 0) return
inizio= db%f(i)%FieldAddress
fine  = db%f(i)%FieldAddress+ db%f(i)%FieldLength- 1
if (db%f(i)%FieldType/="L") return
if (x) then
  stri="T"
 else
  stri="F"
endif
db%Record(inizio: fine)= TRANSFER(stri,db%record, db%f(i)%FieldLength)
end subroutine


subroutine Set_Field_data(db, nomeField, giorno, mese, anno)
TYPE(DBFfile), INTENT(inout):: db
CHARACTER(LEN=*), INTENT(IN):: nomeField
integer, INTENT(IN)::  giorno, mese, anno
CHARACTER(LEN= 20):: stri
INTEGER:: i, inizio, fine
i= iField(db, nomeField)
if (i== 0) return
inizio= db%f(i)%FieldAddress
fine  = db%f(i)%FieldAddress+ db%f(i)%FieldLength- 1
if (db%f(i)%FieldType/="D") return
WRITE(stri,"(I4.4,I2.2,I2.2)") anno, mese, giorno
db%Record(inizio: fine)= TRANSFER(stri,db%record, db%f(i)%FieldLength)
end subroutine

function DATA_DBF()
CHARACTER(LEN=3) data_dbf
CHARACTER(LEN= 8) date
INTEGER:: yy, mm, dd

call DATE_AND_TIME(date)
READ(date, "(I4,I2,I2)") yy, mm, dd
yy= yy-1900
DATA_DBF(1:1)= TRANSFER(yy, DATA_DBF)
DATA_DBF(2:2)= TRANSFER(mm, DATA_DBF)
DATA_DBF(3:3)= TRANSFER(dd, DATA_DBF)
end function


subroutine CreaIndice(db, nomeField)
TYPE(DBFfile), INTENT(inout):: db
CHARACTER(LEN=*), INTENT(IN):: nomeField
INTEGER:: n, i, j
CHARACTER(LEN=40):: aa
if (ASSOCIATED(db%indice)  ) DEALLOCATE(db%indice)
if (db%scrittura) return
n= NumRecords(db)
if (n<1) return
j= iField(db, nomeField )
if (j==0) return
ALLOCATE(db%indice(n) )
do i=1, n
 call LoadRecord(db, i)
 aa= GetField(db, j)
 READ(aa,"(I40)") db%indice(i)
end do
end subroutine


function FindIndice(db, indi)
INTEGER                     :: FindIndice
TYPE(DBFfile), INTENT(inout):: db
INTEGER , INTENT(IN)        :: indi
INTEGER:: i
FindIndice= 0
if (.not. ASSOCIATED(db%indice)  ) return
do i=1, NumRecords(db)
  if (db%indice(i)== indi ) then
       FindIndice= i
       return
  endif
end do
end function

end module


module GestoreCampi
USE DBF
implicit none
private
PUBLIC:: ResetCampi, SetCampo, CreaFileDb
INTEGER,parameter :: maxFields = 128
CHARACTER(LEN=10), DIMENSION(maxFields):: NomiCampi
CHARACTER(LEN=1), DIMENSION(maxFields):: TipoCampi   ! C ascii D date yyyymmdd F numeric N numeric (int?) L logical "YyNnTtFf ?"
integer         , DIMENSION(maxFields):: lungCampi
integer         , DIMENSION(maxFields):: decimCampi

contains

 subroutine ResetCampi()
     NomiCampi = ""
     TipoCampi = ""
     lungCampi =  0
     decimCampi=  0
 end subroutine

 subroutine SetCampo(nome, tipo, lunghezza, decimali)
  INTEGER, INTENT(IN)::  lunghezza, decimali
  CHARACTER(LEN=*), INTENT(IN):: nome, tipo
  INTEGER:: i
  do i=1, maxFields
   if (NomiCampi(i)=="" .or.  NomiCampi(i)== nome) then
     NomiCampi (i)= nome
     TipoCampi (i)= tipo
     lungCampi (i)= lunghezza
     decimCampi(i)= decimali
     return
   end if
  end do
 end subroutine

 function ContaCampi()
 INTEGER:: ContaCampi
 INTEGER:: i
  do i=1, maxFields
   if (NomiCampi(i)=="" ) exit
  end do
  ContaCampi= i- 1
 end function

 function CreaFileDb(nome)
   LOGICAL::  CreaFileDb
   CHARACTER(LEN=*), INTENT(IN):: nome
   INTEGER:: nfields
   nfields= ContaCampi()
   if (nfields<1) return
   call CreaDB(TRIM(nome)//".dbf", nfields, NomiCampi, TipoCampi, lungCampi, decimCampi, esito= CreaFileDb)
 end function CreaFileDb

end module

