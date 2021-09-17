!     Last change:  MC    3 Sep 2008    8:02 am


module AB_Normal
! -------------------------------------------------------------------------------------------------!
! contiene parti non standard                                                                      !
! in particolare OPEN (... , ACCESS="transparent", ...)                                            !
!                WRITE(... , pos= posizione, ...)                                                  !
!                READ (... , pos= posizione, ...)                                                  !
! per Lahey LF95 occorre decommentare le linee marcate con ##LF95## e commentare  ##G95##          !
! per G95        occorre decommentare le linee marcate con ##G95##  e commentare  ##LF95##         !
! -------------------------------------------------------------------------------------------------!
implicit none
contains
subroutine ByteWrite(uni, a, ii, posizione)
  INTEGER,          INTENT(IN):: uni
  CHARACTER(LEN=1), INTENT(IN):: a
  INTEGER,OPTIONAL,  INTENT(OUT):: ii
  INTEGER,OPTIONAL, INTENT(IN):: posizione
  INTEGER:: i
  if (PRESENT(posizione) ) then
!    write(uni, pos= posizione, IOSTAT=i) a  !##G95##
    write(uni, rec= posizione, IOSTAT=i) a  !##LF95##
   else
    write(uni, IOSTAT=i) a
  endif
  if (PRESENT(ii) ) ii= i
end subroutine ByteWrite

subroutine ByteRead(uni, a, ii, posizione)
  INTEGER,          INTENT(IN):: uni
  CHARACTER(LEN=1), INTENT(OUT):: a
  INTEGER,OPTIONAL,  INTENT(OUT):: ii
  INTEGER,OPTIONAL, INTENT(IN):: posizione
  INTEGER:: i
  if (PRESENT(posizione) ) then
!     read(uni, pos= posizione, IOSTAT=i) a   !##G95##
     read(uni, rec= posizione, IOSTAT=i) a   !##LF95##
   else
     read(uni,  IOSTAT=i) a
  endif
  if (PRESENT(ii) ) ii= i
end subroutine  ByteRead

function Byte_Pos(uni)
  INTEGER:: Byte_Pos
  INTEGER, INTENT(IN):: uni
  INTEGER:: i
!  INQUIRE(uni, pos= i )                    !##G95##
  INQUIRE(uni, NEXTREC= i )               !##LF95##
  Byte_Pos = i
end function Byte_Pos

subroutine Byte_Open(uni, filename, action, status, esito)
  INTEGER,          INTENT(IN):: uni  !fortran unit#
  CHARACTER(LEN=*), INTENT(IN):: fileName
  CHARACTER(LEN=*), INTENT(IN):: action
  CHARACTER(LEN=*), INTENT(IN):: status
  LOGICAL,         intent(OUT):: esito
  INTEGER:: ii
!  OPEN(UNIT=uni, FILE= fileName, ACTION= action, STATUS=status, ACCESS="stream", IOSTAT=ii)         !##G95##
  OPEN(UNIT=uni, FILE= fileName, ACTION= action, STATUS=status, ACCESS="transparent", IOSTAT=ii)   !##LF95##  
  esito= (ii==0)
end subroutine Byte_Open

end module AB_Normal



module ByteStreamer
use AB_Normal
implicit none
private
!-------------------------------------------------------------------------------------------------------!
! PARTE PUBBLICA DEL MODULO ByteStreamer                                                                !
!----------------------------------------------------------------!--------------------------------------!
                                                                 !                                      !
PUBLIC:: Stream_Open                                             ! SUB apertura dello stream            !
PUBLIC:: Stream_Close                                            ! SUB chiusura dello stream            !
PUBLIC:: Stream_Pos                                              ! FUN posizione prox accesso           !
PUBLIC:: Stream_Write                                            ! SUB scrittura di dati                !
PUBLIC:: Stream_Read_Integer, Stream_Read_Real, Stream_Read_Char ! FUN lettura di dati                  !
                                                                 !                                      !
!----------------------------------------------------------------!--------------------------------------!
!
! function Stream_Open(uni, fileName, scrittura, stato)    APRE LO STREAM
!  LOGICAL:: Stream_Open                      .true.=> apertura riuscita   .false.=> fallita
!  INTEGER, INTENT(IN):: uni                  fortran unit dello stream
!  CHARACTER(LEN=*), INTENT(IN):: fileName    nome del file (compresa estensione)
!  logical, INTENT(IN), optional:: scrittura  .true.=> consenti scrittura   .false.=>nega scrittura
!  integer, INTENT(IN), optional:: stato      default "UNKNOWN"  1 "OLD" 2 "NEW" 3 "SCRATCH"
!
!
! subroutine Stream_Close(uni)                CHIUDE LO STREAM
!  INTEGER, INTENT(IN):: uni                  fortran unit dello stream
!
! function Stream_Pos(uni)
!    INTEGER            :: Stream_Pos      posizione prossimo accesso (primo byte=1)
!    INTEGER, INTENT(IN):: uni             fortran unit dello stream
!
!
! Stream_Write(uni, aa, nbytes, posizione, unsigned)
!   INTEGER, INTENT(IN):: uni                   fortran unit dello stream
!    ......., INTENT(IN)::  aa                  dati da scrivere (REAL, INTEGER o CHAR) SCALARE
!   INTEGER, INTENT(IN)::            nbytes     numero di bytes da scrivere
!   INTEGER, INTENT(IN), OPTIONAL::  posizione  posizione in cui leggere (primo byte=1) se assente sequenziale
!   logical, INTENT(IN), OPTIONAL::  unsigned   scrive come unsigned (SOLO per aa INTEGER)
!
!function Stream_Read_Integer(uni, nbytes, posizione, unsigned, endfile)
!  INTEGER(KIND=8):: Stream_Read_Integer
!  INTEGER, INTENT(IN):: uni                     fortran unit dello stream
!  INTEGER, INTENT(IN)::             nbytes      numero di bytes da leggere
!  INTEGER, INTENT(IN) , OPTIONAL::  posizione   posizione in cui leggere (primo byte=1) se assente sequenziale
!  logical, INTENT(IN) , OPTIONAL::  unsigned    scrive come unsigned (solo per nbytes= 1 o  2 )
!  LOGICAL, INTENT(OUT), OPTIONAL::  endfile     .true. segnala l'incontrata fine del file
!
!function Stream_Read_Real(uni, nbytes, posizione, unsigned, endfile)
!  REAL(KIND=8)       :: Stream_Read_Real
!  INTEGER, INTENT(IN):: uni                     fortran unit dello stream
!  INTEGER, INTENT(IN)::             nbytes      numero di bytes da leggere
!  INTEGER, INTENT(IN) , OPTIONAL::  posizione   posizione in cui leggere (primo byte=1) se assente sequenziale
!  logical, INTENT(IN) , OPTIONAL::  unsigned    scrive come unsigned (solo per nbytes= 1 o  2 )
!  LOGICAL, INTENT(OUT), OPTIONAL::  endfile     .true. segnala l'incontrata fine del file
!
!function Stream_Read_Char(uni, nbytes, posizione, unsigned, endfile)
!  character(LEN=nbytes)       :: Stream_Read_Char
!  INTEGER, INTENT(IN):: uni                     fortran unit dello stream
!  INTEGER, INTENT(IN)::             nbytes      numero di bytes da leggere
!  INTEGER, INTENT(IN) , OPTIONAL::  posizione   posizione in cui leggere (primo byte=1) se assente sequenziale
!  logical, INTENT(IN) , OPTIONAL::  unsigned    scrive come unsigned (solo per nbytes= 1 o  2 )
!  LOGICAL, INTENT(OUT), OPTIONAL::  endfile     .true. segnala l'incontrata fine del file
!


interface Stream_Write
 module procedure Integer_8_WriteStream
 module procedure Integer_4_WriteStream
 module procedure Integer_2_WriteStream
 module procedure Integer_1_WriteStream
 module procedure REAL_4_WriteStream
 module procedure REAL_8_WriteStream
 module procedure CHAR_WriteStream
end interface


contains

function Stream_Pos(uni)
INTEGER:: Stream_Pos
INTEGER, INTENT(IN):: uni  !fortran unit#
Stream_Pos = Byte_Pos(uni)
end function Stream_Pos


function Stream_Open(uni, fileName, scrittura, stato) RESULT( esito)
  LOGICAL:: esito
  INTEGER, INTENT(IN):: uni  !fortran unit#
  CHARACTER(LEN=*), INTENT(IN):: fileName
  logical, INTENT(IN), optional:: scrittura   ! scrittura consentita
  integer, INTENT(IN), optional:: stato    !OLD    UNKNOWN    (vedi standard F95)
  INTEGER:: ii
  LOGICAL:: flag
  CHARACTER(LEN=32):: ACTION, status
  if (PRESENT(scrittura) ) then
    flag= scrittura
   else
    flag= .false.
  endif
  if (flag) then
    action = "READWRITE"
   else
    action = "READ"
  end if

  if (PRESENT(stato) ) then
    ii= stato
   else
    ii= 0
  endif
  select case (ii)
      case (1)
        status = "OLD"
      case (2)
        status = "NEW"
      case (3)
        status = "SCRATCH"
      case default
        status = "UNKNOWN"
  end select
  call Byte_Open(uni, filename, action, status, esito)
end function  Stream_Open





subroutine Stream_Close(uni)
  INTEGER, INTENT(IN):: uni  !fortran unit#
  CLOSE(uni)
end subroutine Stream_Close



subroutine ScriviStream(uni, aa, nbytes, posizione,  bytesScritti)
  INTEGER, INTENT(IN)         :: uni           ! fortran unit#
  CHARACTER(LEN=*), INTENT(IN):: aa            ! buffer di dati
  INTEGER, INTENT(IN)         :: nbytes        ! numero di bytes da scrivere (max len(a) )
  INTEGER, INTENT(IN)         :: posizione     ! posizione in cui scrivere (primo byte=1)
  INTEGER, INTENT(OUT)        :: bytesScritti  ! restituisce il numero di bytes effettivamente scritti
  INTEGER:: i, ii, j
  bytesScritti = 0
  call ByteWrite(uni, aa(1:1), ii,  posizione)
  if (ii==0) then
    do i=2, MIN(nbytes,  LEN(aa) )
      call ByteWrite(uni, aa(i:i), ii )
      write(uni, IOSTAT=ii)
      if (ii/=0) exit
    enddo
    do j=i, nbytes
      call ByteWrite(uni, CHAR(0))
    end do	
    bytesScritti = j-1
  endif
end subroutine ScriviStream


subroutine CHAR_WriteStream(uni, aa, nbytes, posizione)
  INTEGER, INTENT(IN)::            uni           ! fortran unit#
  CHARACTER(LEN=*), INTENT(IN)::   aa            ! buffer di dati
  INTEGER, INTENT(IN)::            nbytes        ! numero di bytes da scrivere (max len(a) )
  INTEGER, INTENT(IN), OPTIONAL::  posizione     ! posizione in cui scrivere (primo byte=1) se assente sequenziale
  INTEGER:: bytesScritti, iPosizione

  if (PRESENT(posizione) ) then
    iPosizione= posizione
   else
    iPosizione= Byte_Pos(uni)
  endif
  call ScriviStream(uni, aa, nbytes, iPosizione,  bytesScritti)
  if (nbytes/= bytesScritti) then
    stop "non ho scritto tutti i char"
  endif
end subroutine CHAR_WriteStream


subroutine LeggiStream(uni, aa, nbytes, posizione,  bytesLetti)
  INTEGER, INTENT(IN)::            uni          ! fortran unit#
  CHARACTER(LEN=*), INTENT(OUT)::  aa           ! buffer di dati
  INTEGER, INTENT(IN)::            nbytes       ! numero di bytes da leggere (max len(a) )
  INTEGER, INTENT(IN)          ::  posizione    ! posizione in cui leggere (primo byte=1) se assente sequenziale
  INTEGER, INTENT(OUT)         ::  bytesLetti   ! restituisce il numero di bytes effettivamente letti
  INTEGER:: i, ii

  bytesLetti = 0
  aa= REPEAT(CHAR(0), LEN(aa))
  call ByteRead(uni, aa(1:1), ii, posizione)
  if (ii==0) then
    do i=2, MIN(nbytes,  LEN(aa) )
     call ByteRead(uni, aa(i:i), ii)
     if (ii/=0) exit
    enddo
    bytesLetti = i-1
  endif
end subroutine LeggiStream



function Stream_Read_Char(uni, nbytes, posizione, endfile)
  INTEGER, INTENT(IN)::            uni          ! fortran unit#
  INTEGER, INTENT(IN)::            nbytes       ! numero di bytes da scrivere
  INTEGER, INTENT(IN), OPTIONAL::  posizione    ! posizione in cui leggere (primo byte=1) se assente sequenziale
  LOGICAL, INTENT(OUT), OPTIONAL:: endfile      ! segnala l'incontrata fine del file
  character(LEN=nbytes)       :: Stream_Read_Char
  INTEGER:: iPosizione, bytesLetti

  if (PRESENT(posizione) ) then
    iPosizione= posizione
   else
    iPosizione= Byte_Pos(uni)
  endif

  call LeggiStream(uni, Stream_Read_Char, nbytes, iPosizione,  bytesLetti)
  if (PRESENT( endfile ) ) endfile = (bytesLetti/= nbytes)
end function Stream_Read_Char



subroutine Integer_8_WriteStream(uni, aa, nbytes, posizione, unsigned)
  INTEGER, INTENT(IN)::            uni          ! fortran unit#
  INTEGER(KIND=8), INTENT(IN)  ::  aa           ! buffer di dati
  INTEGER, INTENT(IN)::            nbytes       ! numero di bytes da scrivere
  INTEGER, INTENT(IN), OPTIONAL::  posizione    ! posizione in cui leggere (primo byte=1) se assente sequenziale
  logical, INTENT(IN), OPTIONAL::  unsigned     ! scrive come unsigned
  LOGICAL:: unsflag
  INTEGER:: iPosizione
  if (PRESENT(unsigned) ) then
    unsflag= unsigned
   else
    unsflag= .false.
  endif

  if (PRESENT(posizione) ) then
    iPosizione= posizione
   else
    iPosizione= Byte_Pos(uni)
  endif
  call WriteStreamInteger_INTERNAL(uni, INT(aa,8), nbytes, iPosizione, unsflag)
end subroutine Integer_8_WriteStream





subroutine Integer_4_WriteStream(uni, aa, nbytes, posizione, unsigned)
  INTEGER, INTENT(IN)::            uni          ! fortran unit#
  INTEGER(KIND=4), INTENT(IN)  ::  aa           ! buffer di dati
  INTEGER, INTENT(IN)::            nbytes       ! numero di bytes da scrivere
  INTEGER, INTENT(IN), OPTIONAL::  posizione    ! posizione in cui leggere (primo byte=1) se assente sequenziale
  logical, INTENT(IN), OPTIONAL::  unsigned     ! scrive come unsigned (solo per nbytes= 1 o  2 )

  LOGICAL:: unsflag
  INTEGER:: iPosizione
  if (PRESENT(unsigned) ) then
    unsflag= unsigned
   else
    unsflag= .false.
  endif

  if (PRESENT(posizione) ) then
    iPosizione= posizione
   else
    iPosizione= Byte_Pos(uni)
  endif
  call WriteStreamInteger_INTERNAL(uni, INT(aa,8), nbytes, iPosizione, unsflag)
end subroutine Integer_4_WriteStream

subroutine Integer_2_WriteStream(uni, aa, nbytes, posizione, unsigned)
  INTEGER, INTENT(IN)::            uni          ! fortran unit#
  INTEGER(KIND=2), INTENT(IN)  ::  aa           ! buffer di dati
  INTEGER, INTENT(IN)::            nbytes       ! numero di bytes da scrivere
  INTEGER, INTENT(IN), OPTIONAL::  posizione    ! posizione in cui leggere (primo byte=1) se assente sequenziale
  logical, INTENT(IN), OPTIONAL::  unsigned     ! scrive come unsigned (solo per nbytes= 1 o  2 )

  LOGICAL:: unsflag
  INTEGER:: iPosizione
  if (PRESENT(unsigned) ) then
    unsflag= unsigned
   else
    unsflag= .false.
  endif

  if (PRESENT(posizione) ) then
    iPosizione= posizione
   else
    iPosizione= Byte_Pos(uni)
  endif
  call WriteStreamInteger_INTERNAL(uni, INT(aa,8), nbytes, iPosizione, unsflag)
end subroutine Integer_2_WriteStream


subroutine Integer_1_WriteStream(uni, aa, nbytes, posizione, unsigned)
  INTEGER, INTENT(IN)::            uni          ! fortran unit#
  INTEGER(KIND=1), INTENT(IN)  ::  aa           ! buffer di dati
  INTEGER, INTENT(IN)::            nbytes       ! numero di bytes da scrivere
  INTEGER, INTENT(IN), OPTIONAL::  posizione    ! posizione in cui leggere (primo byte=1) se assente sequenziale
  logical, INTENT(IN), OPTIONAL::  unsigned     ! scrive come unsigned (solo per nbytes= 1 o  2 )

  LOGICAL:: unsflag
  INTEGER:: iPosizione
  if (PRESENT(unsigned) ) then
    unsflag= unsigned
   else
    unsflag= .false.
  endif

  if (PRESENT(posizione) ) then
    iPosizione= posizione
   else
    iPosizione= Byte_Pos(uni)
  endif
  call WriteStreamInteger_INTERNAL(uni, INT(aa,8), nbytes, iPosizione, unsflag)
end subroutine Integer_1_WriteStream





subroutine WriteStreamInteger_INTERNAL(uni, aa, nbytes, posizione, unsigned)
  INTEGER, INTENT(IN)        ::  uni          ! fortran unit#
  INTEGER(KIND=8), INTENT(IN)::  aa           ! buffer di dati
  INTEGER, INTENT(IN)        ::  nbytes       ! numero di bytes da scrivere
  INTEGER, INTENT(IN)        ::  posizione    ! posizione in cui leggere (primo byte=1) se assente sequenziale
  logical, INTENT(IN)        ::  unsigned     ! scrive come unsigned (solo per nbytes= 1 o  2 )
  CHARACTER(LEN=8):: buf

  INTEGER:: numeroByte, bytesScritti
  if (unsigned) then
    numeroByte= 8
   else
    numeroByte= nbytes
  endif

  buf= REPEAT(CHAR(0), 8)
  select case (numeroByte )
      case (1)
        buf= TRANSFER(INT(aa, 1), buf)
      case (2)
        buf= TRANSFER(INT(aa, 2), buf)
      case (4)
        buf= TRANSFER(INT(aa, 4), buf)
      case default
        buf= TRANSFER(INT(aa, 8), buf)
  end select

  call ScriviStream(uni, buf, nbytes, posizione,  bytesScritti)
  if (bytesScritti/= nbytes) stop "WriteStreamInteger_INTERNAL"
end subroutine WriteStreamInteger_INTERNAL



subroutine WriteStreamReal_INTERNAL(uni, aa, nbytes, posizione)
  INTEGER, INTENT(IN)        ::  uni          ! fortran unit#
  REAL(KIND=8), INTENT(IN)   ::  aa           ! buffer di dati
  INTEGER, INTENT(IN)        ::  nbytes       ! numero di bytes da scrivere
  INTEGER, INTENT(IN)        ::  posizione    ! posizione in cui leggere (primo byte=1) se assente sequenziale
  CHARACTER(LEN=8):: buf
  INTEGER:: bytesScritti

  buf= REPEAT(CHAR(0), 8)
  select case (nbytes)
      case (4)
        buf= TRANSFER(REAL(aa, 4), buf)
      case default
        buf= TRANSFER(REAL(aa, 8), buf)
  end select

  call ScriviStream(uni, buf, nbytes, posizione,  bytesScritti)
  if (bytesScritti/= nbytes) stop " WriteStreamReal_INTERNAL"
end subroutine  WriteStreamReal_INTERNAL






subroutine REAL_4_WriteStream(uni, aa, nbytes, posizione)
  INTEGER, INTENT(IN)::            uni          ! fortran unit#
  REAL(KIND=4), INTENT(IN)     ::  aa           ! buffer di dati
  INTEGER, INTENT(IN)::            nbytes       ! numero di bytes da scrivere
  INTEGER, INTENT(IN), OPTIONAL::  posizione    ! posizione in cui leggere (primo byte=1) se assente sequenziale

  INTEGER:: iPosizione

  if (PRESENT(posizione) ) then
    iPosizione= posizione
   else
    iPosizione= Byte_Pos(uni)
  endif
  call WriteStreamREAL_INTERNAL(uni, REAL(aa,8), nbytes, iPosizione)
end subroutine REAL_4_WriteStream






subroutine REAL_8_WriteStream(uni, aa, nbytes, posizione)
  INTEGER, INTENT(IN)::            uni          ! fortran unit#
  REAL(KIND=8), INTENT(IN)     ::  aa           ! buffer di dati
  INTEGER, INTENT(IN)::            nbytes       ! numero di bytes da scrivere
  INTEGER, INTENT(IN), OPTIONAL::  posizione    ! posizione in cui leggere (primo byte=1) se assente sequenziale

  INTEGER:: iPosizione

  if (PRESENT(posizione) ) then
    iPosizione= posizione
   else
    iPosizione= Byte_Pos(uni)
  endif
  call WriteStreamREAL_INTERNAL(uni, REAL(aa,8), nbytes, iPosizione)
end subroutine REAL_8_WriteStream



function Stream_Read_Integer(uni, nbytes, posizione, unsigned, endfile)
  INTEGER(KIND=8):: Stream_Read_Integer
  INTEGER, INTENT(IN)::            uni          ! fortran unit#
  INTEGER, INTENT(IN)::            nbytes       ! numero di bytes da scrivere
  INTEGER, INTENT(IN), OPTIONAL::  posizione    ! posizione in cui leggere (primo byte=1) se assente sequenziale
  logical, INTENT(IN), OPTIONAL::  unsigned     ! scrive come unsigned (solo per nbytes= 1 o  2 )
  LOGICAL, INTENT(OUT), OPTIONAL:: endfile      ! segnala l'incontrata fine del file
  CHARACTER(LEN=8):: buf
  LOGICAL:: unsflag
  INTEGER:: numeroByte, iPosizione, bytesLetti
  INTEGER(KIND=1):: i1
  INTEGER(KIND=2):: i2
  INTEGER(KIND=4):: i4
  INTEGER(KIND=8):: i8
  if (PRESENT(unsigned) ) then
    unsflag= unsigned
   else
    unsflag= .false.
  endif
  if (unsflag) then
    numeroByte= 8
   else
    numeroByte= nbytes
  endif

  if (PRESENT(posizione) ) then
    iPosizione= posizione
   else
    iPosizione= Byte_Pos(uni)
  endif

  call LeggiStream(uni, buf, nbytes, iPosizione,  bytesLetti)
  if (PRESENT( endfile ) ) endfile = (bytesLetti/= nbytes)

  select case (numeroByte )
      case (1)
        Stream_Read_Integer= TRANSFER(buf(:1), i1)
      case (2)
        Stream_Read_Integer= TRANSFER(buf(:2), i2)
      case (4)
        Stream_Read_Integer= TRANSFER(buf    , i4)
      case default
        Stream_Read_Integer= TRANSFER(buf    , i8)
  end select
end function Stream_Read_Integer

function Stream_Read_Real(uni, nbytes, posizione, endfile )
  REAL(KIND=8)       :: Stream_Read_Real
  INTEGER, INTENT(IN)::            uni          ! fortran unit#
  INTEGER, INTENT(IN)::            nbytes       ! numero di bytes da leggere
  INTEGER, INTENT(IN), OPTIONAL::  posizione    ! posizione in cui leggere (primo byte=1) se assente sequenziale
  LOGICAL, INTENT(OUT), OPTIONAL:: endfile      ! segnala l'incontrata fine del file
  CHARACTER(LEN=8):: buf
  INTEGER:: iPosizione, bytesLetti
  REAL(KIND=4):: R4
  REAL(KIND=8):: R8

  if (PRESENT(posizione) ) then
    iPosizione= posizione
   else
    iPosizione= Byte_Pos(uni)
  endif

  call LeggiStream(uni, buf, nbytes, iPosizione,  bytesLetti)
  if (PRESENT( endfile ) ) endfile = (bytesLetti/= nbytes)

  select case (nBytes)
      case (4)
        Stream_Read_Real = TRANSFER(buf    , R4)
      case default
        Stream_Read_Real = TRANSFER(buf    , R8)
  end select
end function Stream_Read_Real


end module  ! ----------------------------------------------------------------------------------------





!program test
!use ByteStreamer
!LOGICAL:: eof
!WRITE(*,*)
!if ( Stream_Open(11, "temp.bin", scrittura=.true.) ) then
!    call Stream_Write(11,  4294967296_8, 8, posizione= 1)
!    call Stream_Write(11,           255, 4)
!    call Stream_Write(11,           255, 1, unsigned=.true.)
!    call Stream_Write(11,           255, 2, unsigned=.true.)
!    call Stream_Write(11,            -1, 1)
!    call Stream_Write(11,      1.0E+000, 4)
!    call Stream_Write(11,      1.0E+000, 8)
!    call Stream_Write(11,      1.0D+256, 8)
!    call Stream_Write(11,      1.0D+256, 4)
!    call Stream_Write(11,    "tutta colpa di Alfredo" , 25)


!    WRITE(*,"(I25,A)"  ) Stream_Read_Integer(uni=11, nbytes=8, posizione= 1)      , "   4294967296"
!    WRITE(*,"(I25,A)"  ) Stream_Read_Integer(uni=11, nbytes=4, posizione= 9)      , "          255"
!    WRITE(*,"(I25,A)"  ) Stream_Read_Integer(uni=11, nbytes=1, unsigned=.true. )  , "          255"
!    WRITE(*,"(I25,A)"  ) Stream_Read_Integer(uni=11, nbytes=2, unsigned=.true. )  , "          255"
!    WRITE(*,"(I25,A)"  ) Stream_Read_Integer(uni=11, nbytes=1)                    , "           -1"
!    WRITE(*,"(E25.1,A)") Stream_Read_Real   (uni=11, nbytes=4)                    , "     1.0E+000"
!    WRITE(*,"(E25.1,A)") Stream_Read_Real   (uni=11, nbytes=8)                    , "     1.0E+000"
!    WRITE(*,"(E25.1,A)") Stream_Read_Real   (uni=11, nbytes=8)                    , "     1.0D+256"
!    WRITE(*,"(E25.1,A)") Stream_Read_Real   (uni=11, nbytes=4)                    , "     1.0D+256 (overflow)"
!    WRITE(*,"(A25  ,A)") Stream_Read_Char   (uni=11, nbytes=25)                   , " tutta colpa di Alfredo"
!    WRITE(*,"(I25,A)"  ) Stream_Read_Integer(uni=11, nbytes=4, posizione= 9)      , "          255"
!   WRITE(*,"(A25  ,A)") Stream_Read_Char   (uni=11, nbytes=30,posizione= 41, ENDFILE=eof)
!   if (eof) WRITE(*,*) "incontrato EoF"
!endif
!call Stream_Close(11)
!end program
