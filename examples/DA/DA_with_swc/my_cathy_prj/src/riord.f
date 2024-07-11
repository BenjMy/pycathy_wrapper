C
C**************************  RIORD  ************************************
C
C  sort into ascending order the nodes connecting to each triangular
C  element (case NOD=3) or each tetrahedral element (case NOD=4).
C  NELT is the number of elements and ELTCON gives the element
C  connectivities.
C
C***********************************************************************
C
      SUBROUTINE RIORD(NELT,ELTCON,NOD)
C
      IMPLICIT  NONE
      INTEGER   I,J,K,IIE
      INTEGER   NELT,NOD
      INTEGER   ELTCON(NOD+1,*),E(4)
C
      DO I=1,NELT
         DO J=1,NOD
            E(J)=ELTCON(J,I)
         END DO  
         DO K=1,NOD-1
            DO J=K+1,NOD
               IF(E(K).GT.E(J)) THEN
                  IIE=E(K)
                  E(K)=E(J)
                  E(J)=IIE
               END IF
            END DO  
         END DO  
         DO J=1,NOD
            ELTCON(J,I)=E(J)
         END DO  
      END DO  
C
      RETURN
      END
