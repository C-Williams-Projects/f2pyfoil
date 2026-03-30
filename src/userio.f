C***********************************************************************
C    Module:  userio.f
C
C    Copyright (C) 2000 Mark Drela
C
C    This program is free software; you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation; either version 2 of the License, or
C    (at your option) any later version.
C
C    This program is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with this program; if not, write to the Free Software
C    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
C***********************************************************************
C
C
C==== pure string-parsing and utility routines
C     (interactive terminal I/O routines removed for library use)
C
C
      SUBROUTINE LC2UC(INPUT)
      CHARACTER*(*) INPUT
C
      CHARACTER*26 LCASE, UCASE
      DATA LCASE / 'abcdefghijklmnopqrstuvwxyz' /
      DATA UCASE / 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' /
C
      N = LEN(INPUT)
C
      DO 10 I=1, N
        K = INDEX( LCASE , INPUT(I:I) )
        IF(K.GT.0) INPUT(I:I) = UCASE(K:K)
 10   CONTINUE
C
      RETURN
      END ! LC2UC



      SUBROUTINE GETINT(INPUT,A,N,ERROR)
      CHARACTER*(*) INPUT
      INTEGER A(*)
      LOGICAL ERROR
C----------------------------------------------------------
C     Parses character string INPUT into an array
C     of integer numbers returned in A(1...N)
C
C     Will attempt to extract no more than N numbers,
C     unless N = 0, in which case all numbers present
C     in INPUT will be extracted.
C
C     N returns how many numbers were actually extracted.
C----------------------------------------------------------
      CHARACTER*130 REC
C
C---- only first 128 characters in INPUT will be parsed
      ILEN = MIN( LEN(INPUT) , 128 )
      ILENP = ILEN + 2
C
C---- put input into local work string (which will be munched)
      REC(1:ILENP) = INPUT(1:ILEN) // ' ,'
C
C---- ignore everything after a "!" character
      K = INDEX(REC,'!')
      IF(K.GT.0) REC(1:ILEN) = REC(1:K-1)
C
      NINP = N
C
C---- count up how many numbers are to be extracted
      N = 0
      K = 1
      DO 10 IPASS=1, ILEN
C------ search for next space or comma starting with current index K
        KSPACE = INDEX(REC(K:ILENP),' ') + K - 1
        KCOMMA = INDEX(REC(K:ILENP),',') + K - 1
C
        IF(K.EQ.KSPACE) THEN
C------- just skip this space
         K = K+1
         GO TO 9
        ENDIF
C
        IF(K.EQ.KCOMMA) THEN
C------- comma found.. increment number count and keep looking
         N = N+1
         K = K+1
         GO TO 9
        ENDIF
C
C------ neither space nor comma found, so we ran into a number...
C-    ...increment number counter and keep looking after next space or comma
        N = N+1
        K = MIN(KSPACE,KCOMMA) + 1
C
  9     IF(K.GE.ILEN) GO TO 11
 10   CONTINUE
C
C---- decide on how many numbers to read, and go ahead and read them
 11   IF(NINP.GT.0) N = MIN( N, NINP )
      READ(REC(1:ILEN),*,ERR=20) (A(I),I=1,N)
      ERROR = .FALSE.
      RETURN
C
C---- bzzzt !!!
 20   CONTINUE
      N = 0
      ERROR = .TRUE.
      RETURN
      END


      SUBROUTINE GETFLT(INPUT,A,N,ERROR)
      CHARACTER*(*) INPUT
      REAL A(*)
      LOGICAL ERROR
C----------------------------------------------------------
C     Parses character string INPUT into an array
C     of real numbers returned in A(1...N)
C
C     Will attempt to extract no more than N numbers,
C     unless N = 0, in which case all numbers present
C     in INPUT will be extracted.
C
C     N returns how many numbers were actually extracted.
C----------------------------------------------------------
      CHARACTER*130 REC
C
C---- only first 128 characters in INPUT will be parsed
      ILEN = MIN( LEN(INPUT) , 128 )
      ILENP = ILEN + 2
C
C---- put input into local work string (which will be munched)
      REC(1:ILENP) = INPUT(1:ILEN) // ' ,'
C
C---- ignore everything after a "!" character
      K = INDEX(REC,'!')
      IF(K.GT.0) REC(1:ILEN) = REC(1:K-1)
C
      NINP = N
C
C---- count up how many numbers are to be extracted
      N = 0
      K = 1
      DO 10 IPASS=1, ILEN
C------ search for next space or comma starting with current index K
        KSPACE = INDEX(REC(K:ILENP),' ') + K - 1
        KCOMMA = INDEX(REC(K:ILENP),',') + K - 1
C
        IF(K.EQ.KSPACE) THEN
C------- just skip this space
         K = K+1
         GO TO 9
        ENDIF
C
        IF(K.EQ.KCOMMA) THEN
C------- comma found.. increment number count and keep looking
         N = N+1
         K = K+1
         GO TO 9
        ENDIF
C
C------ neither space nor comma found, so we ran into a number...
C-    ...increment number counter and keep looking after next space or comma
        N = N+1
        K = MIN(KSPACE,KCOMMA) + 1
C
  9     IF(K.GE.ILEN) GO TO 11
 10   CONTINUE
C
C---- decide on how many numbers to read, and go ahead and read them
 11   IF(NINP.GT.0) N = MIN( N, NINP )
      READ(REC(1:ILEN),*,ERR=20) (A(I),I=1,N)
      ERROR = .FALSE.
      RETURN
C
C---- bzzzt !!!
 20   CONTINUE
      N = 0
      ERROR = .TRUE.
      RETURN
      END



      SUBROUTINE STRIP(STRING,NS)
      CHARACTER*(*) STRING
C-------------------------------------------
C     Strips leading blanks off string
C     and returns length of non-blank part.
C-------------------------------------------
      N = LEN(STRING)
C
C---- find last non-blank character
      DO 10 K2=N, 1, -1
        IF(STRING(K2:K2).NE.' ') GO TO 11
   10 CONTINUE
      K2 = 0
   11 CONTINUE
C
C---- find first non-blank character
      DO 20 K1=1, K2
        IF(STRING(K1:K1).NE.' ') GO TO 21
   20 CONTINUE
   21 CONTINUE
C
C---- number of non-blank characters
      NS = K2 - K1 + 1
      IF(NS.EQ.0) RETURN
C
C---- shift STRING so first character is non-blank
      STRING(1:NS) = STRING(K1:K2)
C
C---- pad tail of STRING with blanks
      DO 30 K=NS+1, N
        STRING(K:K) = ' '
   30 CONTINUE
C
      RETURN
      END
