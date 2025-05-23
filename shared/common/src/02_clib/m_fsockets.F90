!F90 ISO_C_BINDING wrapper for socket communication.
!MG got it from https://github.com/i-pi/i-pi/tree/master/drivers

!Copyright (C) 2013, Michele Ceriotti

!Permission is hereby granted, free of charge, to any person obtaining
!a copy of this software and associated documentation files (the
!"Software"), to deal in the Software without restriction, including
!without limitation the rights to use, copy, modify, merge, publish,
!distribute, sublicense, and/or sell copies of the Software, and to
!permit persons to whom the Software is furnished to do so, subject to
!the following conditions:

!The above copyright notice and this permission notice shall be included
!in all copies or substantial portions of the Software.

!THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
!EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
!MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
!IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
!CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
!TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
!SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


!Contains both the functions that transmit data to the socket and read the data
!back out again once finished, and the function which opens the socket initially.

!Functions:
!   open_socket: Opens a socket with the required host server, socket type and port number.
!   write_buffer: Writes a string to the socket.
!   read_buffer: Reads data from the socket.

MODULE m_fsockets

  USE ISO_C_BINDING

  IMPLICIT NONE

  INTERFACE writebuffer
      MODULE PROCEDURE writebuffer_s, writebuffer_d, writebuffer_dv, writebuffer_i
  END INTERFACE writebuffer

  INTERFACE readbuffer
      MODULE PROCEDURE readbuffer_s, readbuffer_dv, readbuffer_d, readbuffer_i
  END INTERFACE readbuffer

  INTERFACE
    SUBROUTINE open_csocket(psockfd, inet, port, host) BIND(C, name="open_socket")
      USE ISO_C_BINDING
      INTEGER(KIND=C_INT)                      :: psockfd, inet, port
      CHARACTER(KIND=C_CHAR), DIMENSION(*)     :: host
    END SUBROUTINE open_csocket

    SUBROUTINE writebuffer_csocket(psockfd, pdata, plen) BIND(C, name="writebuffer")
      USE ISO_C_BINDING
      INTEGER(KIND=C_INT)                      :: psockfd
      TYPE(C_PTR), VALUE                       :: pdata
      INTEGER(KIND=C_INT)                      :: plen
    END SUBROUTINE writebuffer_csocket

    SUBROUTINE readbuffer_csocket(psockfd, pdata, plen) BIND(C, name="readbuffer")
      USE ISO_C_BINDING
      INTEGER(KIND=C_INT)                      :: psockfd
      TYPE(C_PTR), VALUE                       :: pdata
      INTEGER(KIND=C_INT)                      :: plen
    END SUBROUTINE readbuffer_csocket
  END INTERFACE

CONTAINS

   SUBROUTINE open_socket(psockfd, inet, port, host)
      INTEGER, INTENT(IN) :: inet, port
      INTEGER, INTENT(OUT) :: psockfd
      CHARACTER(LEN=1024), INTENT(IN) :: host
      CHARACTER(LEN=1,KIND=C_CHAR) :: chost(1024)

      CALL fstr2cstr(host, chost)
      CALL open_csocket(psockfd, inet, port, host)
   END SUBROUTINE open_socket

   SUBROUTINE fstr2cstr(fstr, cstr, plen)
      CHARACTER(LEN=*), INTENT(IN) :: fstr
      CHARACTER(LEN=1,KIND=C_CHAR), INTENT(OUT) :: cstr(:)
      INTEGER, INTENT(IN), OPTIONAL :: plen

      INTEGER i,n
      IF (PRESENT(plen)) THEN
         n = plen
         DO i=1,n
            cstr(i) = fstr(i:i)
         ENDDO
      ELSE
         n = LEN_TRIM(fstr)
         DO i=1,n
            cstr(i) = fstr(i:i)
         ENDDO
         cstr(n+1) = C_NULL_CHAR
      END IF
   END SUBROUTINE fstr2cstr

  SUBROUTINE writebuffer_d (psockfd, fdata)
    USE ISO_C_BINDING
    INTEGER, INTENT(IN)                      :: psockfd
    REAL(KIND=8), INTENT(IN)                :: fdata

    REAL(KIND=C_DOUBLE), TARGET              :: cdata

    cdata = fdata
    CALL writebuffer_csocket(psockfd, c_loc(cdata), 8)
  END SUBROUTINE writebuffer_d

  SUBROUTINE writebuffer_i (psockfd, fdata)
    USE ISO_C_BINDING
    INTEGER, INTENT(IN)                      :: psockfd, fdata

    INTEGER(KIND=C_INT), TARGET              :: cdata

    cdata = fdata
    CALL writebuffer_csocket(psockfd, c_loc(cdata), 4)
  END SUBROUTINE writebuffer_i

  SUBROUTINE writebuffer_s (psockfd, fstring, plen)
    USE ISO_C_BINDING
    INTEGER, INTENT(IN)                      :: psockfd
    CHARACTER(LEN=*), INTENT(IN)             :: fstring
    INTEGER, INTENT(IN)                      :: plen

    INTEGER                                  :: i
    CHARACTER(LEN=1, KIND=C_CHAR), TARGET    :: cstring(plen)

    DO i = 1,plen
      cstring(i) = fstring(i:i)
    ENDDO
    CALL writebuffer_csocket(psockfd, c_loc(cstring(1)), plen)
  END SUBROUTINE writebuffer_s

  SUBROUTINE writebuffer_dv(psockfd, fdata, plen)
    USE ISO_C_BINDING
    INTEGER, INTENT(IN)                      :: psockfd, plen
    REAL(KIND=8), INTENT(IN), TARGET        :: fdata(plen)

    CALL writebuffer_csocket(psockfd, c_loc(fdata(1)), 8*plen)
  END SUBROUTINE writebuffer_dv

  SUBROUTINE readbuffer_d (psockfd, fdata)
    USE ISO_C_BINDING
    INTEGER, INTENT(IN)                     :: psockfd
    REAL(KIND=8), INTENT(OUT)               :: fdata

    REAL(KIND=C_DOUBLE), TARGET              :: cdata

    CALL readbuffer_csocket(psockfd, c_loc(cdata), 8)
    fdata=cdata
  END SUBROUTINE readbuffer_d

  SUBROUTINE readbuffer_i (psockfd, fdata)
    USE ISO_C_BINDING
    INTEGER, INTENT(IN)                      :: psockfd
    INTEGER, INTENT(OUT)                     :: fdata

    INTEGER(KIND=C_INT), TARGET              :: cdata

    CALL readbuffer_csocket(psockfd, c_loc(cdata), 4)
    fdata = cdata
  END SUBROUTINE readbuffer_i

  SUBROUTINE readbuffer_s (psockfd, fstring, plen)
    USE ISO_C_BINDING
    INTEGER, INTENT(IN)                      :: psockfd
    CHARACTER(LEN=*), INTENT(OUT)            :: fstring
    INTEGER, INTENT(IN)                      :: plen

    INTEGER                                  :: i
    CHARACTER(LEN=1, KIND=C_CHAR), TARGET    :: cstring(plen)

    CALL readbuffer_csocket(psockfd, c_loc(cstring(1)), plen)
    fstring=""
    DO i = 1,plen
       fstring(i:i) = cstring(i)
    ENDDO
  END SUBROUTINE readbuffer_s

  SUBROUTINE readbuffer_dv(psockfd, fdata, plen)
    USE ISO_C_BINDING
    INTEGER, INTENT(IN)                      :: psockfd, plen
    REAL(KIND=8), INTENT(OUT), TARGET       :: fdata(plen)

    CALL readbuffer_csocket(psockfd, c_loc(fdata(1)), 8*plen)
  END SUBROUTINE readbuffer_dv

  SUBROUTINE socket_from_string (srvaddress, socket)
    CHARACTER(len=*), INTENT(IN)  :: srvaddress
    integer, INTENT(out)  :: socket
    CHARACTER(len=len_trim(srvaddress)) :: address
    INTEGER :: port, inet, field_sep_pos

    ! Parses host name, port and socket type
    field_sep_pos = INDEX(srvaddress, ':', back=.true.)
    address = srvaddress(1:field_sep_pos-1)

    ! Check if UNIX type socket
    IF (trim(srvaddress(field_sep_pos+1 :)) == 'UNIX') then
       port = 1234 ! just a place-holder
       inet = 0
       !write(*,*) " Connecting to `", trim(address), "` using UNIX socket"
    ELSE
       read ( srvaddress ( field_sep_pos+1 : ), * ) port
       inet = 1
       !write(*,*) " Connecting to `", trim(address), ":", srvaddress (field_sep_pos+1:), " using INET socket"
    END IF

    ! Create the socket
    CALL open_socket (socket, inet, port, trim(address)//achar(0))

  END SUBROUTINE socket_from_string

END MODULE m_fsockets
