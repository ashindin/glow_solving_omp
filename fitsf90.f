C------------------------------------------------------------------------------
C       The following routines are part of the FITSIO library
C       and are specific to Fortran-90 compilers.
C------------------------------------------------------------------------------
C   This software was prepared by High Energy Astrophysic Science Archive
C   Research Center (HEASARC) at the NASA Goddard Space Flight Center. Users
C   shall not, without prior written permission of the U.S. Government,
C   establish a claim to statutory copyright.  The Government and others acting
C   on its behalf, shall have a royalty-free, non-exclusive, irrevocable,
C   worldwide license for Government purposes to publish, distribute,
C   translate, copy, exhibit, and perform such material. 
C--------------------------------------------------------------------------
        subroutine ftopnf(funit,fname,oldnew,rwmode,block,size,status)

C       low-level, machine-dependent routine to create or open a new file 
C
C       funit   i  Fortran I/O unit number
C       fname   c  name of file to be opened
C       oldnew  i  file status: 0 = open old/existing file; else open new file
C       rwmode  i  file access mode: 0 = readonly; else = read/write
C       block   i  FITS record blocking factor 
C       size    i  min size of file, in bytes
C       status  i  returned error status (0=ok)
C
C       written by Wm Pence, HEASARC/GSFC, December 1996

        integer funit,oldnew,rwmode,block,size,status,lenrec
        character*(*) fname
        character fstat*3, fmode*9
        character dummy*2880

C       some compilers require that dummy be initialized!
        dummy = ' '
        if (oldnew .eq. 0)then
                fstat='OLD'
C               simply return the default block size
                block=1
C               no easy way to determine the size of the file, 
C               so just set value to a huge number
                size=1000000000
        else
C               create a new file
                fstat='NEW'
                size = 0
        end if

        if(rwmode .eq. 0) then
           fmode = 'READ'
        else
           fmode = 'READWRITE'
        end if

C       determine the record length needed to read or write 2880-byte records
C       (could be 720 words or 2880 bytes, depending on the compiler)

        inquire(iolength=lenrec) dummy
        open(unit=funit,file=fname,status=fstat,err=900, action=fmode,
     &       recl=lenrec,form='UNFORMATTED',access='DIRECT')
        return

C       error opening file:
900     status=104 + oldnew
        end
C----------------------------------------------------------------------
        subroutine ftgsdt(dd,mm,yy,status)

C       get the current date from the system

C       dd      i  day of the month (1-31)
C       mm      i  month of the year (1-12)
C       yy      i  last 2 digits of the year (1992 = 92, 2001 = 01)

C     This routine is not valid after the year 1999; must be modified
C     before the year 2000.

        integer dd,mm,yy,status
        integer iarray(8)

        if (status .gt. 0)return

        call date_and_time(values=iarray)

        dd=iarray(3)
        mm=iarray(2)
        yy=mod(iarray(1),100)
        end
C----------------------------------------------------------------------
        subroutine ftpbyt(ounit,nbytes,array,status)

C       write string of data bytes to output buffer.

C       ounit   i  fortran unit number
C       nbytes  i  number of bytes
C       array   i  integer array
C       status  i  output error status
C
C       written by Wm Pence, HEASARC/GSFC, June 1991

        integer nbytes,ounit,status
        integer array((nbytes-1) / 4 + 1)
        character temp

C       Trick the compiler into passing the integer array to a routine
C       that is expecting a character string argument.
C       This is very ugly, but I/O efficiency is of primary importance.

        call ftpxbf(ounit,nbytes,array,temp,status)
        end
C--------------------------------------------------------------------------
        subroutine ftpxbf(ounit,nbytes,array,temp,status)
        integer temp
        integer nbytes,ounit,status
        character*(nbytes) array

C       array is now declared as a string!
        call ftpcbf(ounit,nbytes,array,status) 
        end
C--------------------------------------------------------------------------
        subroutine ftpbyo(ounit,gsize,ngroup,offset,array,status)
 
C       "Put Bytes with Offsets"
C       copy input buffer of bytes to the output character buffer.
 
C       ounit   i  Fortran output unit number
C       gsize   i  size of each group of bytes
C       ngroup  i  number of groups to write
C       offset  i  size of gap between groups
C       array   i  input array of bytes
C       status  i  output error status (0 = ok)
C
C       written by Wm Pence, HEASARC/GSFC, Dec 1996
 
        integer ounit,gsize,ngroup,offset,status
        integer array(*)
        character temp

C       Trick the compiler into passing the integer array to a routine
C       that is expecting a character string argument.
C       This is very ugly, but I/O efficiency is of primary importance.

        call ftpxbo(ounit,gsize,ngroup,offset,array,temp,status)
        end
C--------------------------------------------------------------------------
        subroutine ftpxbo(ounit,gsize,ngroup,offset,array,temp,status)
        integer temp
        integer ounit,gsize,ngroup,offset,status
        character*(gsize*ngroup) array

C       array is now declared as a string!
        call ftpcbo(ounit,gsize,ngroup,offset,array,status)
        end
C----------------------------------------------------------------------
        subroutine ftgbyt(iunit,nbytes,array,status)

C       read string of data bytes from input buffer.

C       iunit   i  fortran unit number
C       nbytes  i  number of bytes
C       array   i  integer array
C       status  i  output error status
C
C       written by Wm Pence, HEASARC/GSFC, June, 1991

        integer nbytes,iunit,status
        integer array(*)
        character temp

C       Trick the compiler into passing the integer array to a routine
C       that is expecting a character string argument.
C       This is very ugly, but I/O efficiency is of primary importance.

        call ftgxbf(iunit,nbytes,array,temp,status)

        end
C--------------------------------------------------------------------------
        subroutine ftgxbf(iunit,nbytes,array,temp,status)
        integer temp
        integer nbytes,iunit,status
        character*(nbytes) array

C       array is now declared as a string!
        call ftgcbf(iunit,nbytes,array,status) 
        end
C--------------------------------------------------------------------------
        subroutine ftgbyo(iunit,gsize,ngroup,offset,array,status)

C       "Get BYtes with Offsets"
C       read bytes from the character buffer.

C       iunit   i  Fortran output unit number
C       gsize   i  size of each group of bytes
C       ngroup  i  number of groups to read
C       offset  i  size of gap between groups
C       array   i  output array of bytes
C       status  i  output error status (0 = ok)
C
C       written by Wm Pence, HEASARC/GSFC, Dec 1996

        integer iunit,gsize,ngroup,offset,status
        integer array(*)
        character temp

C       Trick the compiler into passing the integer array to a routine
C       that is expecting a character string argument.
C       This is very ugly, but I/O efficiency is of primary importance.

        call ftgxbo(iunit,gsize,ngroup,offset,array,temp,status)
        end
C--------------------------------------------------------------------------
        subroutine ftgxbo(iunit,gsize,ngroup,offset,array,temp,status)
        integer temp
        integer iunit,gsize,ngroup,offset,status
        character*(gsize*ngroup) array

C       array is now declared as a string!
        call ftgcbo(iunit,gsize,ngroup,offset,array,status)
        end
C------------------------------------------------------------------------
        subroutine ieevud(dbl1, dbl2, int)
        double precision dbl1,dbl2
        integer int
C       dummy routine; only used on Vax VMS computers
        end
C------------------------------------------------------------------------
        subroutine ieevpd(dbl1, dbl2, int)
        double precision dbl1,dbl2
        integer int
C       dummy routine; only used on Vax VMS computers
        end
C------------------------------------------------------------------------
        integer function cray2ieg(i1,i2,i3,i4,i5,i6,s1)
        integer i1,i2,i3,i4,i5,i6
        character s1
C       dummy routine; only used on Cray supercomputers
        cray2ieg=0
        end
C------------------------------------------------------------------------
        integer function ieg2cray(i1,i2,i3,i4,i5,i6,s1)
        integer i1,i2,i3,i4,i5,i6
        character s1
C       dummy routine; only used on Cray supercomputers
        ieg2cray=0
        end
