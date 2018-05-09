C----------------------------------------------------------------------
        subroutine ftopnx(funit,fname,oldnew,rwmode,block,status)

C       low-level, machine-dependent routine to create or open a new file 
C
C       funit   i  Fortran I/O unit number
C       fname   c  name of file to be opened
C       oldnew  i  file status: 0 = open old/existing file; else open new file
C       rwmode  i  file access mode: 0 = readonly; else = read/write
C       block   i  FITS record blocking factor 
C       status  i  returned error status (0=ok)
C
C       written by Wm Pence, HEASARC/GSFC, June 1991
C       modified Feb 1995

        integer funit,oldnew,rwmode,block,status,i,ibuff,inital,size
        character*(*) fname
        logical igneof,found

C       COMMON BLOCK DEFINITIONS:--------------------------------------------
        integer nb,ne,nf
        parameter (nb = 20)
        parameter (ne = 512)
        parameter (nf = 3000)

        integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
        integer nxtfld
        logical wrmode
        common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb),
     &  wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld

        integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
        integer theap
        double precision tscale,tzero
        common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb),
     &  tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb)
     &  ,theap(nb)

        integer buflun,currnt,reclen,bytnum,maxrec
        common/ftlbuf/buflun(nb),currnt(nb),reclen(nb),
     &  bytnum(nb),maxrec(nb)

        integer pb
        parameter (pb = 20)
        integer maxbuf,logbuf,recnum,pindex
        logical modify
        common/ftpbuf/maxbuf,logbuf(pb),recnum(pb),modify(pb),
     &  pindex(pb)

        integer compid
        common/ftcpid/compid
C       END OF COMMON BLOCK DEFINITIONS-----------------------------------

        real rword
        double precision dword

        save inital
        data inital/0/

        if (status .gt. 0)return

        if (inital .eq. 0)then
C           first time through need to initialize pointers
            nxtfld=0
            maxbuf=pb
            do 2 i=1,nb
                buflun(i)=0
2           continue
            do 4 i=1,pb
                logbuf(i)=0
                recnum(i)=0
                modify(i)=.false.
                pindex(i)=i
4           continue
            inital=1

C           Determine at run time what type of machine we are running on.
C           Initialize a real and double value to arbitrary values.
            rword=1.1111111111
            dword=1.1111111111D+00

C           ftarch looks at the equivalent integer value 
            call ftarch(rword,dword,compid)
        end if

C       check for valid unit number
        if (funit .lt. 1 .or. funit .gt. 199)then
                status=101
                return
        end if

C       find available logical buffer slot for this file
        do 10 i=1,nb
                if (buflun(i) .eq. 0)then
                        ibuff=i
                        go to 20
                end if
10      continue

C       error: no vacant logical buffer slots left
        status=102
        return

20      continue

        if (oldnew .eq. 0)then
            igneof = .false.
C           test if file exists
            inquire(file=fname,exist=found)
            if (.not. found)then
C               error: file doesn't exist??
                status=103
                return
            end if
        else
            igneof = .true.
        end if

        call ftopnf(funit,fname,oldnew,rwmode,block,size,status)

C       initialize the HDU parameters
        maxrec(ibuff)=size

        if (oldnew .eq. 1 .or. block .le. 1)then
C           new files always have a record length of 2880 bytes
            reclen(ibuff)=2880
        else
            reclen(ibuff)=block
        end if

        bufnum(funit)=ibuff
        chdu(ibuff)=1
        hdutyp(ibuff)=0
        maxhdu(ibuff)=1
        hdstrt(ibuff,1)=0
        hdend(ibuff)=0
        nxthdr(ibuff)=0
C       data start location is undefined
        dtstrt(ibuff)=-1000000000

        heapsz(ibuff)=0
        theap(ibuff)=0
        tfield(ibuff)=0
        rowlen(ibuff)=0

C       initialize the logical buffer parameters
        buflun(ibuff)=funit
        currnt(ibuff)=0

        if (rwmode .eq. 0)then
                wrmode(ibuff)=.false.
        else
                wrmode(ibuff)=.true.
        end if

C       load the first record of the file
        call ftldrc(funit,1,igneof,status)
        end
C--------------------------------------------------------------------------
        subroutine ftclsx(iunit,keep,status)

C       low level routine to close a file
C
C       iunit   i  Fortran I/O unit number
C       keep    l  keep the file? (else delete it)
C       status  i  returned error status (0=ok)
C
C       written by Wm Pence, HEASARC/GSFC, Aug 1992

        integer iunit,status
        logical keep

C       COMMON BLOCK DEFINITIONS:--------------------------------------------
        integer nb,ne
        parameter (nb = 20)
        parameter (ne = 512)

        integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
        integer nxtfld
        logical wrmode
        common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb),
     &  wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld

        integer buflun,currnt,reclen,bytnum,maxrec
        common/ftlbuf/buflun(nb),currnt(nb),reclen(nb),
     &  bytnum(nb),maxrec(nb)
C       END OF COMMON BLOCK DEFINITIONS-----------------------------------

        integer ibuff

        ibuff=bufnum(iunit)

        if (ibuff .eq. 0)return

C       reset file common block parameters
        bufnum(iunit)=0
        buflun(ibuff)=0
        wrmode(ibuff)=.false.
        currnt(ibuff)=0
        reclen(ibuff)=0
        bytnum(ibuff)=0

        if (keep)then
                close(iunit,err=900)
        else
                close(iunit,status='DELETE',err=900)
        end if
        return

900     continue
C       set error code, if it has not previous been set
        if (status .le. 0)status=110
        end
C----------------------------------------------------------------------
        subroutine ftmbyt(iunit,bytno,igneof,status)

C       move i/o pointer so that it is pointing to the byte number BYTNUM
C       in the FITS file.  Subsequent read or write operations will begin
C       at this point.

C       iunit   i  fortran unit number
C       bytno   i  number of the byte to point to.
C       igneof  l  ignore end-of-file (107) error?
C       status  i  output error status
C
C       written by Wm Pence, HEASARC/GSFC, June, 1991
C       rewritten Feb, 1995

        integer iunit,bytno,status
        logical igneof

C       COMMON BLOCK DEFINITIONS:--------------------------------------------
        integer nb,ne,pb 
        parameter (nb = 20)
        parameter (ne = 512)
        parameter (pb = 20)
        integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
        integer nxtfld
        logical wrmode
        common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb),
     &  wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld

        integer buflun,currnt,reclen,bytnum,maxrec
        common/ftlbuf/buflun(nb),currnt(nb),reclen(nb),
     &  bytnum(nb),maxrec(nb)

        integer maxbuf,logbuf,recnum,pindex
        logical modify
        common/ftpbuf/maxbuf,logbuf(pb),recnum(pb),modify(pb),
     &  pindex(pb)
C       END OF COMMON BLOCK DEFINITIONS-----------------------------------

        integer lbuff,record,byten

        if (status .gt. 0)then
                return
        else if (bytno .lt. 0)then
C               error: negative byte number
                status=304
        else
                lbuff=bufnum(iunit)

C               calculate the record number and byte offset to move to
                record=bytno/reclen(lbuff)+1
                byten=mod(bytno,reclen(lbuff))

                if (record .ne. recnum(currnt(lbuff)))then
C                       not the current record, so load the new record; 
                        call ftldrc(iunit,record,igneof,status)
                end if
                bytnum(lbuff)=byten
        end if
        end
C--------------------------------------------------------------------------
        subroutine ftpcbf(ounit,nchar,cbuff,status)

C       "Put Character BuFfer"
C       copy input buffer of characters to the output character buffer.

C       ounit   i  Fortran output unit number
C       nchar   i  number of characters in the string
C       cbuff   c  input character string
C       status  i  output error status (0 = ok)
C
C       written by Wm Pence, HEASARC/GSFC, June 1991
C       modified Feb 1995

        integer ounit,nchar,status
        character*(*) cbuff

C       COMMON BLOCK DEFINITIONS:--------------------------------------------
        integer nb,ne,pb
        parameter (nb = 20)
        parameter (ne = 512)
        parameter (pb = 20)

        integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
        integer nxtfld
        logical wrmode
        common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb),
     &  wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld

        integer buflun,currnt,reclen,bytnum,maxrec
        common/ftlbuf/buflun(nb),currnt(nb),reclen(nb),
     &  bytnum(nb),maxrec(nb)

        integer maxbuf,logbuf,recnum,pindex
        logical modify
        common/ftpbuf/maxbuf,logbuf(pb),recnum(pb),modify(pb),
     &  pindex(pb)

C       have to use separate character arrays because of compiler limitations
        character*2880 b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,
     &  b15,b16,b17,b18,b19,b20
        common /ftbuff/b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,
     &  b15,b16,b17,b18,b19,b20
C       END OF COMMON BLOCK DEFINITIONS-----------------------------------

        integer lbuff,pbuff,buflen,lastb,nleft,in1,nbyt,nrec

        if (status .gt. 0)return

        lbuff=bufnum(ounit)
        buflen=reclen(lbuff)

        if (nchar .lt. 0)then
C               error: negative number of bytes to write
                status=306
                return
        else if (.not. wrmode(lbuff))then
C           don't have write access to this file
            status=112
            return
        end if

C       lastb   = position of last byte read from input buffer
C       nleft   = number of bytes left in the input buffer
C       in1     = position of first byte remaining in the input buffer
C       nbyt    = number of bytes to transfer from input to output

        nleft=nchar
        in1=1

C       find the number of bytes that will fit in output buffer
200     pbuff=currnt(lbuff)
        lastb=bytnum(lbuff)
        nbyt=min(nleft,buflen-lastb)
        if (nbyt .gt. 0)then
C           append the input buffer to the output physical buffer
            go to (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,
     &      19,20)pbuff

C               if got here, then pbuff is out of range
                status=101
                return

1               b1(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
                go to 100
2               b2(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
                go to 100
3               b3(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
                go to 100
4               b4(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
                go to 100
5               b5(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
                go to 100
6               b6(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
                go to 100
7               b7(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
                go to 100
8               b8(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
                go to 100
9               b9(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
                go to 100
10              b10(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
                go to 100
11              b11(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
                go to 100
12              b12(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
                go to 100
13              b13(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
                go to 100
14              b14(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
                go to 100
15              b15(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
                go to 100
16              b16(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
                go to 100
17              b17(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
                go to 100
18              b18(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
                go to 100
19              b19(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
                go to 100
20              b20(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)

100         modify(pbuff)=.true.
            bytnum(lbuff)=bytnum(lbuff)+nbyt
            in1=in1+nbyt
            nleft=nleft-nbyt
        end if

C       process more bytes, if any
        if (nleft .gt. 0)then
          nrec=recnum(pbuff)+1

          if (nleft .gt. buflen)then
C           first, flush any current buffers to disk
            call ftflsh(lbuff,status)

C           write whole blocks directly to the FITS file by-passing buffers
150         write(ounit,rec=nrec,err=900)cbuff(in1:in1+buflen-1)
            in1=in1+buflen
            nleft=nleft-buflen
            bytnum(lbuff)=bytnum(lbuff)+buflen
            nrec=nrec+1
            if (nleft .gt. buflen)go to 150

C           Save maximum record written, for comparison in ftread
            maxrec(lbuff) = max(maxrec(lbuff), nrec-1)
          end if

C         load the next file record into a physical buffer
          call ftldrc(ounit,nrec,.true.,status)
          if (status .gt. 0)return
          go to 200
        end if
        return

C       come here if there was a disk write error of some sort
900     status=106
        end
C--------------------------------------------------------------------------
        subroutine ftpcbo(ounit,gsize,ngroup,offset,cbuff,status)
 
C       "Put Character BuFfer with Offsets"
C       copy input buffer of characters to the output character buffer.
 
C       ounit   i  Fortran output unit number
C       gsize   i  size of each group of bytes
C       ngroup  i  number of groups to write
C       offset  i  size of gap between groups
C       cbuff   c  input character string
C       status  i  output error status (0 = ok)
C
C       written by Wm Pence, HEASARC/GSFC, Dec 1996
 
        integer ounit,gsize,ngroup,offset,status
        character cbuff*(*)
 
C       COMMON BLOCK DEFINITIONS:--------------------------------------------
        integer nb,ne,pb
        parameter (nb = 20)
        parameter (ne = 512)
        parameter (pb = 20)
 
        integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
        integer nxtfld
        logical wrmode
        common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb),
     &  wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
 
        integer buflun,currnt,reclen,bytnum,maxrec
        common/ftlbuf/buflun(nb),currnt(nb),reclen(nb),
     &  bytnum(nb),maxrec(nb)
 
        integer maxbuf,logbuf,recnum,pindex
        logical modify
        common/ftpbuf/maxbuf,logbuf(pb),recnum(pb),modify(pb),
     &  pindex(pb)
 
C       have to use separate character arrays because of compiler limitations
        character*2880 b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,
     &  b15,b16,b17,b18,b19,b20
        common /ftbuff/b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,
     &  b15,b16,b17,b18,b19,b20
C       END OF COMMON BLOCK DEFINITIONS-----------------------------------
 
        integer lbuff,pbuff,buflen,lastb,nleft,in1,nbyt
        integer i,bytno,record,oldrec,incre
 
        if (status .gt. 0)return
 
        lbuff=bufnum(ounit)

        if (.not. wrmode(lbuff))then
C           don't have write access to this file
            status=112
            return
        end if

        buflen=reclen(lbuff)
        pbuff =currnt(lbuff)
        oldrec=recnum(pbuff)
C       lastb = position of last byte read or written in FITS buffer
        lastb =bytnum(lbuff)
        bytno =(oldrec-1) * buflen + lastb
C       in1   = position of first byte remaining in the input buffer
        in1   =1
        incre =gsize+offset
        nbyt  = 0

        do 500 i = 1,ngroup

C           nleft   = number of bytes left in the input buffer 
            nleft=gsize
C           nbyt    = number of bytes to transfer from input to output
            nbyt=min(nleft,buflen-lastb)
            if (nbyt .eq. 0)go to 300

200         continue
C           append the input buffer to the output physical buffer
            go to (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,
     &         19,20)pbuff
 
C               if got here, then pbuff is out of range
                status=101
                return
 
1               b1(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
                go to 100
2               b2(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
                go to 100
3               b3(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
                go to 100
4               b4(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
                go to 100
5               b5(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
                go to 100
6               b6(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
                go to 100
7               b7(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
                go to 100
8               b8(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
                go to 100
9               b9(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
                go to 100
10              b10(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
                go to 100
11              b11(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
                go to 100
12              b12(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
                go to 100
13              b13(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
                go to 100
14              b14(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
                go to 100
15              b15(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
                go to 100
16              b16(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
                go to 100
17              b17(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
                go to 100
18              b18(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
                go to 100
19              b19(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)
                go to 100
20              b20(lastb+1:lastb+nbyt)=cbuff(in1:in1+nbyt-1)

100         in1=in1+nbyt
            nleft=nleft-nbyt
 
C           process more bytes, if any
300         continue
            if (nleft .gt. 0)then
C               entire group did not fit in the buffer
C               load the next file record into a physical buffer
                oldrec=oldrec+1
                modify(pbuff)=.true.
                call ftldrc(ounit,oldrec,.true.,status)
                if (status .gt. 0)return
                pbuff=currnt(lbuff)
                lastb=0
                nbyt=nleft
                go to 200
            end if
 
            if (i .ne. ngroup)then
C               move to the position of the next group
                bytno=bytno+incre
                record=bytno/buflen+1
                lastb=mod(bytno,buflen)
 
                if (record .ne. oldrec)then
C                   not the current record, so load the new record;
                    modify(pbuff)=.true.
                    call ftldrc(ounit,record,.true.,status)
                    if (status .gt. 0)return
                    oldrec=record
                    pbuff=currnt(lbuff)
                end if
            end if
500     continue

        modify(pbuff)=.true.
        bytnum(lbuff)=lastb+nbyt
        end
C--------------------------------------------------------------------------
        subroutine ftgcbf(iunit,nchar,array,status)

C       "Get Character BuFfer"
C       read NCHAR characters from the character buffer.

C       iunit   i  Fortran unit number for reading from disk
C       nchar   i  number of characters to read
C       array   c  output character string
C       status  i  output error status (0 = ok)
C
C       written by Wm Pence, HEASARC/GSFC, June 1991
C       modified Feb 1995

        integer iunit,nchar,status
        character*(*) array

C       COMMON BLOCK DEFINITIONS:--------------------------------------------
        integer nb,ne,pb
        parameter (nb = 20)
        parameter (ne = 512)
        parameter (pb = 20)

        integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
        integer nxtfld
        logical wrmode
        common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb),
     &  wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld

        integer buflun,currnt,reclen,bytnum,maxrec
        common/ftlbuf/buflun(nb),currnt(nb),reclen(nb),
     &  bytnum(nb),maxrec(nb)

        integer maxbuf,logbuf,recnum,pindex
        logical modify
        common/ftpbuf/maxbuf,logbuf(pb),recnum(pb),modify(pb),
     &  pindex(pb)

C       have to use separate character arrays because of compiler limitations
        character*2880 b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,
     &  b15,b16,b17,b18,b19,b20
        common /ftbuff/b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,
     &  b15,b16,b17,b18,b19,b20
C       END OF COMMON BLOCK DEFINITIONS-----------------------------------

        integer nleft,nbyt,lastb,in1,lbuff,pbuff,buflen,nrec,ios,i

        if (status .gt. 0)return

        if (nchar .lt. 0)then
C               error: negative number of bytes to read
                status=306
                return
        end if

        lbuff=bufnum(iunit)
        buflen=reclen(lbuff)

C       lastb   = position of last byte read from input buffer
C       nleft   = number of bytes left in the input buffer
C       in1     = position of first byte remaining in the input buffer
C       nbyt    = number of bytes to transfer from input to output

        nleft=nchar
        in1=1

C       find the number of remaining bytes that can be read from buffer
200     pbuff=currnt(lbuff)
        lastb=bytnum(lbuff)
        nbyt=min(nleft,buflen-lastb)

C       get characters from the physical buffer to the output string
        if (nbyt .gt. 0)then

            go to (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,
     &      19,20)pbuff

C               if got here, then pbuff is out of range
                status=101
                return

1               array(in1:in1+nbyt-1)=b1(lastb+1:lastb+nbyt)
                go to 100
2               array(in1:in1+nbyt-1)=b2(lastb+1:lastb+nbyt)
                go to 100
3               array(in1:in1+nbyt-1)=b3(lastb+1:lastb+nbyt)
                go to 100
4               array(in1:in1+nbyt-1)=b4(lastb+1:lastb+nbyt)
                go to 100

C  The SUN F90 compiler gives a segmentation fault on the following
C  statement when executing testprog, while reading a complex (C) column
C  when using the Linux F90 routines (fitsf90_nag.f).
5               array(in1:in1+nbyt-1)=b5(lastb+1:lastb+nbyt)
                go to 100
6               array(in1:in1+nbyt-1)=b6(lastb+1:lastb+nbyt)
                go to 100
7               array(in1:in1+nbyt-1)=b7(lastb+1:lastb+nbyt)
                go to 100
8               array(in1:in1+nbyt-1)=b8(lastb+1:lastb+nbyt)
                go to 100
9               array(in1:in1+nbyt-1)=b9(lastb+1:lastb+nbyt)
                go to 100
10              array(in1:in1+nbyt-1)=b10(lastb+1:lastb+nbyt)
                go to 100
11              array(in1:in1+nbyt-1)=b11(lastb+1:lastb+nbyt)
                go to 100
12              array(in1:in1+nbyt-1)=b12(lastb+1:lastb+nbyt)
                go to 100
13              array(in1:in1+nbyt-1)=b13(lastb+1:lastb+nbyt)
                go to 100
14              array(in1:in1+nbyt-1)=b14(lastb+1:lastb+nbyt)
                go to 100
15              array(in1:in1+nbyt-1)=b15(lastb+1:lastb+nbyt)
                go to 100
16              array(in1:in1+nbyt-1)=b16(lastb+1:lastb+nbyt)
                go to 100
17              array(in1:in1+nbyt-1)=b17(lastb+1:lastb+nbyt)
                go to 100
18              array(in1:in1+nbyt-1)=b18(lastb+1:lastb+nbyt)
                go to 100
19              array(in1:in1+nbyt-1)=b19(lastb+1:lastb+nbyt)
                go to 100
20              array(in1:in1+nbyt-1)=b20(lastb+1:lastb+nbyt)

100         bytnum(lbuff)=bytnum(lbuff)+nbyt
            in1=in1+nbyt
            nleft=nleft-nbyt
        end if

C       process more bytes, if any
        if (nleft .gt. 0)then
          nrec=recnum(pbuff)+1

150       continue

          if (nleft .gt. buflen)then
C           read whole blocks directly from the FITS file by-passing buffers

C           test if desired record exists before trying to read it
            if (nrec + nleft/buflen - 1 .gt. maxrec(lbuff)) then
C               record doesn't exist, so return EOF error
                status=107
                return
            end if

C           check if record is already loaded in one of the physical buffers
C           must read it from buffer since it may have been modified
            do 120 i=1,maxbuf
               if (logbuf(i) .eq. lbuff .and. recnum(i) .eq. nrec)then
C                 found the desired record; don't have to read it
                  go to 170
               end if
120         continue

C           record not already loaded in buffer, so read it from disk
            read(iunit,rec=nrec,iostat=ios)array(in1:in1+buflen-1)

            if (ios .ne. 0)then
C               assume that this error indicates an end of file condition
                status=107
                return
            end if

            bytnum(lbuff)=bytnum(lbuff)+buflen
            in1=in1+buflen
            nleft=nleft-buflen
            nrec=nrec+1
            go to 150
          end if

C         load the next file record into a physical buffer
170       call ftldrc(iunit,nrec,.false.,status)
          if (status .gt. 0)return
          go to 200
        end if
        end
C--------------------------------------------------------------------------
        subroutine ftgcbo(iunit,gsize,ngroup,offset,array,status)

C       "Get Character BuFfer with Offsets"
C       read characters from the character buffer.

C       iunit   i  Fortran output unit number
C       gsize   i  size of each group of bytes
C       ngroup  i  number of groups to read
C       offset  i  size of gap between groups
C       array   c  output character string
C       status  i  output error status (0 = ok)
C
C       written by Wm Pence, HEASARC/GSFC, Dec 1996

        integer iunit,gsize,ngroup,offset,status
        character*(*) array

C       COMMON BLOCK DEFINITIONS:--------------------------------------------
        integer nb,ne,pb
        parameter (nb = 20)
        parameter (ne = 512)
        parameter (pb = 20)

        integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
        integer nxtfld
        logical wrmode
        common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb),
     &  wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld

        integer buflun,currnt,reclen,bytnum,maxrec
        common/ftlbuf/buflun(nb),currnt(nb),reclen(nb),
     &  bytnum(nb),maxrec(nb)

        integer maxbuf,logbuf,recnum,pindex
        logical modify
        common/ftpbuf/maxbuf,logbuf(pb),recnum(pb),modify(pb),
     &  pindex(pb)

C       have to use separate character arrays because of compiler limitations
        character*2880 b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,
     &  b15,b16,b17,b18,b19,b20
        common /ftbuff/b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,
     &  b15,b16,b17,b18,b19,b20
C       END OF COMMON BLOCK DEFINITIONS-----------------------------------

        integer lbuff,pbuff,buflen,lastb,nleft,in1,nbyt
        integer i,bytno,record,oldrec,incre

        if (status .gt. 0)return

        lbuff =bufnum(iunit)
        buflen=reclen(lbuff)
        pbuff =currnt(lbuff)
        oldrec=recnum(pbuff)
C       lastb = position of last byte read from input buffer
        lastb =bytnum(lbuff)
        bytno =(oldrec-1) * buflen + lastb
C       in1   = position of first byte remaining in the input buffer
        in1   =1
        nbyt  =0
        incre =gsize+offset

        do 500 i=1,ngroup

C           nleft   = number of bytes left in the input buffer
            nleft=gsize
C           nbyt    = number of bytes to transfer from input to output
            nbyt=min(nleft,buflen-lastb)
            if (nbyt .eq. 0)go to 300

200         continue
C           get characters from the physical buffer to the output string
            go to (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,
     &         19,20)pbuff

C               if got here, then pbuff is out of range
                status=101
                return

1               array(in1:in1+nbyt-1)=b1(lastb+1:lastb+nbyt)
                go to 100
2               array(in1:in1+nbyt-1)=b2(lastb+1:lastb+nbyt)
                go to 100
3               array(in1:in1+nbyt-1)=b3(lastb+1:lastb+nbyt)
                go to 100
4               array(in1:in1+nbyt-1)=b4(lastb+1:lastb+nbyt)
                go to 100
5               array(in1:in1+nbyt-1)=b5(lastb+1:lastb+nbyt)
                go to 100
6               array(in1:in1+nbyt-1)=b6(lastb+1:lastb+nbyt)
                go to 100
7               array(in1:in1+nbyt-1)=b7(lastb+1:lastb+nbyt)
                go to 100
8               array(in1:in1+nbyt-1)=b8(lastb+1:lastb+nbyt)
                go to 100
9               array(in1:in1+nbyt-1)=b9(lastb+1:lastb+nbyt)
                go to 100
10              array(in1:in1+nbyt-1)=b10(lastb+1:lastb+nbyt)
                go to 100
11              array(in1:in1+nbyt-1)=b11(lastb+1:lastb+nbyt)
                go to 100
12              array(in1:in1+nbyt-1)=b12(lastb+1:lastb+nbyt)
                go to 100
13              array(in1:in1+nbyt-1)=b13(lastb+1:lastb+nbyt)
                go to 100
14              array(in1:in1+nbyt-1)=b14(lastb+1:lastb+nbyt)
                go to 100
15              array(in1:in1+nbyt-1)=b15(lastb+1:lastb+nbyt)
                go to 100
16              array(in1:in1+nbyt-1)=b16(lastb+1:lastb+nbyt)
                go to 100
17              array(in1:in1+nbyt-1)=b17(lastb+1:lastb+nbyt)
                go to 100
18              array(in1:in1+nbyt-1)=b18(lastb+1:lastb+nbyt)
                go to 100
19              array(in1:in1+nbyt-1)=b19(lastb+1:lastb+nbyt)
                go to 100
20              array(in1:in1+nbyt-1)=b20(lastb+1:lastb+nbyt)

100         in1=in1+nbyt
            nleft=nleft-nbyt

C           process more bytes, if any
300         continue
            if (nleft .gt. 0)then
C               load the next file record into a physical buffer
                oldrec=oldrec+1
                call ftldrc(iunit,oldrec,.false.,status)
                if (status .gt. 0)return
                pbuff=currnt(lbuff)
                lastb=0
                nbyt=nleft
                go to 200
            end if

            if (i .ne. ngroup)then
C               move to the position of the next group
                bytno=bytno+incre
                record=bytno/buflen+1
                lastb=mod(bytno,buflen)

                if (record .ne. oldrec)then
C                   not the current record, so load the new record; 
                    call ftldrc(iunit,record,.false.,status)
                    if (status .gt. 0)return
                    oldrec=record
                    pbuff=currnt(lbuff)
                end if
            end if
500     continue

        bytnum(lbuff)=lastb+nbyt
        end
C----------------------------------------------------------------------
        subroutine ftgrsz(iunit,nrows,status)

C       Returns an optimal value for the number of rows that should be 
C       read or written at one time in a binary table for maximum efficiency.
C       Accessing more rows than this may cause excessive flushing and 
C       rereading of buffers to/from disk.

C       iunit   i  fortran unit number
C       nrows   i  optimal number of rows to access 
C       status  i  output error status
C
C       written by Wm Pence, HEASARC/GSFC, December, 1996

        integer iunit,nrows,status

C       COMMON BLOCK DEFINITIONS:--------------------------------------------
        integer nb,nf,ne,pb
        parameter (nb = 20)
        parameter (nf = 3000)
        parameter (ne = 512)
        parameter (pb = 20)

        integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
        integer nxtfld
        logical wrmode
        common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb),
     &  wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
        integer tfield,tstart,tbcol,rowlen,tdtype,trept,tnull,heapsz
        integer theap
        double precision tscale,tzero
        common/ft0002/tfield(nb),tstart(nb),tbcol(nf),rowlen(nb),
     &  tdtype(nf),trept(nf),tscale(nf),tzero(nf),tnull(nf),heapsz(nb)
     &  ,theap(nb)

        integer buflun,currnt,reclen,bytnum,maxrec
        common/ftlbuf/buflun(nb),currnt(nb),reclen(nb),
     &  bytnum(nb),maxrec(nb)

C       END OF COMMON BLOCK DEFINITIONS-----------------------------------

        integer ibuff, ii, jj, unique, nfiles

C       There are pb internal buffers available each reclen(nb) bytes long
        ibuff=bufnum(iunit)

C       if HDU structure is not defined then scan the header keywords
        if (dtstrt(ibuff) .lt. 0)call ftrdef(iunit,status)
        if (status .gt. 0)return

C       determine how many different FITS files are currently open
        nfiles = 0
        do 20 ii = 1,nb
           unique = 1
           do 10 jj = 1, ii-1
             if (buflun(ii) .le. 0 .or. buflun(ii) .eq. buflun(jj))then
               unique = 0
               go to 15
             end if
10         continue
15         continue

           if (unique .eq. 1)nfiles=nfiles+1
20      continue

C       one buffer (at least) is always allocated to each open file.
C       assume record size is 2880 bytes (not necessarily true on Vax)

        nrows = ((pb - nfiles) * 2880) / max(1,rowlen(ibuff))
        nrows = max(1, nrows)
        end
C------------------------------------------------------------------------
        subroutine ftflsh(lbuff,status)

C       flush any modified buffers associated with lbuff to disk.
C       Make the contents of the buffers undefined.
        
C       lbuff   i  logical buffer assocaiated with this file
C       status  i  output error status

        integer lbuff,status

C       COMMON BLOCK DEFINITIONS:--------------------------------------------
        integer nb,pb
        parameter (nb = 20)
        parameter (pb = 20)

        integer buflun,currnt,reclen,bytnum,maxrec
        common/ftlbuf/buflun(nb),currnt(nb),reclen(nb),
     &  bytnum(nb),maxrec(nb)

        integer maxbuf,logbuf,recnum,pindex
        logical modify
        common/ftpbuf/maxbuf,logbuf(pb),recnum(pb),modify(pb),
     &  pindex(pb)
C       END OF COMMON BLOCK DEFINITIONS-----------------------------------

        integer ounit,rlen,i

C       ignore input status and flush buffers regardless of status value

        ounit=buflun(lbuff)
        rlen=reclen(lbuff)

C       find any buffer associated with this file
        do 10 i=1,maxbuf
            if (logbuf(i) .eq. lbuff)then
                if (modify(i))then
C                   write the modified buffer to disk
                    call ftwrit(ounit,recnum(i),rlen,i,status)
                    modify(i)=.false.
                end if

C               erase the association of this buffer with the file
                logbuf(i)=0
                recnum(i)=0
            end if
10      continue
        end            
C----------------------------------------------------------------------
        subroutine ftldrc(iunit,nrec,igneof,status)

C       low-level routine to load a specified record from a file into
C       a physical buffer, if it is not already loaded.  Reset all
C       pointers to make this the new current record for that file.
C       Update ages of all the physical buffers.

C       iunit   i  fortran unit number
C       nrec    i  direct access file record number to be loaded
C       igneof  l  ignore end of file error (107)?
C       status  i  output error status

        integer iunit,nrec,status
        logical igneof

C       COMMON BLOCK DEFINITIONS:--------------------------------------------
        integer nb,ne,pb
        parameter (nb = 20)
        parameter (ne = 512)
        parameter (pb = 20)

        integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
        integer nxtfld
        logical wrmode
        common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb),
     &  wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld
        integer buflun,currnt,reclen,bytnum,maxrec
        common/ftlbuf/buflun(nb),currnt(nb),reclen(nb),
     &  bytnum(nb),maxrec(nb)
        integer maxbuf,logbuf,recnum,pindex
        logical modify
        common/ftpbuf/maxbuf,logbuf(pb),recnum(pb),modify(pb),
     &  pindex(pb)
        integer compid
        common/ftcpid/compid
C       END OF COMMON BLOCK DEFINITIONS-----------------------------------

        integer i,lbuff,pbuff,ounit,olen,orec,a1,tstat

        if (status .gt. 0)return

        lbuff=bufnum(iunit)

C       check if record is already loaded in one of the physical buffers
        do 10 i=1,maxbuf
            if (logbuf(i) .eq. lbuff .and. recnum(i) .eq. nrec)then
C               found the desired record; don't have to read it
                pbuff=i
                go to 20
            end if
10      continue

C       the record is not already loaded, so we have to read it from disk.
C       First, decide which physical buffer into which to read it.
        call ftwhbf(lbuff,pbuff)

        if (modify(pbuff))then
C           old buffer has been modified, so we have to flush it to disk
            ounit=buflun(logbuf(pbuff))
            olen=reclen(logbuf(pbuff))
            orec=recnum(pbuff)
            call ftwrit(ounit,orec,olen,pbuff,status)
            modify(pbuff)=.false.
        end if

C       now read the record into the physical buffer
        olen=reclen(lbuff)
        tstat=0
        call ftread(iunit,nrec,olen,pbuff,tstat)

        if (.not. igneof .and. tstat .eq. 107)then
C           return if hit EOF and told not to ignore it
            status=107
            return
        else if (tstat .eq. 107)then
C           apparently hit end of file

            if (.not. wrmode(lbuff))then
C               just return if we don't have write access to the file
                return
            else 
C               fill the new buffer with the desired value
                if (hdutyp(lbuff) .eq. 1)then
C                   ASCII table: fill buffer with blanks
                    call ftflbl(pbuff)
                else if (compid .ge. -1)then
C                   initialize buffer = 0 (except on Cray machines)
                    call ftflzr(pbuff)
                else
C                   call special routine for Cray machines, since words
C                   are twice as long (integers are 8-bytes long)
                    call ftzrcr(pbuff)
                end if

C               mark the new record as having been modified
                modify(pbuff)=.true.
            end if
        end if

C       define log. buffer and the record number contained in the phys. buffer
        logbuf(pbuff)=lbuff
        recnum(pbuff)=nrec

20      continue
C       this is now the current buffer for this logical buffer
        currnt(lbuff)=pbuff
        bytnum(lbuff)=0

C       find the current position of the buffer in the age index
        do 30 i=1,maxbuf
            if (pindex(i) .eq. pbuff)then
               a1=i
               go to 35
            end if
30      continue

35      continue
C       rebuild the indices giving the chronological ordering of the buffers
        do 40 i=a1,maxbuf-1
                pindex(i)=pindex(i+1)
40      continue
C       this buffer is now the youngest (= last in the index)
        pindex(maxbuf)=pbuff
        end
C------------------------------------------------------------------------
        subroutine ftflzr(pbuff)

C       initalize the common block buffer as efficiently as possible
C       with zeros.  This routine should not be used on Cray computers.

C       pbuff  i  number of the physical buffer to initialize

        integer pbuff,i
        integer pb
        parameter (pb = 20)
        double precision buff
        common /ftbuff/buff(360,pb)

        do 10 i=1,360
                buff(i,pbuff)=0.
10      continue
        end
C------------------------------------------------------------------------
        subroutine ftzrcr(pbuff)

C       initalize the common block buffer with zeros. 
C       This routine is reserved for Cray computers.

C       pbuff  i  number of the physical buffer to initialize

        integer pbuff,i
        integer pb
        parameter (pb = 20)
        integer buff,dummy
        common /ftbuff/buff(360,pb),dummy(360,pb)

C       The dummy array was added to the common block to eliminate
C       compiler warnings on other platforms besides Cray.  On all
C       other machines, an integer is 4 bytes long, hence the dummy
C       array is needed to increase the size of the common block to
C       the same length as other routines where the array is declared
C       as character*2880 buff(pb).

C       This will cause compiler warnings on the CRAY, but these
C       can be ignored since the dummy array is never used.

        do 10 i=1,360
                buff(i,pbuff)=0
10      continue
        end
C------------------------------------------------------------------------
        subroutine ftflbl(pbuff)

C       initalize the common block buffer as efficiently as possible
C       with blanks (ASCII 32) 

C       pbuff  i  number of the physical buffer to initialize

        integer pbuff
        integer pb
        parameter (pb = 20)
        character*2880 cbuff
        common /ftbuff/cbuff(pb)

        cbuff(pbuff) = ' '

        end
C----------------------------------------------------------------------
        subroutine ftwhbf(ilbuff,pbuff)

C       decide which physical buffer to use to load in a new record

C       ilbuff  i  the logical buffer number of the record to be loaded (input)
C       pbuff   i  the physical buffer that should be used (output)

        integer ilbuff,pbuff

C       COMMON BLOCK DEFINITIONS:--------------------------------------------
        integer nb,pb
        parameter (nb = 20)
        parameter (pb = 20)

        integer buflun,currnt,reclen,bytnum,maxrec
        common/ftlbuf/buflun(nb),currnt(nb),reclen(nb),
     &  bytnum(nb),maxrec(nb)

        integer maxbuf,logbuf,recnum,pindex
        logical modify
        common/ftpbuf/maxbuf,logbuf(pb),recnum(pb),modify(pb),
     &  pindex(pb)
C       END OF COMMON BLOCK DEFINITIONS-----------------------------------

        integer i,num,lbuff

C       search from the oldest to the youngest for an unlocked record
        do 10 i=1,maxbuf
C           num is the number of the next oldest buffer
            num=pindex(i)

C           lbuff is the logical buffer associated with this physical buffer
            lbuff=logbuf(num)

C           is the physical buffer a current buffer (i.e., locked)?
            if (lbuff .eq. 0 .or. currnt(lbuff) .ne. num)then
                pbuff=num
                return
            end if
10      continue

C       all the buffers are locked, so we have to reuse the current one
        pbuff=currnt(ilbuff)
        end
C----------------------------------------------------------------------
        subroutine ftwrit(ounit,nrec,length,pbuff,status)

C       lowest-level routine to write a physical buffer to the disk file

C       ounit   i  Fortran unit number to write to 
C       nrec    i  number of the file record to write
C       length  i  number of bytes to write
C       pbuff   i  number of the physical buffer to write from
C       status  i  output error status

        integer ounit,nrec,length,pbuff,status

C       COMMON BLOCK DEFINITIONS:--------------------------------------------
        integer nb,ne
        parameter (nb = 20)
        parameter (ne = 512)

        integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
        integer nxtfld
        logical wrmode
        common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb),
     &  wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld

        integer buflun,currnt,reclen,bytnum,maxrec
        common/ftlbuf/buflun(nb),currnt(nb),reclen(nb),
     &  bytnum(nb),maxrec(nb)

C       have to use separate character arrays because of compiler limitations
        character*2880 b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,
     &  b15,b16,b17,b18,b19,b20
        common /ftbuff/b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,
     &  b15,b16,b17,b18,b19,b20
C       END OF COMMON BLOCK DEFINITIONS-----------------------------------

        integer ibuff

C       Note: performance testing on a SUN showed that writing a
C       c*2888 string was MUCH (11x) faster than writing a C*1(2880) array
C       with write(...)(b1(i),i=1,2880).  It was also 2-3 times faster
C       than if the array was declared as a double and written with
C       write(...)(darray(i),i=1,360).  The VAX took about the same
C       time for all 3 different ways to write the bytes.

        ibuff=bufnum(ounit)
C       Save maximum record written, for comparison in ftread
        maxrec(ibuff) = max(maxrec(ibuff), nrec)

        go to (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,
     &  19,20)pbuff

C       if got here, then pbuff is out of range
        status=101
        return

1       write(ounit,rec=nrec,err=900)b1(1:length)
        return
2       write(ounit,rec=nrec,err=900)b2(1:length)
        return
3       write(ounit,rec=nrec,err=900)b3(1:length)
        return
4       write(ounit,rec=nrec,err=900)b4(1:length)
        return
5       write(ounit,rec=nrec,err=900)b5(1:length)
        return
6       write(ounit,rec=nrec,err=900)b6(1:length)
        return
7       write(ounit,rec=nrec,err=900)b7(1:length)
        return
8       write(ounit,rec=nrec,err=900)b8(1:length)
        return
9       write(ounit,rec=nrec,err=900)b9(1:length)
        return
10      write(ounit,rec=nrec,err=900)b10(1:length)
        return
11      write(ounit,rec=nrec,err=900)b11(1:length)
        return
12      write(ounit,rec=nrec,err=900)b12(1:length)
        return
13      write(ounit,rec=nrec,err=900)b13(1:length)
        return
14      write(ounit,rec=nrec,err=900)b14(1:length)
        return
15      write(ounit,rec=nrec,err=900)b15(1:length)
        return
16      write(ounit,rec=nrec,err=900)b16(1:length)
        return
17      write(ounit,rec=nrec,err=900)b17(1:length)
        return
18      write(ounit,rec=nrec,err=900)b18(1:length)
        return
19      write(ounit,rec=nrec,err=900)b19(1:length)
        return
20      write(ounit,rec=nrec,err=900)b20(1:length)
        return

C       come here if there was a disk write error of some sort
900     status=106
        end
C----------------------------------------------------------------------
        subroutine ftread(iunit,nrec,length,pbuff,status)

C       lowest-level routine to read a disk file record into a physical buffer

C       iunit   i  Fortran unit number to read from
C       nrec    i  number of the file record to read
C       length  i  number of bytes to read
C       pbuff   i  number of the physical buffer to read into
C       status  i  output error status

        integer iunit,nrec,length,pbuff,status

C       COMMON BLOCK DEFINITIONS:--------------------------------------------
        integer nb,ne
        parameter (nb = 20)
        parameter (ne = 512)

        integer bufnum,chdu,hdutyp,maxhdu,hdstrt,hdend,nxthdr,dtstrt
        integer nxtfld
        logical wrmode
        common/ft0001/bufnum(199),chdu(nb),hdutyp(nb),maxhdu(nb),
     &  wrmode(nb),hdstrt(nb,ne),hdend(nb),nxthdr(nb),dtstrt(nb),nxtfld

        integer buflun,currnt,reclen,bytnum,maxrec
        common/ftlbuf/buflun(nb),currnt(nb),reclen(nb),
     &  bytnum(nb),maxrec(nb)

C       have to use separate character arrays because of compiler limitations
        character*2880 b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,
     &  b15,b16,b17,b18,b19,b20
        common /ftbuff/b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,
     &  b15,b16,b17,b18,b19,b20
C       END OF COMMON BLOCK DEFINITIONS-----------------------------------

        integer ibuff,ios

C       test if desired record exists before trying to read it
        ibuff=bufnum(iunit)
        if (nrec .gt. maxrec(ibuff)) then
C             record doesn't exist, so return EOF error
              status=107
              return
        end if

        go to (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,
     &  19,20)pbuff

C       if got here, then pbuff is out of range
        status=101
        return

1       read(iunit,rec=nrec,iostat=ios)b1(1:length)
        go to 100
2       read(iunit,rec=nrec,iostat=ios)b2(1:length)
        go to 100
3       read(iunit,rec=nrec,iostat=ios)b3(1:length)
        go to 100
4       read(iunit,rec=nrec,iostat=ios)b4(1:length)
        go to 100
5       read(iunit,rec=nrec,iostat=ios)b5(1:length)
        go to 100
6       read(iunit,rec=nrec,iostat=ios)b6(1:length)
        go to 100
7       read(iunit,rec=nrec,iostat=ios)b7(1:length)
        go to 100
8       read(iunit,rec=nrec,iostat=ios)b8(1:length)
        go to 100
9       read(iunit,rec=nrec,iostat=ios)b9(1:length)
        go to 100
10      read(iunit,rec=nrec,iostat=ios)b10(1:length)
        go to 100
11      read(iunit,rec=nrec,iostat=ios)b11(1:length)
        go to 100
12      read(iunit,rec=nrec,iostat=ios)b12(1:length)
        go to 100
13      read(iunit,rec=nrec,iostat=ios)b13(1:length)
        go to 100
14      read(iunit,rec=nrec,iostat=ios)b14(1:length)
        go to 100
15      read(iunit,rec=nrec,iostat=ios)b15(1:length)
        go to 100
16      read(iunit,rec=nrec,iostat=ios)b16(1:length)
        go to 100
17      read(iunit,rec=nrec,iostat=ios)b17(1:length)
        go to 100
18      read(iunit,rec=nrec,iostat=ios)b18(1:length)
        go to 100
19      read(iunit,rec=nrec,iostat=ios)b19(1:length)
        go to 100
20      read(iunit,rec=nrec,iostat=ios)b20(1:length)

100     continue
        if (ios .ne. 0)then
C               assume that this error indicates an end of file condition
                status=107
        end if
        end
C----------------------------------------------------------------------
        subroutine ftpi2b(ounit,nvals,incre,i2vals,status)

C       Write an array of Integer*2 bytes to the output FITS file.
C       Does any required translation from internal machine format to FITS.

        integer nvals,incre,ounit,status,offset
        integer*2 i2vals(nvals)

C       ounit   i  fortran unit number
C       nvals   i  number of pixels in the i2vals array
C       incre   i  byte increment between values
C       i2vals  i*2 array of input integer*2 values
C       status  i  output error status

        integer*2 temp(4)
        integer ierr,cray2ieg,remain
        integer compid
        common/ftcpid/compid
        character ctemp*1

        if (compid .eq. 0)then
C           big endian machine (e.g., SUN) doesn't need byte swapping
        else if (compid .eq. -1)then
C           SUN F90 compiler maps I*2 -> I*4; have to pack bytes
            call ftpki2(i2vals,nvals,ctemp)
        else if (compid .ge. 1)then
C           little endian machine (e.g. DEC, VAX, or PC) must be byte swapped
            call ftswby(i2vals,nvals)
        else
C           must be a CRAY
C           convert from cray I*8 to IEEE I*2
C           there is a bug in cray2ieg if the number of values to convert
C           is  1 less than a  multiple of 4 2-byte words.  (3, 7, 11, etc)
            remain=nvals-nvals/4*4
            if (remain .eq. 3)then
              if (nvals .gt. 3)then
                ierr= cray2ieg(7,nvals-3,i2vals,0,i2vals,1,' ')
              end if
              temp(3)=i2vals(nvals)
              temp(2)=i2vals(nvals-1)
              temp(1)=i2vals(nvals-2)
              ierr=cray2ieg(7,4,i2vals(nvals/4+1),0,temp,1,' ')
            else
              ierr=cray2ieg(7,nvals,i2vals,0,i2vals,1,' ')
            end if
        end if

        if (incre .le. 2)then
                call ftpbyt(ounit,nvals*2,i2vals,status)
        else
C               offset is the number of bytes to move between each value
                offset=incre-2
                call ftpbyo(ounit,2,nvals,offset,i2vals,status)
        end if
        end
C----------------------------------------------------------------------
        subroutine ftpi4b(ounit,nvals,incre,i4vals,status)

C       Write an array of Integer*4 bytes to the output FITS file.
C       Does any required translation from internal machine format to FITS.

        integer nvals,incre,ounit,status,offset
        integer i4vals(nvals)

C       ounit   i  fortran unit number
C       nvals   i  number of pixels in the i4vals array
C       incre   i  byte increment between values
C       i4vals  i  array of input integer*4 values
C       status  i  output error status

        integer compid
        common/ftcpid/compid

        integer cray2ieg,neven,ierr

        if (compid .eq. 0 .or. compid .eq. -1)then
C           big endian machine (e.g., SUN) doesn't need byte swapping
        else if (compid .ge. 1)then
C           little endian machine (e.g. DEC, VAX, or PC) must be byte swapped
            call ftswi4(i4vals,nvals)
        else
C         must be a CRAY
C         there is a bug in cray2ieg if the number of values to convert
C         is not a multiple of 8 bytes.  
          neven=nvals/2*2
          if (neven .gt. 0)then
              ierr= cray2ieg(1,neven,i4vals,0,i4vals,1,' ')
          end if

          if (neven .ne. nvals)then
C           have to do the remaining odd word separately
            ierr= cray2ieg(1,1,i4vals(nvals/2+1),0,i4vals(nvals),1,' ')
          end if
        end if

        if (incre .le. 4)then
                call ftpbyt(ounit,nvals*4,i4vals,status)
        else
C               offset is the number of bytes to move between each value
                offset=incre-4
                call ftpbyo(ounit,4,nvals,offset,i4vals,status)
        end if
        end
C----------------------------------------------------------------------
        subroutine ftpr4b(ounit,nvals,incre,r4vals,status)

C       Write an array of Real*4 bytes to the output FITS file.
C       Does any required translation from internal machine format to FITS.

        integer nvals,incre,ounit,status,offset
        real r4vals(nvals)

C       ounit   i  fortran unit number
C       nvals   i  number of pixels in the r4vals array
C       incre   i  byte increment between values
C       r4vals  r  array of input real*4 values
C       status  i  output error status

        integer compid
        common/ftcpid/compid

        integer i,neven,ierr,cray2ieg

        if (compid .eq. 0 .or. compid .eq. -1)then
C           big endian machine (e.g., SUN) doesn't need byte swapping
        else if (compid .eq. 1)then
C           little endian machine (e.g. DEC or PC) must be byte swapped
            call ftswi4(r4vals,nvals)
        else if (compid .ge. 2)then
C           convert from VAX format to IEEE
            do 5 i=1,nvals
                    r4vals(i)=r4vals(i)*0.25
5           continue
            call ftswby(r4vals,nvals*2)
        else
C         must be a CRAY
C         there is a bug in cray2ieg if the number of values to convert
C         is not a multiple of 8 bytes.
          neven=nvals/2*2  
          ierr= cray2ieg(2,neven,r4vals,0,r4vals,1,' ')
          if (neven .ne. nvals)then
C           have to do the remaining odd word separately
            ierr= cray2ieg(2,1,r4vals(nvals/2+1),0,r4vals(nvals),1,' ')
          end if
        end if

        if (incre .le. 4)then
                call ftpbyt(ounit,nvals*4,r4vals,status)
        else
C               offset is the number of bytes to move between each value
                offset=incre-4
                call ftpbyo(ounit,4,nvals,offset,r4vals,status)
        end if
        end
C----------------------------------------------------------------------
        subroutine ftpr8b(ounit,nvals,incre,r8vals,status)

C       Write an array of Real*8 bytes to the output FITS file.
C       Does any required translation from internal machine format to FITS.

        integer nvals,incre,ounit,status,offset
        double precision r8vals(nvals)

C       ounit   i  fortran unit number
C       nvals   i  number of pixels in the r4vals array
C       incre   i  byte increment between values
C       r8vals  d  array of input real*8 values
C       status  i  output error status

        integer compid
        common/ftcpid/compid

        integer i,ierr,cray2ieg

        if (compid .eq. 0 .or. compid .eq. -1)then
C           big endian machine (e.g., SUN) doesn't need byte swapping
        else if (compid .eq. 1)then
C           little endian machine (e.g. DEC or PC) must be byte swapped
            call ftswi8(r8vals,nvals)
        else if (compid .eq. 2)then
C           convert from VAX format to IEEE
            call ieevpd(r8vals,r8vals,nvals)
        else if (compid .eq. 3)then
C           convert from Alpha VMS format to IEEE
            do 5 i=1,nvals
                    r8vals(i)=r8vals(i)*0.25
5           continue
            call ftswby(r8vals,nvals*4)
        else
C           must be a CRAY
            ierr= cray2ieg(3,nvals,r8vals,0,r8vals,1,' ')
        end if

        if (incre .le. 8)then
            call ftpbyt(ounit,nvals*8,r8vals,status)
        else
C           offset is the number of bytes to move between each value
            offset=incre-8
            call ftpbyo(ounit,8,nvals,offset,r8vals,status)
        end if
        end
C----------------------------------------------------------------------
        subroutine ftgi2b(iunit,nvals,incre,i2vals,status)

C       Read an array of Integer*2 bytes from the input FITS file.
C       Does any required translation from FITS to internal machine format 

        integer nvals,iunit,incre,status,offset
        integer*2 i2vals(nvals)

C       iunit   i  fortran unit number
C       nvals   i  number of pixels to read
C       incre   i  byte increment between values
C       i2vals  i*2 output array of integer*2 values
C       status  i  output error status

        integer compid
        common/ftcpid/compid

        integer ierr,ieg2cray,i,nloop,fpixel,ntodo
        integer*2 temp(4)
        character ctemp*1

        if (incre .le. 2)then
                call ftgbyt(iunit,nvals*2,i2vals,status)
        else
C               offset is the number of bytes to move between each value
                offset=incre-2
                call ftgbyo(iunit,2,nvals,offset,i2vals,status)
        end if

        if (compid .eq. 0)then
C           big endian machine (e.g., SUN) doesn't need byte swapping
        else if (compid .eq. -1)then
C           SUN F90 compiler maps I*2 -> I*4; have to unpack bytes
            call ftupi2(i2vals,nvals,ctemp)
        else if (compid .ge. 1)then
C           little endian machine (e.g. DEC, VAX, or PC) must be byte swapped
            call ftswby(i2vals,nvals)
        else
C           must be a CRAY
C           convert from IEEE I*2 to cray I*8

C           have to use a temporary array if nvals = 2 or 3
            if (nvals .le. 3)then
              ierr=ieg2cray(7,nvals,i2vals,0,temp,1,' ')
              do 5 i=1,nvals
                  i2vals(i)=temp(i)
5             continue
            else

C             have to work backwards, so as to not overwrite the input data
              nloop=(nvals-1)/4+1
              fpixel = (nloop*4)-3
              ntodo=nvals-(nloop-1)*4
              do 10 i=nloop,1,-1
                ierr=ieg2cray(7,ntodo,i2vals(i),0,i2vals(fpixel),1,' ')
                fpixel=fpixel-4
                ntodo=4
10            continue
            end if
        end if
        end
C----------------------------------------------------------------------
        subroutine ftgi4b(iunit,nvals,incre,i4vals,status)

C       Read an array of Integer*4 bytes from the input FITS file.
C       Does any required translation from FITS to internal machine format 

        integer nvals,iunit,incre,status,offset
        integer i4vals(nvals)

C       iunit   i  fortran unit number
C       nvals   i  number of pixels to read
C       incre   i  byte increment between values
C       i4vals  i  output array of integer values
C       status  i  output error status

        integer ierr,ieg2cray,nloop,fpixel,ntodo,i

        integer compid
        common/ftcpid/compid

        if (incre .le. 4)then
                call ftgbyt(iunit,nvals*4,i4vals,status)
        else
C               offset is the number of bytes to move between each value
                offset=incre-4
                call ftgbyo(iunit,4,nvals,offset,i4vals,status)
        end if

        if (compid .eq. 0 .or. compid .eq. -1)then
C           big endian machine (e.g., SUN) doesn't need byte swapping
        else if (compid .ge. 1)then
C           little endian machine (e.g. DEC, VAX, or PC) must be byte swapped
            call ftswi4(i4vals,nvals)
        else
C           must be a CRAY
C           convert from IEEE I*4 to cray I*8
C           have to work backwards, so as to not overwrite the input data

            nloop=(nvals+1)/2
            fpixel = (nloop*2)-1
            ntodo=nvals-(nloop-1)*2
            do 10 i=nloop,1,-1
                ierr=ieg2cray(1,ntodo,i4vals(i),0,i4vals(fpixel),1,' ')
                fpixel=fpixel-2
                ntodo=2
10          continue
        end if
        end
C----------------------------------------------------------------------
        subroutine ftgr4b(iunit,nvals,incre,r4vals,status)

C       Read an array of Real*4 bytes from the input FITS file.
C       Does any required translation from FITS to internal machine format.

        integer nvals,iunit,incre,status,offset
        real r4vals(nvals)

C       iunit   i  fortran unit number
C       nvals   i  number of pixels to read
C       incre   i  byte increment between values
C       r4vals  r  output array of real*4 values
C       status  i  output error status

        integer ierr,ieg2cray,nloop,fpixel,ntodo,i

        integer compid
        common/ftcpid/compid

        if (incre .le. 4)then
                call ftgbyt(iunit,nvals*4,r4vals,status)
        else
C               offset is the number of bytes to move between each value
                offset=incre-4
                call ftgbyo(iunit,4,nvals,offset,r4vals,status)
        end if

        if (compid .eq. 0 .or. compid .eq. -1)then
C           big endian machine (e.g., SUN) doesn't need byte swapping
        else if (compid .eq. 1)then
C           little endian machine (e.g. DEC or PC) must be byte swapped
            call ftswi4(r4vals,nvals)
        else if (compid .ge. 2)then
C           convert the values from IEEE format to VAX floating point format
C           first, swap the bytes
            call ftswby(r4vals,nvals*2)
C           then test for IEEE special values and multiply value by 4.0
            call ftr4vx(r4vals,r4vals,nvals)
        else
C           must be a CRAY
C           convert from IEEE R*4 to cray R*8
C           have to work backwards, so as to not overwrite the input data

            nloop=(nvals+1)/2
            fpixel = (nloop*2)-1
            ntodo=nvals-(nloop-1)*2
            do 10 i=nloop,1,-1
                ierr=ieg2cray(2,ntodo,r4vals(i),0,r4vals(fpixel),1,' ')
                fpixel=fpixel-2
                ntodo=2
10          continue
        end if
        end
C----------------------------------------------------------------------
        subroutine ftgr8b(iunit,nvals,incre,r8vals,status)

C       Read an array of Real*8 bytes from the input FITS file.
C       Does any required translation from FITS to internal machine format.

        integer nvals,iunit,incre,status,offset
        double precision r8vals(nvals)

C       iunit   i  fortran unit number
C       nvals   i  number of pixels to read
C       incre   i  byte increment between values
C       r8vals  d  output array of real*8 values
C       status  i  output error status

        integer compid
        common/ftcpid/compid

        integer ierr,ieg2cray,nloop,fpixel,ntodo,i

        if (incre .le. 8)then
                call ftgbyt(iunit,nvals*8,r8vals,status)
        else
C               offset is the number of bytes to move between each value
                offset=incre-8
                call ftgbyo(iunit,8,nvals,offset,r8vals,status)
        end if

        if (compid .eq. 0 .or. compid .eq. -1)then
C           big endian machine (e.g., SUN) doesn't need byte swapping
        else if (compid .eq. 1)then
C           little endian machine (e.g. DEC or PC) must be byte swapped
            call ftswi8(r8vals,nvals)
        else if (compid .eq. 2)then
C           convert the values from IEEE format to VAX double (D) format
            call ieevud(r8vals,r8vals,nvals)
        else if (compid .eq. 3)then
C           convert the values from IEEE format to VMS double (G) format
C           first, swap the bytes
            call ftswby(r8vals,nvals*4)
C           then test for IEEE special values and multiply value by 4.0
            call ftr8vx(r8vals,r8vals,r8vals,nvals)
        else
C           must be a CRAY
C           convert from IEEE R*8 to cray R*16
C           have to work backwards, so as to not overwrite the input data

            nloop=(nvals+1)/2
            fpixel = (nloop*2)-1
            ntodo=nvals-(nloop-1)*2
            do 10 i=nloop,1,-1
                ierr=ieg2cray(3,ntodo,r8vals(i),0,r8vals(fpixel),1,' ')
                fpixel=fpixel-2
                ntodo=2
10          continue
        end if
        end
C----------------------------------------------------------------------
        subroutine ftswby(buffer,npix)
C
C       swap pairs of bytes to convert I*2 values to/from IEEE
C
C       buffer  i  input buffer of 16-bit words to be swapped
C       npix    i  number of 16-bit words to be swapped
C
C       written by Wm Pence, HEASARC/GSFC, February 1991

        integer npix,k
        logical*1 buffer(npix*2),temp

        do 10 k=2,npix*2,2
           temp = buffer(k-1)
           buffer(k-1) = buffer(k)
           buffer(k) = temp
10      continue
        end
C----------------------------------------------------------------------
        subroutine ftswi4(buffer,npix)
C
C       swap 4 successive bytes in values to convert to/from IEEE.
C               if the original byte order is: 1 2 3 4
C               the output order is:           4 3 2 1
C
C       buffer  i  input buffer of bytes to be swapped
C       npix    i  number of 4 bytes-words to be swapped
C
C       written by Wm Pence, HEASARC/GSFC, Dec 1996

C  The equivalence statements in this algorithm cause errors with
C  the Linux NAG F90 compiler.
C        integer npix, i
C        integer buffer(npix)
C        integer*4 temp1,temp2
C        logical*1 l1(4), l2(4)
C        equivalence (temp1,l1)
C        equivalence (temp2,l2)
C
C        do 10 i=1,npix
C            temp1=buffer(i)
C            l2(4)=l1(1)
C            l2(3)=l1(2)
C            l2(2)=l1(3)
C            l2(1)=l1(4)
C            buffer(i)=temp2
C10      continue

C       alternate swapping algorithm, Wm Pence, April 97
        integer npix,i4
        logical*1 buffer(npix*4),temp

        do 10 i4=4,npix*4,4

C          swap 1st and 4th bytes
           temp=buffer(i4-3)
           buffer(i4-3)=buffer(i4)
           buffer(i4)=temp
 
C          swap 2nd and 3rd bytes
           temp=buffer(i4-2)
           buffer(i4-2)=buffer(i4-1)
           buffer(i4-1)=temp

10      continue

        end
C----------------------------------------------------------------------
        subroutine ftswi8(buffer,npix)
C
C       swap 8 successive bytes in values to convert to/from IEEE.
C               if the original byte order is: 1 2 3 4 5 6 7 8
C               the output order is:           8 7 6 5 4 3 2 1
C
C       buffer  i  input buffer of bytes to be swapped
C       npix    i  number of 8 bytes-words to be swapped
C
C       written by Wm Pence, HEASARC/GSFC, Dec 1996

C  The equivalence statements in this algorithm cause errors with
C  the Linux NAG F90 compiler.
C        integer npix, i
C        double precision buffer(npix),temp1,temp2
C        logical*1 l1(8), l2(8)
C        equivalence (temp1,l1)
C        equivalence (temp2,l2)
C
C        do 10 i=1,npix
C            temp1=buffer(i)
C            l2(8)=l1(1)
C            l2(7)=l1(2)
C            l2(6)=l1(3)
C            l2(5)=l1(4)
C            l2(4)=l1(5)
C            l2(3)=l1(6)
C            l2(2)=l1(7)
C            l2(1)=l1(8)
C            buffer(i)=temp2
C10      continue

C       alternate swapping algorithm, Wm Pence, April 97
        integer npix,i8
        logical*1 buffer(npix*8),temp

        do 10 i8=8,npix*8,8

C          swap 1st and 8th bytes
           temp=buffer(i8-7)
           buffer(i8-7)=buffer(i8)
           buffer(i8)=temp
 
C          swap 2nd and 7th bytes
           temp=buffer(i8-6)
           buffer(i8-6)=buffer(i8-1)
           buffer(i8-1)=temp

C          swap 3rd and 6th bytes
           temp=buffer(i8-5)
           buffer(i8-5)=buffer(i8-2)
           buffer(i8-2)=temp

C          swap 4th and 5th bytes
           temp=buffer(i8-4)
           buffer(i8-4)=buffer(i8-3)
           buffer(i8-3)=temp
10      continue
        end
C----------------------------------------------------------------------
        subroutine ftr4vx(r4vals,i2vals,nvals)

C       convert IEEE 32-bit floating point numbers to VAX floating point
C       This routine is only called on Vax and Alpha VMS systems.

        real r4vals(*)
        integer*2 i2vals(*)
        integer nvals,i,j

        j=1
        do 10 i=1,nvals
C           test for NaNs (treat +/- infinity the same as a NaN)
            if (i2vals(j) .ge. 32640 .or. (i2vals(j) .lt. 0 .and.
     &           i2vals(j) .ge. -128))then
                 i2vals(j)=-1
                 i2vals(j+1)=-1

C           set underflows and -0 (8000000 hex) = to zero
            else if (i2vals(j) .le. -32641 .or. (i2vals(j) .ge. 0 .and.
     &           i2vals(j) .le. 127))then
                 r4vals(i)=0.0            
            else
C                Must be a real number, so multiply by 4.0 to convert to Vax
                 r4vals(i)=r4vals(i)*4.0
            end if
            j=j+2
10      continue
        end
C----------------------------------------------------------------------
        subroutine ftr8vx(r8vals,i4vals,i2vals,nvals)

C       convert IEEE 32-bit floating point numbers to VAX floating point
C       This routine is only called on VAX computers.

        double precision r8vals(*)
        integer*2 i2vals(*)
        integer i4vals(*)
        integer nvals,i,j,k

        j=1
        k=1
        do 10 i=1,nvals
C           test for NaNs (treat +/- infinity the same as a NaN)
            if (i2vals(j) .ge. 32752 .or. (i2vals(j) .lt. 0 .and.
     &           i2vals(j) .ge. -16))then
                 i4vals(k)  =-1
                 i4vals(k+1)=-1

C           set underflows and -0 (8000000 hex) = to zero
            else if (i2vals(j) .le. -32753 .or. (i2vals(j) .ge. 0 .and.
     &           i2vals(j) .le. 15))then
                 r8vals(i)=0.0
            else
C                Must be a real number, so multiply by 4.0 to convert to Vax
                 r8vals(i)=r8vals(i)*4.0
            end if
            j=j+4
            k=k+2
10      continue
        end
C------------------------------------------------------------------------
        subroutine ftpki2(i2vals,nvals,temp)

C       pack array of 4-byte integers into sequence of 2-byte integers
C       This routine is only currently used on the SUN Solaris F90
C       which does not directly support integer*2 variables and instead
C       maps them into integer*4 variables.

        integer nvals,ii,jj
        integer*2 temp
        character*1 i2vals(nvals*4)

        jj = 2
        do 10 ii = 4,nvals*4,4
             i2vals(jj-1) = i2vals(ii - 1)
             i2vals(jj)   = i2vals(ii)
             jj = jj +2
10      continue

        end
C------------------------------------------------------------------------
        subroutine ftupi2(i2vals,nvals,temp)

C       unpack array of 2-byte integers into 4-integers
C       This routine is only currently used on the SUN Solaris F90
C       which does not directly support integer*2 variables and instead
C       maps them into integer*4 variables.

        integer nvals,ii,jj
        integer*2 temp
        character*1 i2vals(nvals*4),zero,neg,pos

        zero = char(0)
C       neg is use to extend the sign bit if the value is negative
        neg  = char(255)
        pos  = char(127)

        jj=nvals*4
        do 10 ii = nvals*2,2,-2
             i2vals(jj)=i2vals(ii)
             i2vals(jj-1)=i2vals(ii - 1)

C            fill in the 2 most-significant bytes
             if (i2vals(jj-1) .le. pos)then
                 i2vals(jj-2)=zero
                 i2vals(jj-3)=zero
             else
                 i2vals(jj-2)=neg
                 i2vals(jj-3)=neg
             end if
             jj=jj-4
10      continue
        end
