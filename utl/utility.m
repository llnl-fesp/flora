# $Id: utility.m,v 1.5 2004/07/15 16:33:57 bulmer Exp $

# FLORA utility routines

      subroutine usrmain(argv0, cmdline) # Process command-line args

      implicit none

      Use(CodeInfo)
      Use(CommandLine)

      character*(*) argv0, cmdline
      integer i

      # Process FLORA "-" options...
      codename = argv0
      cmdline0 = cmdline
      i = index(cmdline,"-")
      while (i > 0)
          cmdline = cmdline(i+1:128)
          if (cmdline(1:1) == "d") then
              debug_ = YES
              cmdline = cmdline(i+1:128)
          elseif (cmdline(1:1) == "e") then
              echo_ = YES
              cmdline = cmdline(i+1:128)
          elseif (cmdline(1:1) == "h") then
              write(STDOUT,10)
              write(STDOUT,11)
              write(STDOUT,12)
              call exit(0)
          elseif (cmdline(1:1) == "p") then
              i = index(cmdline," ")
              cmdline = cmdline(i+1:128)
              i = index(cmdline," ")
              probname_ = cmdline(1:i-1)
              cmdline = cmdline(i:128)
          elseif (cmdline(1:1) == "v") then
              verbose_ = YES
              i = index(cmdline," ")
              cmdline = cmdline(i+1:128)
          else
              write(STDOUT,13) char(7), trim(adjustl(cmdline))
              call exit(1)
          endif
          i = index(cmdline,"-")
      endwhile

      # Pass any remaining non "-" options to Basis (never returns)...
      call basmain(argv0,cmdline)

   10 format(/,
     ."FLORA 2D equilibrium and stability code ",
     ."for axisymmetric mirror configurations.",//,
     ."Usage:",//,
     ."        flora [-d] [-e] [-h] [-p string] [-v] [Basis command(s)]",//,
     ."where:",/,
     ."        -d : turn on Basis debug mode",/,
     ."        -e : turn on Basis echo mode",/,
     ."        -h : display this message",/,
     ."        -p : set Basis probname to string",/,
     ."        -v : turn on Basis verbose mode",//,
     ."Note: Options cannot be globbed; multiple Basis commands must be delimited",/,
     ."      with a semicolon and the entire command-string enclosed in quotes.")
   11 format(/,
     .'If the version being executed has been "pickled", the single command-line',/,
     .'argument "-M" will display the manifest of source files used to build it:',//,
     .'        florayymmdd -M')
   12 format(/,"References:",//,
     ."   (1) User's Manual for the FLORA Equilibrium and Stability Code",/,
     ."       Robert P. Freis and Bruce I. Cohen",/,
     ."       LLNL UCID-20400, April 1, 1985",/,
     ."   (2) <modified for kinetic stabilizer mirror systems>",/,
     ."       Jack Byers, 2003",/)
   13 format(a1,"Undefined option: `",a,"', do: flora -help")

      end # usrmain

      subroutine utlinit # Initialize options and session log

      implicit none

      Use(CodeInfo)
      Use(CommandLine)

      integer basopen, io
      character*(500) basgetenv
      character*32 logname, sdir
      character*8 runmach

      # Fill CodeInfo variables...
      call codsysinf(runtime,rundate,runmach,runnode,8)
      call glbtmdat(codetime,codedate)
      call basisexe("date > run_ts")
      call freeus(io)
      open(io,FILE="run_ts")
      read(io,10) run_ts
      close(io)
      call basisexe("rm -f run_ts")

      # Install command-line options...
      call parsestr("debug=debug_")
      call parsestr("echo=echo_")
      call parsestr("probname=probname_")
      call parsestr("verbose=verbose_")

      # Build session ID string...
      call parsestr("code_ts = trim(codedate)")
      sdir = basgetenv("FLORA_SCRIPTS")
      if (index(sdir,"include") == 0) then # Not a pickled code
          call parsestr("code_ts = trim(code_ts)//"":""//substr(codetime,1,2)")
          call parsestr("code_ts = trim(code_ts)//substr(codetime,4,2)")
      endif
      call basisexe("basename "//trim(codename)//" > basename")
      call parsestr("$i = basopen(""basename"",""r"")")
      call parsestr("noisy=yes; $i >> codename; noisy=no")
      call parsestr("call basclose($i)")
      call basisexe("rm basename")
      call parsestr("session = trim(codename)//"" [""")
      call parsestr("session = trim(session)//trim(code_ts)//""]""")
      call parsestr("session = trim(session)//"" on ""//trim(runnode)//"" [""")
      call parsestr("session = trim(session)//trim(run_ts)//""]""")
      call parsestr("session = trim(session)//"" ""//probname")

      # Write code info to graphics file...
      call parsestr("call ezcstxqu(1)")
      if (cmdline0 .ne. " ") then
          call parsestr( \
          "stdplot << ""! ""//trim(codename)//"" ""//trim(cmdline0) << return")
      endif
      call parsestr("stdplot << trim(session) << return")
      call parsestr("output graphics")
      call parsestr("stdplot << ""> version""")
      call parsestr("version")
      call parsestr("stdplot << return << ""> list packages""")
      call parsestr("list packages")
      if (index(sdir,"include") .ne. 0) then # include manifest
          call parsestr("stdplot << return << ""> manifest""")
          call parsestr("call manifest(stdplot)")
      endif
      call parsestr("output tty")

      # Create session log file...
      write(logname,11) trim(adjustl(probname_))
      logunit = basopen(logname,"w")
      call baspecho(logunit)
      write(logunit,12) trim(adjustl(codename)), trim(adjustl(cmdline0))

      # Session ID to screen...
      call parsestr("<< return << session << return")

      return

   10 format(a)
   11 format(a,".log")
   12 format("! ",a," ",a)
      end # utlimit

      function a2i(str) # Convert string to integer

      Use(CodeInfo)

      integer a2i
      character*(*) str
      character*80 msg
      msg = "a2i: cannot convert string to integer"
      read(str,*,err=1) a2i
      return
    1 write(STDOUT,10) char(7), msg  
      call baswline(logunit,msg)
      call kaboom(0)
   10 format(a1,a)

      end # a2i

      subroutine pausef(t) # Interruptable sleep

      implicit none
      
      integer t

      call ruthere
      call sleep(t)
      call ruthere

      return
      end # pausef

      subroutine warning(msg) # Send a warning message
      implicit none
      Use(CodeInfo)
      character*(*) msg

      write(STDOUT,'(a1)') char(7)
      call baswline(STDOUT,msg)

      return
      end # warning

      subroutine fatal(msg) # Send a fatal message and exit
      implicit none
      character*(*) msg

      call warning(msg)
      call kaboom(0)

      end # fatal

      subroutine wmm(a,m,n,astr,fname) # Write a in Matrix Marker format
      implicit none

      integer i, j, m, n, nz, io, basopen
      real a(m,n), nonzero
      character*(*) astr, fname
      character*(80) msg

      # Write matrix a in Matrix Market format for use by MatView:
      #   http://www.csm.ornl.gov/~kohl/MatView
      # Arguments:
      #   a(m,n)  :  real 2D matrix
      #   m       :  integer number of rows
      #   n       :  integer number of columns
      #   astr    :  informational string describing matrix a
      #   fname   :  output file name

      # Count non-zero entries
      nz = 0
      do j = 1,n
          do i = 1,m
              if (a(i,j) .ne. 0) then
                  nz = nz + 1
              endif
          enddo
      enddo
      nonzero = 100*float(nz)/float(m*n)
      write(msg,103) astr, nonzero
      call baswline(STDOUT,msg)

      # Write the file in Matrix Market format:
      #   http://math.nist.gov/MatrixMarket/formats.html
      io = basopen(fname,"w")
      write(io,101)
      write(io,102)
      write(io,103) astr, nonzero
      write(io,104) m, n, nz
      do j = 1,n
          do i = 1,m
              if (a(i,j) .ne. 0) then
                  write(io,105) i, j, a(i,j)
              endif
          enddo
      enddo
      call basclose(io)

  101 format("%%MatrixMarket matrix coordinate real general")
  102 format("% Ref: http://math.nist.gov/MatrixMarket/formats.html")
  103 format("% ",a,", ",F4.1,"% non-zero elements")
  104 format(2I8,I14)
  105 format(2I8,1PE18.8)
      end # wmm
