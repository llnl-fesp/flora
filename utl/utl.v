utl # FLORA utility package

{

}

**** UtilityRoutines:

a2i(string:string) integer function # Convert string to integer

pausef(t:integer) subroutine # Interruptable sleep for t seconds

utlinit subroutine # Initialization routine

wmm(a:real,m:integer,n:integer,astr:string,fname:string) subroutine
	# Write matrix a in Matrix Market format for MatView

**** CodeInfo:
logunit         integer            # Unit number for session log file
cmdline0        character*128      # Original command line
codename        character*128      # Name of executable
codedate        character*8        # Date of build
codetime        character*8        # Time of build
rundate         character*8        # Date of execution
runtime         character*8        # Time of execution
runnode         character*8        # Note of execution
probid          character*80 /" "/ # Problem identification string
note            character*80 /" "/ # Problem comment string
code_ts         character*13 /" "/ # Code time stamp
run_ts          character*28       # Run time stamp
session         character*80 /" "/ # Session identification

**** CommandLine:
debug_          integer /0/
echo_           integer /0/
probname_       character*32 /"untitled"/
verbose_        integer /0/
