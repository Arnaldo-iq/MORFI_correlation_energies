!				PROGRAM CORRECT 
!	AUTHOR: Arnaldo Filho
!	DATE: JULY 2018
! This code corrects the energies obtained by MORFI, as the MP4(SDQ) energies are correct for AA,
! but must be multiplied by a factor of 2 for AB (2) . Note that this is true only for MP3 and MP4(SDQ)
! calculations. For MP2 calculations AA correlatio must be divided by a factor of 2 (1), but no operation
! must be made over AB.
!
! In order to run the program you need three input files:
!     [1] correlation: Contains the uncorrected energies from MORFI
!     [2] lables: Contains the numeric lables for the atoms interacting
!     [3] G09: Contains the MP2, MP3 and MP4 energies, obtanied from Gaussian09 
!
! Make sure to have the module "Parameters.f90" complied along with program "main.f90", the script
! "collect", the output from Gaussian09 and the "test.mout" output from MORFI all in the same folder
!
! FUTURE DEVELOPMENTS: Segregation of the energies computed between bonded and non-bonded terms
! 
PROGRAM MAIN
!
USE Parameters
!
! Variable declaration
!
CHARACTER (len=3) :: mplvl
CHARACTER (len=2) :: version
INTEGER :: natom, i, isize
REAL :: Tot_corr, G09_corr, Err, EMP2, EMP3, EMP4, Tot_int, Tot_self
REAL, DIMENSION(:), ALLOCATABLE :: energy, AA, AB
INTEGER, DIMENSION(:), ALLOCATABLE :: one, two 
TYPE(atom_label), DIMENSION(:), ALLOCATABLE :: inte
!
! Opening units
!
OPEN(UNIT=1, FILE='correlation',STATUS='OLD')
OPEN(UNIT=3, FILE='g09',STATUS='OLD')
OPEN(unit=4, file='output',STATUS='NEW')
OPEN(unit=5, file='numer',STATUS='OLD')
OPEN(unit=6, file='level',STATUS='OLD')
OPEN(unit=7, file='version',STATUS='OLD')
!
! Reading in the information not provided in the input
!
READ(5,*) natom
!
12 READ(6,*) mplvl
!
IF (mplvl.NE."MP2".AND.mplvl.NE."MP3".AND.mplvl.NE."MP4") THEN
PRINT*, "The acceptable values for MPn are MP2, MP3 or MP4, try again"
GO TO 12
END IF
!
! Choosing MPn level
!
	SELECT CASE (mplvl)
		CASE ('MP2')
		READ(3,*) EMP2
		CASE ('MP3')
		READ(3,*) EMP2, EMP3
		CASE ('MP4')
		READ(3,*) EMP2, EMP3, EMP4
		CASE DEFAULT
	END SELECT
!
!PRINT *, "Enter the version of MORFI you are using (AA or AB)"
!
13 READ *, version
!
IF (version.NE."AA".AND.version.NE."AB") THEN
PRINT*, "The acceptable values for version are AA or AB, try again"
GO TO 13
END IF

!
IF (version.EQ."AB") THEN 
!
OPEN(UNIT=2, FILE='labels',STATUS='OLD')
!
!Number of all possible combination of atoms, including self interactions 
!
isize=(natom*(natom+1))/2
!
! Dynamical memory allocation
!
ALLOCATE(energy(isize))
ALLOCATE(one(isize))
ALLOCATE(two(isize))
ALLOCATE(inte(isize))
ALLOCATE(AA(isize))
ALLOCATE(AB(isize))
!
! Reading the uncorrected correlation energies
!
	DO i=1, isize
	READ(1,*) energy(i)
	END DO
!
! Reading atomic lables
!
	DO i=1, isize
	READ(2,*) one(i), two(i)
	END DO
!
	DO i=1, isize
	inte(i)%atom_1 = one(i)
	inte(i)%atom_2 = two(i)
	END DO
!
! Observation (2)
!
	IF (mplvl.EQ."MP2") THEN
		DO i=1, isize
			IF  (inte(i)%atom_1==inte(i)%atom_2) THEN
			inte(i)%corr=(energy(i))*0.5
			ELSE
			inte(i)%corr=(energy(i))
			END IF
		END DO
!
! Observation (1)
!
	ELSE
		DO i=1, isize 
			IF  (inte(i)%atom_1==inte(i)%atom_2) THEN
			inte(i)%corr=energy(i)
			ELSE 
			inte(i)%corr=(energy(i))*2
			END IF
		END DO
	END IF
!
 ELSE
!
isize=natom
!
ALLOCATE(energy(isize))
ALLOCATE(inte(isize))
!
	DO i=1, isize
	READ(1,*) energy(i)
	END DO
!
	DO i=1, isize
	inte(i)%corr=(energy(i))
	END DO
!
END IF
!
	SELECT CASE (mplvl)
		CASE ('MP2')
		G09_corr=EMP2
		CASE ('MP3')
		G09_corr=EMP2+EMP3
		CASE ('MP4')
		G09_corr=EMP2+EMP3+EMP4
	END SELECT
!
! Calculating integration error
!
Total_corr=SUM(inte%corr)
Err=(G09_corr-Total_corr)*2625.5
!
! Writing output
!
WRITE(4,*) "This is a ", mplvl, " calculation. ", "Number of atoms is ", natom 
WRITE(4,*) "The corrected energies are"
!
IF (version.EQ."AB") THEN
!
	DO i=1, isize
		IF  (inte(i)%atom_1==inte(i)%atom_2) THEN
		WRITE(4,*) "AA interaction ", "Atom 1 label = ", inte(i)%atom_1, "Atom 2 label = ", inte(i)%atom_2, "Correlation energy (a.u.)", inte(i)%corr 
		ELSE
		WRITE(4,*) "AB interaction ", "Atom 1 label = ", inte(i)%atom_1, "Atom 2 label = ", inte(i)%atom_2, "Correlation energy (a.u.)", inte(i)%corr
		END IF
	END DO
	
	DO i=1, isize
	IF  (inte(i)%atom_1==inte(i)%atom_2) THEN
		AA(i)=inte(i)%corr
		ELSE
		AB(i)=inte(i)%corr
	END IF
	END DO
!
Tot_self=sum(AA)
Tot_int=sum(AB)
!
WRITE(4,*) "The Total atomic correlation is: ", Tot_self
WRITE(4,*) "The Total bond/int correlation is: ", Tot_int
!
ELSE
	DO i=1, isize
	WRITE(4,*) "The correlation energy (a.u.) for atom ", i, "is", inte(i)%corr
	END DO
!
END IF
!
WRITE(4,*) "MORFI's is energy:(a.u.)", Total_corr, "GAUSSIAN09's is energy (a.u.):", G09_corr, "The integration error is", Err, "kJ.mol-1"
!
END PROGRAM MAIN
