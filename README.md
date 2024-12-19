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
