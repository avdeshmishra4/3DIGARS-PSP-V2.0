module myGAfuns_dat

!      INTEGER, PARAMETER:: ChromosomeLength = 66 !33! 30 !27 !33 ! ChromosomeLength = 33 for 3 variable of length 11 bit (55 => 5 x 11)
	INTEGER, PARAMETER:: ChromosomeLength = 55 !33! 30 !27 !33 ! ChromosomeLength = 33 for 3 variable of length 11 bit (55 => 5 x 11)
       INTEGER, PARAMETER:: KGA_K=200, KGA_P = 2
      TYPE genotype    
          REAL*8:: Fitness  ! Now, Fitness (Total Fitness) should be the sum of all the target. Fitness = F_PCC + (-Fit_AvgTMLowEMdl)
          REAL*8:: PCC      ! Pearson's Corelation Cofficient, negative quantity
          REAL*8:: TMLowEMdl! Avg. TMscore of the lowest energy model, positive quantity
		  REAL*8:: Native_Count
		  REAL*8:: ZSCORE
          REAL*8:: PCC1      
          REAL*8:: TMLowEMdl1
		  REAL*8:: Native_Count1
		  REAL*8:: ZSCORE1
          REAL*8:: PCC2      
          REAL*8:: TMLowEMdl2
		  REAL*8:: Native_Count2
		  REAL*8:: ZSCORE2
		  REAL*8:: PCC3      
          REAL*8:: TMLowEMdl3
		  REAL*8:: Native_Count3
		  REAL*8:: ZSCORE3
          Integer, DIMENSION(ChromosomeLength):: Chromo 
        END TYPE genotype
 
  INTEGER, PARAMETER:: Max_generation =20000, PopSize =200   
  REAL*8, PARAMETER:: ElitRate = 0.05, CrossOverRate =0.9, MutationRate = 0.5                      
  REAL*8:: TotalPopFitness, WorstPopFitness
 
  type (genotype), DIMENSION(PopSize+10), target ::popA, popB  
  type (genotype), DIMENSION(ChromosomeLength+10):: popAMBottom, popAMUpper
  type (genotype), DIMENSION(10):: tmpstruct

  type (genotype), pointer :: ptrPop1(:), ptrPop2(:) ! '(:)', since array type pointer  

 Real*8:: Ks(1000)   ! 12/Dec/2012 -> changed from 10 to 100: this was important but often missed.

  
  Integer:: TOGGLE_SWAP_POP=0, generation, CurrentPosNewPopulation, CrossOverCount, myFileCount, TotRec, CurrLoc,myFileCount1, TotRec1, CurrLoc1, myFileCount2, TotRec2, CurrLoc2, myFileCount3, TotRec3, CurrLoc3   

  TYPE genotypeFile  
           Integer:: RecCount
           Integer:: StartIdx
           Integer:: EndIdx  
           Real*8::  PCC          ! Pearson's Corelation Cofficient  
           Real*8::  TMLowEMdl    ! TMscore of the lowest energy model
		   Real*8::	 Z_score
  END TYPE genotypeFile
 
type (genotypeFile), allocatable::FileList(:), FileList1(:),FileList2(:),FileList3(:)


  TYPE genotypeRec  
           Real*8:: TMScore
           Real*8:: ThrDIGARS_Egy
           Real*8:: ASA_Egy
           Real*8:: uPhi_Egy
           Real*8:: uPsi_Egy
           Real*8:: ThrDIGARSThr_Egy
           Real*8:: ASA_Egy_3DIGARS_Data_Ram
           Real*8:: ASA_Egy_SSD_Data_Ram
           Real*8:: ASA_Egy_Spinex_Data_Ram
		   Real*8:: ASA_Egy_3DIGARS_Data_Hoque
		   Real*8:: ASA_Egy_SSD_Data_Hoque
		   Real*8:: ASA_Egy_Spinex_Data_Hoque
		   Real*8:: ASA_SS_Egy_3DIGARS_Data_Ram
		   Real*8:: ASA_SS_Egy_SSD_Data_Ram
		   Real*8:: ASA_SS_Egy_Spinex_Data_Ram
		   Real*8:: ASA_SS_Egy_3DIGARS_Data_Hoque
		   Real*8:: ASA_SS_Egy_SSD_Data_Hoque
		   Real*8:: ASA_SS_Egy_Spinex_Data_Hoque
		   Real*8:: Phi_Egy_3DIGARS_Data_Ram
		   Real*8:: Phi_Egy_SSD_Data_Ram
		   Real*8:: Phi_Egy_Spinex_Data_Ram
		   Real*8:: Phi_Egy_3DIGARS_Data_Hoque
		   Real*8:: Phi_Egy_SSD_Data_Hoque
		   Real*8:: Phi_Egy_Spinex_Data_Hoque
		   Real*8:: Psi_Egy_3DIGARS_Data_Ram
		   Real*8:: Psi_Egy_SSD_Data_Ram
		   Real*8:: Psi_Egy_Spinex_Data_Ram
		   Real*8:: Psi_Egy_3DIGARS_Data_Hoque
		   Real*8:: Psi_Egy_SSD_Data_Hoque
		   Real*8:: Psi_Egy_Spinex_Data_Hoque
		   Real*8:: ASAr_Egy_3DIGARS_Data_Ram
		   Real*8:: ASAr_Egy_SSD_Data_Ram
		   Real*8:: ASAr_Egy_Spinex_Data_Ram
		   Real*8:: ASAr_Egy_3DIGARS_Data_Hoque
		   Real*8:: ASAr_Egy_SSD_Data_Hoque
		   Real*8:: ASAr_Egy_Spinex_Data_Hoque
		   Real*8:: ASAr_SS_Egy_3DIGARS_Data_Ram
		   Real*8:: ASAr_SS_Egy_SSD_Data_Ram
		   Real*8:: ASAr_SS_Egy_Spinex_Data_Ram
		   Real*8:: ASAr_SS_Egy_3DIGARS_Data_Hoque
		   Real*8:: ASAr_SS_Egy_SSD_Data_Hoque
		   Real*8:: ASAr_SS_Egy_Spinex_Data_Hoque
		   Real*8:: Phir_Egy_3DIGARS_Data_Ram
		   Real*8:: Phir_Egy_SSD_Data_Ram
		   Real*8:: Phir_Egy_Spinex_Data_Ram
		   Real*8:: Phir_Egy_3DIGARS_Data_Hoque
		   Real*8:: Phir_Egy_SSD_Data_Hoque
		   Real*8:: Phir_Egy_Spinex_Data_Hoque
		   Real*8:: Psir_Egy_3DIGARS_Data_Ram
		   Real*8:: Psir_Egy_SSD_Data_Ram
		   Real*8:: Psir_Egy_Spinex_Data_Ram
		   Real*8:: Psir_Egy_3DIGARS_Data_Hoque
		   Real*8:: Psir_Egy_SSD_Data_Hoque
		   Real*8:: Psir_Egy_Spinex_Data_Hoque
		   Real*8:: ASA_Egy_Comb_Data_Ram
		   Real*8:: ASA_SS_Egy_Comb_Data_Ram
		   Real*8:: ASA_Egy_Comb_Data_Hoque
		   Real*8:: ASA_SS_Egy_Comb_Data_Hoque
		   Real*8:: Phi_Egy_Comb_Data_Ram
		   Real*8:: Phi_Egy_Comb_Data_Hoque
		   Real*8:: Psi_Egy_Comb_Data_Ram
		   Real*8:: Psi_Egy_Comb_Data_Hoque
		   Real*8:: ASAr_Egy_Comb_Data_Ram
		   Real*8:: ASAr_SS_Egy_Comb_Data_Ram
		   Real*8:: ASAr_Egy_Comb_Data_Hoque
		   Real*8:: ASAr_SS_Egy_Comb_Data_Hoque
		   Real*8:: Phir_Egy_Comb_Data_Ram
		   Real*8:: Phir_Egy_Comb_Data_Hoque
		   Real*8:: Psir_Egy_Comb_Data_Ram
		   Real*8:: Psir_Egy_Comb_Data_Hoque
		   Real*8:: ASA_Egy_Triplet_Hoque
		   Real*8:: Phi_Egy_Triplet_Hoque
		   Real*8:: Psi_Egy_Triplet_Hoque
		   Real*8:: uPhi_Egy_Hoque
		   Real*8:: uPsi_Egy_Hoque
		   Real*8:: o3DIGARS
		   
  END TYPE genotypeRec

type (genotypeRec), allocatable::RecList(:),RecList1(:), RecList2(:), RecList3(:)
type (genotypeRec)::tmpRecList(2)


contains  !! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   NOTE: CONTAINs the following >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!=========================== START of Swap_Population_Pointer ============================================
 ! Swappopulation using pointers, in Fortran pointer does not return the address of the pointer, so had to      
 subroutine Swap_Population_Pointer                        
      IF(TOGGLE_SWAP_POP==0) THEN
          ptrPop1=> popB
          ptrPop2=> popA
          TOGGLE_SWAP_POP=1
          if (myid==0) THEN
              ! write (*,*) '========= Swapped population and popB is now pointed by ptrPop1 ========'  
          END IF
      ELSE 
          TOGGLE_SWAP_POP=0
          ptrPop1=> popA
          ptrPop2=> popB
          if (myid==0) THEN
            !write (*,*) '========= Swapped population and popA is now pointed by ptrPop1 ========'
          END IF
 
      END IF
 end subroutine Swap_Population_Pointer
!=========================== END of Swap_Population_Pointer ============================================


!------------------------------------------FreeUnit ----------
integer function freeunit() result(ip)
   logical:: used
   do ip = 10, 99
      inquire(unit=ip, opened=used)
      if(.not. used) return
   end do
   write(0,*) "No Available Fortran I/O Units between 10 -- 99"
end function
!----------------------------------------- Binary2Decimal --------
Integer FUNCTION Binary2Decimal(i, PopType, St_i, End_i)  
    Integer:: sum, m, n
    Integer, INTENT(IN):: i, PopType, St_i, End_i
    
    sum = 0
    n = 1

 IF (PopType ==1) THEN   
    DO m = End_i, St_i, -1
        IF (ptrPop1(i)%Chromo(m) == 1) THEN
             sum = sum + n
        END IF
        n = n * 2
    END DO 
 END IF   
 


 IF (PopType ==2) THEN   
    DO m = End_i, St_i, -1
        IF (ptrPop2(i)%Chromo(m) == 1) THEN
             sum = sum + n
        END IF
        n = n * 2
    END DO
 END IF

 IF (PopType ==3) THEN   
    DO m = End_i, St_i, -1
        IF (tmpstruct(i)%Chromo(m) == 1) THEN
             sum = sum + n
        END IF
        n = n * 2
    END DO
 END IF

    Binary2Decimal = sum

END FUNCTION Binary2Decimal
!---------------------------------------------------------------------------------
SUBROUTINE ChromoToVarValue(i, PopType) ! this will read the binary chromosome and will return the real values ranges: -5.12 to +5.12
 INTEGER, INTENT(IN):: i, PopType
 INTEGER:: j

  DO j = 1,5
    ! Ks(j) = ((Binary2Decimal(i, PopType, ((j-1) * 10) + 1, ((j-1) * 10) + 10)) - 512.00) / 100.00    !for ranges: -5.12 to +5.12, 10 bit
     Ks(j) = (Binary2Decimal(i, PopType, ((j-1) * 11) + 1, ((j-1) * 11) + 11)) / 1000.00             ! for 0 to 2.048,             11 bit
   ! Ks(j) = (Binary2Decimal(i, PopType, ((j-1) * 9) + 1, ((j-1) * 9) + 9)) / 100.00                    !for ranges: 0 to +5.12,     9 bit
 
  END DO
 
END SUBROUTINE ChromoToVarValue
!-----------------------------------------
!===========================START of RandNumTam =================================
    real*8 FUNCTION RandNumTam(myI)  ! this function extract the milisecond (3 digit) of the current cpu time and square and feed it for the random seed.
      real*8:: myrand                      ! incrementing, i (1,2,3,4) can be used to call sequentially CPUs to allocate random#s,   
      INTEGER,  DIMENSION(1):: RandSeed  ! The function will return a number inbetween 0 to 1. 
      INTEGER :: K, myint
      character(10):: T 
      INTEGER, INTENT(IN) :: myI
     
       CALL RANDOM_SEED(SIZE=K)
       call DATE_AND_TIME(time=T)
       read (T(8:10),'(I10)') myint  ! convert the last 3 digit character into integer
       myint = myint * myI * myI * myI
       RandSeed(1) = myint*myint  
	CALL RANDOM_SEED(PUT=RandSeed(1:K))
       call random_number(myrand) 
      RandNumTam=myrand
    END FUNCTION RandNumTam
!===========================END of RandNumTam =======================================


!===========================START: Random Chromosomes ===============================
SUBROUTINE Random_Chromosome (i, PopType)
  INTEGER::j, i, PopType ! i = location in pop, PopType:: 1 = pointed by ptrPop1 and 2 = pointed by ptrPop2

  Integer:: n  
  real *8:: random1 

  IF (PopType == 1) THEN
       DO j = 1, ChromosomeLength
              !n = (RandNumTam((j+PopType)*i) * 2)
              n= random1() * 2
              !write(*,*) 'n=',n  
              ptrPop1(i)%chromo(j)=n!  (RandNumTam((j+PopType)*i) * 2)
       END DO
      ptrPop1(i)%Fitness = CalculateFitness(i,1) ! for => CalculateFitness(ptrPop1(i)) 
      
  END IF
       
  
IF (PopType == 2) THEN
       DO j = 1, ChromosomeLength
              !n = (RandNumTam((j+PopType)*i) * 2)
              !write(*,*) 'n=',n  
              n= random1() * 2
              ptrPop2(i)%chromo(j)= n! (RandNumTam((j+PopType)*i) * 2)
       END DO
      ptrPop2(i)%Fitness = CalculateFitness(i,2) ! for => CalculateFitness(ptrPop2(i)) 
 
  END IF


END SUBROUTINE random_chromosome
!=========================== END: Random Chromosomes ===============================

!================================= START: initilise the population + AM (Associ. Memory) + compute fitness ===========
SUBROUTINE PopInitialization
    INTEGER:: i, j
     
    ptrPop1=> popA
    ptrPop2=> popB
    TOGGLE_SWAP_POP=0
    DO i = 1, PopSize 
      call Random_Chromosome (i, 1)         
   END DO

   Write(*,*) 'initialization is done successfully ...'

END SUBROUTINE PopInitialization
!================================= END: initilise the population + AM (Associ. Memory) + compute fitness ===========

!================================= START: sort the population pointed by ptrPop1 ===========
SUBROUTINE SortPopulation(PopType)
    INTEGER:: i, j, PopType
    

IF (PopType ==1) THEN 
     DO i = 1, PopSize-1 
       DO j = (i+1), PopSize
         IF ((ptrPop1(i)%Fitness) > (ptrPop1(j)%Fitness)) THEN  ! the more negative pearson_Co_co the better, so should be front of the population
              tmpstruct(1) = ptrPop1(i)
              ptrPop1(i) = ptrPop1(j)
              ptrPop1(j) = tmpstruct(1)
          END IF      
       END DO
    END DO
END IF

IF (PopType ==2) THEN 
     DO i = 1, PopSize-1 
       DO j = (i+1), PopSize
         IF ((ptrPop2(i)%Fitness) > (ptrPop2(j)%Fitness)) THEN  ! the more negative pearson_Co_co the better, so should be front of the population
              tmpstruct(1) = ptrPop2(i)
              ptrPop2(i) = ptrPop2(j)
              ptrPop2(j) = tmpstruct(1)
          END IF      
       END DO
    END DO
END IF


END SUBROUTINE  SortPopulation
!================================= END: sort the population pointed by ptrPop1 ===========


!================================= START: Initilize the AM-Memory POP (Upper, Bottom) with the best chromosome pointed by ptrPop1 ===========
SUBROUTINE Initialize_AM_Pop      ! It should have been upper-traingle and lower triangular, but as it is one time and to make it simple 
                                  ! I have copied the full chromosome structure in the AM_Pop while initilizing them
    INTEGER:: i, j
    
    DO i = 1, ChromosomeLength    ! Note: not the PopSize but the ChromosomeLength
         popAMBottom(i) = ptrPop1(1)
         popAMUpper(i)= ptrPop1(1)   
    END DO

END SUBROUTINE Initialize_AM_Pop
!================================= END: sort population pointed by ptrPop1 ===========

!===========================START: Elite FN ============
Integer FUNCTION Elitist()
 INTEGER:: k, i
 k = ElitRate * PopSize
   DO i = 1,k  
      ptrPop2(i) = ptrPop1(i)
   END DO
  Elitist = k
End Function Elitist
!!===========================END: Elite FN ============

!===========================START: ComputeFitnessSum_WorstFitness ============
subroutine ComputeFitnessSum_WorstFitness ! Fitness could be -ve, 0, +ve and -ve is better
  INTEGER :: I 
  TotalPopFitness = ptrPop1(1)%Fitness
  WorstPopFitness = ptrPop1(1)%Fitness
  
  DO I = 2, PopSize
     TotalPopFitness =TotalPopFitness + ptrPop1(I)%Fitness
     If (WorstPopFitness < ptrPop1(I)%Fitness) THEN
        WorstPopFitness = ptrPop1(I)%Fitness  
     END IF 
  END DO 

 TotalPopFitness = ABS(TotalPopFitness - (WorstPopFitness * PopSize))+ PopSize  ! so the worst fitness is always ABS('-1') (checked), 
                                                                                ! irrespective of +ve or -ve
 
 !The worst one will be (abs(-1)) or finally a +1 after the translation 
 !Any relative fitness should be calculated as: [ABS(ptrPop1(I)%Fitness - WorstPopFitness)+1]
 !And thus the proportionate selection should be <[ABS(ptrPop1(I)%Fitness - WorstPopFitness)+1] / (TotalPopFitness)>

End subroutine ComputeFitnessSum_WorstFitness
!!===========================END: ComputeFitnessSum_WorstFitness ============

! ===============  START: Proportionate Selection ============================
Integer FUNCTION SelectChromosome(k)

real*8:: random1     ! for random1()

Integer, INTENT(IN) :: k   ! the k will help vary the random seed if it is called too fast sequentially 
Integer:: i
real*8:: TheValue, RunningSum 

!TheValue = (RandNumTam(k) * TotalPopFitness)+1       !   +1 becasue it will produce from 0 to (TotalPopFitness-1) [in case an integer range]
TheValue = (random1()* TotalPopFitness)+1

RunningSum = 0
i = 0

DO While ((RunningSum <= TheValue) .AND. (i < PopSize))
   i = i + 1
   RunningSum = RunningSum + Abs(ptrPop1(i)%Fitness - WorstPopFitness) + 1 !See "ComputeFitnessSum_WorstFitness" subroutine for comments and understanding 
END DO

SelectChromosome = i 


End FUNCTION SelectChromosome
! ==========================  END: Proportionate Selection ============================

!===============  START: Associative Memory Based Crossover ============================

subroutine My_AM_Crossover(i, j, n)

 INTEGER:: i, j, n, NewPoploc
 ! popAMBottom, popAMUpper
 
  NewPoploc = CurrentPosNewPopulation
  
  ptrPop2(NewPoploc +1) = ptrPop1(i)
  ptrPop2(NewPoploc +2) = ptrPop1(j)


                          !! PART (1): I(upper) with J(Bottom) or AM(Bottom)
                           !  ================================
 tmpstruct(1)%chromo(1:n) = ptrPop1(i)%chromo(1:n)
 tmpstruct(1)%chromo(n+1:ChromosomeLength)= ptrPop1(j)%chromo(n+1:ChromosomeLength)
 
 tmpstruct(2)%chromo(1:n) = ptrPop1(i)%chromo(1:n)
 tmpstruct(2)%chromo(n+1:ChromosomeLength)= popAMBottom(n)%chromo(n+1:ChromosomeLength)
 
 tmpstruct(1)%Fitness = CalculateFitness(1,3) ! CalculateFitness(tmpstruct(1)), here PopType=3=> tmpstruct
 tmpstruct(2)%Fitness = Calculatefitness(2,3) ! CalculateFitness(tmpstruct(2))

 
 IF (tmpstruct(1)%Fitness < tmpstruct(2)%Fitness) THEN
      IF (tmpstruct(1)%Fitness < ptrPop2((NewPopLoc+1))%Fitness) THEN
        popAMBottom(n)%chromo(n+1:ChromosomeLength) = ptrPop1(j)%chromo(n+1:ChromosomeLength)
        ptrPop2(NewPopLoc+1) = tmpstruct(1)   ! So, I with J was better and thus fragment of J is updated in AM(Bottom)
      END IF
  ELSE
        IF (tmpstruct(2)%Fitness < ptrPop2((NewPopLoc+1))%Fitness) THEN
           ptrPop2((NewPopLoc+1)) = tmpstruct(2)   !So, I(Upper) with AM(Bottom(n)) is  better
        END IF
 END IF 
                           !! PART (2): J(Upper) or AM(Upper) with I(Bottom) 
                           !  ===============================================
   
 tmpstruct(1)%chromo(1:n) = ptrPop1(j)%chromo(1:n)
 tmpstruct(1)%chromo(n+1:ChromosomeLength) = ptrPop1(i)%chromo(n+1:ChromosomeLength)
 tmpstruct(2)%chromo(1:n) = popAMUpper(n)%chromo(1:n)
 tmpstruct(2)%chromo(n+1:ChromosomeLength) = ptrPop1(i)%chromo(n+1:ChromosomeLength)

 tmpstruct(1)%Fitness = CalculateFitness(1,3) ! CalculateFitness(tmpstruct(1))
 tmpstruct(2)%Fitness = CalculateFitness(2,3) !CalculateFitness(tmpstruct(2))

 
 If (tmpstruct(1)%Fitness < tmpstruct(2)%Fitness) THEN
      IF (tmpstruct(1)%Fitness < ptrPop2((NewPopLoc+2))%Fitness) THEN
        popAMUpper(n)%chromo(1:n) = ptrPop1(j)%chromo(1:n)
        ptrPop2((NewPopLoc+2)) = tmpstruct(1)   ! So, I with J was better and thus fragment of J is updated in AM(Upper)
      END IF
 ELSE
      IF (tmpstruct(2)%Fitness < ptrPop2((NewPopLoc+2))%Fitness) THEN
        ptrPop2((NewPopLoc+2)) = tmpstruct(2)   ! So, AM(Upper(n)) with I(Bottom) is better
      END IF
 END IF 

END subroutine 
!===============  END: Associative Memory Based Crossover ============================

! ================= START: Mutation =================================================
SUBROUTINE Mutation(i, n)
   
INTEGER, INTENT (INOUT):: i
   INTEGER, INTENT (IN)::n 
   Real*8:: m

   m=random1()*2
   tmpstruct(1)=ptrPop2(i)   
   tmpstruct(1)%chromo(n) = m !(RandNumTam(i*n) * 2)
   tmpstruct(1)%Fitness=CalculateFitness(1,3)
    
   IF  (tmpstruct(1)%Fitness < ptrPop2(i)%Fitness) THEN
           ptrPop2(i)%chromo(n) = tmpstruct(1)%chromo(n)          ! change will be done IFF the mutated one is better.
   END IF 
  
  ! Mutation = 1 ! return 1 to indicate that the job is done rather than success

END SUBROUTINE Mutation
! ================= END: Mutation =====================================

!----------------------------------------------------------------------------
real*8 FUNCTION CalculateFitness(i,PopType)
Integer, INTENT(IN):: i, PopType

real*8:: TMScore_avg, o3DIGARS_Avg, Upper1, Lower1, Lower2, PCC_Avg, LowEMdl, TMLowEMdl_Avg, native_count_Avg, nat_count

Integer:: j, k, idx_lowE, native_index

  Call ChromoToVarValue(i, PopType) ! this will calculate and store the values in Ks(1), Ks(2), Ks(3), ...

  
! Ks(1)=1; Ks(2)=1; Ks(3)=1; Ks(4)=1; Ks(5)=1; Ks(6)=1    ! Was forced these values to test the correctness


  !Compute xDFIRE for all
  DO j =1,TotRec
	 RecList(j)%o3DIGARS = RecList(j)%ASA_Egy_Spinex_Data_Hoque + Ks(1)*RecList(j)%ThrDIGARS_Egy+Ks(2)*RecList(j)%ASA_Egy+Ks(3)*RecList(j)%Psi_Egy_Comb_Data_Hoque+Ks(4)*RecList(j)%Psi_Egy_SSD_Data_Hoque+Ks(5)*RecList(j)%Psi_Egy_Triplet_Hoque
	 !write(*,*) RecList(j)%o3DIGARS
  END DO

	
  Do j=1, myFileCount  
       TMScore_avg=0; o3DIGARS_Avg=0; !xDFIRE_Avg
       LowEMdl=RecList(FileList(j)%StartIdx)%o3DIGARS; idx_lowE=FileList(j)%StartIdx ! NewLine
     
 !write(*,*) j, FileList(j)%StartIdx, FileList(j)%EndIdx

    Do k=FileList(j)%StartIdx, (FileList(j)%EndIdx-1)
           TMScore_avg= TMScore_avg + RecList(k)%TMScore
           o3DIGARS_Avg = o3DIGARS_Avg  + RecList(k)%o3DIGARS 
			
          IF ((RecList(k)%o3DIGARS < LowEMdl) .AND. (k/=FileList(j)%EndIdx)) THEN     ! NewLines: 2nd term ensures that the native is not included
               LowEMdl=RecList(k)%o3DIGARS; idx_lowE=k
!                write(*,*) 'k=',k, LowEMdl
           END IF 
      END Do  
			
			native_index = FileList(j)%EndIdx;
			IF(RecList(native_index)%o3DIGARS < LowEMdl) THEN
				nat_count=nat_count+1.0;
			END IF
			
		!	write(*,*) 'casp8lowEmdl=',LowEMdl
		!	write(*,*) 'casp8native_egy=',RecList(native_index)%o3DIGARS
		!	write(*,*) 'casp8native_count=',nat_count
			
           TMScore_avg = TMScore_avg / (FileList(j)%RecCount-1)
           o3DIGARS_Avg = o3DIGARS_Avg / (FileList(j)%RecCount-1)
		   
           FileList(j)%TMLowEMdl = RecList(idx_lowE)%TMScore ! NewLine for minEofModel(TMscore) 
           		   		   
           Upper1=0; Lower1=0; Lower2=0  ! Upper and Lowers are the U and L of the formula to compute PCC
		   Zscore=0;
		   std_dev=0;

      Do k=FileList(j)%StartIdx, (FileList(j)%EndIdx-1)
           Upper1 =  Upper1 + (RecList(k)%o3DIGARS - o3DIGARS_Avg) * (RecList(k)%TMScore - TMscore_Avg)            
           Lower1 = Lower1+ (RecList(k)%o3DIGARS - o3DIGARS_Avg)**2
           Lower2 = Lower2 + (RecList(k)%TMScore - TMscore_Avg)**2
      END Do
		  std_dev = Lower1/(FileList(j)%RecCount-1)
		  std_dev = std_dev ** 0.5;
          Lower1 = Lower1 ** 0.5; Lower2 = Lower2 ** 0.5;
          FileList(j)%PCC= (Upper1) /(Lower1 * Lower2)
		  
		  Zscore = (RecList(native_index)%o3DIGARS - o3DIGARS_Avg)/std_dev
		  FileList(j)%Z_score = Zscore
        
         ! write (*,*) 'FileList(j)%PCC =', FileList(j)%PCC
         ! write (*,*) 'RecList(idx_lowE)%TMScore=', RecList(idx_lowE)%TMScore 
!STOP
 
  END DO


        PCC_Avg=0; TMLowEMdl_Avg=0; Z_score_Avg=0;
		
     Do j=1, myFileCount
        PCC_Avg=PCC_Avg+FileList(j)%PCC
        TMLowEMdl_Avg = TMLowEMdl_Avg + FileList(j)%TMLowEMdl
		Z_score_Avg = Z_score_Avg+FileList(j)%Z_score
     END DO 
        PCC_Avg = PCC_Avg / myFileCount
        TMLowEMdl_Avg = TMLowEMdl_Avg / myFileCount
		Z_score_Avg = Z_score_Avg/myFileCount
		native_count_Avg = nat_count/myFileCount
	!	write(*,*) 'casp8nat_count=',REAL(nat_count)
	!	write(*,*) 'casp8myFileCount=',REAL(myFileCount)
	!	write(*,*) 'casp8native_count_Avg=',native_count_Avg
	!	write (*,*) 'PCC_Avg =', PCC_Avg
	!	write (*,*) 'TMLowEMdl_Avg =', TMLowEMdl_Avg
    !    CalculateFitness =  (1.0)*(PCC_Avg +(-1 * TMLowEMdl_Avg)) + (1.0) * (CalculateFitness1(i,PopType)) + (1.0)* (CalculateFitness2(i,PopType))
	CalculateFitness =  (1.0)*(PCC_Avg +(-1 * TMLowEMdl_Avg)+(-1*native_count_Avg)+Z_score_Avg) + (1.0) * (CalculateFitness1(i,PopType)) + (1.0) * (CalculateFitness2(i,PopType)) + (1.0) * (CalculateFitness3(i,PopType))
		

 IF (PopType==1) THEN
   ptrPop1(i)%PCC=PCC_Avg
   ptrPop1(i)%TMLowEMdl=TMLowEMdl_Avg 
   ptrPop1(i)%Native_Count=native_count_Avg*myFileCount
   ptrPop1(i)%ZSCORE=Z_score_Avg
 else IF (PopType==2) THEN
   ptrPop2(i)%PCC=PCC_Avg
   ptrPop2(i)%TMLowEMdl=TMLowEMdl_Avg 
   ptrPop2(i)%Native_Count=native_count_Avg*myFileCount
   ptrPop2(i)%ZSCORE=Z_score_Avg
 else
   tmpstruct(i)%PCC=PCC_Avg
   tmpstruct(i)%TMLowEMdl=TMLowEMdl_Avg
   tmpstruct(i)%Native_Count=native_count_Avg*myFileCount
   tmpstruct(i)%ZSCORE=Z_score_Avg
 END IF 


END FUNCTION CalculateFitness
!-----------------------------------------
!----------------------------------------------------------------------------
real*8 FUNCTION CalculateFitness1(i,PopType)
Integer, INTENT(IN):: i, PopType

real*8:: TMScore_avg, o3DIGARS_Avg, Upper1, Lower1, Lower2, PCC_Avg, LowEMdl, TMLowEMdl_Avg, native_count_Avg, nat_count

Integer:: j, k, idx_lowE, native_index

  Call ChromoToVarValue(i, PopType) ! this will calculate and store the values in Ks(1), Ks(2), Ks(3), ...

  
! Ks(1)=1; Ks(2)=1; Ks(3)=1    ! Was forced these values to test the correctness


  !Compute xDFIRE for all
  DO j =1,TotRec1
	RecList1(j)%o3DIGARS = RecList1(j)%ASA_Egy_Spinex_Data_Hoque + Ks(1)*RecList1(j)%ThrDIGARS_Egy+Ks(2)*RecList1(j)%ASA_Egy+Ks(3)*RecList1(j)%Psi_Egy_Comb_Data_Hoque+Ks(4)*RecList1(j)%Psi_Egy_SSD_Data_Hoque+Ks(5)*RecList1(j)%Psi_Egy_Triplet_Hoque
  END DO

	
  Do j=1, myFileCount1  
       TMScore_avg=0; o3DIGARS_Avg=0; 
       LowEMdl=RecList1(FileList1(j)%StartIdx)%o3DIGARS; idx_lowE=FileList1(j)%StartIdx ! NewLine
     
 !write(*,*) j, FileList1(j)%StartIdx, FileList1(j)%EndIdx

    Do k=FileList1(j)%StartIdx, (FileList1(j)%EndIdx-1)
           TMScore_avg= TMScore_avg + RecList1(k)%TMScore
           o3DIGARS_Avg = o3DIGARS_Avg  + RecList1(k)%o3DIGARS 
       
          IF ((RecList1(k)%o3DIGARS < LowEMdl) .AND. (k/=FileList1(j)%EndIdx)) THEN     ! NewLines: 2nd term ensures that the native is not included
               LowEMdl=RecList1(k)%o3DIGARS; idx_lowE=k
!                write(*,*) 'k=',k, LowEMdl
           END IF 
      END Do  

			native_index = FileList1(j)%EndIdx;
			IF(RecList1(native_index)%o3DIGARS < LowEMdl) THEN
				nat_count=nat_count+1;
			END IF
		!	write(*,*) 'casp9lowEmdl=',LowEMdl
		!	write(*,*) 'casp9native_egy=',RecList1(native_index)%o3DIGARS
           TMScore_avg = TMScore_avg / (FileList1(j)%RecCount-1)
           o3DIGARS_Avg = o3DIGARS_Avg / (FileList1(j)%RecCount-1)
          
           FileList1(j)%TMLowEMdl = RecList1(idx_lowE)%TMScore ! NewLine for minEofModel(TMscore) 
                  
           Upper1=0; Lower1=0; Lower2=0  ! Upper and Lowers are the U and L of the formula to compute PCC
			 Zscore=0;
			 std_dev=0;
      Do k=FileList1(j)%StartIdx, (FileList1(j)%EndIdx-1)
           Upper1 =  Upper1 + (RecList1(k)%o3DIGARS - o3DIGARS_Avg) * (RecList1(k)%TMScore - TMscore_Avg)            
           Lower1 = Lower1+ (RecList1(k)%o3DIGARS - o3DIGARS_Avg)**2
           Lower2 = Lower2 + (RecList1(k)%TMScore - TMscore_Avg)**2
      END Do
		  std_dev = Lower1/(FileList1(j)%RecCount-1)
		  std_dev = std_dev ** 0.5;
          Lower1 = Lower1 ** 0.5; Lower2 = Lower2 ** 0.5; 
          FileList1(j)%PCC= (Upper1) /(Lower1 * Lower2)
		  
		  Zscore = (RecList1(native_index)%o3DIGARS - o3DIGARS_Avg)/std_dev
		  FileList1(j)%Z_score = Zscore
			
        !  write (*,*) 'FileList1(j)%PCC =', FileList1(j)%PCC
       !   write (*,*) 'RecList1(idx_lowE)%TMScore=', RecList1(idx_lowE)%TMScore 
!STOP
 
  END DO


        PCC_Avg=0; TMLowEMdl_Avg=0; Z_score_Avg=0; native_count_Avg = 0;
     Do j=1, myFileCount1
        PCC_Avg=PCC_Avg+FileList1(j)%PCC
        TMLowEMdl_Avg = TMLowEMdl_Avg + FileList1(j)%TMLowEMdl
		Z_score_Avg = Z_score_Avg+FileList1(j)%Z_score
     END DO 
        PCC_Avg = PCC_Avg / myFileCount1
        TMLowEMdl_Avg = TMLowEMdl_Avg / myFileCount1
		Z_score_Avg = Z_score_Avg/myFileCount1
		native_count_Avg = nat_count/myFileCount1
		
        CalculateFitness1 =  PCC_Avg +(-1 * TMLowEMdl_Avg)+(-1 * native_count_Avg)+Z_score_Avg

 IF (PopType==1) THEN
   ptrPop1(i)%PCC1=PCC_Avg
   ptrPop1(i)%TMLowEMdl1=TMLowEMdl_Avg
   ptrPop1(i)%Native_Count1=native_count_Avg*myFileCount1
   ptrPop1(i)%ZSCORE1=Z_score_Avg
 else IF (PopType==2) THEN
   ptrPop2(i)%PCC1=PCC_Avg
   ptrPop2(i)%TMLowEMdl1=TMLowEMdl_Avg
   ptrPop2(i)%Native_Count1=native_count_Avg*myFileCount1
   ptrPop2(i)%ZSCORE1=Z_score_Avg
 else
   tmpstruct(i)%PCC1=PCC_Avg
   tmpstruct(i)%TMLowEMdl1=TMLowEMdl_Avg
   tmpstruct(i)%Native_Count1=native_count_Avg*myFileCount1
   tmpstruct(i)%ZSCORE1=Z_score_Avg
 END IF 


END FUNCTION CalculateFitness1
!-----------------------------------------

!----------------------------------------------------------------------------
real*8 FUNCTION CalculateFitness2(i,PopType)
Integer, INTENT(IN):: i, PopType

real*8:: TMScore_avg, o3DIGARS_Avg, Upper1, Lower1, Lower2, PCC_Avg, LowEMdl, TMLowEMdl_Avg, native_count_Avg, nat_count

Integer:: j, k, idx_lowE, native_index

  Call ChromoToVarValue(i, PopType) ! this will calculate and store the values in Ks(1), Ks(2), Ks(3), ...

  
! Ks(1)=1; Ks(2)=1; Ks(3)=1    ! Was forced these values to test the correctness


  !Compute xDFIRE for all
  DO j =1,TotRec2
	RecList2(j)%o3DIGARS = RecList2(j)%ASA_Egy_Spinex_Data_Hoque + Ks(1)*RecList2(j)%ThrDIGARS_Egy+Ks(2)*RecList2(j)%ASA_Egy+Ks(3)*RecList2(j)%Psi_Egy_Comb_Data_Hoque+Ks(4)*RecList2(j)%Psi_Egy_SSD_Data_Hoque+Ks(5)*RecList2(j)%Psi_Egy_Triplet_Hoque
  END DO


  Do j=1, myFileCount2  
       TMScore_avg=0; o3DIGARS_Avg=0; 
       LowEMdl=RecList2(FileList2(j)%StartIdx)%o3DIGARS; idx_lowE=FileList2(j)%StartIdx ! NewLine
     
 !write(*,*) j, FileList2(j)%StartIdx, FileList2(j)%EndIdx

    Do k=FileList2(j)%StartIdx, (FileList2(j)%EndIdx-1)
           TMScore_avg= TMScore_avg + RecList2(k)%TMScore
           o3DIGARS_Avg = o3DIGARS_Avg  + RecList2(k)%o3DIGARS 
       
          IF ((RecList2(k)%o3DIGARS < LowEMdl) .AND. (k/=FileList2(j)%EndIdx)) THEN     ! NewLines: 2nd term ensures that the native is not included
               LowEMdl=RecList2(k)%o3DIGARS; idx_lowE=k
!                write(*,*) 'k=',k, LowEMdl
           END IF 
      END Do  

			native_index = FileList2(j)%EndIdx;
			IF(RecList2(native_index)%o3DIGARS < LowEMdl) THEN
				nat_count=nat_count+1;
			END IF
		!	write(*,*) 'casp10lowEmdl=',LowEMdl
		!	write(*,*) 'casp10native_egy=',RecList2(native_index)%o3DIGARS
			
           TMScore_avg = TMScore_avg / (FileList2(j)%RecCount-1)
           o3DIGARS_Avg = o3DIGARS_Avg / (FileList2(j)%RecCount-1)
          
           FileList2(j)%TMLowEMdl = RecList2(idx_lowE)%TMScore ! NewLine for minEofModel(TMscore) 
                  
           Upper1=0; Lower1=0; Lower2=0  ! Upper and Lowers are the U and L of the formula to compute PCC
			Zscore=0;
			std_dev=0;
      Do k=FileList2(j)%StartIdx, (FileList2(j)%EndIdx-1)
           Upper1 =  Upper1 + (RecList2(k)%o3DIGARS - o3DIGARS_Avg) * (RecList2(k)%TMScore - TMscore_Avg)            
           Lower1 = Lower1+ (RecList2(k)%o3DIGARS - o3DIGARS_Avg)**2
           Lower2 = Lower2 + (RecList2(k)%TMScore - TMscore_Avg)**2
      END Do
          std_dev = Lower1/(FileList2(j)%RecCount-1)
		  std_dev = std_dev ** 0.5;
		  Lower1 = Lower1 ** 0.5; Lower2 = Lower2 ** 0.5; 
          FileList2(j)%PCC= (Upper1) /(Lower1 * Lower2)
		  

			Zscore = (RecList2(native_index)%o3DIGARS - o3DIGARS_Avg)/std_dev
		  FileList2(j)%Z_score = Zscore		  
        
        !  write (*,*) 'FileList2(j)%PCC =', FileList2(j)%PCC
       !   write (*,*) 'RecList2(idx_lowE)%TMScore=', RecList2(idx_lowE)%TMScore 
!STOP
 
  END DO


        PCC_Avg=0; TMLowEMdl_Avg=0
     Do j=1, myFileCount2
        PCC_Avg=PCC_Avg+FileList2(j)%PCC
        TMLowEMdl_Avg = TMLowEMdl_Avg + FileList2(j)%TMLowEMdl
		Z_score_Avg = Z_score_Avg+FileList2(j)%Z_score
     END DO 
        PCC_Avg = PCC_Avg / myFileCount2
        TMLowEMdl_Avg = TMLowEMdl_Avg / myFileCount2
		Z_score_Avg = Z_score_Avg/myFileCount2
		native_count_Avg = nat_count/myFileCount2
		
        CalculateFitness2 =  PCC_Avg +(-1 * TMLowEMdl_Avg)+(-1 * native_count_Avg)+Z_score_Avg

 IF (PopType==1) THEN
   ptrPop1(i)%PCC2=PCC_Avg
   ptrPop1(i)%TMLowEMdl2=TMLowEMdl_Avg 
   ptrPop1(i)%Native_Count2=native_count_Avg*myFileCount2
   ptrPop1(i)%ZSCORE2=Z_score_Avg
 else IF (PopType==2) THEN
   ptrPop2(i)%PCC2=PCC_Avg
   ptrPop2(i)%TMLowEMdl2=TMLowEMdl_Avg
   ptrPop2(i)%Native_Count2=native_count_Avg*myFileCount2
   ptrPop2(i)%ZSCORE2=Z_score_Avg
 else
   tmpstruct(i)%PCC2=PCC_Avg
   tmpstruct(i)%TMLowEMdl2=TMLowEMdl_Avg
   tmpstruct(i)%Native_Count2=native_count_Avg*myFileCount2
   tmpstruct(i)%ZSCORE2=Z_score_Avg
 END IF 


END FUNCTION CalculateFitness2
!-----------------------------------------


real*8 FUNCTION CalculateFitness3(i,PopType)
Integer, INTENT(IN):: i, PopType

real*8:: TMScore_avg, o3DIGARS_Avg, Upper1, Lower1, Lower2, PCC_Avg, LowEMdl, TMLowEMdl_Avg, native_count_Avg, nat_count

Integer:: j, k, idx_lowE, native_index

  Call ChromoToVarValue(i, PopType) ! this will calculate and store the values in Ks(1), Ks(2), Ks(3), ...

  
! Ks(1)=1; Ks(2)=1; Ks(3)=1    ! Was forced these values to test the correctness


  !Compute xDFIRE for all
  DO j =1,TotRec3
	RecList3(j)%o3DIGARS = RecList3(j)%ASA_Egy_Spinex_Data_Hoque + Ks(1)*RecList3(j)%ThrDIGARS_Egy+Ks(2)*RecList3(j)%ASA_Egy+Ks(3)*RecList3(j)%Psi_Egy_Comb_Data_Hoque+Ks(4)*RecList3(j)%Psi_Egy_SSD_Data_Hoque+Ks(5)*RecList3(j)%Psi_Egy_Triplet_Hoque
  END DO


  Do j=1, myFileCount3  
       TMScore_avg=0; o3DIGARS_Avg=0; 
       LowEMdl=RecList3(FileList3(j)%StartIdx)%o3DIGARS; idx_lowE=FileList3(j)%StartIdx ! NewLine
     
 !write(*,*) j, FileList3(j)%StartIdx, FileList3(j)%EndIdx

    Do k=FileList3(j)%StartIdx, (FileList3(j)%EndIdx-1)
           TMScore_avg= TMScore_avg + RecList3(k)%TMScore
           o3DIGARS_Avg = o3DIGARS_Avg  + RecList3(k)%o3DIGARS 
       
          IF ((RecList3(k)%o3DIGARS < LowEMdl) .AND. (k/=FileList3(j)%EndIdx)) THEN     ! NewLines: 2nd term ensures that the native is not included
               LowEMdl=RecList3(k)%o3DIGARS; idx_lowE=k
!                write(*,*) 'k=',k, LowEMdl
           END IF 
      END Do  

			native_index = FileList3(j)%EndIdx;
			IF(RecList3(native_index)%o3DIGARS < LowEMdl) THEN
				nat_count=nat_count+1;
			END IF
		!	write(*,*) 'casp11lowEmdl=',LowEMdl
		!	write(*,*) 'casp11native_egy=',RecList3(native_index)%o3DIGARS
           TMScore_avg = TMScore_avg / (FileList3(j)%RecCount-1)
           o3DIGARS_Avg = o3DIGARS_Avg / (FileList3(j)%RecCount-1)
          
           FileList3(j)%TMLowEMdl = RecList3(idx_lowE)%TMScore ! NewLine for minEofModel(TMscore) 
                  
           Upper1=0; Lower1=0; Lower2=0  ! Upper and Lowers are the U and L of the formula to compute PCC
			   Zscore=0;
			   std_dev=0;
      Do k=FileList3(j)%StartIdx, (FileList3(j)%EndIdx-1)
           Upper1 =  Upper1 + (RecList3(k)%o3DIGARS - o3DIGARS_Avg) * (RecList3(k)%TMScore - TMscore_Avg)            
           Lower1 = Lower1+ (RecList3(k)%o3DIGARS - o3DIGARS_Avg)**2
           Lower2 = Lower2 + (RecList3(k)%TMScore - TMscore_Avg)**2
      END Do
          std_dev = Lower1/(FileList3(j)%RecCount-1)
		  std_dev = std_dev ** 0.5;
		  Lower1 = Lower1 ** 0.5; Lower2 = Lower2 ** 0.5; 
          FileList3(j)%PCC= (Upper1) /(Lower1 * Lower2)    
        
		Zscore = (RecList3(native_index)%o3DIGARS - o3DIGARS_Avg)/std_dev
		  FileList3(j)%Z_score = Zscore
		
        !  write (*,*) 'FileList3(j)%PCC =', FileList3(j)%PCC
       !   write (*,*) 'RecList3(idx_lowE)%TMScore=', RecList3(idx_lowE)%TMScore 
!STOP
 
  END DO


        PCC_Avg=0; TMLowEMdl_Avg=0; Z_score_Avg=0; native_count_Avg = 0;
     Do j=1, myFileCount3
        PCC_Avg=PCC_Avg+FileList3(j)%PCC
        TMLowEMdl_Avg = TMLowEMdl_Avg + FileList3(j)%TMLowEMdl
		Z_score_Avg = Z_score_Avg+FileList3(j)%Z_score
     END DO 
        PCC_Avg = PCC_Avg / myFileCount3
        TMLowEMdl_Avg = TMLowEMdl_Avg / myFileCount3
		Z_score_Avg = Z_score_Avg/myFileCount3
		native_count_Avg = nat_count/myFileCount3
		
        CalculateFitness3 =  PCC_Avg +(-1 * TMLowEMdl_Avg) +(-1 * native_count_Avg) + Z_score_Avg

 IF (PopType==1) THEN
   ptrPop1(i)%PCC3=PCC_Avg
   ptrPop1(i)%TMLowEMdl3=TMLowEMdl_Avg 
   ptrPop1(i)%Native_Count3=native_count_Avg*myFileCount3
   ptrPop1(i)%ZSCORE3=Z_score_Avg
 else IF (PopType==2) THEN
   ptrPop2(i)%PCC3=PCC_Avg
   ptrPop2(i)%TMLowEMdl3=TMLowEMdl_Avg 
   ptrPop2(i)%Native_Count3=native_count_Avg*myFileCount3
   ptrPop2(i)%ZSCORE3=Z_score_Avg
 else
   tmpstruct(i)%PCC3=PCC_Avg
   tmpstruct(i)%TMLowEMdl3=TMLowEMdl_Avg
   tmpstruct(i)%Native_Count3=native_count_Avg*myFileCount3
   tmpstruct(i)%ZSCORE3=Z_score_Avg
 END IF 


END FUNCTION CalculateFitness3
!-----------------------------------------


!------------------------- START: RealToStr --------------------------------------------
character(100) FUNCTION RealToStr(MyRealNumber)
REAL*8, INTENT(IN):: MyRealNumber
character(len=100) :: CR
  write(CR,*) MyRealNumber
  CR=adjustl(CR)
  RealToStr=trim(CR)
END FUNCTION RealToStr
!-------------------------- END: RealToStr ---------------------------------------------

!------------------------- START: IntegerToStr --------------------------------------------
character(100) FUNCTION IntegerToStr (MyIntNumber)
INTEGER, INTENT(IN):: MyIntNumber
character(len=100) :: CI
  write(CI,*) MyIntNumber
  CI=adjustl(CI)
  IntegerToStr=trim(CI)
END FUNCTION IntegerToStr 
!-------------------------- END: RealToStr ---------------------------------------------

END MODULE

!=====================================================================================================================================================
!///////////////////////////////////////////////////////// MODULE BOUNDARY [ENDs] ////////////////////////////////////////////////////////////////////
!=====================================================================================================================================================

!-------------------------------------------------------------
SUBROUTINE LoadTheDataFromInput1

use myGAfuns_dat

Integer:: f1, f2, ierr, IOstatus, i
Logical:: myEOF
character (len=30):: RdLine

call system('ls input1/*.csv > TableList1.txt')
f1=freeunit(); open (f1, file='TableList1.txt', action ='read');Call Fseek(f1,0,0,ierr)


myEOF=.false.; myFileCount1=0
DO while (myEOF == .false.)        
  read(f1, '(a)', IOSTAT=IOstatus) RdLine
  IF ( IOstatus < 0) Then  
       myEOF = .true.
   ELSE
      myFileCount1=myFileCount1+1           
   END IF 
END DO !While END

 
IF(ALLOCATED(FileList1)) DEALLOCATE(FileList1)
allocate (FileList1(myFileCount1+10))


Call Fseek(f1,0,0,ierr); CurrLoc1=0
DO i=1, myFileCount1
  read(f1, '(a)', IOSTAT=IOstatus) RdLine 
  Call CountTableDetailForEachFile1(i,RdLine)
END DO



TotRec1=0
DO i=1, myFileCount1
   TotRec1=TotRec1 + FileList1(i)%RecCount
END DO

!write(*,*) 'TotRec1=', TotRec1



IF(ALLOCATED(RecList1)) DEALLOCATE(RecList1)
allocate (RecList1(TotRec1+10))



CurrLoc1=0
Call Fseek(f1,0,0,ierr)
DO i=1, myFileCount1
  read(f1, '(a)', IOSTAT=IOstatus) RdLine  
  Call LoadTableDetailForEachFile1(i, RdLine)
END DO

 Write (*,*) 'Load-Table-Detail-For-Each-File is done successfully 1...'

close(f1)

END SUBROUTINE LoadTheDataFromInput1

!-------------------------------------------------

SUBROUTINE LoadTheDataFromInput2

use myGAfuns_dat

Integer:: f1, f2, ierr, IOstatus, i
Logical:: myEOF
character (len=30):: RdLine

call system('ls input2/*.csv > TableList2.txt')
f1=freeunit(); open (f1, file='TableList2.txt', action ='read');Call Fseek(f1,0,0,ierr)

myEOF=.false.; myFileCount2=0
DO while (myEOF == .false.)        
  read(f1, '(a)', IOSTAT=IOstatus) RdLine
  IF ( IOstatus < 0) Then  
       myEOF = .true.
   ELSE
      myFileCount2=myFileCount2+1           
   END IF 
END DO !While END

 
IF(ALLOCATED(FileList2)) DEALLOCATE(FileList2)
allocate (FileList2(myFileCount2+10))


Call Fseek(f1,0,0,ierr); CurrLoc2=0
DO i=1, myFileCount2
  read(f1, '(a)', IOSTAT=IOstatus) RdLine 
  Call CountTableDetailForEachFile2(i,RdLine)
END DO



TotRec2=0
DO i=1, myFileCount2
   TotRec2=TotRec2 + FileList2(i)%RecCount
END DO

!write(*,*) 'TotRec2=', TotRec2



IF(ALLOCATED(RecList2)) DEALLOCATE(RecList2)
allocate (RecList2(TotRec2+10))


CurrLoc2=0
Call Fseek(f1,0,0,ierr)
DO i=1, myFileCount2
  read(f1, '(a)', IOSTAT=IOstatus) RdLine  
  Call LoadTableDetailForEachFile2(i, RdLine)
END DO
 Write (*,*) 'Load-Table-Detail-For-Each-File is done successfully 2...'

close(f1)

END SUBROUTINE LoadTheDataFromInput2

!-------------------------------------------------------------

!-------------------------------------------------

SUBROUTINE LoadTheDataFromInput3

use myGAfuns_dat

Integer:: f1, f2, ierr, IOstatus, i
Logical:: myEOF
character (len=30):: RdLine

call system('ls input2/*.csv > TableList2.txt')
f1=freeunit(); open (f1, file='TableList3.txt', action ='read');Call Fseek(f1,0,0,ierr)

myEOF=.false.; myFileCount3=0
DO while (myEOF == .false.)        
  read(f1, '(a)', IOSTAT=IOstatus) RdLine
  IF ( IOstatus < 0) Then  
       myEOF = .true.
   ELSE
      myFileCount3=myFileCount3+1           
   END IF 
END DO !While END

 
IF(ALLOCATED(FileList3)) DEALLOCATE(FileList3)
allocate (FileList3(myFileCount3+10))


Call Fseek(f1,0,0,ierr); CurrLoc3=0
DO i=1, myFileCount3
  read(f1, '(a)', IOSTAT=IOstatus) RdLine 
  Call CountTableDetailForEachFile3(i,RdLine)
END DO



TotRec3=0
DO i=1, myFileCount3
   TotRec3=TotRec3 + FileList3(i)%RecCount
END DO

!write(*,*) 'TotRec3=', TotRec3



IF(ALLOCATED(RecList3)) DEALLOCATE(RecList3)
allocate (RecList3(TotRec3+10))


CurrLoc3=0
Call Fseek(f1,0,0,ierr)
DO i=1, myFileCount3
  read(f1, '(a)', IOSTAT=IOstatus) RdLine  
  Call LoadTableDetailForEachFile3(i, RdLine)
END DO
 Write (*,*) 'Load-Table-Detail-For-Each-File is done successfully 3...'

close(f1)

END SUBROUTINE LoadTheDataFromInput3

!-------------------------------------------------------------

!-------------------------------------------------
SUBROUTINE LoadTheDataFromInput

use myGAfuns_dat

Integer:: f1, f2, ierr, IOstatus, i
Logical:: myEOF
character (len=30):: RdLine

call system('ls input/*.csv > TableList.txt')
f1=freeunit(); open (f1, file='TableList.txt', action ='read');Call Fseek(f1,0,0,ierr)

myEOF=.false.; myFileCount=0
DO while (myEOF == .false.)        
  read(f1, '(a)', IOSTAT=IOstatus) RdLine
  IF ( IOstatus < 0) Then  
       myEOF = .true.
   ELSE
      myFileCount=myFileCount+1           
   END IF 
END DO !While END

 
IF(ALLOCATED(FileList)) DEALLOCATE(FileList)
allocate (FileList(myFileCount+10))

Call Fseek(f1,0,0,ierr); CurrLoc=0
DO i=1, myFileCount
  read(f1, '(a)', IOSTAT=IOstatus) RdLine 
  Call CountTableDetailForEachFile(i,RdLine)
END DO



TotRec=0
DO i=1, myFileCount
   TotRec=TotRec + FileList(i)%RecCount
END DO

!write(*,*) 'TotRec=', TotRec

IF(ALLOCATED(RecList)) DEALLOCATE(RecList)
allocate (RecList(TotRec+10))



CurrLoc=0
Call Fseek(f1,0,0,ierr)
DO i=1, myFileCount
  read(f1, '(a)', IOSTAT=IOstatus) RdLine  
  Call LoadTableDetailForEachFile(i, RdLine)
END DO
 Write (*,*) 'Load-Table-Detail-For-Each-File is done successfully ...'

close(f1)

END SUBROUTINE LoadTheDataFromInput
!-------------------------------------------------------------

SUBROUTINE CountTableDetailForEachFile(i, fname)

use myGAfuns_dat

 Integer, INTENT(IN):: i 
 character (len=30), INTENT(IN):: fname
 Logical:: myEOF
 Integer:: ierr, IOstatus, myCount
 character (len=30):: RdLine


 Integer:: f2
 

 f2=freeunit(); open (f2, file=trim(fname), action ='read');Call Fseek(f2,0,0,ierr)
 read(f2, '(a)', IOSTAT=IOstatus) RdLine  ! Read off the header
 
 !CurrLoc = CurrLoc + 1
 myEOF=.false.; myCount=0
 DO while (myEOF == .false.)        
    read(f2, '(a)', IOSTAT=IOstatus) RdLine
   IF ( IOstatus < 0) Then  
       myEOF = .true.
   ELSE
      myCount=myCount+1           
   END IF 
 END DO !While END

 FileList(i)%RecCount = myCount
 FileList(i)%Startidx = CurrLoc+1
 FileList(i)%Endidx = CurrLoc+ myCount

 CurrLoc= FileList(i)%Endidx 

 close(f2) 

END SUBROUTINE CountTableDetailForEachFile
!-------------------------------------------------------------
SUBROUTINE CountTableDetailForEachFile1(i, fname)

use myGAfuns_dat

 Integer, INTENT(IN):: i 
 character (len=30), INTENT(IN):: fname
 Logical:: myEOF
 Integer:: ierr, IOstatus, myCount1
 character (len=30):: RdLine


 Integer:: f2
 

 f2=freeunit(); open (f2, file=trim(fname), action ='read');Call Fseek(f2,0,0,ierr)
 read(f2, '(a)', IOSTAT=IOstatus) RdLine  ! Read off the header
 
 !CurrLoc = CurrLoc + 1
 myEOF=.false.; myCount1=0
 DO while (myEOF == .false.)        
    read(f2, '(a)', IOSTAT=IOstatus) RdLine
   IF ( IOstatus < 0) Then  
       myEOF = .true.
   ELSE
      myCount1=myCount1+1           
   END IF 
 END DO !While END

 FileList1(i)%RecCount = myCount1
 FileList1(i)%Startidx = CurrLoc1+1
 FileList1(i)%Endidx = CurrLoc1+ myCount1

 CurrLoc1= FileList1(i)%Endidx 

 close(f2) 

END SUBROUTINE CountTableDetailForEachFile1
!-------------------------------------------------------------
SUBROUTINE CountTableDetailForEachFile2(i, fname)

use myGAfuns_dat

 Integer, INTENT(IN):: i 
 character (len=30), INTENT(IN):: fname
 Logical:: myEOF
 Integer:: ierr, IOstatus, myCount2
 character (len=30):: RdLine


 Integer:: f2
 

 f2=freeunit(); open (f2, file=trim(fname), action ='read');Call Fseek(f2,0,0,ierr)
 read(f2, '(a)', IOSTAT=IOstatus) RdLine  ! Read off the header
 
 !CurrLoc = CurrLoc + 1
 myEOF=.false.; myCount2=0
 DO while (myEOF == .false.)        
    read(f2, '(a)', IOSTAT=IOstatus) RdLine
   IF ( IOstatus < 0) Then  
       myEOF = .true.
   ELSE
      myCount2=myCount2+1           
   END IF 
 END DO !While END

 FileList2(i)%RecCount = myCount2
 FileList2(i)%Startidx = CurrLoc2+1
 FileList2(i)%Endidx = CurrLoc2+ myCount2

 CurrLoc2= FileList2(i)%Endidx 

 close(f2) 

END SUBROUTINE CountTableDetailForEachFile2
!-------------------------------------------------------------

!-------------------------------------------------------------
SUBROUTINE CountTableDetailForEachFile3(i, fname)

use myGAfuns_dat

 Integer, INTENT(IN):: i 
 character (len=30), INTENT(IN):: fname
 Logical:: myEOF
 Integer:: ierr, IOstatus, myCount3
 character (len=30):: RdLine


 Integer:: f2
 

 f2=freeunit(); open (f2, file=trim(fname), action ='read');Call Fseek(f2,0,0,ierr)
 read(f2, '(a)', IOSTAT=IOstatus) RdLine  ! Read off the header
 
 !CurrLoc = CurrLoc + 1
 myEOF=.false.; myCount3=0
 DO while (myEOF == .false.)        
    read(f2, '(a)', IOSTAT=IOstatus) RdLine
   IF ( IOstatus < 0) Then  
       myEOF = .true.
   ELSE
      myCount3=myCount3+1           
   END IF 
 END DO !While END

 FileList3(i)%RecCount = myCount3
 FileList3(i)%Startidx = CurrLoc3+1
 FileList3(i)%Endidx = CurrLoc3+ myCount3

 CurrLoc3= FileList3(i)%Endidx 

 close(f2) 

END SUBROUTINE CountTableDetailForEachFile3
!-------------------------------------------------------------

!-------------------------------------------------------------
SUBROUTINE LoadTableDetailForEachFile(i, fname)

use myGAfuns_dat

 Integer, INTENT(IN):: i
 character (len=30), INTENT(IN):: fname
 Logical:: myEOF
 Integer:: ierr, IOstatus, myCount, j, k, k2, k3
 character (len=400):: RdLine, tmpstr 
 real*8:: dDFIRE_read

 Integer:: f2
 
 f2=freeunit(); open (f2, file=trim(fname), action ='read'); Call Fseek(f2,0,0,ierr)
 read(f2, '(a)', IOSTAT=IOstatus) RdLine  ! Read off the header
 
k2=0
DO j=FileList(i)%StartIdx,FileList(i)%EndIdx 
   read(f2,*) k, tmpstr, k2, RecList(j)%TMscore,RecList(j)%ThrDIGARS_Egy,RecList(j)%ASA_Egy,RecList(j)%uPhi_Egy, RecList(j)%uPsi_Egy, RecList(j)%ThrDIGARSThr_Egy,RecList(j)%ASA_Egy_3DIGARS_Data_Ram,RecList(j)%ASA_Egy_SSD_Data_Ram,RecList(j)%ASA_Egy_Spinex_Data_Ram,RecList(j)%ASA_Egy_3DIGARS_Data_Hoque,RecList(j)%ASA_Egy_SSD_Data_Hoque,RecList(j)%ASA_Egy_Spinex_Data_Hoque,RecList(j)%ASA_SS_Egy_3DIGARS_Data_Ram,RecList(j)%ASA_SS_Egy_SSD_Data_Ram,RecList(j)%ASA_SS_Egy_Spinex_Data_Ram,RecList(j)%ASA_SS_Egy_3DIGARS_Data_Hoque,RecList(j)%ASA_SS_Egy_SSD_Data_Hoque,RecList(j)%ASA_SS_Egy_Spinex_Data_Hoque,RecList(j)%Phi_Egy_3DIGARS_Data_Ram,RecList(j)%Phi_Egy_SSD_Data_Ram,RecList(j)%Phi_Egy_Spinex_Data_Ram,RecList(j)%Phi_Egy_3DIGARS_Data_Hoque,RecList(j)%Phi_Egy_SSD_Data_Hoque,RecList(j)%Phi_Egy_Spinex_Data_Hoque,RecList(j)%Psi_Egy_3DIGARS_Data_Ram,RecList(j)%Psi_Egy_SSD_Data_Ram,RecList(j)%Psi_Egy_Spinex_Data_Ram,RecList(j)%Psi_Egy_3DIGARS_Data_Hoque,RecList(j)%Psi_Egy_SSD_Data_Hoque,RecList(j)%Psi_Egy_Spinex_Data_Hoque,RecList(j)%ASAr_Egy_3DIGARS_Data_Ram,RecList(j)%ASAr_Egy_SSD_Data_Ram,RecList(j)%ASAr_Egy_Spinex_Data_Ram,RecList(j)%ASAr_Egy_3DIGARS_Data_Hoque,RecList(j)%ASAr_Egy_SSD_Data_Hoque,RecList(j)%ASAr_Egy_Spinex_Data_Hoque,RecList(j)%ASAr_SS_Egy_3DIGARS_Data_Ram,RecList(j)%ASAr_SS_Egy_SSD_Data_Ram,RecList(j)%ASAr_SS_Egy_Spinex_Data_Ram,RecList(j)%ASAr_SS_Egy_3DIGARS_Data_Hoque,RecList(j)%ASAr_SS_Egy_SSD_Data_Hoque,RecList(j)%ASAr_SS_Egy_Spinex_Data_Hoque,RecList(j)%Phir_Egy_3DIGARS_Data_Ram,RecList(j)%Phir_Egy_SSD_Data_Ram,RecList(j)%Phir_Egy_Spinex_Data_Ram,RecList(j)%Phir_Egy_3DIGARS_Data_Hoque,RecList(j)%Phir_Egy_SSD_Data_Hoque,RecList(j)%Phir_Egy_Spinex_Data_Hoque,RecList(j)%Psir_Egy_3DIGARS_Data_Ram,RecList(j)%Psir_Egy_SSD_Data_Ram,RecList(j)%Psir_Egy_Spinex_Data_Ram,RecList(j)%Psir_Egy_3DIGARS_Data_Hoque,RecList(j)%Psir_Egy_SSD_Data_Hoque,RecList(j)%Psir_Egy_Spinex_Data_Hoque,RecList(j)%ASA_Egy_Comb_Data_Ram,RecList(j)%ASA_SS_Egy_Comb_Data_Ram,RecList(j)%ASA_Egy_Comb_Data_Hoque,RecList(j)%ASA_SS_Egy_Comb_Data_Hoque,RecList(j)%Phi_Egy_Comb_Data_Ram,RecList(j)%Phi_Egy_Comb_Data_Hoque,RecList(j)%Psi_Egy_Comb_Data_Ram,RecList(j)%Psi_Egy_Comb_Data_Hoque,RecList(j)%ASAr_Egy_Comb_Data_Ram,RecList(j)%ASAr_SS_Egy_Comb_Data_Ram,RecList(j)%ASAr_Egy_Comb_Data_Hoque,RecList(j)%ASAr_SS_Egy_Comb_Data_Hoque,RecList(j)%Phir_Egy_Comb_Data_Ram,RecList(j)%Phir_Egy_Comb_Data_Hoque,RecList(j)%Psir_Egy_Comb_Data_Ram,RecList(j)%Psir_Egy_Comb_Data_Hoque,RecList(j)%ASA_Egy_Triplet_Hoque,RecList(j)%Phi_Egy_Triplet_Hoque,RecList(j)%Psi_Egy_Triplet_Hoque,RecList(j)%uPhi_Egy_Hoque,RecList(j)%uPsi_Egy_Hoque
  ! write(*,*) RecList(j)%TMscore,RecList(j)%ThrDIGARS_Egy,RecList(j)%ASA_Egy,RecList(j)%uPhi_Egy,RecList(j)%uPsi_Egy,RecList(j)%ThrDIGARSThr_Egy,RecList(j)%ASA_Egy_3DIGARS_Data_Ram,RecList(j)%ASA_Egy_SSD_Data_Ram,RecList(j)%ASA_Egy_Spinex_Data_Ram,RecList(j)%ASA_Egy_3DIGARS_Data_Hoque,RecList(j)%ASA_Egy_SSD_Data_Hoque,RecList(j)%ASA_Egy_Spinex_Data_Hoque,RecList(j)%ASA_SS_Egy_3DIGARS_Data_Ram,RecList(j)%ASA_SS_Egy_SSD_Data_Ram,RecList(j)%ASA_SS_Egy_Spinex_Data_Ram,RecList(j)%ASA_SS_Egy_3DIGARS_Data_Hoque,RecList(j)%ASA_SS_Egy_SSD_Data_Hoque,RecList(j)%ASA_SS_Egy_Spinex_Data_Hoque,RecList(j)%Phi_Egy_3DIGARS_Data_Ram,RecList(j)%Phi_Egy_SSD_Data_Ram,RecList(j)%Phi_Egy_Spinex_Data_Ram,RecList(j)%Phi_Egy_3DIGARS_Data_Hoque,RecList(j)%Phi_Egy_SSD_Data_Hoque,RecList(j)%Phi_Egy_Spinex_Data_Hoque,RecList(j)%Psi_Egy_3DIGARS_Data_Ram,RecList(j)%Psi_Egy_SSD_Data_Ram,RecList(j)%Psi_Egy_Spinex_Data_Ram,RecList(j)%Psi_Egy_3DIGARS_Data_Hoque,RecList(j)%Psi_Egy_SSD_Data_Hoque,RecList(j)%Psi_Egy_Spinex_Data_Hoque
    
	 IF (k2==1) THEN 
          k3=j
     END IF  
END DO

  !Swap the native to the last 
  IF (k3/= FileList(i)%EndIdx) THEN    ! if K3 is not already last of the list then swap the native to put at the end of the list
   tmpRecList(1) = RecList(k3)
   RecList(k3) = RecList(FileList(i)%EndIdx)
   RecList(FileList(i)%EndIdx)=tmpRecList(1)
  END if 


 close(f2) 

END SUBROUTINE LoadTableDetailForEachFile
!--------------------------------------------------------
!-------------------------------------------------------------
SUBROUTINE LoadTableDetailForEachFile1(i, fname)

use myGAfuns_dat

 Integer, INTENT(IN):: i
 character (len=30), INTENT(IN):: fname
 Logical:: myEOF
 Integer:: ierr, IOstatus, myCount1, j, k, k2, k3
 character (len=400):: RdLine, tmpstr 
 real*8:: dDFIRE_read

 Integer:: f2
 
 f2=freeunit(); open (f2, file=trim(fname), action ='read'); Call Fseek(f2,0,0,ierr)
 read(f2, '(a)', IOSTAT=IOstatus) RdLine  ! Read off the header
 


k2=0
DO j=FileList1(i)%StartIdx,FileList1(i)%EndIdx 
   !read(f2,*) k, tmpstr, k2, RecList1(j)%TMscore,RecList1(j)%ThrDIGARS_Egy,RecList1(j)%ASA_Egy,RecList1(j)%uPhi_Egy, RecList1(j)%uPsi_Egy, RecList1(j)%ThrDIGARSThr_Egy,RecList1(j)%ASA_Egy_3DIGARS_Data_Ram,RecList1(j)%ASA_Egy_SSD_Data_Ram,RecList1(j)%ASA_Egy_Spinex_Data_Ram,RecList1(j)%ASA_Egy_3DIGARS_Data_Hoque,RecList1(j)%ASA_Egy_SSD_Data_Hoque,RecList1(j)%ASA_Egy_Spinex_Data_Hoque,RecList1(j)%ASA_SS_Egy_3DIGARS_Data_Ram,RecList1(j)%ASA_SS_Egy_SSD_Data_Ram,RecList1(j)%ASA_SS_Egy_Spinex_Data_Ram,RecList1(j)%ASA_SS_Egy_3DIGARS_Data_Hoque,RecList1(j)%ASA_SS_Egy_SSD_Data_Hoque,RecList1(j)%ASA_SS_Egy_Spinex_Data_Hoque,RecList1(j)%Phi_Egy_3DIGARS_Data_Ram,RecList1(j)%Phi_Egy_SSD_Data_Ram,RecList1(j)%Phi_Egy_Spinex_Data_Ram,RecList1(j)%Phi_Egy_3DIGARS_Data_Hoque,RecList1(j)%Phi_Egy_SSD_Data_Hoque,RecList1(j)%Phi_Egy_Spinex_Data_Hoque,RecList1(j)%Psi_Egy_3DIGARS_Data_Ram,RecList1(j)%Psi_Egy_SSD_Data_Ram,RecList1(j)%Psi_Egy_Spinex_Data_Ram,RecList1(j)%Psi_Egy_3DIGARS_Data_Hoque,RecList1(j)%Psi_Egy_SSD_Data_Hoque,RecList1(j)%Psi_Egy_Spinex_Data_Hoque
    read(f2,*) k, tmpstr, k2, RecList1(j)%TMscore,RecList1(j)%ThrDIGARS_Egy,RecList1(j)%ASA_Egy,RecList1(j)%uPhi_Egy, RecList1(j)%uPsi_Egy, RecList1(j)%ThrDIGARSThr_Egy,RecList1(j)%ASA_Egy_3DIGARS_Data_Ram,RecList1(j)%ASA_Egy_SSD_Data_Ram,RecList1(j)%ASA_Egy_Spinex_Data_Ram,RecList1(j)%ASA_Egy_3DIGARS_Data_Hoque,RecList1(j)%ASA_Egy_SSD_Data_Hoque,RecList1(j)%ASA_Egy_Spinex_Data_Hoque,RecList1(j)%ASA_SS_Egy_3DIGARS_Data_Ram,RecList1(j)%ASA_SS_Egy_SSD_Data_Ram,RecList1(j)%ASA_SS_Egy_Spinex_Data_Ram,RecList1(j)%ASA_SS_Egy_3DIGARS_Data_Hoque,RecList1(j)%ASA_SS_Egy_SSD_Data_Hoque,RecList1(j)%ASA_SS_Egy_Spinex_Data_Hoque,RecList1(j)%Phi_Egy_3DIGARS_Data_Ram,RecList1(j)%Phi_Egy_SSD_Data_Ram,RecList1(j)%Phi_Egy_Spinex_Data_Ram,RecList1(j)%Phi_Egy_3DIGARS_Data_Hoque,RecList1(j)%Phi_Egy_SSD_Data_Hoque,RecList1(j)%Phi_Egy_Spinex_Data_Hoque,RecList1(j)%Psi_Egy_3DIGARS_Data_Ram,RecList1(j)%Psi_Egy_SSD_Data_Ram,RecList1(j)%Psi_Egy_Spinex_Data_Ram,RecList1(j)%Psi_Egy_3DIGARS_Data_Hoque,RecList1(j)%Psi_Egy_SSD_Data_Hoque,RecList1(j)%Psi_Egy_Spinex_Data_Hoque,RecList1(j)%ASAr_Egy_3DIGARS_Data_Ram,RecList1(j)%ASAr_Egy_SSD_Data_Ram,RecList1(j)%ASAr_Egy_Spinex_Data_Ram,RecList1(j)%ASAr_Egy_3DIGARS_Data_Hoque,RecList1(j)%ASAr_Egy_SSD_Data_Hoque,RecList1(j)%ASAr_Egy_Spinex_Data_Hoque,RecList1(j)%ASAr_SS_Egy_3DIGARS_Data_Ram,RecList1(j)%ASAr_SS_Egy_SSD_Data_Ram,RecList1(j)%ASAr_SS_Egy_Spinex_Data_Ram,RecList1(j)%ASAr_SS_Egy_3DIGARS_Data_Hoque,RecList1(j)%ASAr_SS_Egy_SSD_Data_Hoque,RecList1(j)%ASAr_SS_Egy_Spinex_Data_Hoque,RecList1(j)%Phir_Egy_3DIGARS_Data_Ram,RecList1(j)%Phir_Egy_SSD_Data_Ram,RecList1(j)%Phir_Egy_Spinex_Data_Ram,RecList1(j)%Phir_Egy_3DIGARS_Data_Hoque,RecList1(j)%Phir_Egy_SSD_Data_Hoque,RecList1(j)%Phir_Egy_Spinex_Data_Hoque,RecList1(j)%Psir_Egy_3DIGARS_Data_Ram,RecList1(j)%Psir_Egy_SSD_Data_Ram,RecList1(j)%Psir_Egy_Spinex_Data_Ram,RecList1(j)%Psir_Egy_3DIGARS_Data_Hoque,RecList1(j)%Psir_Egy_SSD_Data_Hoque,RecList1(j)%Psir_Egy_Spinex_Data_Hoque,RecList1(j)%ASA_Egy_Comb_Data_Ram,RecList1(j)%ASA_SS_Egy_Comb_Data_Ram,RecList1(j)%ASA_Egy_Comb_Data_Hoque,RecList1(j)%ASA_SS_Egy_Comb_Data_Hoque,RecList1(j)%Phi_Egy_Comb_Data_Ram,RecList1(j)%Phi_Egy_Comb_Data_Hoque,RecList1(j)%Psi_Egy_Comb_Data_Ram,RecList1(j)%Psi_Egy_Comb_Data_Hoque,RecList1(j)%ASAr_Egy_Comb_Data_Ram,RecList1(j)%ASAr_SS_Egy_Comb_Data_Ram,RecList1(j)%ASAr_Egy_Comb_Data_Hoque,RecList1(j)%ASAr_SS_Egy_Comb_Data_Hoque,RecList1(j)%Phir_Egy_Comb_Data_Ram,RecList1(j)%Phir_Egy_Comb_Data_Hoque,RecList1(j)%Psir_Egy_Comb_Data_Ram,RecList1(j)%Psir_Egy_Comb_Data_Hoque,RecList1(j)%ASA_Egy_Triplet_Hoque,RecList1(j)%Phi_Egy_Triplet_Hoque,RecList1(j)%Psi_Egy_Triplet_Hoque,RecList1(j)%uPhi_Egy_Hoque,RecList1(j)%uPsi_Egy_Hoque
    !write(*,*) RecList1(j)%TMscore,RecList1(j)%ePNP,RecList1(j)%eNPNP   


     IF (k2==1) THEN 
          k3=j
     END IF  
END DO

  !Swap the native to the last 
  IF (k3/= FileList1(i)%EndIdx) THEN    ! if K3 is not already last of the list then swap the native to put at the end of the list
   tmpRecList(1) = RecList1(k3)
   
   !write (*,*) k3, FileList1(i)%EndIdx
   RecList1(k3) = RecList1(FileList1(i)%EndIdx)    
   RecList1(FileList1(i)%EndIdx)=tmpRecList(1)
  END if 
 close(f2) 


END SUBROUTINE LoadTableDetailForEachFile1
!------------------------------------------------------------------
!-------------------------------------------------------------
SUBROUTINE LoadTableDetailForEachFile2(i, fname)

use myGAfuns_dat

 Integer, INTENT(IN):: i
 character (len=30), INTENT(IN):: fname
 Logical:: myEOF
 Integer:: ierr, IOstatus, myCount2, j, k, k2, k3
 character (len=300):: RdLine, tmpstr 
 real*8:: dDFIRE_read

 Integer:: f2
 
 f2=freeunit(); open (f2, file=trim(fname), action ='read'); Call Fseek(f2,0,0,ierr)
 read(f2, '(a)', IOSTAT=IOstatus) RdLine  ! Read off the header
 
k2=0
DO j=FileList2(i)%StartIdx,FileList2(i)%EndIdx 
   read(f2,*) k, tmpstr, k2, RecList2(j)%TMscore,RecList2(j)%ThrDIGARS_Egy,RecList2(j)%ASA_Egy,RecList2(j)%uPhi_Egy, RecList2(j)%uPsi_Egy, RecList2(j)%ThrDIGARSThr_Egy,RecList2(j)%ASA_Egy_3DIGARS_Data_Ram,RecList2(j)%ASA_Egy_SSD_Data_Ram,RecList2(j)%ASA_Egy_Spinex_Data_Ram,RecList2(j)%ASA_Egy_3DIGARS_Data_Hoque,RecList2(j)%ASA_Egy_SSD_Data_Hoque,RecList2(j)%ASA_Egy_Spinex_Data_Hoque,RecList2(j)%ASA_SS_Egy_3DIGARS_Data_Ram,RecList2(j)%ASA_SS_Egy_SSD_Data_Ram,RecList2(j)%ASA_SS_Egy_Spinex_Data_Ram,RecList2(j)%ASA_SS_Egy_3DIGARS_Data_Hoque,RecList2(j)%ASA_SS_Egy_SSD_Data_Hoque,RecList2(j)%ASA_SS_Egy_Spinex_Data_Hoque,RecList2(j)%Phi_Egy_3DIGARS_Data_Ram,RecList2(j)%Phi_Egy_SSD_Data_Ram,RecList2(j)%Phi_Egy_Spinex_Data_Ram,RecList2(j)%Phi_Egy_3DIGARS_Data_Hoque,RecList2(j)%Phi_Egy_SSD_Data_Hoque,RecList2(j)%Phi_Egy_Spinex_Data_Hoque,RecList2(j)%Psi_Egy_3DIGARS_Data_Ram,RecList2(j)%Psi_Egy_SSD_Data_Ram,RecList2(j)%Psi_Egy_Spinex_Data_Ram,RecList2(j)%Psi_Egy_3DIGARS_Data_Hoque,RecList2(j)%Psi_Egy_SSD_Data_Hoque,RecList2(j)%Psi_Egy_Spinex_Data_Hoque,RecList2(j)%ASAr_Egy_3DIGARS_Data_Ram,RecList2(j)%ASAr_Egy_SSD_Data_Ram,RecList2(j)%ASAr_Egy_Spinex_Data_Ram,RecList2(j)%ASAr_Egy_3DIGARS_Data_Hoque,RecList2(j)%ASAr_Egy_SSD_Data_Hoque,RecList2(j)%ASAr_Egy_Spinex_Data_Hoque,RecList2(j)%ASAr_SS_Egy_3DIGARS_Data_Ram,RecList2(j)%ASAr_SS_Egy_SSD_Data_Ram,RecList2(j)%ASAr_SS_Egy_Spinex_Data_Ram,RecList2(j)%ASAr_SS_Egy_3DIGARS_Data_Hoque,RecList2(j)%ASAr_SS_Egy_SSD_Data_Hoque,RecList2(j)%ASAr_SS_Egy_Spinex_Data_Hoque,RecList2(j)%Phir_Egy_3DIGARS_Data_Ram,RecList2(j)%Phir_Egy_SSD_Data_Ram,RecList2(j)%Phir_Egy_Spinex_Data_Ram,RecList2(j)%Phir_Egy_3DIGARS_Data_Hoque,RecList2(j)%Phir_Egy_SSD_Data_Hoque,RecList2(j)%Phir_Egy_Spinex_Data_Hoque,RecList2(j)%Psir_Egy_3DIGARS_Data_Ram,RecList2(j)%Psir_Egy_SSD_Data_Ram,RecList2(j)%Psir_Egy_Spinex_Data_Ram,RecList2(j)%Psir_Egy_3DIGARS_Data_Hoque,RecList2(j)%Psir_Egy_SSD_Data_Hoque,RecList2(j)%Psir_Egy_Spinex_Data_Hoque,RecList2(j)%ASA_Egy_Comb_Data_Ram,RecList2(j)%ASA_SS_Egy_Comb_Data_Ram,RecList2(j)%ASA_Egy_Comb_Data_Hoque,RecList2(j)%ASA_SS_Egy_Comb_Data_Hoque,RecList2(j)%Phi_Egy_Comb_Data_Ram,RecList2(j)%Phi_Egy_Comb_Data_Hoque,RecList2(j)%Psi_Egy_Comb_Data_Ram,RecList2(j)%Psi_Egy_Comb_Data_Hoque,RecList2(j)%ASAr_Egy_Comb_Data_Ram,RecList2(j)%ASAr_SS_Egy_Comb_Data_Ram,RecList2(j)%ASAr_Egy_Comb_Data_Hoque,RecList2(j)%ASAr_SS_Egy_Comb_Data_Hoque,RecList2(j)%Phir_Egy_Comb_Data_Ram,RecList2(j)%Phir_Egy_Comb_Data_Hoque,RecList2(j)%Psir_Egy_Comb_Data_Ram,RecList2(j)%Psir_Egy_Comb_Data_Hoque,RecList2(j)%ASA_Egy_Triplet_Hoque,RecList2(j)%Phi_Egy_Triplet_Hoque,RecList2(j)%Psi_Egy_Triplet_Hoque,RecList2(j)%uPhi_Egy_Hoque,RecList2(j)%uPsi_Egy_Hoque
   !write(*,*) RecList2(j)%TMscore,RecList2(j)%ePNP,RecList2(j)%eNPNP   
     IF (k2==1) THEN 
          k3=j
     END IF  
END DO

  !Swap the native to the last 
  IF (k3/= FileList2(i)%EndIdx) THEN    ! if K3 is not already last of the list then swap the native to put at the end of the list
   tmpRecList(1) = RecList2(k3)
   RecList2(k3) = RecList2(FileList2(i)%EndIdx)
   RecList2(FileList2(i)%EndIdx)=tmpRecList(1)
  END if 
 close(f2) 

END SUBROUTINE LoadTableDetailForEachFile2

!-------------------------------------------------------------
SUBROUTINE LoadTableDetailForEachFile3(i, fname)

use myGAfuns_dat

 Integer, INTENT(IN):: i
 character (len=30), INTENT(IN):: fname
 Logical:: myEOF
 Integer:: ierr, IOstatus, myCount3, j, k, k2, k3
 character (len=300):: RdLine, tmpstr 
 real*8:: dDFIRE_read

 Integer:: f2
 
 f2=freeunit(); open (f2, file=trim(fname), action ='read'); Call Fseek(f2,0,0,ierr)
 read(f2, '(a)', IOSTAT=IOstatus) RdLine  ! Read off the header
 
k2=0
DO j=FileList3(i)%StartIdx,FileList3(i)%EndIdx 
   read(f2,*) k, tmpstr, k2, RecList3(j)%TMscore,RecList3(j)%ThrDIGARS_Egy,RecList3(j)%ASA_Egy,RecList3(j)%uPhi_Egy, RecList3(j)%uPsi_Egy, RecList3(j)%ThrDIGARSThr_Egy,RecList3(j)%ASA_Egy_3DIGARS_Data_Ram,RecList3(j)%ASA_Egy_SSD_Data_Ram,RecList3(j)%ASA_Egy_Spinex_Data_Ram,RecList3(j)%ASA_Egy_3DIGARS_Data_Hoque,RecList3(j)%ASA_Egy_SSD_Data_Hoque,RecList3(j)%ASA_Egy_Spinex_Data_Hoque,RecList3(j)%ASA_SS_Egy_3DIGARS_Data_Ram,RecList3(j)%ASA_SS_Egy_SSD_Data_Ram,RecList3(j)%ASA_SS_Egy_Spinex_Data_Ram,RecList3(j)%ASA_SS_Egy_3DIGARS_Data_Hoque,RecList3(j)%ASA_SS_Egy_SSD_Data_Hoque,RecList3(j)%ASA_SS_Egy_Spinex_Data_Hoque,RecList3(j)%Phi_Egy_3DIGARS_Data_Ram,RecList3(j)%Phi_Egy_SSD_Data_Ram,RecList3(j)%Phi_Egy_Spinex_Data_Ram,RecList3(j)%Phi_Egy_3DIGARS_Data_Hoque,RecList3(j)%Phi_Egy_SSD_Data_Hoque,RecList3(j)%Phi_Egy_Spinex_Data_Hoque,RecList3(j)%Psi_Egy_3DIGARS_Data_Ram,RecList3(j)%Psi_Egy_SSD_Data_Ram,RecList3(j)%Psi_Egy_Spinex_Data_Ram,RecList3(j)%Psi_Egy_3DIGARS_Data_Hoque,RecList3(j)%Psi_Egy_SSD_Data_Hoque,RecList3(j)%Psi_Egy_Spinex_Data_Hoque,RecList3(j)%ASAr_Egy_3DIGARS_Data_Ram,RecList3(j)%ASAr_Egy_SSD_Data_Ram,RecList3(j)%ASAr_Egy_Spinex_Data_Ram,RecList3(j)%ASAr_Egy_3DIGARS_Data_Hoque,RecList3(j)%ASAr_Egy_SSD_Data_Hoque,RecList3(j)%ASAr_Egy_Spinex_Data_Hoque,RecList3(j)%ASAr_SS_Egy_3DIGARS_Data_Ram,RecList3(j)%ASAr_SS_Egy_SSD_Data_Ram,RecList3(j)%ASAr_SS_Egy_Spinex_Data_Ram,RecList3(j)%ASAr_SS_Egy_3DIGARS_Data_Hoque,RecList3(j)%ASAr_SS_Egy_SSD_Data_Hoque,RecList3(j)%ASAr_SS_Egy_Spinex_Data_Hoque,RecList3(j)%Phir_Egy_3DIGARS_Data_Ram,RecList3(j)%Phir_Egy_SSD_Data_Ram,RecList3(j)%Phir_Egy_Spinex_Data_Ram,RecList3(j)%Phir_Egy_3DIGARS_Data_Hoque,RecList3(j)%Phir_Egy_SSD_Data_Hoque,RecList3(j)%Phir_Egy_Spinex_Data_Hoque,RecList3(j)%Psir_Egy_3DIGARS_Data_Ram,RecList3(j)%Psir_Egy_SSD_Data_Ram,RecList3(j)%Psir_Egy_Spinex_Data_Ram,RecList3(j)%Psir_Egy_3DIGARS_Data_Hoque,RecList3(j)%Psir_Egy_SSD_Data_Hoque,RecList3(j)%Psir_Egy_Spinex_Data_Hoque,RecList3(j)%ASA_Egy_Comb_Data_Ram,RecList3(j)%ASA_SS_Egy_Comb_Data_Ram,RecList3(j)%ASA_Egy_Comb_Data_Hoque,RecList3(j)%ASA_SS_Egy_Comb_Data_Hoque,RecList3(j)%Phi_Egy_Comb_Data_Ram,RecList3(j)%Phi_Egy_Comb_Data_Hoque,RecList3(j)%Psi_Egy_Comb_Data_Ram,RecList3(j)%Psi_Egy_Comb_Data_Hoque,RecList3(j)%ASAr_Egy_Comb_Data_Ram,RecList3(j)%ASAr_SS_Egy_Comb_Data_Ram,RecList3(j)%ASAr_Egy_Comb_Data_Hoque,RecList3(j)%ASAr_SS_Egy_Comb_Data_Hoque,RecList3(j)%Phir_Egy_Comb_Data_Ram,RecList3(j)%Phir_Egy_Comb_Data_Hoque,RecList3(j)%Psir_Egy_Comb_Data_Ram,RecList3(j)%Psir_Egy_Comb_Data_Hoque,RecList3(j)%ASA_Egy_Triplet_Hoque,RecList3(j)%Phi_Egy_Triplet_Hoque,RecList3(j)%Psi_Egy_Triplet_Hoque,RecList3(j)%uPhi_Egy_Hoque,RecList3(j)%uPsi_Egy_Hoque
   !write(*,*) RecList3(j)%TMscore,RecList3(j)%ePNP,RecList3(j)%eNPNP   
     IF (k2==1) THEN 
          k3=j
     END IF  
END DO

  !Swap the native to the last 
  IF (k3/= FileList3(i)%EndIdx) THEN    ! if K3 is not already last of the list then swap the native to put at the end of the list
   tmpRecList(1) = RecList3(k3)
   RecList3(k3) = RecList3(FileList3(i)%EndIdx)
   RecList3(FileList3(i)%EndIdx)=tmpRecList(1)
  END if 
 close(f2) 

END SUBROUTINE LoadTableDetailForEachFile3


!------------ main sub-routine -------------------------
SUBROUTINE ga_run
 use myGAfuns_dat

 Integer:: i, j, f100
 real*8:: random1
 character(len=5000)::SaveStr  
  
   ! load the data from set of .CVS to calculate the fitness
   ! function_xyz
  
f100=freeunit(); open (f100, file='output/Results.csv', action ='write'); Call Fseek(f100,0,0,ierr) 
write(f100, *) 'Gen, TotalFitness, PCC_Avg, TMScoreLowEnMdl_Avg,PCC_Avg1, TMScoreLowEnMdl_Avg1 ... , K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, K24, K25, K26, K27 ... k67'


    Call PopInitialization 
    Call SortPopulation(1)
    Call Initialize_AM_Pop

    generation = 0 
    CurrentPosNewPopulation = 0

    DO WHILE (generation < Max_generation)
         generation = generation + 1            
         CurrentPosNewPopulation = Elitist() 
         Call ComputeFitnessSum_WorstFitness       
 
         !! == START of Crossover related issues =======================================
          CrossOverCount = PopSize - (ElitRate * PopSize)
          CrossOverCount = CrossOverCount * CrossOverRate
          CrossOverCount = CurrentPosNewPopulation + CrossOverCount
        
          ! Crossover Loop next:
          DO WHILE ( CurrentPosNewPopulation <= CrossOverCount)
               i = SelectChromosome(generation*CrossOverCount)  
               j = SelectChromosome(i)

               !n = ((RandNumTam(j) * (Chromosomelength - 1))+1)
               n = ((Random1() * (Chromosomelength - 1))+1) 
               Call My_AM_Crossover(i, j, n)
               CurrentPosNewPopulation =CurrentPosNewPopulation+2

          END DO  
          !! == END of Crossover related issues =======================================
          
          !! Fill in the gap by random structure if there is a gap, basically Initialization (from 'CurrentPosNewPopulation +1' to 'PopSize')   
          DO while (CurrentPosNewPopulation <= PopSize)
                   CurrentPosNewPopulation = CurrentPosNewPopulation + 1
                   Call random_chromosome (CurrentPosNewPopulation, 2)    !2 => ptrPop2 (i.e the new population)   
          END DO

          !! Mutation should be on the >>NEW POP<< at a lower rate (since otherwise will have to replace, which is meaning less)
          m_mut = 0
             MutationCount = PopSize - (ElitRate * PopSize)
             MutationCount = MutationCount * MutationRate 
           
         !! Mutation Loop NEXT:   
         DO While (m_mut <= (MutationCount))
       
             i = SelectChromosome(m_mut*generation)    ! (...) => just to vary the random-seed incase it is a quick call
             !n = (RandNumTam(i) * Chromosomelength)+1  ! For Len =30, it is [0..29]+1 or, [1..30]
             n = (Random1() * Chromosomelength)+1  

             Call Mutation(i, n) 
             m_mut=m_mut+1   
           END DO
          
           !Instead of copying the new_pop as old_pop I am swaping the pointers only, which is basically equivalent to => CopyNewPopToPop
           Call Swap_Population_Pointer ! This ptrPop1 now having the current population 
          
           DO j = 1, PopSize
             ptrPop1(j)%Fitness = CalculateFitness(j,1)   ! If needed this loop
           END DO

           Call SortPopulation(1)
		   
		   ptrPop1(1)%Fitness = CalculateFitness(1,1)
           
     ! SaveStr = trim(IntegertoStr(generation))//','//trim(RealtoStr(ptrPop1(1)%Fitness))//','//trim(RealtoStr(ptrPop1(1)%PCC))//','//trim(RealtoStr(ptrPop1(1)%TMLowEMdl))//','//trim(RealtoStr(ptrPop1(1)%PCC1))//','//trim(RealtoStr(ptrPop1(1)%TMLowEMdl1))//','//trim(RealToStr(Ks(1)))//','//trim(RealToStr(Ks(2)))//','//trim(RealToStr(Ks(3)))//','//trim(RealToStr(Ks(4)))//','//trim(RealToStr(Ks(5)))//','//trim(RealToStr(Ks(6)))//','//trim(RealToStr(Ks(7)))//','//trim(RealToStr(Ks(8)))//','//trim(RealToStr(Ks(9)))//','//trim(RealToStr(Ks(10)))//','//trim(RealToStr(Ks(11)))//','//trim(RealToStr(Ks(12)))//','//trim(RealToStr(Ks(13)))//','//trim(RealToStr(Ks(14)))//','//trim(RealToStr(Ks(15)))//','//trim(RealToStr(Ks(16)))//','//trim(RealToStr(Ks(17)))//','//trim(RealToStr(Ks(18)))//','//trim(RealToStr(Ks(19)))//','//trim(RealToStr(Ks(20)))//','//trim(RealToStr(Ks(21)))//','//trim(RealToStr(Ks(22)))//','//trim(RealToStr(Ks(23)))//','//trim(RealToStr(Ks(24)))//','//trim(RealToStr(Ks(25)))//','//trim(RealToStr(Ks(26)))//','//trim(RealToStr(Ks(27)))
	   SaveStr = trim(IntegertoStr(generation))//','//trim(RealtoStr(ptrPop1(1)%Fitness))//','//trim(RealtoStr(ptrPop1(1)%PCC))//','//trim(RealtoStr(ptrPop1(1)%TMLowEMdl))//','//trim(RealtoStr(ptrPop1(1)%Native_Count))//','//trim(RealtoStr(ptrPop1(1)%ZSCORE))//','//trim(RealtoStr(ptrPop1(1)%PCC1))//','//trim(RealtoStr(ptrPop1(1)%TMLowEMdl1))//','//trim(RealtoStr(ptrPop1(1)%Native_Count1))//','//trim(RealtoStr(ptrPop1(1)%ZSCORE1))//','//trim(RealtoStr(ptrPop1(1)%PCC2))//','//trim(RealtoStr(ptrPop1(1)%TMLowEMdl2))//','//trim(RealtoStr(ptrPop1(1)%Native_Count2))//','//trim(RealtoStr(ptrPop1(1)%ZSCORE2))//','//trim(RealtoStr(ptrPop1(1)%PCC3))//','//trim(RealtoStr(ptrPop1(1)%TMLowEMdl3))//','//trim(RealtoStr(ptrPop1(1)%Native_Count3))//','//trim(RealtoStr(ptrPop1(1)%ZSCORE3))//','//trim(RealToStr(Ks(1)))//','//trim(RealToStr(Ks(2)))//','//trim(RealToStr(Ks(3)))//','//trim(RealToStr(Ks(4)))//','//trim(RealToStr(Ks(5)))
          
      !     write (*,"(A)",advance="no") trim(SaveStr) 
      !     write (*,"(A)",advance="yes") ''

           write (f100,"(A)",advance="no") trim(SaveStr)
           write (f100,"(A)",advance="yes") ''
           call flush(f100)   

    END DO ! While END

  close(f100)

END SUBROUTINE ga_run



