program main
   
   use myGAfuns_dat

   Integer:: i, n

   real *8:: random1
  
   write(*,*) 'Starting GA to optimize parameters'  

   Call LoadTheDataFromInput    ! CASP8 data
   Call LoadTheDataFromInput1   ! CASP9 data
   Call LoadTheDataFromInput2   ! CASP10 data
   Call LoadTheDataFromInput3   ! CASP11 data

   Call ga_run
   
   write(*,*) 'GA based parameter optimization is done ... '

end program

