program main
  use iso_fortran_env, only: pr => real64
  use hard_spheres_mix, only: alloc
  implicit none

  integer   :: nc
  integer   :: i, j

  integer   :: funit_input, funit_output

  call read_file("test.txt")
  call run
 
 contains
  subroutine read_file(filename)
    use hard_spheres_mix, only: read_eos_parameters, alloc
    character(len=*), intent(in) :: filename

    open (newunit=funit_input, file=trim(filename))
    read (funit_input, *) model, nc
    call alloc(nc)
    call ReadParameters (funit_input, funit_output, model, NC)
    close (unit=funit_input)

  end subroutine

  subroutine run
 
    integer, parameter :: funit_data = 10, L = 2, V = 1
 
    ! Input variables
    real(pr) :: beta = 0, feed(NC), K(NC), Pcalc, Pexp, t, v_ig, w(NC,2), x, y, Z(2)


    integer :: comp1, comp2, iostat_i, iostat_j, ndata, point
    logical, dimension(NC) :: gas = .FALSE. , heavy = .FALSE., trace = .TRUE.
    character     :: calc_type = "T", TAB = achar(9)
    character(2)  :: itext
    character(20) :: data_author(:), data_publ(:), data_year(:), dat_vol(:), data_pag(:)
    
    open (unit = funit_data, file="data.tsv")
    read (funit_data, *) number_sets
    i = 0
    do while (iostat_i == 0)
    
      i = i + 1
      write (itext, '(I2)') i
      open (unit=funit_set, file="set_"//itext//".tsv")
      read (funit_data, '(2I)', iostat=iostat_i) comp1, comp2
	if (iostat_i /= 0) exit      
      read (funit_data, *) data_authors(ndata), data_publ(ndata), data_year(ndata), data_vol(ndata), data_pag(ndata)
      feed(:) = 0
      feed(comp1) = x(j)
      trace(comp1) = .FALSE.
      feed(comp2) = 1 - x(j)
      trace(comp2) = .FALSE.
        
      do while (iostat_j == 0)
        
        read (funit_data, '(4G)', iostat=iostat_j) T, Pexp, x, y
        
        if (iostat_j /= 0) exit
        
        call Flash (model, NC, calc_type, T, Pcalc, feed, trace, gas, heavy, beta, error_thermo, K, w, Z, var_spec, dXi_dS, lnPhi_out)
        v_ig = Rgas*T/Pcalc
        write (funit_set, *) Pexp, TAB, x, TAB, y, TAB, Pcalc, TAB, w(comp1,L), TAB, w(comp2,L), TAB, w(comp1,V), TAB, w(comp2,V), TAB, Z(L)*v_ig, TAB, Z(V)*v_IG
          
      enddo
      close (unit=funit_set)
  
    enddo
  end subroutine
end program
