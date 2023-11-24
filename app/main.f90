program main
  use iso_fortran_env, only: pr => real64
  use hard_spheres_mix, only: fugacity, alloc
  implicit none

  integer   :: nc
  integer   :: i, n_deriv, t_deriv, p_deriv, phase_type

  integer   :: funit_input, funit_output
  character :: TAB = achar(9)   

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
    call run
  end subroutine

  subroutine run
 
    integer, parameter :: funit_data = 10, L = 2, V = 1
 
    ! Input variables
    real(pr) :: beta = 0, p, t
    real(pr) :: n(nc)

    ! Output variables
    real(pr) :: lnphi(nc), dlnphi_dn(nc, nc), dlnphi_dt(nc), dlnphi_dp(nc)

    integer :: comp1, comp2, i, ndata, point
    logical :: gas = .FALSE. , heavy = .FALSE., trace = .TRUE.
    character(20) :: data_author(:), data_publ(:), data_year(:), dat_vol(:), data_pag(:)

    open (unit = funit_data, file="data.tsv")
    read (funit_data, *) number_sets
    i = 0
    do while (iostat_i == 0)
    
      i = i + 1
      write (itext, '(I2)') i
      open (unit=funit_set, file="set_"//itext//".tsv")
      read (funit_data, *, i) comp1, comp2
      read (funit_data, *) data_authors(ndata), data_publ(ndata), data_year(ndata), data_vol(ndata), data_pag(ndata)
      z(:) = 0
      feed(comp1) = x(j)
      trace(comp1) = .FALSE.
      feed(comp2) = 1 - x(j)
      trace(comp2) = .FALSE.
        
      do while (iostat_j == 0)
        
        read (funit_data, '(4G)', iostat=ios) T(j), P(j), x(j), y(j)
        
        if (iostat_j /= 0) exit
        
        call Flash (model, NC, calc_type, T, Pcalc, feed, trace, gas, heavy, beta, error_thermo, K, w, Z, var_spec, dXi_dS, lnPhi_out)
        v_ig = Rgas*T/P
        write (funit_set, *) Pexp, x, y, Pcalc, TAB, w(comp1,L), w(comp2,L), w(comp1,V), w(comp2,V), Z(L)*v_ig, Z(V)*v_IG
          
      enddo
      close (unit=funit_set)
  
    enddo
  end subroutine
end program
