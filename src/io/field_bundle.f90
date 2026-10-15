module m_field_bundle
    !! Reverse C-order dimensions once, preserving precisely the ctypes buffer.
    use, intrinsic :: iso_fortran_env, only: real64, int64
    use m_batch_support
    use m_npy, only: read_npy
    use m_namelist, only: read_namelist
    use m_allcom, only: nx, ny, nz, nspec, qm
    implicit none
    private
    public :: load_field_bundle
contains
    subroutine load_field_bundle(directory, ispec, use_electric_field, use_magnetic_field, fields, dimensions)
        character(*), intent(in) :: directory
        integer, intent(in) :: ispec
        logical, intent(in) :: use_electric_field, use_magnetic_field
        real(real64), allocatable, intent(out) :: fields(:, :, :, :)
        integer, intent(out) :: dimensions(3)
        character(:), allocatable :: metadata, species
        real(real64), allocatable :: flat(:)
        integer(int64), allocatable :: npy_shape(:)
        integer :: p, first, value, ios, i
        logical :: found
        metadata = read_text(path_join(directory, 'bundle.json'))
        if (json_string(metadata, 'schema') /= 'vdsolverf.field_bundle' .or. &
            json_string(metadata, 'version') /= '0.1.0') call batch_fail('unsupported field bundle schema/version')
        dimensions = [json_integer(metadata, 'nx'), json_integer(metadata, 'ny'), json_integer(metadata, 'nz')]
        if (any(dimensions < 1) .or. any(dimensions == huge(1))) call batch_fail('invalid bundle dimensions')
        if (use_electric_field .neqv. json_logical(metadata, 'use_electric_field')) &
            call batch_fail('use_electric_field must match the field bundle')
        if (use_magnetic_field .neqv. json_logical(metadata, 'use_magnetic_field')) &
            call batch_fail('use_magnetic_field must match the field bundle')
        species = json_member(metadata, 'ispecs')
        if (species(1:1) /= '[' .or. species(len(species):) /= ']') call batch_fail('invalid bundle ispecs')
        p = 2
        found = .false.
        do
            do while (p < len(species))
                if (index(' '//achar(9)//achar(10)//achar(13), species(p:p)) == 0) exit
                p = p + 1
            end do
            if (p == len(species)) exit
            first = p
            do while (p < len(species))
                if (index('0123456789', species(p:p)) == 0) exit
                p = p + 1
            end do
            if (first == p) call batch_fail('invalid bundle species index')
            read(species(first:p-1), *, iostat=ios) value
            if (ios /= 0) call batch_fail('invalid bundle species index')
            if (value == ispec) found = .true.
            do while (p < len(species))
                if (index(' '//achar(9)//achar(10)//achar(13), species(p:p)) == 0) exit
                p = p + 1
            end do
            if (p == len(species)) exit
            if (species(p:p) /= ',') call batch_fail('invalid bundle species separator')
            p = p + 1
        end do
        if (.not. found) call batch_fail('species '//number(ispec)//' is missing from the field bundle')
        call read_namelist(path_join(directory, 'plasma-vdsolverf.inp'))
        if (any(dimensions /= [nx, ny, nz])) call batch_fail('bundle and namelist dimensions differ')
        if (ispec < 0 .or. ispec >= nspec .or. ispec >= size(qm)) &
            call batch_fail('ispec is outside the namelist species range')
        call read_npy(path_join(directory, 'ebvalues_ispec'//number(ispec)//'.npy'), flat, npy_shape)
        if (size(npy_shape) /= 4) call batch_fail('field npy must have shape (nz+1,ny+1,nx+1,9)')
        if (any(npy_shape /= [int(dimensions(3:1:-1), int64)+1, 9_int64])) &
            call batch_fail('field npy shape disagrees with bundle dimensions')
        allocate(fields(9, nx+1, ny+1, nz+1), stat=i)
        if (i /= 0) call batch_fail('cannot allocate field bundle')
        fields = reshape(flat, shape(fields))
    end subroutine
end module
