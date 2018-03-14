	program China
	implicit none

	Integer i, j
	Integer, parameter :: nla=10801 ! latitude
	Integer, parameter :: nlo=21601 ! longtitude
	Real,	 dimension(nla, nlo) :: a, b, height
	
! Open and read
	open(51, file='ETOPO1_Ice_g_int.xyz')
	write(*,*) "...... reading all data ......"
	do j=1, nla
		do i=1, nlo
			read(51,*) a(j,i),b(j,i),height(j,i)
		enddo
	enddo 
	close(51)

	write(*,*) "...... all data are ready ......"

!Open and write
	open(52,file='China.xyz')
	write(*,*) "...... Save China data ......"
	do j=2101, 5401
		do i=15001, 19201
	!the value of i, j are shown in Draft.nb file.
			write(52,104) a(j,i),b(j,i),height(j,i)
		enddo
	enddo 
	close(52)



104 format(f8.4,f8.4,f9.2)
	
	write(*,*) "...... Well done!!! ......"
	end program China