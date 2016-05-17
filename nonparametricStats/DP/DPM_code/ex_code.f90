subroutine gibbs(npost, nthin, ps)
	integer,intent(in) :: npost, nthin
	real*8,intent(out) :: ps(npost,2)
	integer :: i,j
	real*8 :: x,y,gamrnd,normrnd

	call rndstart() !GetRNGstate()
	do i=1,npost
		do j=1,nthin
			x=gamrnd(3.d0,1.d0/(y*y+4))
			y=normrnd(1.d0/(x+1),1.d0/dsqrt(2*x+2))
		end do
		ps(i,1)=x; ps(i,2)=y
	end do
	call rndend() !PutRNGstate()
end subroutine gibbs