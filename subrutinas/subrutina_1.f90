subroutine algo(a,b,c)
	integer,parameter::pcs=8
	real(kind=pcs),intent(in)::a,b
	real(kind=pcs),intent(out)::c
	
	c=a+b
	
end subroutine
