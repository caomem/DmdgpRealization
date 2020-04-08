using LinearAlgebra

Rx(θ::Real) = [1 0 0 0; 
        0 cos(θ) -sin(θ) 0; 
        0 sin(θ) cos(θ) 0;
        0 0 0 1]

Ry(θ::Real) = [cos(θ) 0 -sin(θ) 0; 
        0 1 0 0; 
        sin(θ) 0 cos(θ) 0;
        0 0 0 1]

Rz(θ::Real) = [cos(θ) -sin(θ) 0 0;
        sin(θ) cos(θ) 0 0;
        0 0 1 0;
        0 0 0 1]

Tx(d::Real) = [1 0 0 d;
        0 1 0 0;
        0 0 1 0;
        0 0 0 1]

π = pi
d = 2
angToRad(ang::Real) = (π/180)*ang
θ(ang::Real = 150) = angToRad(ang) #70
ω(ang::Real = 160) = angToRad(ang)

function bondangle(i,D)#i=3,...,n
	c = (-D[i-2,i]^2+ D[i-1,i]^2+ D[i-2,i-1]^2)/(2.0*D[i-1,i]*D[i-2,i-1])
	if c<-1.0
		c=-1.0
	end
	if c>1.0
		c=1.0
	end
	s = sqrt(1.0-c^2)
	return c,s
end

function torsionangle(i,D)#i=4,...,n
	d12=D[i-3,i-2]
	d13=D[i-3,i-1]
	d14=D[i-3,i]
	d23=D[i-2,i-1]
	d24=D[i-2,i]
	d34=D[i-1,i]
	a = d12*d12 + d24*d24 - d14*d14
	a = a/(2.0*d12*d24)
	b = d24*d24 + d23*d23 - d34*d34        
	b = b / (2.0*d24*d23)
	c = d12*d12 + d23*d23 - d13*d13
	c = c / (2.0*d12*d23)
	e = 1.0 - b^2;
	f = 1.0 - c^2;
	if (e < 0.0 || f < 0.0)  
		return -2
	end
	e = sqrt(e)
	f = sqrt(f)
	valc = (a - b*c)/(e*f)
	if (valc < -1.0)  
		valc = -1.0
	end
	if (valc >  1.0)  
		valc =  1.0
	end
	vals=sqrt(1.0-valc^2)
	return valc,vals
end

D = [0 1.526 2.4918 2.99897 0 4.44352;1.526 0 1.525426 2.49207 1.526 3.828471; 2.4918 1.525426 0 1.52604 2.4918 2.49212; 2.998972 2.49207 1.52604 0 2.998972 1.525877;0 1.526 2.4918 2.99897 0 4.44352;4.44352 3.82847 2.49212 1.52587 4.44352 0]

B1 = I(4) #nada a fazer

B2 =  Ry(π)*Tx(d)

B3 = Rx(π)*Rz(θ())*Ry(π)*Tx(d)

Bi(ω = ω(), θ = θ(), d = d) = Rx(ω)*Rx(π)*Rz(θ)*Ry(π)*Tx(d)

Bi′(B) = B*inv(Tx(d))

unk = [1 0 0 0;
        0 1 0 0;
        0 0 1 0;
        1 1 1 1]

function sysGenerate(B, α = 1)
        return ["x/|x|" "y/|y|" "z/|z|" "i"; B*unk*α]
end

function flatProjection(A, x̂ = [cos(π/4) -cos(π/4)]'*0.5)
	return [A[2:3,1]+x̂*A[1,1] A[2:3,2]+x̂*A[1,2] A[2:3,3]+x̂*A[1,3] A[2:3,4]+x̂*A[1,4]]
end

#rotationFactor = Rz(π/2)*Rx(π/2)*Rz(-θ(30))*Rx(π)

print("DmdgpRealization included\n")
#sysGenerate(flatProjection(rotationFactor*B2*B3*Bi(ω(230),θ(120))), 20)
#Bo = B2*B3*Bi(ω(230),θ(120))
#B = Bi(ω(230),θ(120))
#sysGenerate(flatProjection(rotationFactor*Bo*B*B*B), 20)