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

B1 = I(4) #nada a fazer

#B2 =  Ry(π)*Tx(d)

B2 = zeros(4,4)
B2[1,1] = -1.0
B2[2,2] = 1.0
B2[3,3] = -1.0
B2[4,4] = 1.0
B2[1,4] = -d

#B3 = Rx(π)*Rz(θ())*Ry(π)*Tx(d)

# tird atom
D12 = d
D13 = 3.8637
D23 = 2
D14 = 0.0
D24 = 0.0
D34 = 0.0
cθ,sθ = (cos(150),sin(150))
cω,sω = (0.0,0.0)
B3 = zeros(4,4)
B3[1,1] = -cθ
B3[1,2] = -sθ
B3[1,4] = -D23*cθ
B3[2,1] = sθ
B3[2,2] = -cθ
B3[2,4] = D23*sθ
B3[3,3] = 1.0
B3[4,4] = 1.0

function torsionmatrix(cosθ,sinθ,cosω,sinω,d34,sign::Bool)
	if sign == true
		
		B=zeros(4,4)
		B[1,1] = -cosθ
		B[1,2] = -sinθ
		B[1,4] = -d34*cosθ
		B[2,1] = sinθ*cosω
		B[2,2] = -cosθ*cosω
		B[2,3] = -sinω
		B[2,4] = d34*sinθ*cosω
		B[3,1] = sinθ*sinω
		B[3,2] = -cosθ*sinω
		B[3,3] = cosω
		B[3,4] = d34*sinθ*sinω 
		B[4,4] = 1
	else
		
		B=zeros(4,4)
		B[1,1] = -cosθ
		B[1,2] = -sinθ
		B[1,4] = -d34*cosθ
		B[2,1] = sinθ*cosω
		B[2,2] = -cosθ*cosω
		B[2,3] = sinω
		B[2,4] = d34*sinθ*cosω
		B[3,1] = -sinθ*sinω
		B[3,2] = cosθ*sinω
		B[3,3] = cosω
		B[3,4] = -d34*sinθ*sinω 
		B[4,4] = 1
		
	end
	return B
end

Bi(ω = ω(), θ = θ(), d = d) = torsionmatrix(cos(θ), sin(θ), cos(ω), sin(ω), d, true) #Bi(ω = ω(), θ = θ(), d = d) = Rx(ω)*Rx(π)*Rz(θ)*Ry(π)*Tx(d)

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