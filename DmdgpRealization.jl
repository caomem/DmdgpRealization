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