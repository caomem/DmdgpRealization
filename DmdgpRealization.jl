using LinearAlgebra, Plots

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
cθ,sθ = (cos(θ()),sin(θ()))
B3 = zeros(4,4)
B3[1,1] = -cθ
B3[1,2] = -sθ
B3[1,4] = -d*cθ
B3[2,1] = sθ
B3[2,2] = -cθ
B3[2,4] = d*sθ
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

Bi(ω = ω(), θ = θ(), d = d) = torsionmatrix(cos(angToRad(θ)), sin(angToRad(θ)), cos(angToRad(ω)), sin(angToRad(ω)), d, false) #Bi(ω = ω(), θ = θ(), d = d) = Rx(ω)*Rx(π)*Rz(θ)*Ry(π)*Tx(d)

Bi′(B) = B*inv(Tx(d))

unk = [1 0 0 0;
        0 1 0 0;
        0 0 1 0;
        1 1 1 1]

function sysGenerate(B, α = 1; title = true)
    if title
		return ["x/|x|" "y/|y|" "z/|z|" "i"; B*unk*α]		
	else
		return B*unk*α
	end
end

function flatProjection(A, x̂ = [cos(π/4) -cos(π/4)]'*0.5; x :: Bool = true)
	if x
		return [A[2:3,1]+x̂*A[1,1] A[2:3,2]+x̂*A[1,2] A[2:3,3]+x̂*A[1,3] A[2:3,4]+x̂*A[1,4]]
	else
		return [A[1:2,1]+x̂*A[3,1] A[1:2,2]+x̂*A[3,2] A[1:2,3]+x̂*A[3,3] A[1:2,4]+x̂*A[3,4]]
	end
end	

function plotFlatProjection(B; sprain=1.01, static = false)
	param = [B[1,4]+B[2,4]im, B[1,4]+B[2,4]im+B[1,1]+B[2,1]im, Inf, B[1,4]+B[2,4]im, B[1,4]+B[2,4]im+B[1,2]+B[2,2]im, Inf, B[1,4]+B[2,4]im, B[1,4]+B[2,4]im+B[1,3]+B[2,3]im]
	if static
		plot!(param, arrow = 2)
	else
		plot(param, arrow = 2)
	end
	annotate!(B[1,4]*sprain, B[2,4]*sprain, text("v", :black, :right, 10))
	annotate!(B[1,4]+B[1,1]*sprain, B[2,4]+B[2,1]*sprain, text("x", :black, :right, 10))
	annotate!(B[1,4]+B[1,2]*sprain, B[2,4]+B[2,2]*sprain, text("y", :black, :right, 10))
	annotate!(B[1,4]+B[1,3]*sprain, B[2,4]+B[2,3]*sprain, text("z", :black, :right, 10))
	# ylims!((-1.5,2)); xlims!((-1.5,1.5))
end

print("DmdgpRealization included\n")

#rotationFactor = Rz(π/2)*Rx(π/2)*Rz(-θ(30))*Rx(π)
#sysGenerate(flatProjection(rotationFactor*B2*B3*Bi(ω(230),θ(120))), 20)
#Bo = B2*B3*Bi(ω(230),θ(120))
#B = Bi(ω(230),θ(120))
#sysGenerate(flatProjection(rotationFactor*Bo*B*B*B), 20)